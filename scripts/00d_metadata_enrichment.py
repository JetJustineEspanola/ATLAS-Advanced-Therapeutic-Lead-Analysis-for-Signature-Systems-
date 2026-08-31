#!/usr/bin/env python3
"""
ATLAS — Stage 00D
Metadata Enrichment + Download Planning

Purpose
-------
Take 00A-00C discovery candidates and enrich them at study/sample/file level
before downloading large data.

Current adapters
----------------
GEO:
  - download one GEO family SOFT file per GSE
  - parse Series and all GSM sample metadata locally
  - infer sample phenotype fields conservatively

SRA:
  - resolve SRR/SRX/SRP accessions through NCBI E-utilities
  - fetch RunInfo CSV
  - populate run/file planning fields

Outputs
-------
DuckDB:
  samples
  files
  datasets (enriched flags / scores)

Files:
  data/enriched/sample_metadata.csv
  data/enriched/download_manifest.csv
  data/enriched/dataset_enrichment_summary.csv

Important
---------
Phenotype inference is conservative. Any ambiguous resistant/sensitive,
treatment, clone, or replicate field is marked unresolved and sent to manual
review instead of being guessed.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import csv
import gzip
import io
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

import duckdb
import pandas as pd
import requests


CATALOG = PROJECT_ROOT / "data" / "catalog" / "atlas_catalog.duckdb"
OUT_DIR = PROJECT_ROOT / "data" / "enriched"
CACHE_DIR = PROJECT_ROOT / "data" / "cache" / "metadata"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def clean(x: Any) -> str:
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()


def get(url: str, params=None, timeout=(5, 30), retries=2) -> requests.Response:
    last = None
    for attempt in range(retries + 1):
        try:
            r = requests.get(
                url,
                params=params,
                timeout=timeout,
                headers={"User-Agent": "ATLAS-00D/0.1"},
            )
            r.raise_for_status()
            return r
        except Exception as e:
            last = e
            if attempt < retries:
                time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET failed: {url}: {last}")


def cached_text(key: str, url: str, params=None) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", key)
    path = CACHE_DIR / f"{safe}.txt"
    if path.exists() and path.stat().st_size > 20:
        return path.read_text(encoding="utf-8", errors="replace")
    text = get(url, params=params).text
    path.write_text(text, encoding="utf-8")
    return text


def parse_soft(text: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for line in text.splitlines():
        if not line.startswith("!"):
            continue
        if " = " not in line:
            continue
        key, value = line[1:].split(" = ", 1)
        out.setdefault(key.strip(), []).append(value.strip())
    return out


def first(meta: dict[str, list[str]], *keys: str) -> str:
    for key in keys:
        vals = meta.get(key, [])
        if vals:
            return clean(vals[0])
    return ""


def joined(meta: dict[str, list[str]], *keys: str) -> str:
    vals = []
    for key in keys:
        vals.extend(meta.get(key, []))
    return " | ".join(clean(x) for x in vals if clean(x))


def phenotype_from_text(text: str) -> dict[str, Any]:
    t = clean(text)
    low = t.lower()

    resistant = bool(re.search(
        r"\b(resistant|resistance|refractory|adapted)\b", low
    ))
    sensitive = bool(re.search(
        r"\b(sensitive|parental|wild[- ]?type|control)\b", low
    ))

    if resistant and not sensitive:
        resistance_status = "RESISTANT"
    elif sensitive and not resistant:
        resistance_status = "SENSITIVE_OR_PARENTAL"
    elif resistant and sensitive:
        resistance_status = "AMBIGUOUS"
    else:
        resistance_status = "UNRESOLVED"

    trastuzumab = any(x in low for x in ["trastuzumab", "herceptin"])

    # Common HER2+ cell lines.
    cell_lines = []
    for label, patterns in {
        "BT474": ["bt474", "bt-474", "bt 474"],
        "SKBR3": ["skbr3", "sk-br-3", "sk br 3"],
        "HCC1954": ["hcc1954"],
        "AU565": ["au565", "au-565"],
        "MDA-MB-453": ["mda-mb-453", "mda mb 453"],
    }.items():
        if any(p in low for p in patterns):
            cell_lines.append(label)

    cell_line = " | ".join(cell_lines)

    clone_match = re.search(
        r"\b(?:clone|cl)[\s_-]*([a-z0-9.-]+)\b",
        low,
        flags=re.I,
    )
    clone_id = clone_match.group(1) if clone_match else ""

    # Technical replicate language is only marked when explicitly stated.
    if re.search(r"\btechnical replicate\b|\brna prep\b", low):
        replicate_type = "TECHNICAL"
    elif re.search(r"\bbiological replicate\b", low):
        replicate_type = "BIOLOGICAL"
    else:
        replicate_type = "UNRESOLVED"

    # Generic replicate number if explicitly present.
    rep_match = re.search(
        r"\b(?:replicate|rep|r)[\s_-]*([0-9]+)\b",
        low,
        flags=re.I,
    )
    replicate_id = rep_match.group(1) if rep_match else ""

    return {
        "resistance_status": resistance_status,
        "trastuzumab_exposure_text_flag": trastuzumab,
        "cell_line": cell_line,
        "clone_id": clone_id,
        "replicate_type": replicate_type,
        "replicate_id": replicate_id,
    }


def geo_series_bucket(accession: str) -> str:
    """
    GEO FTP bucket convention:
      GSE132055 -> GSE132nnn
      GSE89216  -> GSE89nnn
    """
    m = re.fullmatch(r"(GSE)(\d+)", accession.upper())
    if not m:
        raise ValueError(f"Invalid GEO Series accession: {accession}")

    digits = m.group(2)
    if len(digits) <= 3:
        prefix = ""
    else:
        prefix = digits[:-3]

    return f"GSE{prefix}nnn"


def fetch_geo_family_soft(accession: str) -> str:
    """
    Download the GEO family SOFT file.

    Using the family SOFT is more robust than the accession viewer for sample
    enumeration because it contains the Series section and all Sample sections
    in one machine-readable file.
    """
    bucket = geo_series_bucket(accession)
    cache_path = CACHE_DIR / f"GEO_{accession}_family.soft"

    if cache_path.exists() and cache_path.stat().st_size > 100:
        return cache_path.read_text(
            encoding="utf-8",
            errors="replace",
        )

    url = (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{bucket}/{accession}/soft/{accession}_family.soft.gz"
    )

    r = get(url, timeout=(5, 60), retries=2)

    try:
        text = gzip.decompress(r.content).decode(
            "utf-8",
            errors="replace",
        )
    except Exception as e:
        raise RuntimeError(
            f"Could not decompress GEO family SOFT for {accession}: {e}"
        )

    cache_path.write_text(
        text,
        encoding="utf-8",
    )

    return text


def split_geo_family_soft(
    text: str,
) -> tuple[dict[str, list[str]], dict[str, dict[str, list[str]]]]:
    """
    Parse one GEO family SOFT document into:
      series metadata
      {GSM: sample metadata}
    """
    series_meta: dict[str, list[str]] = {}
    samples: dict[str, dict[str, list[str]]] = {}

    current_kind = ""
    current_id = ""
    current_meta: dict[str, list[str]] | None = None

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")

        if line.startswith("^SERIES"):
            current_kind = "SERIES"
            current_id = line.split("=", 1)[1].strip()
            current_meta = series_meta
            continue

        if line.startswith("^SAMPLE"):
            current_kind = "SAMPLE"
            current_id = line.split("=", 1)[1].strip()
            current_meta = samples.setdefault(
                current_id,
                {},
            )
            continue

        # Skip platform/table sections until the next SAMPLE/SERIES section.
        if line.startswith("^PLATFORM") or line.startswith("^DATABASE"):
            current_kind = "OTHER"
            current_id = ""
            current_meta = None
            continue

        if (
            current_meta is not None
            and line.startswith("!")
            and " = " in line
        ):
            key, value = line[1:].split(" = ", 1)
            current_meta.setdefault(
                key.strip(),
                [],
            ).append(
                value.strip()
            )

    return series_meta, samples


def geo_sample_from_meta(
    dataset_id: str,
    gsm: str,
    meta: dict[str, list[str]],
) -> dict[str, Any]:
    title = first(meta, "Sample_title")
    source_name = first(meta, "Sample_source_name_ch1")
    characteristics = joined(meta, "Sample_characteristics_ch1")
    treatment = joined(meta, "Sample_treatment_protocol_ch1")
    growth = joined(meta, "Sample_growth_protocol_ch1")
    description = joined(meta, "Sample_description")
    organism = first(meta, "Sample_organism_ch1")
    library_strategy = first(meta, "Sample_library_strategy")
    library_source = first(meta, "Sample_library_source")
    platform_id = first(meta, "Sample_platform_id")

    combined = " | ".join(
        x for x in [
            title,
            source_name,
            characteristics,
            treatment,
            growth,
            description,
        ]
        if x
    )

    pheno = phenotype_from_text(combined)

    relation = joined(meta, "Sample_relation")

    # GEO Sample relations may reference SRA Experiment, BioSample, or SRA run.
    sra_rel = ""
    m = re.search(
        r"\b(SR[AXRP]\d+)\b",
        relation,
    )
    if m:
        sra_rel = m.group(1)

    return {
        "dataset_id": dataset_id,
        "sample_id": gsm,
        "title": title,
        "source_name": source_name,
        "characteristics": characteristics,
        "treatment": treatment,
        "description": description,
        "organism": organism,
        "platform_id": platform_id,
        "library_strategy": library_strategy,
        "library_source": library_source,
        "biological_group": pheno["resistance_status"],
        "resistance_status": pheno["resistance_status"],
        "trastuzumab_exposure_text_flag": pheno[
            "trastuzumab_exposure_text_flag"
        ],
        "cell_line": pheno["cell_line"],
        "patient_id": "",
        "clone_id": pheno["clone_id"],
        "replicate_id": pheno["replicate_id"],
        "replicate_type": pheno["replicate_type"],
        "sra_relation": sra_rel,
        "qc_status": "NOT_RUN",
        "qc_reason": "",
        "metadata_review_required": (
            pheno["resistance_status"]
            in {"UNRESOLVED", "AMBIGUOUS"}
            or pheno["replicate_type"] == "UNRESOLVED"
        ),
    }


def fetch_geo_gse(
    accession: str,
) -> tuple[dict[str, Any], dict[str, dict[str, list[str]]]]:
    text = fetch_geo_family_soft(accession)
    series_meta, samples = split_geo_family_soft(text)

    study = {
        "title": first(series_meta, "Series_title"),
        "summary": joined(
            series_meta,
            "Series_summary",
            "Series_overall_design",
        ),
        "overall_design": joined(
            series_meta,
            "Series_overall_design",
        ),
        "sample_n": len(samples),
    }

    return study, samples


def sra_uid(accession: str) -> str:
    js = get(
        f"{EUTILS}/esearch.fcgi",
        params={
            "db": "sra",
            "term": f"{accession}[Accession]",
            "retmode": "json",
            "retmax": 5,
            "tool": "ATLAS",
        },
    ).json()
    ids = js.get("esearchresult", {}).get("idlist", [])
    return ids[0] if ids else ""


def fetch_sra_runinfo(dataset_id: str, accession: str) -> list[dict[str, Any]]:
    uid = sra_uid(accession)
    if not uid:
        return []

    text = cached_text(
        f"SRA_{accession}_runinfo",
        f"{EUTILS}/efetch.fcgi",
        params={
            "db": "sra",
            "id": uid,
            "rettype": "runinfo",
            "retmode": "text",
            "tool": "ATLAS",
        },
    )

    df = pd.read_csv(io.StringIO(text))
    rows = []

    for _, r in df.iterrows():
        run = clean(r.get("Run"))
        sample = clean(r.get("Sample"))
        experiment = clean(r.get("Experiment"))
        biosample = clean(r.get("BioSample"))
        library = clean(r.get("LibraryStrategy"))
        platform = clean(r.get("Platform"))
        scientific_name = clean(r.get("ScientificName"))
        bases = pd.to_numeric(
            pd.Series([r.get("bases")]), errors="coerce"
        ).iloc[0]
        size = pd.to_numeric(
            pd.Series([r.get("size_MB")]), errors="coerce"
        ).iloc[0]

        sample_text = " | ".join(
            x for x in [
                clean(r.get("SampleName")),
                clean(r.get("LibraryName")),
                clean(r.get("CenterName")),
            ] if x
        )
        pheno = phenotype_from_text(sample_text)

        rows.append({
            "dataset_id": dataset_id,
            "sample_id": sample or run,
            "title": sample_text,
            "source_name": "",
            "characteristics": "",
            "treatment": "",
            "description": "",
            "organism": scientific_name,
            "platform_id": platform,
            "library_strategy": library,
            "library_source": clean(r.get("LibrarySource")),
            "biological_group": pheno["resistance_status"],
            "resistance_status": pheno["resistance_status"],
            "trastuzumab_exposure_text_flag": pheno[
                "trastuzumab_exposure_text_flag"
            ],
            "cell_line": pheno["cell_line"],
            "patient_id": "",
            "clone_id": pheno["clone_id"],
            "replicate_id": pheno["replicate_id"],
            "replicate_type": pheno["replicate_type"],
            "sra_relation": run,
            "qc_status": "NOT_RUN",
            "qc_reason": "",
            "metadata_review_required": True,
            "_run": run,
            "_experiment": experiment,
            "_biosample": biosample,
            "_bases": None if pd.isna(bases) else int(bases),
            "_size_mb": None if pd.isna(size) else float(size),
        })

    return rows


def build_file_rows(sample_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    files = []
    seen = set()

    for r in sample_rows:
        run = clean(r.get("_run")) or clean(r.get("sra_relation"))
        if not re.fullmatch(r"SRR\d+", run):
            continue

        key = (r["dataset_id"], run)
        if key in seen:
            continue
        seen.add(key)

        size_mb = r.get("_size_mb")
        size_bytes = (
            int(float(size_mb) * 1024 * 1024)
            if size_mb not in [None, ""]
            else None
        )

        files.append({
            "dataset_id": r["dataset_id"],
            "file_id": run,
            "file_name": f"{run}.sra",
            "file_type": "SRA_RUN",
            "local_path": str(
                PROJECT_ROOT / "data" / "raw" / "sra" / run
            ),
            "remote_url": (
                "SRA_TOOLKIT:" + run
            ),
            "size_bytes": size_bytes,
            "md5": "",
            "sha256": "",
            "download_status": "PLANNED",
            "source_release": "",
        })

    return files


def ensure_schema(con: duckdb.DuckDBPyConnection) -> None:
    # Add richer sample columns if the original 00B schema is already present.
    sample_cols = {
        "title": "VARCHAR",
        "source_name": "VARCHAR",
        "characteristics": "VARCHAR",
        "description": "VARCHAR",
        "organism": "VARCHAR",
        "platform_id": "VARCHAR",
        "library_strategy": "VARCHAR",
        "library_source": "VARCHAR",
        "trastuzumab_exposure_text_flag": "BOOLEAN",
        "sra_relation": "VARCHAR",
        "metadata_review_required": "BOOLEAN",
    }

    existing = {
        r[1] for r in con.execute(
            "PRAGMA table_info('samples')"
        ).fetchall()
    }

    for col, typ in sample_cols.items():
        if col not in existing:
            con.execute(
                f'ALTER TABLE samples ADD COLUMN "{col}" {typ}'
            )


def update_dataset_flags(
    con: duckdb.DuckDBPyConnection,
    dataset_id: str,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        return

    rdf = pd.DataFrame(rows)

    statuses = set(
        rdf["resistance_status"].dropna().astype(str)
    )

    groups_defined = (
        "RESISTANT" in statuses
        and "SENSITIVE_OR_PARENTAL" in statuses
    )

    # Biological replication is only considered established when at least
    # two samples per phenotype are present AND they are not explicitly
    # labeled technical replicates.
    bio_rep = False
    if groups_defined:
        ok = rdf[
            ~rdf["replicate_type"].eq("TECHNICAL")
        ]
        counts = ok["resistance_status"].value_counts()
        bio_rep = (
            counts.get("RESISTANT", 0) >= 2
            and counts.get("SENSITIVE_OR_PARENTAL", 0) >= 2
        )

    complete_metadata = bool(
        rdf["title"].fillna("").astype(str).str.len().gt(0).all()
        and rdf["organism"].fillna("").astype(str).str.len().gt(0).all()
    )

    phenotype_confidence = (
        "HIGH"
        if groups_defined and bio_rep
        else "MODERATE"
        if groups_defined
        else "LOW"
    )

    con.execute(
        """
        UPDATE datasets
        SET
            resistant_sensitive_groups_defined = ?,
            biological_replication = ?,
            complete_sample_metadata = ?,
            phenotype_confidence = ?,
            sample_count = ?
        WHERE dataset_id = ?
        """,
        [
            bool(groups_defined),
            bool(bio_rep),
            bool(complete_metadata),
            phenotype_confidence,
            int(len(rdf)),
            str(dataset_id),
        ],
    )


def write_samples(con, sample_df: pd.DataFrame) -> None:
    if sample_df.empty:
        return

    sample_df = sample_df.copy()
    for col in [
        "trastuzumab_exposure_text_flag",
        "metadata_review_required",
    ]:
        if col in sample_df.columns:
            sample_df[col] = sample_df[col].map(
                lambda x: None if pd.isna(x) else bool(x)
            )

    public_cols = [
        c for c in sample_df.columns if not c.startswith("_")
    ]
    sample_df = sample_df[public_cols].copy()

    # Replace dataset rows to keep reruns idempotent.
    for dataset_id in sample_df["dataset_id"].unique():
        con.execute(
            "DELETE FROM samples WHERE dataset_id = ?",
            [dataset_id],
        )

    con.register("new_samples", sample_df)

    table_cols = [
        r[1] for r in con.execute(
            "PRAGMA table_info('samples')"
        ).fetchall()
    ]
    cols = [c for c in table_cols if c in sample_df.columns]

    con.execute(
        f"""
        INSERT INTO samples ({", ".join(cols)})
        SELECT {", ".join('"' + c + '"' for c in cols)}
        FROM new_samples
        """
    )


def write_files(con, files_df: pd.DataFrame) -> None:
    if files_df.empty:
        return

    for dataset_id in files_df["dataset_id"].unique():
        con.execute(
            "DELETE FROM files WHERE dataset_id = ?",
            [dataset_id],
        )

    con.register("new_files", files_df)

    table_cols = [
        r[1] for r in con.execute(
            "PRAGMA table_info('files')"
        ).fetchall()
    ]
    cols = [c for c in table_cols if c in files_df.columns]

    con.execute(
        f"""
        INSERT INTO files ({", ".join(cols)})
        SELECT {", ".join('"' + c + '"' for c in cols)}
        FROM new_files
        """
    )


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--min-score",
        type=int,
        default=40,
        help="Only enrich datasets with 00C score >= this value.",
    )
    p.add_argument(
        "--sources",
        default="GEO,SRA",
    )
    p.add_argument(
        "--accessions",
        default="",
        help="Optional comma-separated source accessions to enrich explicitly, e.g. GSE114575,GSE121105.",
    )
    p.add_argument(
        "--max-datasets",
        type=int,
        default=30,
    )
    p.add_argument(
        "--workers",
        type=int,
        default=6,
    )
    args = p.parse_args()

    sources = [
        s.strip().upper()
        for s in args.sources.split(",")
        if s.strip()
    ]

    if not CATALOG.exists():
        print(f"ERROR: catalog not found: {CATALOG}")
        return 1

    con = duckdb.connect(str(CATALOG))
    ensure_schema(con)

    qmarks = ",".join(["?"] * len(sources))
    explicit_accessions = [
        x.strip().upper()
        for x in args.accessions.split(",")
        if x.strip()
    ]

    if explicit_accessions:
        aq = ",".join(["?"] * len(explicit_accessions))
        datasets = con.execute(
            f"""
            SELECT *
            FROM datasets
            WHERE upper(source) IN ({qmarks})
              AND (
                    eligibility_score >= ?
                    OR upper(source_accession) IN ({aq})
                  )
            ORDER BY
              CASE WHEN upper(source_accession) IN ({aq}) THEN 0 ELSE 1 END,
              eligibility_score DESC,
              source,
              source_accession
            LIMIT ?
            """,
            [*sources, args.min_score, *explicit_accessions, *explicit_accessions, args.max_datasets],
        ).fetchdf()
    else:
        datasets = con.execute(
            f"""
            SELECT *
            FROM datasets
            WHERE eligibility_score >= ?
              AND upper(source) IN ({qmarks})
            ORDER BY eligibility_score DESC, source, source_accession
            LIMIT ?
            """,
            [args.min_score, *sources, args.max_datasets],
        ).fetchdf()

    if datasets.empty:
        print("No datasets matched enrichment criteria.")
        con.close()
        return 0

    print("=" * 78)
    print("ATLAS — Stage 00D Metadata Enrichment + Download Planning")
    print("=" * 78)
    print(f"Datasets selected: {len(datasets)}")
    print(f"Sources: {sources}")
    print(f"Workers: {args.workers}")

    all_samples: list[dict[str, Any]] = []
    enrichment_rows = []

    for _, ds in datasets.iterrows():
        dataset_id = clean(ds["dataset_id"])
        source = clean(ds["source"]).upper()
        accession = clean(ds["source_accession"])

        print(f"\n[{source}] {accession}", flush=True)

        try:
            rows: list[dict[str, Any]] = []

            if source == "GEO":
                study, sample_meta = fetch_geo_gse(
                    accession
                )

                print(
                    f"  family SOFT samples discovered: {len(sample_meta)}",
                    flush=True,
                )

                # All GSM metadata is already present in the family SOFT,
                # so no per-GSM network requests are needed.
                for gsm, meta in sample_meta.items():
                    try:
                        rows.append(
                            geo_sample_from_meta(
                                dataset_id,
                                gsm,
                                meta,
                            )
                        )
                    except Exception as e:
                        print(
                            f"  {gsm}: ERROR {e}",
                            flush=True,
                        )

            elif source == "SRA":
                rows = fetch_sra_runinfo(
                    dataset_id,
                    accession,
                )

            else:
                continue

            all_samples.extend(rows)
            update_dataset_flags(
                con,
                dataset_id,
                rows,
            )

            rdf = pd.DataFrame(rows)

            if rdf.empty:
                resistant_n = sensitive_n = unresolved_n = 0
            else:
                counts = rdf[
                    "resistance_status"
                ].value_counts()
                resistant_n = int(counts.get("RESISTANT", 0))
                sensitive_n = int(
                    counts.get("SENSITIVE_OR_PARENTAL", 0)
                )
                unresolved_n = int(
                    counts.get("UNRESOLVED", 0)
                    + counts.get("AMBIGUOUS", 0)
                )

            enrichment_rows.append({
                "dataset_id": dataset_id,
                "source": source,
                "source_accession": accession,
                "sample_rows": len(rows),
                "resistant_n": resistant_n,
                "sensitive_or_parental_n": sensitive_n,
                "unresolved_n": unresolved_n,
                "enrichment_status": (
                    "ENRICHED" if rows else "NO_SAMPLE_METADATA"
                ),
            })

            print(
                f"  enriched samples: {len(rows)} "
                f"(R={resistant_n}, S/P={sensitive_n}, unresolved={unresolved_n})",
                flush=True,
            )

        except Exception as e:
            enrichment_rows.append({
                "dataset_id": dataset_id,
                "source": source,
                "source_accession": accession,
                "sample_rows": 0,
                "resistant_n": 0,
                "sensitive_or_parental_n": 0,
                "unresolved_n": 0,
                "enrichment_status": "ERROR: " + str(e)[:300],
            })
            print(f"  ERROR: {e}", flush=True)

    sample_df = pd.DataFrame(all_samples)
    file_rows = build_file_rows(all_samples)
    files_df = pd.DataFrame(file_rows)

    write_samples(con, sample_df)
    write_files(con, files_df)

    # Recompute a practical post-enrichment quality view.
    enriched_datasets = con.execute(
        """
        SELECT *
        FROM datasets
        ORDER BY
          resistant_sensitive_groups_defined DESC,
          biological_replication DESC,
          eligibility_score DESC,
          source,
          source_accession
        """
    ).fetchdf()

    con.close()

    sample_out = OUT_DIR / "sample_metadata.csv"
    manifest_out = OUT_DIR / "download_manifest.csv"
    summary_out = OUT_DIR / "dataset_enrichment_summary.csv"

    # IMPORTANT: export the complete catalog state, not only the current run.
    # This preserves previously enriched datasets when later doing targeted
    # enrichment of low-score but highly relevant accessions.
    con_export = duckdb.connect(str(CATALOG))
    all_samples_df = con_export.execute(
        "SELECT * FROM samples ORDER BY dataset_id, sample_id"
    ).fetchdf()
    all_files_df = con_export.execute(
        "SELECT * FROM files ORDER BY dataset_id, file_id"
    ).fetchdf()
    con_export.close()

    all_samples_df.to_csv(sample_out, index=False)
    all_files_df.to_csv(manifest_out, index=False)
    pd.DataFrame(enrichment_rows).to_csv(
        summary_out,
        index=False,
    )

    try:
        all_samples_df.to_parquet(
            OUT_DIR / "sample_metadata.parquet",
            index=False,
        )
        all_files_df.to_parquet(
            OUT_DIR / "download_manifest.parquet",
            index=False,
        )
    except Exception:
        pass

    print("\n" + "=" * 78)
    print("STAGE 00D COMPLETE")
    print("=" * 78)
    print(f"Sample metadata rows (catalog total): {len(all_samples_df)}")
    print(f"Planned SRA run files (catalog total): {len(all_files_df)}")
    print(f"\n{sample_out}")
    print(manifest_out)
    print(summary_out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
