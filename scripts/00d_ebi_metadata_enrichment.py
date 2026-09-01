#!/usr/bin/env python3
"""
ATLAS 00D-EBI — BioStudies / ArrayExpress metadata enrichment

Purpose
-------
Enrich BIOSTUDIES_ARRAYEXPRESS discovery leads without pretending that
cross-database mirrors are independent validation datasets.

What it does
------------
1. Reads data/discovery/ebi_external_candidates.csv.
2. Queries BioStudies study/info/files APIs.
3. Detects GEO mirrors from the cross-reference columns.
4. For E-MTAB-only studies, finds and parses SDRF sample metadata when available.
5. Writes rich EBI metadata to separate CSVs.
6. Inserts only non-GEO-mirror EBI samples into the existing DuckDB `samples` table.
7. Leaves resistance_status as UNRESOLVED. Existing ATLAS phenotype stages must decide it.
8. Updates conservative dataset metadata (sample_count, metadata_complete_flag, etc.).

Outputs
-------
data/enriched/ebi_dataset_enrichment_summary.csv
data/enriched/ebi_sample_metadata.csv
data/enriched/ebi_file_inventory.csv

Scientific guardrails
---------------------
- E-GEOD records linked to an existing GSE are marked CROSS_DATABASE_MIRROR_OF_GEO
  and are NOT inserted as independent sample sets.
- Discovery relevance scores do not become phenotype calls.
- No sample is automatically labeled resistant/sensitive here.
- A successful SDRF parse means metadata was retrieved, not that the dataset is
  suitable for primary validation.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DISC = ROOT / "data" / "discovery"
ENRICHED = ROOT / "data" / "enriched"
CATALOG = ROOT / "data" / "catalog" / "atlas_catalog.duckdb"

CANDIDATES = DISC / "ebi_external_candidates.csv"
OUT_SUMMARY = ENRICHED / "ebi_dataset_enrichment_summary.csv"
OUT_SAMPLES = ENRICHED / "ebi_sample_metadata.csv"
OUT_FILES = ENRICHED / "ebi_file_inventory.csv"

BASE = "https://www.ebi.ac.uk/biostudies/api/v1"
UA = "ATLAS-EBI-metadata-enrichment/1.0"

GEO_RE = re.compile(r"\bGSE\d+\b", re.I)
ENA_RUN_RE = re.compile(r"\b[SED]RR\d+\b", re.I)
ENA_SAMPLE_RE = re.compile(r"\b[SED]RS\d+\b", re.I)
BIOSAMPLE_RE = re.compile(r"\b(?:SAMEA|SAMN|SAMD)\d+\b", re.I)

def get_bytes(url: str, tries: int = 3) -> bytes:
    last = None
    for attempt in range(tries):
        try:
            req = Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
            with urlopen(req, timeout=60) as r:
                return r.read()
        except (HTTPError, URLError, TimeoutError) as e:
            last = e
            if attempt + 1 < tries:
                time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET failed: {url}: {last}")

def get_json(url: str) -> dict:
    return json.loads(get_bytes(url).decode("utf-8", errors="replace"))

def flatten_strings(obj):
    vals = []
    def walk(x):
        if isinstance(x, dict):
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
        elif isinstance(x, (str, int, float)):
            s = str(x).strip()
            if s:
                vals.append(s)
    walk(obj)
    return vals

def collect_file_nodes(obj):
    """Recursively collect dicts that look like BioStudies file records."""
    out = []
    def walk(x):
        if isinstance(x, dict):
            keys = {str(k).lower() for k in x}
            if keys & {"path", "file", "filename", "name"}:
                candidate = (
                    x.get("path") or x.get("file") or x.get("fileName")
                    or x.get("filename") or x.get("name")
                )
                if isinstance(candidate, str) and "." in candidate:
                    out.append(x)
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(obj)
    return out

def file_path(node):
    for k in ("path", "file", "fileName", "filename", "name"):
        v = node.get(k)
        if isinstance(v, str) and v.strip():
            return v.strip()
    return ""

def int_or_none(v):
    try:
        return int(v)
    except Exception:
        return None

def build_download_url(info: dict, rel_file: str) -> str:
    """
    BioStudies info returns ftpLink to study root. Files live under /Files/.
    Prefer HTTPS transport via ftp.ebi.ac.uk.
    """
    ftp = str(info.get("ftpLink") or "").rstrip("/")
    rel_file = rel_file.lstrip("/")
    if ftp:
        root = ftp.replace("ftp://ftp.ebi.ac.uk", "https://ftp.ebi.ac.uk")
        if "/Files/" in rel_file:
            return root + "/" + quote(rel_file, safe="/()[]_.+-")
        return root + "/Files/" + quote(rel_file, safe="/()[]_.+-")

    rel_path = str(info.get("relPath") or "").strip("/")
    if rel_path:
        mode = str(info.get("storageMode") or info.get("storage") or "fire").strip("/")
        return (
            f"https://ftp.ebi.ac.uk/biostudies/{mode}/{rel_path}/Files/"
            + quote(rel_file, safe="/()[]_.+-")
        )
    return ""

def norm_col(c):
    return re.sub(r"\s+", " ", str(c).strip())

def find_col(cols, patterns):
    for c in cols:
        lc = c.lower()
        if any(p in lc for p in patterns):
            return c
    return None

def first_nonempty(row, cols):
    for c in cols:
        if c and c in row:
            v = str(row[c]).strip()
            if v and v.lower() not in {"nan", "none", "na", "n/a"}:
                return v
    return ""

def join_design_fields(row, columns):
    preferred = []
    for c in columns:
        lc = c.lower()
        if (
            "factor value" in lc
            or "characteristic" in lc
            or "treatment" in lc
            or "compound" in lc
            or "stimulus" in lc
            or "disease" in lc
            or "phenotype" in lc
            or "genotype" in lc
        ):
            v = str(row.get(c, "")).strip()
            if v and v.lower() not in {"nan", "none", "na", "n/a"}:
                preferred.append(f"{c}={v}")
    return " | ".join(preferred[:20])

def parse_sdrf(accession: str, text: str):
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str, keep_default_na=False)
    df.columns = [norm_col(c) for c in df.columns]
    cols = list(df.columns)

    source_col = find_col(cols, ["source name"])
    sample_col = find_col(cols, ["sample name"])
    assay_col = find_col(cols, ["assay name"])
    cell_col = find_col(cols, ["characteristics[cell line", "cell line"])
    organism_col = find_col(cols, ["characteristics[organism", "organism"])
    biosample_col = find_col(cols, ["biosample", "bio sample"])
    ena_sample_col = find_col(cols, ["ena_sample", "ena sample"])
    ena_run_col = find_col(cols, ["ena_run", "ena run"])

    rows = []
    seen = set()
    for idx, r in df.iterrows():
        d = r.to_dict()
        sample_id = first_nonempty(
            d, [source_col, sample_col, assay_col]
        ) or f"{accession}_row_{idx+1}"

        # SDRFs often have several assay rows per biological source. Keep one
        # biological row per source/sample/design combination.
        design = join_design_fields(d, cols)
        dedupe_key = (sample_id, design)
        if dedupe_key in seen:
            continue
        seen.add(dedupe_key)

        blob = " ".join(str(v) for v in d.values())
        bios = first_nonempty(d, [biosample_col])
        if not bios:
            m = BIOSAMPLE_RE.search(blob)
            bios = m.group(0).upper() if m else ""
        ena_sample = first_nonempty(d, [ena_sample_col])
        if not ena_sample:
            m = ENA_SAMPLE_RE.search(blob)
            ena_sample = m.group(0).upper() if m else ""
        ena_run = first_nonempty(d, [ena_run_col])
        if not ena_run:
            runs = sorted({x.upper() for x in ENA_RUN_RE.findall(blob)})
            ena_run = "|".join(runs)

        cell_line = first_nonempty(d, [cell_col])
        organism = first_nonempty(d, [organism_col])

        rows.append({
            "dataset_id": f"BIOSTUDIES:{accession}",
            "source_accession": accession,
            "sample_id": sample_id,
            "biological_group": design or sample_id,
            "resistance_status": "UNRESOLVED",
            "treatment": design,
            "cell_line": cell_line,
            "patient_id": "",
            "clone_id": "",
            "replicate_id": "",
            "replicate_type": "",
            "qc_status": "METADATA_ONLY",
            "qc_reason": "EBI SDRF imported; phenotype requires ATLAS audit",
            "organism": organism,
            "biosample_accession": bios,
            "ena_sample_accession": ena_sample,
            "ena_run_accessions": ena_run,
            "raw_design_text": design,
        })
    return rows, len(df)

def choose_sdrf(paths):
    scored = []
    for p in paths:
        lp = p.lower()
        if "sdrf" not in lp:
            continue
        score = 0
        if lp.endswith(".txt"): score += 3
        if lp.endswith(".tsv"): score += 3
        if "files/" in lp: score += 1
        scored.append((score, p))
    return sorted(scored, reverse=True)[0][1] if scored else ""

def existing_geo_accessions():
    if not CATALOG.exists():
        return set()
    con = duckdb.connect(str(CATALOG), read_only=True)
    try:
        vals = con.execute(
            "SELECT upper(source_accession) FROM datasets WHERE upper(source)='GEO'"
        ).fetchall()
        return {r[0] for r in vals if r and r[0]}
    finally:
        con.close()

def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        pd.DataFrame().to_csv(path, index=False)
        return
    pd.DataFrame(rows).to_csv(path, index=False)

def update_duckdb(summary_rows, sample_rows):
    con = duckdb.connect(str(CATALOG))
    try:
        con.execute("BEGIN")

        # E-GEOD records are GEO mirrors and must never persist as
        # independent BioStudies sample rows.
        con.execute("""
            DELETE FROM samples
            WHERE dataset_id LIKE 'BIOSTUDIES:E-GEOD-%'
        """)

        # Idempotent: replace only EBI sample rows processed by this stage.
        ids = sorted({r["dataset_id"] for r in sample_rows})
        for did in ids:
            con.execute("DELETE FROM samples WHERE dataset_id = ?", [did])

        if sample_rows:
            sample_df = pd.DataFrame(sample_rows)
            table_cols = [r[1] for r in con.execute("PRAGMA table_info('samples')").fetchall()]
            common = [c for c in table_cols if c in sample_df.columns]
            con.register("ebi_samples_df", sample_df)
            ins = ", ".join(common)
            sel = ", ".join(f'"{c}"' for c in common)
            con.execute(f"INSERT INTO samples ({ins}) SELECT {sel} FROM ebi_samples_df")

        # Conservative dataset metadata update.
        for s in summary_rows:
            if s["relationship_role"] == "CROSS_DATABASE_MIRROR_OF_GEO":
                con.execute(
                    """UPDATE datasets
                       SET metadata_complete_flag = TRUE,
                           eligibility_reasons = COALESCE(eligibility_reasons,'') ||
                             CASE WHEN COALESCE(eligibility_reasons,'')='' THEN '' ELSE ' | ' END ||
                             ?
                       WHERE dataset_id = ?""",
                    [
                        "EBI cross-database GEO mirror; do not count independently",
                        s["dataset_id"],
                    ],
                )
                continue

            if s["status"] == "ENRICHED_SDRF":
                con.execute(
                    """UPDATE datasets
                       SET sample_count = ?,
                           metadata_complete_flag = TRUE,
                           complete_sample_metadata = TRUE,
                           phenotype_confidence = COALESCE(phenotype_confidence, 'UNRESOLVED')
                       WHERE dataset_id = ?""",
                    [int(s["unique_samples"]), s["dataset_id"]],
                )
            else:
                con.execute(
                    """UPDATE datasets
                       SET metadata_complete_flag = FALSE,
                           complete_sample_metadata = FALSE,
                           manual_review_required = TRUE
                       WHERE dataset_id = ?""",
                    [s["dataset_id"]],
                )

        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    finally:
        con.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", nargs="*", default=None,
                    help="optional E-MTAB/E-GEOD accessions to enrich")
    ap.add_argument("--sleep", type=float, default=0.15)
    args = ap.parse_args()

    if not CANDIDATES.exists():
        print(f"ERROR: missing {CANDIDATES}")
        return 1
    if not CATALOG.exists():
        print(f"ERROR: missing DuckDB catalog: {CATALOG}")
        return 1

    cand = pd.read_csv(CANDIDATES, dtype=str).fillna("")
    if args.accessions:
        wanted = {a.upper() for a in args.accessions}
        cand = cand[cand["source_accession"].str.upper().isin(wanted)].copy()

    geo_catalog = existing_geo_accessions()

    summary_rows = []
    sample_rows = []
    file_rows = []

    print("=" * 88)
    print("ATLAS — 00D-EBI BIOSTUDIES / ARRAYEXPRESS METADATA ENRICHMENT")
    print("=" * 88)
    print(f"Candidates selected: {len(cand)}")
    print("Phenotype policy: NO automatic resistant/sensitive assignment")

    for i, r in cand.iterrows():
        acc = str(r.get("source_accession", "")).strip()
        did = f"BIOSTUDIES:{acc}"
        geo_links = {
            x.upper() for x in str(r.get("geo_accessions", "")).split("|") if x.strip()
        }

        # ArrayExpress E-GEOD accessions are GEO imports/mirrors by definition.
        # Never treat them as independent validation datasets, even if the
        # corresponding GSE was not selected into the current GEO catalog.
        if acc.upper().startswith("E-GEOD-"):
            suffix = acc.upper().replace("E-GEOD-", "", 1)
            inferred_geo = f"GSE{suffix}" if suffix.isdigit() else ""
            if inferred_geo:
                geo_links.add(inferred_geo)

        mirrored = sorted(geo_links)

        print(f"\n[{acc}]", flush=True)

        if acc.upper().startswith("E-GEOD-") or mirrored:
            print(f"  GEO mirror: {', '.join(mirrored)} -> metadata only, no independent samples")
            summary_rows.append({
                "dataset_id": did,
                "source_accession": acc,
                "status": "GEO_MIRROR",
                "relationship_role": "CROSS_DATABASE_MIRROR_OF_GEO",
                "mirrors": "|".join(mirrored),
                "sdrf_file": "",
                "sdrf_rows": 0,
                "unique_samples": 0,
                "error": "",
            })
            continue

        try:
            study = get_json(f"{BASE}/studies/{quote(acc)}")
            info = get_json(f"{BASE}/studies/{quote(acc)}/info")
            try:
                files_json = get_json(f"{BASE}/files/{quote(acc)}")
            except Exception:
                files_json = {}

            nodes = collect_file_nodes({"study": study, "info": info, "files": files_json})
            paths = sorted({file_path(n) for n in nodes if file_path(n)})

            # Also mine arbitrary strings because ArrayExpress migration records may
            # expose file paths in nested link structures rather than file dicts.
            for s in flatten_strings({"study": study, "info": info, "files": files_json}):
                if ("sdrf" in s.lower() or "idf" in s.lower()) and "." in s:
                    paths.append(s)
            paths = sorted(set(paths))

            for p in paths:
                file_rows.append({
                    "dataset_id": did,
                    "source_accession": acc,
                    "file_path": p,
                    "file_kind": (
                        "SDRF" if "sdrf" in p.lower()
                        else "IDF" if "idf" in p.lower()
                        else "OTHER"
                    ),
                })

            sdrf = choose_sdrf(paths)
            if not sdrf:
                print("  No SDRF discovered -> manual review")
                summary_rows.append({
                    "dataset_id": did,
                    "source_accession": acc,
                    "status": "NO_SDRF_MANUAL_REVIEW",
                    "relationship_role": "INDEPENDENCE_NOT_YET_ESTABLISHED",
                    "mirrors": "",
                    "sdrf_file": "",
                    "sdrf_rows": 0,
                    "unique_samples": 0,
                    "error": "",
                })
                continue

            url = build_download_url(info, sdrf)
            if not url:
                raise RuntimeError(f"Could not build download URL for SDRF: {sdrf}")

            print(f"  SDRF: {sdrf}")
            raw = get_bytes(url)
            text = raw.decode("utf-8-sig", errors="replace")
            rows, raw_n = parse_sdrf(acc, text)

            sample_rows.extend(rows)
            print(f"  SDRF rows: {raw_n}; unique metadata samples: {len(rows)}")

            summary_rows.append({
                "dataset_id": did,
                "source_accession": acc,
                "status": "ENRICHED_SDRF",
                "relationship_role": "INDEPENDENCE_NOT_YET_ESTABLISHED",
                "mirrors": "",
                "sdrf_file": sdrf,
                "sdrf_rows": raw_n,
                "unique_samples": len(rows),
                "error": "",
            })

        except Exception as e:
            print(f"  ERROR: {e}")
            summary_rows.append({
                "dataset_id": did,
                "source_accession": acc,
                "status": "ERROR",
                "relationship_role": "INDEPENDENCE_NOT_YET_ESTABLISHED",
                "mirrors": "",
                "sdrf_file": "",
                "sdrf_rows": 0,
                "unique_samples": 0,
                "error": str(e)[:500],
            })

        time.sleep(args.sleep)

    write_csv(OUT_SUMMARY, summary_rows)
    write_csv(OUT_SAMPLES, sample_rows)
    write_csv(OUT_FILES, file_rows)

    update_duckdb(summary_rows, sample_rows)

    # --------------------------------------------------------------
    # Synchronize current EBI samples into the common ATLAS sample
    # metadata file consumed by 00D1/00D2.
    #
    # Remove prior BIOSTUDIES rows first, then append only this run's
    # non-mirror E-MTAB samples. This makes reruns idempotent and
    # prevents stale E-GEOD mirror samples from surviving.
    # --------------------------------------------------------------
    shared_samples = ENRICHED / "sample_metadata.csv"

    if shared_samples.exists():
        base = pd.read_csv(
            shared_samples,
            dtype=str,
            keep_default_na=False,
        )

        # Purge all previous BioStudies rows.
        if "dataset_id" in base.columns:
            base = base[
                ~base["dataset_id"]
                .astype(str)
                .str.startswith("BIOSTUDIES:")
            ].copy()

        ebi = pd.DataFrame(sample_rows).fillna("")

        if not ebi.empty:
            if "dataset_id" not in ebi.columns or "sample_id" not in ebi.columns:
                raise RuntimeError(
                    "EBI metadata is missing dataset_id/sample_id"
                )

            # Match the schema expected by existing ATLAS phenotype stages.
            for c in base.columns:
                if c not in ebi.columns:
                    ebi[c] = ""

            ebi = ebi[
                [c for c in base.columns if c in ebi.columns]
            ].copy()

            combined = pd.concat(
                [base, ebi],
                ignore_index=True,
            )

            combined = combined.drop_duplicates(
                subset=["dataset_id", "sample_id"],
                keep="last",
            )
        else:
            combined = base

        combined.to_csv(
            shared_samples,
            index=False,
        )

        current_ebi_n = int(
            combined["dataset_id"]
            .astype(str)
            .str.startswith("BIOSTUDIES:")
            .sum()
        )

        print(
            f"Shared sample metadata synchronized: "
            f"{len(combined)} total rows; "
            f"{current_ebi_n} BioStudies rows"
        )

    print("\n" + "=" * 88)
    print("00D-EBI COMPLETE")
    print("=" * 88)
    if summary_rows:
        print(pd.DataFrame(summary_rows)["status"].value_counts().to_string())
    print(f"EBI sample rows written: {len(sample_rows)}")
    print(OUT_SUMMARY)
    print(OUT_SAMPLES)
    print(OUT_FILES)
    print("IMPORTANT: EBI samples remain UNRESOLVED until ATLAS phenotype/design validation.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
