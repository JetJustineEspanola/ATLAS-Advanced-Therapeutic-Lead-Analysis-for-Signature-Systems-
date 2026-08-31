#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd

from atlas_data.common import (
    PROJECT_ROOT, DISCOVERY_DIR, MANIFEST_DIR,
    canonical_row, clean, http_get, load_config
)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GDC_FILES = "https://api.gdc.cancer.gov/files"
CBIO = "https://www.cbioportal.org/api"

def ncbi_search(db: str, query: str, retmax: int):
    r = http_get(
        f"{EUTILS}/esearch.fcgi",
        params={
            "db": db, "term": query, "retmode": "json",
            "retmax": retmax, "tool": "ATLAS"
        },
    )
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    s = http_get(
        f"{EUTILS}/esummary.fcgi",
        params={
            "db": db, "id": ",".join(ids), "retmode": "json",
            "tool": "ATLAS"
        },
    ).json()
    result = s.get("result", {})
    return [result[i] for i in ids if i in result]

def discover_geo(query_family: str, query: str, retmax: int):
    rows = []
    for d in ncbi_search("gds", query, retmax):
        accession = clean(d.get("accession")) or clean(d.get("gse"))
        title = clean(d.get("title"))
        summary = clean(d.get("summary"))
        gds_type = clean(d.get("gdstype"))
        n_samples = d.get("n_samples")
        rows.append(canonical_row(
            dataset_id=f"GEO:{accession}",
            source="GEO",
            source_accession=accession,
            title=title,
            summary=summary,
            organism=clean(d.get("taxon")),
            assay_type=gds_type,
            sample_count=n_samples,
            publication_date=clean(d.get("pdat")),
            source_url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
            processed_data_available=True,
            raw_data_available=None,
            query_family=query_family,
            query_text=query,
        ))
    return rows

def discover_sra(query_family: str, query: str, retmax: int):
    rows = []
    for d in ncbi_search("sra", query, retmax):
        accession = clean(d.get("runs")) or clean(d.get("accession"))
        title = clean(d.get("title"))
        expxml = clean(d.get("expxml"))
        summary = clean(d.get("summary")) or expxml
        # Extract first SRX/SRP/SRR-like token if summary lacks a clean accession.
        m = re.search(r"\b(SR[APRX]\d+)\b", accession + " " + expxml)
        acc = m.group(1) if m else accession
        rows.append(canonical_row(
            dataset_id=f"SRA:{acc}",
            source="SRA",
            source_accession=acc,
            title=title,
            summary=summary,
            organism=clean(d.get("organism")),
            assay_type="RAW_SEQUENCING",
            raw_data_available=True,
            processed_data_available=False,
            source_url=f"https://www.ncbi.nlm.nih.gov/sra/?term={acc}",
            query_family=query_family,
            query_text=query,
        ))
    return rows

def discover_gdc(query_family: str, query: str, retmax: int):
    # GDC does not offer free-text study discovery like GEO, so query TCGA-BRCA
    # harmonized RNA-seq files as a translational validation cohort.
    if query_family not in {"trastuzumab_resistance", "tgfb_resistance", "checkpoint_resistance"}:
        return []
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-BRCA"]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        ],
    }
    fields = ",".join([
        "file_id","file_name","data_type","data_format","experimental_strategy",
        "cases.submitter_id","cases.samples.sample_type","cases.project.project_id"
    ])
    r = http_get(
        GDC_FILES,
        params={
            "filters": json.dumps(filters),
            "fields": fields,
            "format": "JSON",
            "size": min(retmax, 100),
        },
    )
    hits = r.json().get("data", {}).get("hits", [])
    rows = []
    for h in hits:
        fid = clean(h.get("file_id"))
        fname = clean(h.get("file_name"))
        cases = h.get("cases", []) or []
        case_ids = ";".join(clean(c.get("submitter_id")) for c in cases if c)
        rows.append(canonical_row(
            dataset_id=f"GDC:{fid}",
            source="GDC",
            source_accession=fid,
            title=f"TCGA-BRCA {fname}",
            summary=f"TCGA-BRCA harmonized transcriptome file; case={case_ids}",
            organism="Homo sapiens",
            assay_type=clean(h.get("data_type")),
            platform="BULK_RNA_SEQ",
            sample_count=1,
            raw_data_available=False,
            processed_data_available=True,
            source_url=f"https://portal.gdc.cancer.gov/files/{fid}",
            query_family=query_family,
            query_text=query,
        ))
    return rows

def discover_cbioportal(query_family: str, query: str, retmax: int):
    # Retrieve studies and locally filter by breast/HER2-related text.
    try:
        r = http_get(
            f"{CBIO}/studies",
            params={"pageSize": min(max(retmax * 4, 100), 1000), "pageNumber": 0},
        )
        studies = r.json()
    except Exception:
        return []
    qwords = {w.lower() for w in re.findall(r"[A-Za-z0-9-]{4,}", query)}
    rows = []
    for st in studies:
        text = " ".join([
            clean(st.get("name")),
            clean(st.get("description")),
            clean(st.get("cancerTypeId")),
        ])
        low = text.lower()
        if "breast" not in low and "brca" not in low:
            continue
        # lightweight relevance filter
        if qwords and not any(w in low for w in qwords if w not in {"and", "with"}):
            continue
        sid = clean(st.get("studyId"))
        rows.append(canonical_row(
            dataset_id=f"CBIOPORTAL:{sid}",
            source="CBIOPORTAL",
            source_accession=sid,
            title=clean(st.get("name")),
            summary=clean(st.get("description")),
            organism="Homo sapiens",
            assay_type="CANCER_GENOMICS_STUDY",
            raw_data_available=False,
            processed_data_available=True,
            source_url=f"https://www.cbioportal.org/study/summary?id={sid}",
            query_family=query_family,
            query_text=query,
        ))
        if len(rows) >= retmax:
            break
    return rows

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sources", default="geo,sra,gdc,cbioportal")
    p.add_argument("--retmax", type=int, default=25)
    p.add_argument("--query-family", default="all")
    args = p.parse_args()

    config = load_config()
    families = config["query_families"]
    if args.query_family != "all":
        families = {args.query_family: families[args.query_family]}

    adapters = {
        "geo": discover_geo,
        "sra": discover_sra,
        "gdc": discover_gdc,
        "cbioportal": discover_cbioportal,
    }
    selected = [x.strip().lower() for x in args.sources.split(",") if x.strip()]

    rows = []
    for family, queries in families.items():
        for query in queries:
            print(f"\n[{family}] {query}", flush=True)
            for source in selected:
                fn = adapters.get(source)
                if fn is None:
                    print(f"  {source}: unsupported", flush=True)
                    continue
                try:
                    found = fn(family, query, args.retmax)
                    rows.extend(found)
                    print(f"  {source}: {len(found)}", flush=True)
                except Exception as e:
                    print(f"  {source}: ERROR {e}", flush=True)

    df = pd.DataFrame(rows)
    if df.empty:
        print("No datasets discovered.")
        return 1

    # Deduplicate by source + accession, retaining all query families/queries.
    agg = {}
    for _, row in df.iterrows():
        key = (row["source"], row["source_accession"])
        d = agg.setdefault(key, row.to_dict())
        if row["query_family"] not in clean(d.get("query_family")).split(" | "):
            d["query_family"] = " | ".join(
                sorted(set(filter(None, clean(d.get("query_family")).split(" | ") + [row["query_family"]])))
            )
        if row["query_text"] not in clean(d.get("query_text")):
            d["query_text"] = clean(d.get("query_text")) + " || " + row["query_text"]

    out = pd.DataFrame(list(agg.values()))
    out = out.sort_values(["source", "source_accession"]).reset_index(drop=True)

    csv_path = DISCOVERY_DIR / "dataset_candidates.csv"
    parquet_path = DISCOVERY_DIR / "dataset_candidates.parquet"
    out.to_csv(csv_path, index=False)
    try:
        out.to_parquet(parquet_path, index=False)
    except Exception as e:
        print(f"Parquet skipped: {e}")

    manifest = {
        "created_utc": pd.Timestamp.utcnow().isoformat(),
        "sources": selected,
        "query_family": args.query_family,
        "retmax": args.retmax,
        "unique_candidates": int(len(out)),
        "output_csv": str(csv_path),
        "output_parquet": str(parquet_path),
    }
    (MANIFEST_DIR / "00a_discovery_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    print("\n=== 00A COMPLETE ===")
    print(f"Unique candidates: {len(out)}")
    print(out["source"].value_counts().to_string())
    print(f"\n{csv_path}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
