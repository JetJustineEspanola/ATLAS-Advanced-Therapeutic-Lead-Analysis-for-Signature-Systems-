#!/usr/bin/env python3
"""
ATLAS — 00F Primary Validation Matrix Inspection

Reads the downloaded primary-validation count matrices, reports shape/columns,
and attempts to map matrix sample columns to GEO sample metadata. This is the
safety checkpoint before running differential expression.

Outputs:
  data/validation_expression/matrix_inspection_summary.csv
  data/validation_expression/GSE121105_column_mapping.csv
  data/validation_expression/GSE237606_column_mapping.csv
"""

from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/validation_expression"
META = ROOT / "data/enriched/sample_metadata.csv"

FILES = {
    "GSE121105": DATA / "GSE121105/GSE121105_geneCount.csv.gz",
    "GSE237606": DATA / "GSE237606/GSE237606_RawCounts.txt.gz",
}


def clean(x):
    return "" if pd.isna(x) else str(x).strip()


def normalize(s):
    return re.sub(r"[^a-z0-9]+", "", clean(s).lower())


def read_matrix(path: Path):
    if path.suffixes[-2:] == [".csv", ".gz"] or path.name.endswith(".csv.gz"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def best_match(col, meta):
    ncol = normalize(col)
    best = None
    score = -1

    for _, r in meta.iterrows():
        candidates = {
            "sample_id": clean(r.get("sample_id")),
            "title": clean(r.get("title")),
            "source_name": clean(r.get("source_name")),
            "description": clean(r.get("description")),
        }

        for field, val in candidates.items():
            nv = normalize(val)
            if not nv:
                continue

            s = 0
            if ncol == nv:
                s = 100
            elif ncol in nv or nv in ncol:
                s = 80
            else:
                # reward shared alphanumeric token fragments
                toks_c = set(re.findall(r"[a-z]+\d*|\d+", clean(col).lower()))
                toks_v = set(re.findall(r"[a-z]+\d*|\d+", val.lower()))
                overlap = len(toks_c & toks_v)
                s = overlap * 10

            if s > score:
                score = s
                best = (r, field, val)

    if best is None:
        return {
            "matrix_column": col,
            "matched_sample_id": "",
            "matched_title": "",
            "matched_resistance_status": "",
            "primary_contrast_included": "",
            "match_field": "",
            "match_score": 0,
        }

    r, field, val = best
    return {
        "matrix_column": col,
        "matched_sample_id": clean(r.get("sample_id")),
        "matched_title": clean(r.get("title")),
        "matched_resistance_status": clean(r.get("resistance_status")),
        "primary_contrast_included": clean(r.get("primary_contrast_included")),
        "match_field": field,
        "match_score": score,
    }


def main():
    meta = pd.read_csv(META)
    summary = []

    print("=" * 78)
    print("ATLAS — 00F PRIMARY VALIDATION MATRIX INSPECTION")
    print("=" * 78)

    for acc, path in FILES.items():
        print(f"\n[{acc}]")

        if not path.exists():
            print(f"  ERROR: missing {path}")
            continue

        df = read_matrix(path)

        print(f"  shape: {df.shape[0]} rows x {df.shape[1]} columns")
        print(f"  first columns: {list(df.columns[:10])}")

        # Assume first column is gene identifier unless clearly numeric sample data.
        gene_col = df.columns[0]
        sample_cols = list(df.columns[1:])

        m = meta[meta["dataset_id"] == f"GEO:{acc}"].copy()
        mappings = [best_match(c, m) for c in sample_cols]
        map_df = pd.DataFrame(mappings)

        out = DATA / f"{acc}_column_mapping.csv"
        map_df.to_csv(out, index=False)

        high = int((map_df["match_score"] >= 80).sum()) if not map_df.empty else 0
        print(f"  candidate gene-id column: {gene_col}")
        print(f"  sample columns: {len(sample_cols)}")
        print(f"  high-confidence metadata matches: {high}/{len(sample_cols)}")

        if not map_df.empty:
            print(map_df.head(20).to_string(index=False))

        summary.append({
            "accession": acc,
            "rows": df.shape[0],
            "columns": df.shape[1],
            "gene_id_column": gene_col,
            "sample_columns": len(sample_cols),
            "high_confidence_matches": high,
            "mapping_file": str(out),
        })

    sout = DATA / "matrix_inspection_summary.csv"
    pd.DataFrame(summary).to_csv(sout, index=False)

    print("\nOutput:")
    print(sout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
