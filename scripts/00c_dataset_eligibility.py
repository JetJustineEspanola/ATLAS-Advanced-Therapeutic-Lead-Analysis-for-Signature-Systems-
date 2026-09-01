#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
from pathlib import Path
import duckdb
import pandas as pd

from atlas_data.common import CATALOG_PATH, DISCOVERY_DIR, load_config, clean

def boolish(v):
    if v is None:
        return False
    try:
        if pd.isna(v):
            return False
    except Exception:
        pass
    if isinstance(v, str):
        return v.strip().lower() in {"true","1","yes","y"}
    return bool(v)

def score_row(row, weights):
    reasons = []
    score = 0

    def add(key, condition, label):
        nonlocal score
        if condition:
            score += int(weights[key])
            reasons.append(f"+{weights[key]} {label}")

    add("direct_trastuzumab_resistance",
        boolish(row.get("direct_trastuzumab_resistance")),
        "direct trastuzumab-resistance wording")

    add("her2_positive_confirmed",
        boolish(row.get("her2_positive_confirmed")),
        "HER2/ERBB2 context")

    add("resistant_sensitive_groups_defined",
        boolish(row.get("resistant_sensitive_groups_defined")),
        "resistant/sensitive groups defined")

    add("biological_replication",
        boolish(row.get("biological_replication")),
        "biological replication documented")

    add("raw_data_available",
        boolish(row.get("raw_data_available")),
        "raw data available")

    add("complete_sample_metadata",
        boolish(row.get("complete_sample_metadata")),
        "complete sample metadata")

    add("independent_model_or_cohort",
        boolish(row.get("independent_model_or_cohort", True)),
        "independent model/cohort")

    add("pd1_pdl1_or_tgfb_relevance",
        boolish(row.get("pd1_pdl1_or_tgfb_relevance")),
        "PD-1/PD-L1 or TGF-beta relevance")

    # Evidence-quality penalties.
    sample_count = pd.to_numeric(
        pd.Series([row.get("sample_count")]), errors="coerce"
    ).iloc[0]

    if pd.notna(sample_count) and sample_count < 2:
        score -= 10
        reasons.append("-10 singleton/no replication")

    if clean(row.get("organism")).lower() not in {
        "", "homo sapiens", "human", "9606"
    }:
        score -= 10
        reasons.append("-10 non-human organism")

    return max(0, min(100, int(score))), " | ".join(reasons)

def category(score, thresholds):
    if score >= thresholds["PRIMARY_VALIDATION"]:
        return "PRIMARY_VALIDATION"
    if score >= thresholds["SUPPORTING_VALIDATION"]:
        return "SUPPORTING_VALIDATION"
    if score >= thresholds["EXPLORATORY"]:
        return "EXPLORATORY"
    return "EXCLUDE_OR_MANUAL_REVIEW"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--catalog", default=str(CATALOG_PATH))
    args = p.parse_args()

    config = load_config()
    weights = config["eligibility_weights"]
    thresholds = config["eligibility_thresholds"]

    con = duckdb.connect(args.catalog)
    df = con.execute("SELECT * FROM datasets").fetchdf()

    if df.empty:
        print("No datasets in catalog.")
        return 1

    scored = []
    for _, row in df.iterrows():
        s, reasons = score_row(row, weights)
        cat = category(s, thresholds)

        # Automatically require manual review unless primary and all critical
        # phenotype metadata are resolved.
        critical = (
            boolish(row.get("direct_trastuzumab_resistance"))
            and boolish(row.get("her2_positive_confirmed"))
            and boolish(row.get("resistant_sensitive_groups_defined"))
        )
        manual = not (cat == "PRIMARY_VALIDATION" and critical)

        scored.append((row["dataset_id"], s, cat, reasons, manual))

    score_df = pd.DataFrame(
        scored,
        columns=[
            "dataset_id","eligibility_score","eligibility_category",
            "eligibility_reasons","manual_regulatory_review_required_dummy"
        ],
    )

    con.register("score_df", score_df)
    con.execute("""
        UPDATE datasets
        SET
            eligibility_score = s.eligibility_score,
            eligibility_category = s.eligibility_category,
            eligibility_reasons = s.eligibility_reasons,
            manual_review_required = s.manual_regulatory_review_required_dummy
        FROM score_df s
        WHERE datasets.dataset_id = s.dataset_id
    """)

    out = con.execute("""
        SELECT *
        FROM datasets
        ORDER BY eligibility_score DESC, source, source_accession
    """).fetchdf()

    con.close()

    csv_path = DISCOVERY_DIR / "dataset_candidates_scored.csv"
    parquet_path = DISCOVERY_DIR / "dataset_candidates_scored.parquet"
    out.to_csv(csv_path, index=False)
    try:
        out.to_parquet(parquet_path, index=False)
    except Exception as e:
        print(f"Parquet skipped: {e}")

    print("=== 00C COMPLETE ===")
    print(
        out["eligibility_category"]
        .value_counts()
        .to_string()
    )
    print("\nTop candidates:")
    cols = [
        "eligibility_score","eligibility_category","source",
        "source_accession","title","manual_review_required"
    ]
    print(out[cols].head(20).to_string(index=False))
    print(f"\n{csv_path}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
