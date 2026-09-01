#!/usr/bin/env python3
"""
ATLAS 00C2 — Re-score datasets after 00D sample-level enrichment.
"""

from __future__ import annotations
from pathlib import Path
import sys
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import duckdb
import pandas as pd

CATALOG = PROJECT_ROOT / "data/catalog/atlas_catalog.duckdb"
CONFIG = PROJECT_ROOT / "config/dataset_queries.json"
OUT = PROJECT_ROOT / "data/enriched/dataset_candidates_rescored.csv"


def b(v):
    if v is None:
        return False
    try:
        if pd.isna(v):
            return False
    except Exception:
        pass
    return bool(v)


def main():
    cfg = json.loads(CONFIG.read_text())
    w = cfg["eligibility_weights"]
    th = cfg["eligibility_thresholds"]

    con = duckdb.connect(str(CATALOG))
    df = con.execute("SELECT * FROM datasets").fetchdf()

    rows = []
    for _, r in df.iterrows():
        score = 0
        reasons = []

        checks = [
            ("direct_trastuzumab_resistance", b(r["direct_trastuzumab_resistance"])),
            ("her2_positive_confirmed", b(r["her2_positive_confirmed"])),
            ("resistant_sensitive_groups_defined", b(r["resistant_sensitive_groups_defined"])),
            ("biological_replication", b(r["biological_replication"])),
            ("raw_data_available", b(r["raw_data_available"])),
            ("complete_sample_metadata", b(r["complete_sample_metadata"])),
            ("independent_model_or_cohort", b(r["independent_model_or_cohort"])),
            ("pd1_pdl1_or_tgfb_relevance", b(r["pd1_pdl1_or_tgfb_relevance"])),
        ]

        for key, flag in checks:
            if flag:
                score += int(w[key])
                reasons.append(f"+{w[key]} {key}")

        score = max(0, min(100, score))

        if score >= th["PRIMARY_VALIDATION"]:
            cat = "PRIMARY_VALIDATION"
        elif score >= th["SUPPORTING_VALIDATION"]:
            cat = "SUPPORTING_VALIDATION"
        elif score >= th["EXPLORATORY"]:
            cat = "EXPLORATORY"
        else:
            cat = "EXCLUDE_OR_MANUAL_REVIEW"

        critical = (
            b(r["direct_trastuzumab_resistance"])
            and b(r["her2_positive_confirmed"])
            and b(r["resistant_sensitive_groups_defined"])
        )

        manual = not (cat == "PRIMARY_VALIDATION" and critical)

        con.execute(
            """
            UPDATE datasets
            SET eligibility_score=?, eligibility_category=?,
                eligibility_reasons=?, manual_review_required=?
            WHERE dataset_id=?
            """,
            [score, cat, " | ".join(reasons), manual, r["dataset_id"]],
        )

        rows.append({
            **r.to_dict(),
            "eligibility_score": score,
            "eligibility_category": cat,
            "eligibility_reasons": " | ".join(reasons),
            "manual_review_required": manual,
        })

    out = pd.DataFrame(rows).sort_values(
        ["eligibility_score", "source", "source_accession"],
        ascending=[False, True, True],
    )
    out.to_csv(OUT, index=False)
    con.close()

    print("=== 00C2 POST-ENRICHMENT RESCORE COMPLETE ===")
    print(out["eligibility_category"].value_counts().to_string())
    print("\nTop enriched candidates:")
    cols = [
        "eligibility_score","eligibility_category","source",
        "source_accession","title","phenotype_confidence",
        "resistant_sensitive_groups_defined","biological_replication"
    ]
    print(out[cols].head(20).to_string(index=False))
    print(f"\n{OUT}")


if __name__ == "__main__":
    main()
