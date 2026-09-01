#!/usr/bin/env python3
"""
ATLAS — Stage 00C3
Phenotype-Gated Dataset Eligibility Re-scoring

Fixes the main problem in 00C2:
a dataset cannot become SUPPORTING_VALIDATION merely by accumulating
metadata/raw-data points while its resistant-vs-sensitive phenotype is unresolved.

Scientific gates
----------------
PRIMARY_VALIDATION requires:
- direct trastuzumab-resistance evidence
- HER2/ERBB2 context
- resistant and sensitive/parental groups defined
- biological replication
- reasonably complete sample metadata

SUPPORTING_VALIDATION requires:
- HER2/ERBB2 context
- resistant and sensitive/parental groups defined
- at least one of:
    biological replication
    direct trastuzumab-resistance evidence

EXPLORATORY:
- potentially relevant but phenotype/group structure is incomplete

EXCLUDE_OR_MANUAL_REVIEW:
- insufficient relevance or incompatible evidence
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


def txt(v):
    if v is None:
        return ""
    try:
        if pd.isna(v):
            return ""
    except Exception:
        pass
    return str(v).strip()


def main():
    cfg = json.loads(CONFIG.read_text())
    w = cfg["eligibility_weights"]

    con = duckdb.connect(str(CATALOG))
    df = con.execute("SELECT * FROM datasets").fetchdf()

    rows = []

    for _, r in df.iterrows():
        score = 0
        reasons = []

        flags = {
            "direct_trastuzumab_resistance": b(r.get("direct_trastuzumab_resistance")),
            "her2_positive_confirmed": b(r.get("her2_positive_confirmed")),
            "resistant_sensitive_groups_defined": b(r.get("resistant_sensitive_groups_defined")),
            "biological_replication": b(r.get("biological_replication")),
            "raw_data_available": b(r.get("raw_data_available")),
            "complete_sample_metadata": b(r.get("complete_sample_metadata")),
            "independent_model_or_cohort": b(r.get("independent_model_or_cohort")),
            "pd1_pdl1_or_tgfb_relevance": b(r.get("pd1_pdl1_or_tgfb_relevance")),
        }

        for key, flag in flags.items():
            if flag:
                score += int(w[key])
                reasons.append(f"+{w[key]} {key}")

        # Basic penalties
        sample_count = pd.to_numeric(
            pd.Series([r.get("sample_count")]),
            errors="coerce",
        ).iloc[0]

        if pd.notna(sample_count) and sample_count < 2:
            score -= 10
            reasons.append("-10 singleton_dataset_unit")

        phenotype_conf = txt(r.get("phenotype_confidence")).upper()
        if phenotype_conf in {"LOW", "UNRESOLVED", ""}:
            score -= 5
            reasons.append("-5 low_or_unresolved_phenotype_confidence")

        score = max(0, min(100, int(score)))

        # -------------------------------------------------------------
        # Phenotype gates
        # -------------------------------------------------------------
        primary_gate = (
            flags["direct_trastuzumab_resistance"]
            and flags["her2_positive_confirmed"]
            and flags["resistant_sensitive_groups_defined"]
            and flags["biological_replication"]
            and flags["complete_sample_metadata"]
        )

        supporting_gate = (
            flags["her2_positive_confirmed"]
            and flags["resistant_sensitive_groups_defined"]
            and (
                flags["biological_replication"]
                or flags["direct_trastuzumab_resistance"]
            )
        )

        exploratory_gate = (
            flags["her2_positive_confirmed"]
            or flags["direct_trastuzumab_resistance"]
            or flags["pd1_pdl1_or_tgfb_relevance"]
        )

        if primary_gate and score >= 80:
            cat = "PRIMARY_VALIDATION"
        elif supporting_gate and score >= 60:
            cat = "SUPPORTING_VALIDATION"
        elif exploratory_gate and score >= 35:
            cat = "EXPLORATORY"
        else:
            cat = "EXCLUDE_OR_MANUAL_REVIEW"

        # Run-level SRA records must never be promoted by themselves.
        source = txt(r.get("source")).upper()
        accession = txt(r.get("source_accession")).upper()

        if source == "SRA" and accession.startswith("SRR"):
            if cat in {"PRIMARY_VALIDATION", "SUPPORTING_VALIDATION"}:
                cat = "EXPLORATORY"
                reasons.append(
                    "CATEGORY_CAP: SRR is a run, not an independent study-level validation dataset"
                )

        manual = cat != "PRIMARY_VALIDATION"

        con.execute(
            """
            UPDATE datasets
            SET eligibility_score=?,
                eligibility_category=?,
                eligibility_reasons=?,
                manual_review_required=?
            WHERE dataset_id=?
            """,
            [
                score,
                cat,
                " | ".join(reasons),
                manual,
                r["dataset_id"],
            ],
        )

        row = r.to_dict()
        row.update({
            "eligibility_score": score,
            "eligibility_category": cat,
            "eligibility_reasons": " | ".join(reasons),
            "manual_review_required": manual,
        })
        rows.append(row)

    out = pd.DataFrame(rows).sort_values(
        [
            "eligibility_category",
            "eligibility_score",
            "source",
            "source_accession",
        ],
        ascending=[True, False, True, True],
    )

    out.to_csv(OUT, index=False)
    try:
        out.to_parquet(
            OUT.with_suffix(".parquet"),
            index=False,
        )
    except Exception:
        pass

    con.close()

    print("=== 00C3 PHENOTYPE-GATED RESCORE COMPLETE ===")
    print(out["eligibility_category"].value_counts().to_string())

    print("\nTop candidates:")
    cols = [
        "eligibility_score",
        "eligibility_category",
        "source",
        "source_accession",
        "title",
        "phenotype_confidence",
        "resistant_sensitive_groups_defined",
        "biological_replication",
    ]
    print(out[cols].head(25).to_string(index=False))
    print(f"\n{OUT}")


if __name__ == "__main__":
    main()
