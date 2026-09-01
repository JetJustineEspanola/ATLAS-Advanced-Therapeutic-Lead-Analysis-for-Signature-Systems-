#!/usr/bin/env python3
"""
ATLAS — Stage 00D2b
Curated Study-Specific Phenotype Mapping

Purpose
-------
Apply explicit, study-aware phenotype rules to high-value GEO validation cohorts
whose metadata cannot be safely handled by a global classifier.

Current curated rules
---------------------
GSE121105:
  PRIMARY comparison:
    parental BT474 untreated -> SENSITIVE_OR_PARENTAL
    BT-TR1 / BT-TR2 -> RESISTANT
  EXCLUDE from the trastuzumab-resistance contrast:
    BT474 + trastuzumab acute-treatment samples
    BT474 + trastuzumab + pertuzumab acute-treatment samples
    BT-TPR1 / BT-TPR2 dual trastuzumab+pertuzumab-resistant samples

GSE123754:
  Trastuzumab-sensitive SK -> SENSITIVE_OR_PARENTAL
  Trastuzumab-resistant SKTR -> RESISTANT
  Because n=1/group, downstream validation should be directional/exploratory only.

GSE114575:
  SKBr3 -> SENSITIVE_OR_PARENTAL
  Her/trastuz resist -> RESISTANT
  Lapatinib-only and sequential multi-drug resistant samples -> EXCLUDE
  Because n=1/group for the clean trastuzumab comparison, downstream validation
  should be directional/exploratory only.

GSE245486:
  NOT auto-mapped here.
  The supplied metadata does not explicitly establish that HR6 means
  trastuzumab-resistant. Keep it for manual review until confirmed from the
  study-level methods/publication.

Behavior
--------
- --dry-run: write proposals only
- default: apply curated phenotype labels to DuckDB samples table
- excluded samples are labeled EXCLUDE_FROM_PRIMARY_CONTRAST

Outputs
-------
data/enriched/curated_phenotype_mapping_proposals.csv
data/enriched/curated_phenotype_dataset_summary.csv
"""

from __future__ import annotations

from pathlib import Path
import argparse
import re
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import duckdb
import pandas as pd

CATALOG = PROJECT_ROOT / "data/catalog/atlas_catalog.duckdb"
INFILE = PROJECT_ROOT / "data/enriched/sample_metadata.csv"
OUTDIR = PROJECT_ROOT / "data/enriched"
OUTDIR.mkdir(parents=True, exist_ok=True)


def clean(x):
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()


def label_text(row):
    return " | ".join(
        clean(row.get(c))
        for c in ["title", "source_name", "characteristics", "treatment", "description"]
        if clean(row.get(c))
    )


def classify(dataset_id: str, row) -> tuple[str, str, str]:
    title = clean(row.get("title"))
    source = clean(row.get("source_name"))
    text = label_text(row)

    # ------------------------------------------------------------
    # GSE121105
    # ------------------------------------------------------------
    if dataset_id == "GEO:GSE121105":
        if re.match(r"^BT474_replicate\d+$", title):
            return (
                "SENSITIVE_OR_PARENTAL",
                "HIGH",
                "Curated GSE121105 rule: untreated parental BT474",
            )

        if re.match(r"^BT-TR[12]_replicate\d+$", title):
            return (
                "RESISTANT",
                "HIGH",
                "Curated GSE121105 rule: BT-TR1/BT-TR2 trastuzumab-resistant clone",
            )

        if (
            title.startswith("BT474+trastuzumab_")
            or title.startswith("BT474+trastuzumab+pertuzumab_")
            or re.match(r"^BT-TPR[12]_replicate\d+$", title)
        ):
            return (
                "EXCLUDE_FROM_PRIMARY_CONTRAST",
                "HIGH",
                "Curated GSE121105 rule: acute-treatment or dual-resistance group excluded from clean trastuzumab-resistance contrast",
            )

        return ("UNRESOLVED", "LOW", "No curated GSE121105 rule matched")

    # ------------------------------------------------------------
    # GSE123754
    # ------------------------------------------------------------
    if dataset_id == "GEO:GSE123754":
        if "Trastuzumab-sensitive model" in title:
            return (
                "SENSITIVE_OR_PARENTAL",
                "HIGH",
                "Curated GSE123754 rule: explicitly trastuzumab-sensitive SK model",
            )

        if "Trastuzumab-resistant model" in title:
            return (
                "RESISTANT",
                "HIGH",
                "Curated GSE123754 rule: explicitly trastuzumab-resistant SKTR model",
            )

        return ("UNRESOLVED", "LOW", "No curated GSE123754 rule matched")

    # ------------------------------------------------------------
    # GSE114575
    # ------------------------------------------------------------
    if dataset_id == "GEO:GSE114575":
        if source == "SKBr3":
            return (
                "SENSITIVE_OR_PARENTAL",
                "HIGH",
                "Curated GSE114575 rule: parental SKBr3",
            )

        if source == "Her/trastuz resist":
            return (
                "RESISTANT",
                "HIGH",
                "Curated GSE114575 rule: explicitly trastuzumab-resistant sample",
            )

        if source in {
            "Lap/lapatinib resist",
            "Her/Lap resist - sequential",
            "Lap/Her resist - sequential",
        }:
            return (
                "EXCLUDE_FROM_PRIMARY_CONTRAST",
                "HIGH",
                "Curated GSE114575 rule: lapatinib or sequential multidrug resistance excluded from trastuzumab-only contrast",
            )

        return ("UNRESOLVED", "LOW", "No curated GSE114575 rule matched")

    # ------------------------------------------------------------
    # GSE245486
    # ------------------------------------------------------------
    if dataset_id == "GEO:GSE245486":
        return (
            "MANUAL_REVIEW_REQUIRED",
            "LOW",
            "HR6 identity is not explicitly defined as trastuzumab-resistant in the current sample metadata; do not infer automatically",
        )

    return ("UNCHANGED", "LOW", "Dataset not in curated rule set")


def ensure_columns(con):
    existing = {
        r[1] for r in con.execute("PRAGMA table_info('samples')").fetchall()
    }
    additions = {
        "phenotype_confidence": "VARCHAR",
        "phenotype_rule": "VARCHAR",
        "primary_contrast_included": "BOOLEAN",
    }
    for col, typ in additions.items():
        if col not in existing:
            con.execute(f'ALTER TABLE samples ADD COLUMN "{col}" {typ}')


def update_dataset_flags(con, dataset_id):
    sdf = con.execute(
        """
        SELECT resistance_status, primary_contrast_included, replicate_type
        FROM samples
        WHERE dataset_id=?
        """,
        [dataset_id],
    ).fetchdf()

    if sdf.empty:
        return

    included = sdf[
        sdf["primary_contrast_included"].fillna(False)
    ].copy()

    counts = included["resistance_status"].value_counts()

    r_n = int(counts.get("RESISTANT", 0))
    s_n = int(counts.get("SENSITIVE_OR_PARENTAL", 0))
    groups_defined = r_n >= 1 and s_n >= 1

    nontech = included[
        included["replicate_type"].fillna("UNRESOLVED") != "TECHNICAL"
    ]
    nt = nontech["resistance_status"].value_counts()

    replication = (
        groups_defined
        and int(nt.get("RESISTANT", 0)) >= 2
        and int(nt.get("SENSITIVE_OR_PARENTAL", 0)) >= 2
    )

    # Important: n=1/group should never become HIGH dataset confidence.
    if groups_defined and replication:
        confidence = "HIGH"
    elif groups_defined:
        confidence = "MODERATE"
    else:
        confidence = "LOW"

    con.execute(
        """
        UPDATE datasets
        SET resistant_sensitive_groups_defined=?,
            biological_replication=?,
            phenotype_confidence=?
        WHERE dataset_id=?
        """,
        [
            bool(groups_defined),
            bool(replication),
            confidence,
            str(dataset_id),
        ],
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    p.add_argument(
        "--accessions",
        default="GSE121105,GSE123754,GSE114575,GSE245486",
    )
    args = p.parse_args()

    wanted = {
        f"GEO:{x.strip().upper()}"
        for x in args.accessions.split(",")
        if x.strip()
    }

    if not INFILE.exists():
        print(f"ERROR: missing {INFILE}")
        return 1

    df = pd.read_csv(INFILE)
    df = df[df["dataset_id"].isin(wanted)].copy()

    proposals = []

    for _, row in df.iterrows():
        status, confidence, reason = classify(row["dataset_id"], row)

        include = status in {"RESISTANT", "SENSITIVE_OR_PARENTAL"}

        proposals.append({
            "dataset_id": row["dataset_id"],
            "sample_id": row["sample_id"],
            "title": clean(row.get("title")),
            "source_name": clean(row.get("source_name")),
            "proposed_status": status,
            "phenotype_confidence": confidence,
            "primary_contrast_included": include,
            "mapping_reason": reason,
        })

    pdf = pd.DataFrame(proposals)

    summaries = []
    for dataset_id, g in pdf.groupby("dataset_id"):
        c = g["proposed_status"].value_counts()
        summaries.append({
            "dataset_id": dataset_id,
            "resistant_n": int(c.get("RESISTANT", 0)),
            "sensitive_parental_n": int(c.get("SENSITIVE_OR_PARENTAL", 0)),
            "excluded_n": int(c.get("EXCLUDE_FROM_PRIMARY_CONTRAST", 0)),
            "manual_review_n": int(c.get("MANUAL_REVIEW_REQUIRED", 0)),
            "unresolved_n": int(c.get("UNRESOLVED", 0)),
            "clean_groups_defined": (
                c.get("RESISTANT", 0) > 0
                and c.get("SENSITIVE_OR_PARENTAL", 0) > 0
            ),
            "replicated_clean_contrast": (
                c.get("RESISTANT", 0) >= 2
                and c.get("SENSITIVE_OR_PARENTAL", 0) >= 2
            ),
        })

    sdf = pd.DataFrame(summaries).sort_values(
        ["replicated_clean_contrast", "resistant_n", "sensitive_parental_n"],
        ascending=[False, False, False],
    )

    p1 = OUTDIR / "curated_phenotype_mapping_proposals.csv"
    p2 = OUTDIR / "curated_phenotype_dataset_summary.csv"

    pdf.to_csv(p1, index=False)
    sdf.to_csv(p2, index=False)

    print("=" * 78)
    print("ATLAS — 00D2b CURATED STUDY-SPECIFIC PHENOTYPE MAPPING")
    print("=" * 78)
    print("\nCurated dataset summary:")
    print(sdf.to_string(index=False))

    if args.dry_run:
        print("\nDRY RUN: DuckDB was not modified.")
    else:
        con = duckdb.connect(str(CATALOG))
        ensure_columns(con)

        for _, r in pdf.iterrows():
            status = r["proposed_status"]

            if status in {"UNCHANGED", "MANUAL_REVIEW_REQUIRED", "UNRESOLVED"}:
                continue

            if status == "EXCLUDE_FROM_PRIMARY_CONTRAST":
                db_status = "EXCLUDED"
                include = False
            else:
                db_status = status
                include = True

            con.execute(
                """
                UPDATE samples
                SET resistance_status=?,
                    biological_group=?,
                    phenotype_confidence=?,
                    phenotype_rule=?,
                    primary_contrast_included=?
                WHERE dataset_id=? AND sample_id=?
                """,
                [
                    db_status,
                    db_status,
                    r["phenotype_confidence"],
                    r["mapping_reason"],
                    bool(include),
                    r["dataset_id"],
                    r["sample_id"],
                ],
            )

        for dataset_id in sorted(wanted):
            update_dataset_flags(con, dataset_id)

        con.close()
        print("\nCurated mappings applied to DuckDB.")

    print("\nOutputs:")
    print(p1)
    print(p2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
