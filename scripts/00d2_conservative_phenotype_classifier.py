#!/usr/bin/env python3
"""
ATLAS — Stage 00D2 v3
Field-Aware Conservative Phenotype Classification

Why v3?
-------
The v2 dry run showed BT474R-containing samples becoming AMBIGUOUS because
sample-specific resistant labels were being mixed with generic parental/sensitive
language elsewhere in the sample metadata.

v3 introduces evidence hierarchy:

Tier 1: sample identity fields
  - title
  - source_name

Tier 2: sample-specific descriptive fields
  - characteristics
  - treatment
  - description

Rules
-----
1. Explicit resistant/sensitive evidence in Tier 1 takes precedence.
2. A Tier 1 resistant alias (e.g. BT474R) is not cancelled by generic
   parental/sensitive wording appearing only in Tier 2.
3. Ambiguity is retained when conflicting evidence occurs within Tier 1.
4. Tier 2 is only used when Tier 1 has no phenotype call.
5. Parental aliases (e.g. BT474) are only used if the same dataset contains
   the paired resistant alias (e.g. BT474R).
6. No low-confidence calls are automatically applied.

Outputs
-------
data/enriched/phenotype_classification_proposals_v3.csv
data/enriched/phenotype_dataset_classification_summary_v3.csv

Use --dry-run first.
"""

from __future__ import annotations

from pathlib import Path
import re
import sys
from collections import Counter

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import duckdb
import pandas as pd

CATALOG = PROJECT_ROOT / "data/catalog/atlas_catalog.duckdb"
INFILE = PROJECT_ROOT / "data/enriched/sample_metadata.csv"
OUTDIR = PROJECT_ROOT / "data/enriched"
OUTDIR.mkdir(parents=True, exist_ok=True)

IDENTITY_FIELDS = ["title", "source_name"]
SECONDARY_FIELDS = ["characteristics", "treatment", "description", "cell_line"]

STRONG_RESISTANCE_PHRASES = [
    r"\btrastuzumab[- ]resistant\b",
    r"\bherceptin[- ]resistant\b",
    r"\bacquired resistance\b",
    r"\bresistant clone\b",
    r"\bresistant cells?\b",
]

STRONG_SENSITIVE_PHRASES = [
    r"\btrastuzumab[- ]sensitive\b",
    r"\bherceptin[- ]sensitive\b",
    r"\bsensitive cells?\b",
    r"\bparental cells?\b",
]

RESISTANT_ALIASES = {
    "bt474r": "BT474",
    "bt-474r": "BT474",
    "skbr3r": "SKBR3",
    "sk-br-3r": "SKBR3",
}

PARENTAL_ALIASES = {
    "bt474": "BT474",
    "bt-474": "BT474",
    "skbr3": "SKBR3",
    "sk-br-3": "SKBR3",
}


def clean(x):
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()


def join_fields(row, fields):
    return " | ".join(clean(row.get(c)) for c in fields if clean(row.get(c)))


def lexical_tokens(text):
    return set(
        re.findall(
            r"[A-Za-z0-9]+(?:[-_/][A-Za-z0-9]+)*",
            text.lower(),
        )
    )


def matched_patterns(text, patterns):
    return [p for p in patterns if re.search(p, text, flags=re.I)]


def detect_aliases(text):
    toks = lexical_tokens(text)

    resistant = [
        (alias, lineage)
        for alias, lineage in RESISTANT_ALIASES.items()
        if alias in toks
    ]

    parental = [
        (alias, lineage)
        for alias, lineage in PARENTAL_ALIASES.items()
        if alias in toks
    ]

    return resistant, parental


def build_dataset_context(df):
    context = {}

    for dataset_id, g in df.groupby("dataset_id"):
        res_lineages = Counter()
        par_lineages = Counter()

        # Use identity fields only to establish paired aliases.
        for _, row in g.iterrows():
            text = join_fields(row, IDENTITY_FIELDS)
            r_aliases, p_aliases = detect_aliases(text)

            for _, lineage in r_aliases:
                res_lineages[lineage] += 1

            for _, lineage in p_aliases:
                par_lineages[lineage] += 1

        paired = {
            lineage
            for lineage in set(res_lineages) | set(par_lineages)
            if res_lineages[lineage] > 0 and par_lineages[lineage] > 0
        }

        context[dataset_id] = {
            "paired_lineages": paired,
            "res_lineages": dict(res_lineages),
            "par_lineages": dict(par_lineages),
        }

    return context


def evidence_for_text(text, paired_lineages):
    res_ev = []
    sen_ev = []

    for p in matched_patterns(text, STRONG_RESISTANCE_PHRASES):
        res_ev.append(f"phrase:{p}")

    for p in matched_patterns(text, STRONG_SENSITIVE_PHRASES):
        sen_ev.append(f"phrase:{p}")

    r_aliases, p_aliases = detect_aliases(text)

    for alias, lineage in r_aliases:
        res_ev.append(f"resistant_alias:{alias}")

    for alias, lineage in p_aliases:
        if lineage in paired_lineages:
            sen_ev.append(f"paired_parental_alias:{alias}")

    return res_ev, sen_ev


def classify(row, context):
    paired = context[row["dataset_id"]]["paired_lineages"]

    identity_text = join_fields(row, IDENTITY_FIELDS)
    secondary_text = join_fields(row, SECONDARY_FIELDS)

    id_res, id_sen = evidence_for_text(identity_text, paired)
    sec_res, sec_sen = evidence_for_text(secondary_text, paired)

    # Tier 1: identity fields decide first.
    if id_res and not id_sen:
        return (
            "RESISTANT",
            "HIGH",
            "IDENTITY_FIELD",
            id_res,
            id_sen,
            sec_res,
            sec_sen,
        )

    if id_sen and not id_res:
        return (
            "SENSITIVE_OR_PARENTAL",
            "HIGH",
            "IDENTITY_FIELD",
            id_res,
            id_sen,
            sec_res,
            sec_sen,
        )

    if id_res and id_sen:
        return (
            "AMBIGUOUS",
            "LOW",
            "IDENTITY_FIELD_CONFLICT",
            id_res,
            id_sen,
            sec_res,
            sec_sen,
        )

    # Tier 2 only if identity fields did not decide.
    if sec_res and not sec_sen:
        return (
            "RESISTANT",
            "MODERATE",
            "SECONDARY_FIELD",
            id_res,
            id_sen,
            sec_res,
            sec_sen,
        )

    if sec_sen and not sec_res:
        return (
            "SENSITIVE_OR_PARENTAL",
            "MODERATE",
            "SECONDARY_FIELD",
            id_res,
            id_sen,
            sec_res,
            sec_sen,
        )

    if sec_res and sec_sen:
        return (
            "AMBIGUOUS",
            "LOW",
            "SECONDARY_FIELD_CONFLICT",
            id_res,
            id_sen,
            sec_res,
            sec_sen,
        )

    return (
        "UNRESOLVED",
        "LOW",
        "NO_DECISIVE_EVIDENCE",
        id_res,
        id_sen,
        sec_res,
        sec_sen,
    )


def ensure_columns(con):
    existing = {
        r[1] for r in con.execute(
            "PRAGMA table_info('samples')"
        ).fetchall()
    }

    additions = {
        "phenotype_confidence": "VARCHAR",
        "phenotype_rule": "VARCHAR",
    }

    for col, typ in additions.items():
        if col not in existing:
            con.execute(
                f'ALTER TABLE samples ADD COLUMN "{col}" {typ}'
            )


def update_dataset_flags(con, dataset_id):
    sdf = con.execute(
        """
        SELECT resistance_status, replicate_type, phenotype_confidence
        FROM samples
        WHERE dataset_id=?
        """,
        [dataset_id],
    ).fetchdf()

    if sdf.empty:
        return

    counts = sdf["resistance_status"].value_counts()
    resistant_n = int(counts.get("RESISTANT", 0))
    sensitive_n = int(counts.get("SENSITIVE_OR_PARENTAL", 0))

    groups_defined = resistant_n >= 1 and sensitive_n >= 1

    nontech = sdf[
        sdf["replicate_type"].fillna("UNRESOLVED") != "TECHNICAL"
    ]
    nt_counts = nontech["resistance_status"].value_counts()

    basic_replication = (
        groups_defined
        and int(nt_counts.get("RESISTANT", 0)) >= 2
        and int(nt_counts.get("SENSITIVE_OR_PARENTAL", 0)) >= 2
    )

    confidence = (
        "HIGH"
        if groups_defined and basic_replication
        else "MODERATE"
        if groups_defined
        else "LOW"
    )

    con.execute(
        """
        UPDATE datasets
        SET resistant_sensitive_groups_defined=?,
            biological_replication=?,
            phenotype_confidence=?
        WHERE dataset_id=?
        """,
        [
            groups_defined,
            basic_replication,
            confidence,
            dataset_id,
        ],
    )


def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Generate proposals without modifying DuckDB.",
    )
    args = p.parse_args()

    if not INFILE.exists():
        print(f"ERROR: missing {INFILE}")
        return 1

    df = pd.read_csv(INFILE)
    context = build_dataset_context(df)

    proposals = []

    for _, row in df.iterrows():
        (
            status,
            confidence,
            rule_tier,
            id_res,
            id_sen,
            sec_res,
            sec_sen,
        ) = classify(row, context)

        # Automatic updates are intentionally restricted to HIGH-confidence
        # identity-field classifications.
        auto = (
            status in {"RESISTANT", "SENSITIVE_OR_PARENTAL"}
            and confidence == "HIGH"
            and rule_tier == "IDENTITY_FIELD"
        )

        proposals.append({
            "dataset_id": row["dataset_id"],
            "sample_id": row["sample_id"],
            "title": clean(row.get("title")),
            "source_name": clean(row.get("source_name")),
            "characteristics": clean(row.get("characteristics")),
            "treatment": clean(row.get("treatment")),
            "proposed_resistance_status": status,
            "phenotype_confidence": confidence,
            "rule_tier": rule_tier,
            "identity_resistance_evidence": " | ".join(id_res),
            "identity_sensitive_evidence": " | ".join(id_sen),
            "secondary_resistance_evidence": " | ".join(sec_res),
            "secondary_sensitive_evidence": " | ".join(sec_sen),
            "apply_automatically": auto,
        })

    pdf = pd.DataFrame(proposals)

    summary_rows = []

    for dataset_id, g in pdf.groupby("dataset_id"):
        counts = g["proposed_resistance_status"].value_counts()

        summary_rows.append({
            "dataset_id": dataset_id,
            "sample_n": len(g),
            "resistant_n": int(counts.get("RESISTANT", 0)),
            "sensitive_parental_n": int(
                counts.get("SENSITIVE_OR_PARENTAL", 0)
            ),
            "ambiguous_n": int(counts.get("AMBIGUOUS", 0)),
            "unresolved_n": int(counts.get("UNRESOLVED", 0)),
            "high_confidence_auto_n": int(g["apply_automatically"].sum()),
            "groups_defined_by_proposal": (
                counts.get("RESISTANT", 0) > 0
                and counts.get("SENSITIVE_OR_PARENTAL", 0) > 0
            ),
        })

    sdf = pd.DataFrame(summary_rows).sort_values(
        [
            "groups_defined_by_proposal",
            "resistant_n",
            "sensitive_parental_n",
            "sample_n",
        ],
        ascending=[False, False, False, False],
    )

    proposal_path = (
        OUTDIR / "phenotype_classification_proposals_v3.csv"
    )
    summary_path = (
        OUTDIR / "phenotype_dataset_classification_summary_v3.csv"
    )

    pdf.to_csv(proposal_path, index=False)
    sdf.to_csv(summary_path, index=False)

    print("=" * 78)
    print("ATLAS — 00D2 v3 FIELD-AWARE PHENOTYPE CLASSIFICATION")
    print("=" * 78)

    print("\nDataset proposals:")
    print(
        sdf[
            [
                "dataset_id",
                "sample_n",
                "resistant_n",
                "sensitive_parental_n",
                "ambiguous_n",
                "unresolved_n",
                "high_confidence_auto_n",
                "groups_defined_by_proposal",
            ]
        ].head(20).to_string(index=False)
    )

    if args.dry_run:
        print("\nDRY RUN: DuckDB was not modified.")
    else:
        con = duckdb.connect(str(CATALOG))
        ensure_columns(con)

        applied = pdf[pdf["apply_automatically"]].copy()

        for _, r in applied.iterrows():
            rule = (
                f"{r['rule_tier']} | "
                f"{r['identity_resistance_evidence']} | "
                f"{r['identity_sensitive_evidence']}"
            ).strip(" |")

            con.execute(
                """
                UPDATE samples
                SET resistance_status=?,
                    biological_group=?,
                    phenotype_confidence=?,
                    phenotype_rule=?
                WHERE dataset_id=? AND sample_id=?
                """,
                [
                    r["proposed_resistance_status"],
                    r["proposed_resistance_status"],
                    r["phenotype_confidence"],
                    rule,
                    r["dataset_id"],
                    r["sample_id"],
                ],
            )

        for dataset_id in applied["dataset_id"].drop_duplicates():
            update_dataset_flags(con, dataset_id)

        con.close()

        print(
            f"\nApplied HIGH-confidence identity-field phenotype calls: "
            f"{len(applied)}"
        )

    print("\nOutputs:")
    print(proposal_path)
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
