#!/usr/bin/env python3
"""
ATLAS — Stage 00D3
Dataset Relationship / Sample-Overlap Audit

Purpose
-------
Prevent double-counting related GEO series as independent validation datasets.

This script compares sample IDs across enriched datasets and identifies:
- exact duplicate sample sets
- subset / superseries relationships
- high sample overlap
- likely independent datasets

It does NOT modify DuckDB by default.

Outputs
-------
data/enriched/dataset_relationship_audit.csv
data/enriched/dataset_overlap_pairs.csv
"""

from __future__ import annotations

from pathlib import Path
import itertools
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd

INFILE = PROJECT_ROOT / "data" / "enriched" / "sample_metadata.csv"
OUTDIR = PROJECT_ROOT / "data" / "enriched"
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


def main():
    if not INFILE.exists():
        print(f"ERROR: missing {INFILE}")
        return 1

    df = pd.read_csv(INFILE)

    sets = {}
    for dataset_id, g in df.groupby("dataset_id"):
        sets[dataset_id] = {
            clean(x)
            for x in g["sample_id"]
            if clean(x)
        }

    pair_rows = []
    related = {d: [] for d in sets}

    for a, b in itertools.combinations(sorted(sets), 2):
        sa, sb = sets[a], sets[b]
        if not sa or not sb:
            continue

        inter = sa & sb
        if not inter:
            continue

        union = sa | sb
        jaccard = len(inter) / len(union)
        frac_a = len(inter) / len(sa)
        frac_b = len(inter) / len(sb)

        if sa == sb:
            relationship = "EXACT_SAMPLE_DUPLICATE"
        elif sa < sb:
            relationship = "A_SUBSET_OF_B"
        elif sb < sa:
            relationship = "B_SUBSET_OF_A"
        elif max(frac_a, frac_b) >= 0.8:
            relationship = "HIGH_OVERLAP"
        else:
            relationship = "PARTIAL_OVERLAP"

        pair_rows.append({
            "dataset_a": a,
            "dataset_b": b,
            "samples_a": len(sa),
            "samples_b": len(sb),
            "shared_samples": len(inter),
            "fraction_a_shared": round(frac_a, 4),
            "fraction_b_shared": round(frac_b, 4),
            "jaccard": round(jaccard, 4),
            "relationship": relationship,
        })

        related[a].append((b, relationship, len(inter), frac_a))
        related[b].append((a, relationship, len(inter), frac_b))

    pair_df = pd.DataFrame(pair_rows)

    audit_rows = []
    for dataset_id, s in sets.items():
        rels = related[dataset_id]

        likely_independent = True
        role = "INDEPENDENT_OR_NO_OVERLAP_DETECTED"
        reason = ""

        for other, rel, shared_n, frac_self in rels:
            if rel == "EXACT_SAMPLE_DUPLICATE":
                likely_independent = False
                role = "DUPLICATE_SAMPLE_SET"
                reason = f"Exact sample set duplicated with {other}"
                break

            if rel == "A_SUBSET_OF_B" and dataset_id < other:
                # lexical order is not meaningful; infer from set relationship below instead
                pass

        # Explicitly determine if this dataset contains or is contained by another.
        containers = []
        subsets = []
        for other, so in sets.items():
            if other == dataset_id or not s or not so:
                continue
            if s < so:
                containers.append(other)
            elif so < s:
                subsets.append(other)

        if containers:
            likely_independent = False
            role = "SUBSERIES_OR_NESTED_DATASET"
            reason = "All samples are contained in: " + ", ".join(sorted(containers))
        elif subsets:
            likely_independent = False
            role = "SUPERSERIES_OR_UMBRELLA_DATASET"
            reason = "Contains complete sample sets from: " + ", ".join(sorted(subsets))
        elif rels:
            max_overlap = max(r[3] for r in rels)
            if max_overlap >= 0.8:
                likely_independent = False
                role = "HIGH_OVERLAP_DATASET"
                reason = ">=80% of samples overlap another dataset"

        audit_rows.append({
            "dataset_id": dataset_id,
            "sample_n": len(s),
            "relationship_role": role,
            "likely_independent_validation_unit": likely_independent,
            "relationship_reason": reason,
        })

    audit_df = pd.DataFrame(audit_rows).sort_values(
        ["likely_independent_validation_unit", "sample_n"],
        ascending=[True, False],
    )

    p1 = OUTDIR / "dataset_relationship_audit.csv"
    p2 = OUTDIR / "dataset_overlap_pairs.csv"

    audit_df.to_csv(p1, index=False)
    pair_df.to_csv(p2, index=False)

    print("=" * 78)
    print("ATLAS — 00D3 DATASET RELATIONSHIP / OVERLAP AUDIT")
    print("=" * 78)

    print("\nPotentially non-independent datasets:")
    flagged = audit_df[
        ~audit_df["likely_independent_validation_unit"]
    ]
    if flagged.empty:
        print("None detected.")
    else:
        print(
            flagged[
                [
                    "dataset_id",
                    "sample_n",
                    "relationship_role",
                    "relationship_reason",
                ]
            ].to_string(index=False)
        )

    if not pair_df.empty:
        print("\nOverlapping dataset pairs:")
        print(
            pair_df.sort_values(
                ["shared_samples", "jaccard"],
                ascending=False,
            ).head(20).to_string(index=False)
        )

    print("\nOutputs:")
    print(p1)
    print(p2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
