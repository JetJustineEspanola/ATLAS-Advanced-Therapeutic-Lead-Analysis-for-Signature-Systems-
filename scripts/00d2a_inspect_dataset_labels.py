#!/usr/bin/env python3
"""
ATLAS — Stage 00D2a
Targeted Phenotype Label Inspector

Prints the exact sample-level metadata labels for selected datasets so that
dataset-specific phenotype rules can be written without guessing.

Usage:
    python -u scripts/00d2a_inspect_dataset_labels.py \
      --accessions GSE121105,GSE123754,GSE114575,GSE245486
"""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd

INFILE = PROJECT_ROOT / "data/enriched/sample_metadata.csv"


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
    p = argparse.ArgumentParser()
    p.add_argument("--accessions", required=True)
    args = p.parse_args()

    if not INFILE.exists():
        print(f"ERROR: missing {INFILE}")
        return 1

    wanted = {
        f"GEO:{x.strip().upper()}"
        for x in args.accessions.split(",")
        if x.strip()
    }

    df = pd.read_csv(INFILE)
    df = df[df["dataset_id"].isin(wanted)].copy()

    if df.empty:
        print("No matching datasets found.")
        return 1

    cols = [
        "dataset_id",
        "sample_id",
        "title",
        "source_name",
        "characteristics",
        "treatment",
        "description",
        "resistance_status",
        "cell_line",
        "replicate_type",
    ]

    for dataset_id, g in df.groupby("dataset_id"):
        print("\n" + "=" * 100)
        print(dataset_id)
        print("=" * 100)

        for _, r in g.iterrows():
            print(f"\n[{clean(r.get('sample_id'))}]")
            for c in cols[2:]:
                val = clean(r.get(c))
                if val:
                    print(f"{c}: {val}")

        print("\n--- Distinct title/source_name combinations ---")
        label_df = (
            g.assign(
                _label=g.apply(
                    lambda r: " | ".join(
                        x for x in [
                            f"title={clean(r.get('title'))}" if clean(r.get('title')) else "",
                            f"source_name={clean(r.get('source_name'))}" if clean(r.get('source_name')) else "",
                        ]
                        if x
                    ),
                    axis=1,
                )
            )
            .groupby("_label", dropna=False)
            .size()
            .reset_index(name="n")
            .sort_values("n", ascending=False)
        )
        print(label_df.to_string(index=False))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
