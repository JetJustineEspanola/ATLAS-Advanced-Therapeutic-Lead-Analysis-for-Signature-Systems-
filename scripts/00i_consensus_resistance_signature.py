#!/usr/bin/env python3
"""
ATLAS — 00I Consensus Resistance Signature

Builds a conservative consensus signature from replicated primary DEGs.

Selection:
  - replicated in both primary cohorts
  - same direction
  - FDR < 0.05 in both
  - minimum absolute log2FC >= 1.0 in BOTH studies

Ranking:
  geometric mean absolute effect size × combined significance

Outputs:
  results/external_validation/consensus_resistance_signature.csv
  results/external_validation/consensus_top100_up.csv
  results/external_validation/consensus_top100_down.csv
  results/external_validation/consensus_signature_summary.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INFILE = ROOT / "results/external_validation/replicated_primary_DEGs.csv"
OUTDIR = ROOT / "results/external_validation"
OUTDIR.mkdir(parents=True, exist_ok=True)

MIN_ABS_LFC = 1.0


def main():
    df = pd.read_csv(INFILE)

    df["abs_lfc_121105"] = df["GSE121105_log2FC"].abs()
    df["abs_lfc_237606"] = df["GSE237606_log2FC"].abs()

    strong = df[
        (df["replicated_primary"] == True)
        & (df["abs_lfc_121105"] >= MIN_ABS_LFC)
        & (df["abs_lfc_237606"] >= MIN_ABS_LFC)
    ].copy()

    eps = 1e-300
    p1 = strong["GSE121105_padj"].clip(lower=eps).fillna(1.0)
    p2 = strong["GSE237606_padj"].clip(lower=eps).fillna(1.0)

    strong["min_abs_log2FC"] = strong[
        ["abs_lfc_121105", "abs_lfc_237606"]
    ].min(axis=1)

    strong["geom_mean_abs_log2FC"] = np.sqrt(
        strong["abs_lfc_121105"] * strong["abs_lfc_237606"]
    )

    strong["combined_significance"] = -np.log10(p1) - np.log10(p2)

    strong["consensus_rank_score"] = (
        strong["geom_mean_abs_log2FC"] * strong["combined_significance"]
    )

    strong = strong.sort_values(
        ["consensus_rank_score", "min_abs_log2FC"],
        ascending=[False, False],
    ).reset_index(drop=True)

    strong["consensus_rank"] = np.arange(1, len(strong) + 1)

    cols = [
        "consensus_rank",
        "gene_symbol",
        "replication_direction",
        "GSE121105_log2FC",
        "GSE121105_padj",
        "GSE237606_log2FC",
        "GSE237606_padj",
        "min_abs_log2FC",
        "geom_mean_abs_log2FC",
        "combined_significance",
        "consensus_rank_score",
    ]
    strong = strong[cols]

    out_all = OUTDIR / "consensus_resistance_signature.csv"
    out_up = OUTDIR / "consensus_top100_up.csv"
    out_down = OUTDIR / "consensus_top100_down.csv"
    out_summary = OUTDIR / "consensus_signature_summary.csv"

    strong.to_csv(out_all, index=False)

    up = strong[
        strong["replication_direction"] == "UP_IN_RESISTANT"
    ].head(100)
    down = strong[
        strong["replication_direction"] == "DOWN_IN_RESISTANT"
    ].head(100)

    up.to_csv(out_up, index=False)
    down.to_csv(out_down, index=False)

    summary = pd.DataFrame([{
        "replicated_input_genes": len(df),
        "consensus_abs_log2fc_threshold_each_study": MIN_ABS_LFC,
        "consensus_genes": len(strong),
        "consensus_up": int((strong["replication_direction"] == "UP_IN_RESISTANT").sum()),
        "consensus_down": int((strong["replication_direction"] == "DOWN_IN_RESISTANT").sum()),
        "top100_up_n": len(up),
        "top100_down_n": len(down),
    }])
    summary.to_csv(out_summary, index=False)

    print("=" * 78)
    print("ATLAS — 00I CONSENSUS RESISTANCE SIGNATURE")
    print("=" * 78)
    print(f"Replicated input genes: {len(df)}")
    print(f"Consensus genes with |log2FC| >= {MIN_ABS_LFC} in both studies: {len(strong)}")
    print(f"  Up in resistant: {(strong['replication_direction'] == 'UP_IN_RESISTANT').sum()}")
    print(f"  Down in resistant: {(strong['replication_direction'] == 'DOWN_IN_RESISTANT').sum()}")

    print("\nTop consensus genes:")
    print(strong.head(25).to_string(index=False))

    print("\nOutputs:")
    print(out_all)
    print(out_up)
    print(out_down)
    print(out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
