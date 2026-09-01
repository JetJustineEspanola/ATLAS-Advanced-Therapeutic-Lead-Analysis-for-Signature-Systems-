#!/usr/bin/env python3
"""
ATLAS — 00O Reproducible TGF-beta Module Validation

Builds a conservative "positive TGF-beta resistance module" from the 16
leading-edge genes shared by discovery and GSE237606, then checks:
- direction and significance in all three datasets
- overlap with the 242-gene strict three-dataset core
- whether GSE121105 reverses the module as a coordinated block

Outputs:
  results/external_validation/pathway_validation/
    tgfb_reproducible_module_gene_table.csv
    tgfb_reproducible_module_summary.csv
"""

from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PV = ROOT / "results/external_validation/pathway_validation"
STRICT = ROOT / "results/external_validation/three_dataset_strict_core_genes.csv"
AUDIT = PV / "tgfb_gene_level_audit.csv"
OVERLAP = PV / "tgfb_leading_edge_overlap.csv"

OUT_GENE = PV / "tgfb_reproducible_module_gene_table.csv"
OUT_SUM = PV / "tgfb_reproducible_module_summary.csv"


def main():
    audit = pd.read_csv(AUDIT)
    overlap = pd.read_csv(OVERLAP)
    strict = pd.read_csv(STRICT)

    module = overlap[
        overlap["leading_edge_pattern"] == "DISCOVERY_GSE237606_ONLY"
    ][["gene_symbol"]].drop_duplicates()

    df = module.merge(audit, on="gene_symbol", how="left")

    strict_set = set(strict["gene_symbol"].astype(str))
    df["in_three_dataset_strict_core"] = df["gene_symbol"].isin(strict_set)

    df["discovery_up"] = df["discovery_log2FC"] > 0
    df["GSE237606_up"] = df["GSE237606_log2FC"] > 0
    df["GSE121105_up"] = df["GSE121105_log2FC"] > 0

    df["GSE121105_reversed_vs_positive_module"] = df["GSE121105_log2FC"] < 0
    df["GSE121105_fdr05"] = (
        df["GSE121105_padj"].notna() & (df["GSE121105_padj"] < 0.05)
    )

    df = df.sort_values(
        ["in_three_dataset_strict_core", "GSE121105_reversed_vs_positive_module", "gene_symbol"],
        ascending=[False, False, True],
    )

    df.to_csv(OUT_GENE, index=False)

    n = len(df)
    reversed_n = int(df["GSE121105_reversed_vs_positive_module"].sum())
    reversed_sig_n = int(
        (
            df["GSE121105_reversed_vs_positive_module"]
            & df["GSE121105_fdr05"]
        ).sum()
    )
    strict_n = int(df["in_three_dataset_strict_core"].sum())

    summary = pd.DataFrame([{
        "reproducible_positive_tgfb_module_genes": n,
        "up_in_discovery": int(df["discovery_up"].sum()),
        "up_in_GSE237606": int(df["GSE237606_up"].sum()),
        "down_in_GSE121105": reversed_n,
        "down_in_GSE121105_and_fdr05": reversed_sig_n,
        "overlap_with_242_gene_strict_core": strict_n,
        "gse121105_reversal_fraction": reversed_n / n if n else np.nan,
    }])
    summary.to_csv(OUT_SUM, index=False)

    print("=" * 78)
    print("ATLAS — 00O REPRODUCIBLE TGF-BETA MODULE VALIDATION")
    print("=" * 78)
    print(summary.to_string(index=False))

    print("\nModule genes:")
    cols = [
        "gene_symbol",
        "discovery_log2FC",
        "GSE121105_log2FC",
        "GSE121105_padj",
        "GSE237606_log2FC",
        "in_three_dataset_strict_core",
        "GSE121105_reversed_vs_positive_module",
    ]
    print(df[cols].to_string(index=False))

    print("\nOutputs:")
    print(OUT_GENE)
    print(OUT_SUM)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
