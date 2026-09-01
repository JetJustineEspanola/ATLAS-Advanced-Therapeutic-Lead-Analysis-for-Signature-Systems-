#!/usr/bin/env python3
"""
ATLAS — 00J v3 Strict Discovery + External Consensus Validation

The external consensus already requires:
- FDR < 0.05 in both external primary cohorts
- same direction in both
- |log2FC| >= 1 in both

This stage now additionally requires discovery significance, avoiding the
overstatement that direction-only agreement equals three-dataset validation.

Outputs:
- discovery_consensus_overlap.csv
- three_dataset_significant_validated_genes.csv
- three_dataset_strict_core_genes.csv
- discovery_consensus_summary.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DISCOVERY = ROOT / "results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv"
CONSENSUS = ROOT / "results/external_validation/consensus_resistance_signature.csv"
OUTDIR = ROOT / "results/external_validation"
OUTDIR.mkdir(parents=True, exist_ok=True)

DISCOVERY_FDR = 0.05
STRICT_DISCOVERY_ABS_LFC = 1.0


def norm(x):
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def main():
    d = pd.read_csv(DISCOVERY)
    c = pd.read_csv(CONSENSUS)

    gene_col = "Gene name"
    lfc_col = "log2FoldChange"
    padj_col = "padj"

    d["gene_symbol_norm"] = d[gene_col].map(norm)
    d["discovery_log2FC"] = pd.to_numeric(d[lfc_col], errors="coerce")
    d["discovery_padj"] = pd.to_numeric(d[padj_col], errors="coerce")

    # One annotated row per symbol; choose the most significant, then largest effect.
    d["abs_discovery_log2FC"] = d["discovery_log2FC"].abs()
    d = (
        d[d["gene_symbol_norm"] != ""]
        .sort_values(
            ["gene_symbol_norm", "discovery_padj", "abs_discovery_log2FC"],
            ascending=[True, True, False],
            na_position="last",
        )
        .drop_duplicates("gene_symbol_norm", keep="first")
    )

    c["gene_symbol_norm"] = c["gene_symbol"].map(norm)

    merged = c.merge(
        d[
            [
                "gene_symbol_norm",
                "discovery_log2FC",
                "discovery_padj",
                "abs_discovery_log2FC",
            ]
        ],
        on="gene_symbol_norm",
        how="left",
    )

    merged["present_in_discovery"] = merged["discovery_log2FC"].notna()

    ext_sign = np.where(
        merged["replication_direction"].eq("UP_IN_RESISTANT"),
        1,
        np.where(
            merged["replication_direction"].eq("DOWN_IN_RESISTANT"), -1, 0
        ),
    )
    disc_sign = np.sign(merged["discovery_log2FC"].fillna(0))

    merged["same_direction_as_discovery"] = (
        merged["present_in_discovery"] & (ext_sign == disc_sign)
    )

    merged["discovery_fdr05"] = (
        merged["discovery_padj"].notna()
        & (merged["discovery_padj"] < DISCOVERY_FDR)
    )

    # Validated in all three datasets by significance + direction.
    merged["three_dataset_significant_validated"] = (
        merged["same_direction_as_discovery"]
        & merged["discovery_fdr05"]
    )

    # Stricter core also requires large effect in discovery.
    merged["three_dataset_strict_core"] = (
        merged["three_dataset_significant_validated"]
        & (merged["abs_discovery_log2FC"] >= STRICT_DISCOVERY_ABS_LFC)
    )

    significant = merged[
        merged["three_dataset_significant_validated"]
    ].copy()

    strict = merged[
        merged["three_dataset_strict_core"]
    ].copy()

    merged = merged.sort_values(
        ["three_dataset_strict_core", "three_dataset_significant_validated", "consensus_rank"],
        ascending=[False, False, True],
    )
    significant = significant.sort_values("consensus_rank")
    strict = strict.sort_values("consensus_rank")

    all_out = OUTDIR / "discovery_consensus_overlap.csv"
    sig_out = OUTDIR / "three_dataset_significant_validated_genes.csv"
    strict_out = OUTDIR / "three_dataset_strict_core_genes.csv"
    sum_out = OUTDIR / "discovery_consensus_summary.csv"

    merged.to_csv(all_out, index=False)
    significant.to_csv(sig_out, index=False)
    strict.to_csv(strict_out, index=False)

    overlap_n = int(merged["present_in_discovery"].sum())
    dir_n = int(merged["same_direction_as_discovery"].sum())

    summary = pd.DataFrame([{
        "external_consensus_genes": len(merged),
        "present_in_discovery": overlap_n,
        "same_direction_as_discovery": dir_n,
        "directional_agreement_among_overlap": (
            dir_n / overlap_n if overlap_n else np.nan
        ),
        "three_dataset_significant_validated": len(significant),
        "three_dataset_strict_core_abs_log2fc_ge_1_all_three": len(strict),
        "strict_core_up": int(
            (strict["replication_direction"] == "UP_IN_RESISTANT").sum()
        ),
        "strict_core_down": int(
            (strict["replication_direction"] == "DOWN_IN_RESISTANT").sum()
        ),
    }])
    summary.to_csv(sum_out, index=False)

    print("=" * 78)
    print("ATLAS — 00J v3 STRICT THREE-DATASET VALIDATION")
    print("=" * 78)
    print(f"External consensus genes: {len(merged)}")
    print(f"Present in discovery: {overlap_n}")
    print(f"Same direction as discovery: {dir_n}")
    if overlap_n:
        print(f"Directional agreement among overlap: {dir_n / overlap_n:.1%}")
    print(
        f"Significant in discovery too (FDR<0.05 + same direction): {len(significant)}"
    )
    print(
        f"Strict core (also |discovery log2FC|>=1): {len(strict)}"
    )
    print(
        f"  Up: {(strict['replication_direction'] == 'UP_IN_RESISTANT').sum()}"
    )
    print(
        f"  Down: {(strict['replication_direction'] == 'DOWN_IN_RESISTANT').sum()}"
    )

    print("\nTop strict-core genes:")
    cols = [
        "consensus_rank",
        "gene_symbol",
        "replication_direction",
        "discovery_log2FC",
        "discovery_padj",
        "GSE121105_log2FC",
        "GSE121105_padj",
        "GSE237606_log2FC",
        "GSE237606_padj",
    ]
    if strict.empty:
        print("None")
    else:
        print(strict[cols].head(25).to_string(index=False))

    print("\nOutputs:")
    print(all_out)
    print(sig_out)
    print(strict_out)
    print(sum_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
