#!/usr/bin/env python3
"""
ATLAS — 00H Cross-Study DEG Concordance

Harmonizes gene symbols from the two primary validation DE tables and measures
replication across studies without pretending the studies are identical.

Primary replication rule:
  - gene present in both studies
  - same log2FC direction
  - FDR < 0.05 in both studies

Also reports:
  - shared-gene log2FC Spearman correlation
  - directional agreement
  - genes significant in either/both studies
  - a ranked consensus table

Duplicate symbol handling:
  If multiple rows map to the same symbol within one study, keep the row with
  the highest baseMean (deterministic and avoids cherry-picking by p-value).

Outputs:
  results/external_validation/cross_study_concordance.csv
  results/external_validation/replicated_primary_DEGs.csv
  results/external_validation/cross_study_concordance_summary.csv
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
INDIR = ROOT / "results/external_validation"
OUTDIR = INDIR

FILES = {
    "GSE121105": INDIR / "GSE121105_DE.csv",
    "GSE237606": INDIR / "GSE237606_DE.csv",
}


def canonical_symbol(gene_id: str, all_ids=None) -> str:
    s = str(gene_id).strip()
    if "|" in s:
        s = s.split("|", 1)[0].strip()

    # Only strip a trailing _number if the unsuffixed base symbol exists
    # somewhere in the same table (e.g. ZYX and ZYX_1).
    if all_ids is not None:
        m = re.match(r"^(.+)_\d+$", s)
        if m and m.group(1) in all_ids:
            s = m.group(1)

    return s


def load_study(accession, path):
    df = pd.read_csv(path)

    raw_ids = []
    for x in df["gene_id"].astype(str):
        raw = x.split("|", 1)[0].strip()
        raw_ids.append(raw)
    raw_set = set(raw_ids)

    df["gene_symbol"] = [
        canonical_symbol(x, raw_set) for x in df["gene_id"].astype(str)
    ]

    for c in ["baseMean", "log2FoldChange", "pvalue", "padj"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Deterministic duplicate resolution: highest mean abundance.
    df = (
        df.sort_values(["gene_symbol", "baseMean"], ascending=[True, False])
          .drop_duplicates("gene_symbol", keep="first")
          .copy()
    )

    keep = ["gene_symbol", "baseMean", "log2FoldChange", "pvalue", "padj"]
    df = df[keep].rename(columns={
        "baseMean": f"{accession}_baseMean",
        "log2FoldChange": f"{accession}_log2FC",
        "pvalue": f"{accession}_pvalue",
        "padj": f"{accession}_padj",
    })

    return df


def main():
    print("=" * 78)
    print("ATLAS — 00H CROSS-STUDY DEG CONCORDANCE")
    print("=" * 78)

    a = load_study("GSE121105", FILES["GSE121105"])
    b = load_study("GSE237606", FILES["GSE237606"])

    merged = a.merge(b, on="gene_symbol", how="inner")

    lfc1 = merged["GSE121105_log2FC"]
    lfc2 = merged["GSE237606_log2FC"]

    merged["same_direction"] = (
        ((lfc1 > 0) & (lfc2 > 0))
        | ((lfc1 < 0) & (lfc2 < 0))
    )

    merged["GSE121105_fdr05"] = merged["GSE121105_padj"].notna() & (
        merged["GSE121105_padj"] < 0.05
    )
    merged["GSE237606_fdr05"] = merged["GSE237606_padj"].notna() & (
        merged["GSE237606_padj"] < 0.05
    )

    merged["significant_both"] = (
        merged["GSE121105_fdr05"] & merged["GSE237606_fdr05"]
    )
    merged["replicated_primary"] = (
        merged["same_direction"] & merged["significant_both"]
    )

    merged["replication_direction"] = np.where(
        merged["replicated_primary"] & (lfc1 > 0),
        "UP_IN_RESISTANT",
        np.where(
            merged["replicated_primary"] & (lfc1 < 0),
            "DOWN_IN_RESISTANT",
            ""
        ),
    )

    # Conservative consensus score for ranking only, not a probability.
    # Rewards significance in both studies and effect-size consistency.
    eps = 1e-300
    p1 = merged["GSE121105_padj"].clip(lower=eps).fillna(1.0)
    p2 = merged["GSE237606_padj"].clip(lower=eps).fillna(1.0)
    merged["consensus_score"] = (
        (-np.log10(p1) - np.log10(p2))
        * np.sqrt(np.abs(lfc1 * lfc2))
        * np.where(merged["same_direction"], 1.0, -1.0)
    )

    # Correlations.
    valid = merged[["GSE121105_log2FC", "GSE237606_log2FC"]].dropna()
    rho_all, p_all = spearmanr(
        valid["GSE121105_log2FC"],
        valid["GSE237606_log2FC"],
    )

    sig_either = merged[
        merged["GSE121105_fdr05"] | merged["GSE237606_fdr05"]
    ][["GSE121105_log2FC", "GSE237606_log2FC"]].dropna()

    if len(sig_either) >= 3:
        rho_sig, p_sig = spearmanr(
            sig_either["GSE121105_log2FC"],
            sig_either["GSE237606_log2FC"],
        )
    else:
        rho_sig, p_sig = np.nan, np.nan

    replicated = merged[merged["replicated_primary"]].copy()
    replicated = replicated.sort_values(
        ["consensus_score"],
        ascending=False,
    )

    merged = merged.sort_values(
        ["replicated_primary", "consensus_score"],
        ascending=[False, False],
    )

    all_out = OUTDIR / "cross_study_concordance.csv"
    rep_out = OUTDIR / "replicated_primary_DEGs.csv"
    summary_out = OUTDIR / "cross_study_concordance_summary.csv"

    merged.to_csv(all_out, index=False)
    replicated.to_csv(rep_out, index=False)

    shared_n = len(merged)
    same_n = int(merged["same_direction"].sum())
    both_sig_n = int(merged["significant_both"].sum())
    rep_n = len(replicated)
    rep_up = int((replicated["replication_direction"] == "UP_IN_RESISTANT").sum())
    rep_down = int((replicated["replication_direction"] == "DOWN_IN_RESISTANT").sum())

    summary = pd.DataFrame([{
        "shared_genes": shared_n,
        "same_direction_genes": same_n,
        "direction_agreement_fraction": same_n / shared_n if shared_n else np.nan,
        "significant_both_fdr05": both_sig_n,
        "replicated_same_direction_fdr05_both": rep_n,
        "replicated_up": rep_up,
        "replicated_down": rep_down,
        "spearman_rho_all_shared": rho_all,
        "spearman_p_all_shared": p_all,
        "spearman_rho_significant_either": rho_sig,
        "spearman_p_significant_either": p_sig,
    }])
    summary.to_csv(summary_out, index=False)

    print(f"\nShared harmonized genes: {shared_n}")
    print(f"Same-direction genes: {same_n} ({same_n/shared_n:.1%})")
    print(f"FDR<0.05 in both studies: {both_sig_n}")
    print(f"Replicated (same direction + FDR<0.05 in both): {rep_n}")
    print(f"  Up in resistant: {rep_up}")
    print(f"  Down in resistant: {rep_down}")
    print(f"Spearman rho, all shared log2FC: {rho_all:.3f} (p={p_all:.3g})")
    print(f"Spearman rho, significant in either: {rho_sig:.3f} (p={p_sig:.3g})")

    print("\nTop replicated genes:")
    cols = [
        "gene_symbol",
        "GSE121105_log2FC",
        "GSE121105_padj",
        "GSE237606_log2FC",
        "GSE237606_padj",
        "replication_direction",
        "consensus_score",
    ]
    print(replicated[cols].head(25).to_string(index=False))

    print("\nOutputs:")
    print(all_out)
    print(rep_out)
    print(summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
