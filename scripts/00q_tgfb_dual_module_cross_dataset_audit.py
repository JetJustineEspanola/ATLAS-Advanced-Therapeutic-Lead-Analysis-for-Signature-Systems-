#!/usr/bin/env python3
"""
ATLAS — 00Q TGF-beta Dual-Module Cross-Dataset Audit

Tests whether trastuzumab resistance is better described as coordinated
TGF-beta remodeling with:
  1) POSITIVE_TGFB_16 activation module
  2) NEGATIVE_TGFB_7 repression module

Uses per-gene DE statistics from:
- discovery
- GSE121105
- GSE237606

Outputs:
  results/external_validation/pathway_validation/
    tgfb_dual_module_gene_effects.csv
    tgfb_dual_module_cross_dataset_summary.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
PV = ROOT / "results/external_validation/pathway_validation"
PV.mkdir(parents=True, exist_ok=True)

DISCOVERY = ROOT / "results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv"
G121 = ROOT / "results/external_validation/GSE121105_DE.csv"
G237 = ROOT / "results/external_validation/GSE237606_DE.csv"

POSITIVE = {
    "ACVR1","ARID4B","FURIN","HDAC1","LTBP2","MAP3K7","NOG","RAB31",
    "RHOA","SMAD1","SMURF2","SPTBN1","TGFB1","TGFBR1","TRIM33","XIAP"
}
NEGATIVE = {"HIPK2","ID1","ID3","SKI","SLC20A1","SMAD7","THBS1"}


def norm(x):
    if x is None or pd.isna(x):
        return ""
    s = str(x).strip()
    if "|" in s:
        s = s.split("|", 1)[0]
    return s.upper()


def load_discovery():
    d = pd.read_csv(DISCOVERY)
    out = pd.DataFrame({
        "gene_symbol": d["Gene name"].map(norm),
        "log2FC": pd.to_numeric(d["log2FoldChange"], errors="coerce"),
        "padj": pd.to_numeric(d["padj"], errors="coerce"),
    })
    out["abs_lfc"] = out["log2FC"].abs()
    return (
        out[out["gene_symbol"] != ""]
        .sort_values(["gene_symbol","padj","abs_lfc"], ascending=[True,True,False], na_position="last")
        .drop_duplicates("gene_symbol")
        .drop(columns="abs_lfc")
    )


def load_external(path):
    d = pd.read_csv(path)
    out = pd.DataFrame({
        "gene_symbol": d["gene_id"].map(norm),
        "log2FC": pd.to_numeric(d["log2FoldChange"], errors="coerce"),
        "padj": pd.to_numeric(d["padj"], errors="coerce"),
    })
    out["abs_lfc"] = out["log2FC"].abs()
    return (
        out[out["gene_symbol"] != ""]
        .sort_values(["gene_symbol","padj","abs_lfc"], ascending=[True,True,False], na_position="last")
        .drop_duplicates("gene_symbol")
        .drop(columns="abs_lfc")
    )


def summarize(dataset, df, module_name, genes):
    sub = df[df["gene_symbol"].isin(genes)].copy()
    vals = sub["log2FC"].dropna()

    if len(vals) >= 2 and not np.allclose(vals, 0):
        try:
            _, p = wilcoxon(vals)
        except Exception:
            p = np.nan
    else:
        p = np.nan

    return sub.assign(dataset=dataset, module=module_name), {
        "dataset": dataset,
        "module": module_name,
        "genes_present": len(sub),
        "mean_log2FC": vals.mean() if len(vals) else np.nan,
        "median_log2FC": vals.median() if len(vals) else np.nan,
        "up_genes": int((vals > 0).sum()),
        "down_genes": int((vals < 0).sum()),
        "fdr05_genes": int((sub["padj"] < 0.05).fillna(False).sum()),
        "wilcoxon_vs_zero_p": p,
    }


def main():
    datasets = {
        "discovery": load_discovery(),
        "GSE121105": load_external(G121),
        "GSE237606": load_external(G237),
    }

    all_rows = []
    summaries = []

    print("=" * 78)
    print("ATLAS — 00Q TGF-BETA DUAL-MODULE CROSS-DATASET AUDIT")
    print("=" * 78)

    for dataset, df in datasets.items():
        for module_name, genes in [
            ("POSITIVE_TGFB_16", POSITIVE),
            ("NEGATIVE_TGFB_7", NEGATIVE),
        ]:
            rows, summary = summarize(dataset, df, module_name, genes)
            all_rows.append(rows)
            summaries.append(summary)

            print(
                f"{dataset:10s} {module_name:18s} "
                f"mean={summary['mean_log2FC']:.3f} "
                f"median={summary['median_log2FC']:.3f} "
                f"up={summary['up_genes']} down={summary['down_genes']} "
                f"FDR<0.05={summary['fdr05_genes']} "
                f"Wilcoxon p={summary['wilcoxon_vs_zero_p']:.4g}"
            )

    genes_out = PV / "tgfb_dual_module_gene_effects.csv"
    summary_out = PV / "tgfb_dual_module_cross_dataset_summary.csv"

    pd.concat(all_rows, ignore_index=True).to_csv(genes_out, index=False)
    sdf = pd.DataFrame(summaries)
    sdf.to_csv(summary_out, index=False)

    print("\nOutputs:")
    print(genes_out)
    print(summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
