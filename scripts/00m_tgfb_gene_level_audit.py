#!/usr/bin/env python3
"""
ATLAS — 00M TGF-beta Gene-Level Audit Across Three Datasets

Purpose
-------
Explain why Hallmark TGF-beta Signaling is:
- positively enriched in discovery
- negatively enriched in GSE121105
- positively enriched in GSE237606

This script:
1. Retrieves the Hallmark TGF-beta gene set from Enrichr/MSigDB via gseapy.
2. Harmonizes gene symbols across the three DE tables.
3. Extracts per-gene log2FC, p-value, and adjusted p-value.
4. Classifies each gene as:
   - SAME_DIRECTION_ALL_THREE
   - GSE121105_FLIP
   - DISCOVERY_FLIP
   - GSE237606_FLIP
   - MIXED_OR_ZERO
   - MISSING
5. Ranks likely drivers of the GSE121105 sign reversal.
6. Summarizes agreement and discordance.

Inputs
------
results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv
results/external_validation/GSE121105_DE.csv
results/external_validation/GSE237606_DE.csv

Outputs
-------
results/external_validation/pathway_validation/
  tgfb_gene_level_audit.csv
  tgfb_gene_level_audit_summary.csv
  tgfb_gse121105_flip_drivers.csv
"""

from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DISCOVERY = ROOT / "results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv"
G121 = ROOT / "results/external_validation/GSE121105_DE.csv"
G237 = ROOT / "results/external_validation/GSE237606_DE.csv"
OUTDIR = ROOT / "results/external_validation/pathway_validation"
OUTDIR.mkdir(parents=True, exist_ok=True)

LIBRARY = "MSigDB_Hallmark_2020"
TARGET_TERM = "TGF-beta Signaling"

try:
    import gseapy as gp
except ImportError:
    print("ERROR: gseapy is not installed.")
    print("Install with: python -m pip install -U gseapy")
    raise SystemExit(2)


def norm_symbol(x):
    if x is None or pd.isna(x):
        return ""
    s = str(x).strip()
    if "|" in s:
        s = s.split("|", 1)[0].strip()
    return s.upper()


def load_discovery():
    df = pd.read_csv(DISCOVERY)
    if "Gene name" not in df.columns:
        raise RuntimeError("Discovery table is missing 'Gene name'.")

    out = pd.DataFrame({
        "gene_symbol": df["Gene name"].map(norm_symbol),
        "discovery_log2FC": pd.to_numeric(df["log2FoldChange"], errors="coerce"),
        "discovery_pvalue": pd.to_numeric(df["pvalue"], errors="coerce"),
        "discovery_padj": pd.to_numeric(df["padj"], errors="coerce"),
        "discovery_stat": pd.to_numeric(df["stat"], errors="coerce")
        if "stat" in df.columns else np.nan,
    })

    out["abs_lfc"] = out["discovery_log2FC"].abs()
    out = (
        out[out["gene_symbol"] != ""]
        .sort_values(
            ["gene_symbol", "discovery_padj", "abs_lfc"],
            ascending=[True, True, False],
            na_position="last",
        )
        .drop_duplicates("gene_symbol", keep="first")
        .drop(columns=["abs_lfc"])
    )
    return out


def load_external(path: Path, prefix: str):
    df = pd.read_csv(path)

    out = pd.DataFrame({
        "gene_symbol": df["gene_id"].map(norm_symbol),
        f"{prefix}_log2FC": pd.to_numeric(df["log2FoldChange"], errors="coerce"),
        f"{prefix}_pvalue": pd.to_numeric(df["pvalue"], errors="coerce"),
        f"{prefix}_padj": pd.to_numeric(df["padj"], errors="coerce"),
        f"{prefix}_stat": pd.to_numeric(df["stat"], errors="coerce")
        if "stat" in df.columns else np.nan,
    })

    out["abs_lfc"] = out[f"{prefix}_log2FC"].abs()
    out = (
        out[out["gene_symbol"] != ""]
        .sort_values(
            ["gene_symbol", f"{prefix}_padj", "abs_lfc"],
            ascending=[True, True, False],
            na_position="last",
        )
        .drop_duplicates("gene_symbol", keep="first")
        .drop(columns=["abs_lfc"])
    )
    return out


def fetch_tgfb_genes():
    libs = gp.get_library(name=LIBRARY, organism="human")
    if TARGET_TERM not in libs:
        # Be tolerant to minor naming changes.
        matches = [k for k in libs if "TGF" in k.upper()]
        if not matches:
            raise RuntimeError(
                f"Could not find TGF-beta pathway in {LIBRARY}. "
                f"Available TGF-like terms: {matches}"
            )
        term = matches[0]
    else:
        term = TARGET_TERM

    genes = sorted({norm_symbol(g) for g in libs[term] if norm_symbol(g)})
    return term, genes


def sign_of(x):
    if pd.isna(x) or x == 0:
        return 0
    return 1 if x > 0 else -1


def pattern(row):
    a = sign_of(row["discovery_log2FC"])
    b = sign_of(row["GSE121105_log2FC"])
    c = sign_of(row["GSE237606_log2FC"])

    if 0 in {a, b, c}:
        return "MIXED_OR_ZERO"

    if a == b == c:
        return "SAME_DIRECTION_ALL_THREE"

    if a == c and b != a:
        return "GSE121105_FLIP"

    if b == c and a != b:
        return "DISCOVERY_FLIP"

    if a == b and c != a:
        return "GSE237606_FLIP"

    return "MIXED_OR_ZERO"


def safe_neglog10(x):
    if pd.isna(x):
        return 0.0
    x = max(float(x), 1e-300)
    return -math.log10(x)


def main():
    print("=" * 78)
    print("ATLAS — 00M TGF-BETA GENE-LEVEL AUDIT")
    print("=" * 78)

    term, genes = fetch_tgfb_genes()
    print(f"Hallmark term: {term}")
    print(f"TGF-beta genes in set: {len(genes)}")

    d = load_discovery()
    a = load_external(G121, "GSE121105")
    b = load_external(G237, "GSE237606")

    base = pd.DataFrame({"gene_symbol": genes})

    merged = (
        base.merge(d, on="gene_symbol", how="left")
            .merge(a, on="gene_symbol", how="left")
            .merge(b, on="gene_symbol", how="left")
    )

    merged["present_all_three"] = merged[
        ["discovery_log2FC", "GSE121105_log2FC", "GSE237606_log2FC"]
    ].notna().all(axis=1)

    merged["direction_pattern"] = merged.apply(pattern, axis=1)

    merged["discovery_sig_fdr05"] = (
        merged["discovery_padj"].notna() & (merged["discovery_padj"] < 0.05)
    )
    merged["GSE121105_sig_fdr05"] = (
        merged["GSE121105_padj"].notna() & (merged["GSE121105_padj"] < 0.05)
    )
    merged["GSE237606_sig_fdr05"] = (
        merged["GSE237606_padj"].notna() & (merged["GSE237606_padj"] < 0.05)
    )

    # Driver score: prioritize genes that flip in GSE121105 while the other two
    # agree, and reward magnitude/significance in all three.
    sig_strength = (
        merged["discovery_padj"].map(safe_neglog10)
        + merged["GSE121105_padj"].map(safe_neglog10)
        + merged["GSE237606_padj"].map(safe_neglog10)
    )

    effect_strength = (
        merged[["discovery_log2FC", "GSE121105_log2FC", "GSE237606_log2FC"]]
        .abs()
        .min(axis=1)
        .fillna(0.0)
    )

    merged["flip_driver_score"] = np.where(
        merged["direction_pattern"].eq("GSE121105_FLIP"),
        sig_strength * effect_strength,
        0.0,
    )

    merged = merged.sort_values(
        ["direction_pattern", "flip_driver_score", "gene_symbol"],
        ascending=[True, False, True],
    )

    out_all = OUTDIR / "tgfb_gene_level_audit.csv"
    merged.to_csv(out_all, index=False)

    flips = merged[
        merged["direction_pattern"] == "GSE121105_FLIP"
    ].sort_values("flip_driver_score", ascending=False)

    out_flips = OUTDIR / "tgfb_gse121105_flip_drivers.csv"
    flips.to_csv(out_flips, index=False)

    present = merged[merged["present_all_three"]].copy()

    summary_rows = []

    if not present.empty:
        counts = present["direction_pattern"].value_counts()

        same = int(counts.get("SAME_DIRECTION_ALL_THREE", 0))
        g121_flip = int(counts.get("GSE121105_FLIP", 0))
        disc_flip = int(counts.get("DISCOVERY_FLIP", 0))
        g237_flip = int(counts.get("GSE237606_FLIP", 0))
        mixed = int(counts.get("MIXED_OR_ZERO", 0))

        summary_rows.append({
            "hallmark_term": term,
            "pathway_genes_total": len(genes),
            "genes_present_all_three": len(present),
            "same_direction_all_three": same,
            "gse121105_flip": g121_flip,
            "discovery_flip": disc_flip,
            "gse237606_flip": g237_flip,
            "mixed_or_zero": mixed,
            "gse121105_flip_fraction_present": g121_flip / len(present),
        })

    summary = pd.DataFrame(summary_rows)
    out_summary = OUTDIR / "tgfb_gene_level_audit_summary.csv"
    summary.to_csv(out_summary, index=False)

    print(f"Genes present in all three: {len(present)}")

    if not present.empty:
        print("\nDirection-pattern counts:")
        print(present["direction_pattern"].value_counts().to_string())

    print("\nTop GSE121105 flip drivers:")
    cols = [
        "gene_symbol",
        "discovery_log2FC",
        "discovery_padj",
        "GSE121105_log2FC",
        "GSE121105_padj",
        "GSE237606_log2FC",
        "GSE237606_padj",
        "flip_driver_score",
    ]
    if flips.empty:
        print("None")
    else:
        print(flips[cols].head(25).to_string(index=False))

    print("\nInterpretation guide:")
    print("- SAME_DIRECTION_ALL_THREE: consistent TGF-beta gene behavior.")
    print("- GSE121105_FLIP: discovery + GSE237606 agree, GSE121105 reverses.")
    print("- A large GSE121105_FLIP group supports context-dependent pathway direction.")
    print("- A small flip group dominated by a few genes suggests pathway NES may be driver-sensitive.")

    print("\nOutputs:")
    print(out_all)
    print(out_flips)
    print(out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
