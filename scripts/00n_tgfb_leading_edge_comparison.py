#!/usr/bin/env python3
"""
ATLAS — 00N TGF-beta Leading-Edge Comparison Across Three Datasets

Purpose
-------
Identify the leading-edge genes driving Hallmark TGF-beta Signaling in:
- discovery
- GSE121105
- GSE237606

This explains which TGF-beta subprogram is consistently associated with
trastuzumab resistance and which component is reversed in GSE121105.

Outputs
-------
results/external_validation/pathway_validation/
  discovery_tgfb_leading_edge.csv
  GSE121105_tgfb_leading_edge.csv
  GSE237606_tgfb_leading_edge.csv
  tgfb_leading_edge_overlap.csv
  tgfb_leading_edge_summary.csv
"""

from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import pandas as pd
import gseapy as gp

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/external_validation/pathway_validation"
OUTDIR.mkdir(parents=True, exist_ok=True)

INPUTS = {
    "discovery": ROOT / "results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv",
    "GSE121105": ROOT / "results/external_validation/GSE121105_DE.csv",
    "GSE237606": ROOT / "results/external_validation/GSE237606_DE.csv",
}

LIBRARY = "MSigDB_Hallmark_2020"
TARGET_TERM = "TGF-beta Signaling"


def norm_symbol(x):
    if x is None or pd.isna(x):
        return ""
    s = str(x).strip()
    if "|" in s:
        s = s.split("|", 1)[0].strip()
    return s.upper()


def prep_rank(label: str, path: Path):
    df = pd.read_csv(path)

    if label == "discovery":
        gene_col = "Gene name"
    else:
        gene_col = "gene_id"

    score_col = "stat" if "stat" in df.columns else "log2FoldChange"

    df["gene_symbol"] = df[gene_col].map(norm_symbol)
    df["rank_score"] = pd.to_numeric(df[score_col], errors="coerce")

    df = df.dropna(subset=["gene_symbol", "rank_score"])
    df = df[df["gene_symbol"] != ""]

    df["abs_rank"] = df["rank_score"].abs()

    df = (
        df.sort_values("abs_rank", ascending=False)
          .drop_duplicates("gene_symbol", keep="first")
    )

    return df[["gene_symbol", "rank_score"]].sort_values(
        "rank_score", ascending=False
    )


def parse_leading_edge(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []

    s = str(value).strip()

    # GSEApy usually returns semicolon-delimited genes.
    if ";" in s:
        parts = s.split(";")
    elif "," in s:
        parts = s.split(",")
    else:
        parts = s.split()

    return [norm_symbol(x) for x in parts if norm_symbol(x)]


def run_dataset(label: str, path: Path):
    print(f"\n[{label}]")

    rnk = prep_rank(label, path)

    pre = gp.prerank(
        rnk=rnk,
        gene_sets=LIBRARY,
        outdir=None,
        min_size=10,
        max_size=500,
        permutation_num=1000,
        seed=42,
        verbose=False,
    )

    res = pre.res2d.copy()
    term_col = "Term" if "Term" in res.columns else res.columns[0]

    hit = res[
        res[term_col].astype(str).str.contains("TGF", case=False, na=False)
    ].copy()

    if hit.empty:
        raise RuntimeError(f"{label}: TGF-beta Hallmark result not found")

    row = hit.iloc[0]

    nes = float(row["NES"])
    p = row.get("NOM p-val", np.nan)
    fdr = row.get("FDR q-val", np.nan)

    ledge_raw = row.get("Lead_genes", row.get("ledge_genes", ""))
    genes = parse_leading_edge(ledge_raw)

    # If gseapy column name differs, inspect any column containing lead.
    if not genes:
        for col in res.columns:
            if "lead" in str(col).lower():
                genes = parse_leading_edge(row.get(col))
                if genes:
                    break

    print(f"  NES={nes:.3f}, FDR={fdr}")
    print(f"  leading-edge genes: {len(genes)}")

    if genes:
        print("  " + ", ".join(genes[:30]))

    # Attach rank scores to leading-edge genes.
    rank_map = dict(zip(rnk["gene_symbol"], rnk["rank_score"]))
    out_df = pd.DataFrame({
        "gene_symbol": genes,
        "rank_score": [rank_map.get(g, np.nan) for g in genes],
        "dataset": label,
        "NES": nes,
        "FDR_q": fdr,
    })

    out = OUTDIR / f"{label}_tgfb_leading_edge.csv"
    out_df.to_csv(out, index=False)

    return {
        "dataset": label,
        "NES": nes,
        "FDR_q": fdr,
        "leading_edge_genes": genes,
        "leading_edge_n": len(genes),
    }


def main():
    print("=" * 78)
    print("ATLAS — 00N TGF-BETA LEADING-EDGE COMPARISON")
    print("=" * 78)

    results = [
        run_dataset(label, path)
        for label, path in INPUTS.items()
    ]

    by_dataset = {
        r["dataset"]: set(r["leading_edge_genes"])
        for r in results
    }

    all_genes = sorted(set().union(*by_dataset.values()))

    rows = []
    for gene in all_genes:
        d = gene in by_dataset["discovery"]
        a = gene in by_dataset["GSE121105"]
        b = gene in by_dataset["GSE237606"]

        if d and a and b:
            pattern = "LEADING_EDGE_ALL_THREE"
        elif d and b and not a:
            pattern = "DISCOVERY_GSE237606_ONLY"
        elif a and not d and not b:
            pattern = "GSE121105_ONLY"
        elif a and d and not b:
            pattern = "DISCOVERY_GSE121105_ONLY"
        elif a and b and not d:
            pattern = "GSE121105_GSE237606_ONLY"
        elif d and not a and not b:
            pattern = "DISCOVERY_ONLY"
        elif b and not d and not a:
            pattern = "GSE237606_ONLY"
        else:
            pattern = "OTHER"

        rows.append({
            "gene_symbol": gene,
            "in_discovery_leading_edge": d,
            "in_GSE121105_leading_edge": a,
            "in_GSE237606_leading_edge": b,
            "leading_edge_pattern": pattern,
        })

    overlap = pd.DataFrame(rows)

    # Add gene-level log2FC audit information if available.
    audit_path = OUTDIR / "tgfb_gene_level_audit.csv"
    if audit_path.exists():
        audit = pd.read_csv(audit_path)
        keep = [
            "gene_symbol",
            "direction_pattern",
            "discovery_log2FC",
            "GSE121105_log2FC",
            "GSE237606_log2FC",
        ]
        keep = [c for c in keep if c in audit.columns]
        overlap = overlap.merge(
            audit[keep],
            on="gene_symbol",
            how="left",
        )

    overlap = overlap.sort_values(
        ["leading_edge_pattern", "gene_symbol"]
    )

    out_overlap = OUTDIR / "tgfb_leading_edge_overlap.csv"
    overlap.to_csv(out_overlap, index=False)

    # Summary metrics.
    d = by_dataset["discovery"]
    a = by_dataset["GSE121105"]
    b = by_dataset["GSE237606"]

    summary = pd.DataFrame([{
        "discovery_NES": results[0]["NES"],
        "discovery_FDR": results[0]["FDR_q"],
        "discovery_leading_edge_n": len(d),

        "GSE121105_NES": results[1]["NES"],
        "GSE121105_FDR": results[1]["FDR_q"],
        "GSE121105_leading_edge_n": len(a),

        "GSE237606_NES": results[2]["NES"],
        "GSE237606_FDR": results[2]["FDR_q"],
        "GSE237606_leading_edge_n": len(b),

        "leading_edge_all_three": len(d & a & b),
        "discovery_GSE237606_overlap": len(d & b),
        "discovery_GSE237606_only": len((d & b) - a),
        "GSE121105_only": len(a - d - b),
        "GSE121105_overlap_discovery": len(a & d),
        "GSE121105_overlap_GSE237606": len(a & b),
    }])

    out_summary = OUTDIR / "tgfb_leading_edge_summary.csv"
    summary.to_csv(out_summary, index=False)

    print("\nLeading-edge overlap summary:")
    print(summary.to_string(index=False))

    print("\nKey categories:")

    cats = [
        "LEADING_EDGE_ALL_THREE",
        "DISCOVERY_GSE237606_ONLY",
        "GSE121105_ONLY",
    ]

    for cat in cats:
        genes = overlap.loc[
            overlap["leading_edge_pattern"] == cat,
            "gene_symbol",
        ].tolist()

        print(f"\n{cat} ({len(genes)}):")
        if genes:
            print(", ".join(genes[:50]))
        else:
            print("None")

    print("\nInterpretation:")
    print("- DISCOVERY_GSE237606_ONLY genes likely drive the reproducible positive TGF-beta program.")
    print("- GSE121105_ONLY genes likely drive the negative/context-specific TGF-beta program.")
    print("- LEADING_EDGE_ALL_THREE genes are shared pathway components, but their expression direction must be checked in the gene-level audit.")

    print("\nOutputs:")
    for r in results:
        print(OUTDIR / f"{r['dataset']}_tgfb_leading_edge.csv")
    print(out_overlap)
    print(out_summary)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
