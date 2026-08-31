#!/usr/bin/env python3
"""
ATLAS — 00L Targeted TGF-beta Ranked Validation

Runs preranked GSEA separately in:
- original discovery DE
- GSE121105
- GSE237606

This is more appropriate than strict-core ORA for testing whether the
TGF-beta pathway replicates as a coordinated genome-wide signal.

Outputs:
  results/external_validation/pathway_validation/
    discovery_hallmark_prerank.csv
    GSE121105_hallmark_prerank.csv
    GSE237606_hallmark_prerank.csv
    tgfb_ranked_validation_summary.csv
"""

from pathlib import Path
import re
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


def symbol_from_external(x):
    s = str(x).strip()
    if "|" in s:
        s = s.split("|", 1)[0].strip()
    return s


def prepare_rank(label, path):
    df = pd.read_csv(path)

    if label == "discovery":
        gene_col = "Gene name"
    else:
        gene_col = "gene_id"

    score_col = "stat" if "stat" in df.columns else "log2FoldChange"

    if label == "discovery":
        df["gene_symbol"] = df[gene_col].astype(str).str.strip()
    else:
        df["gene_symbol"] = df[gene_col].map(symbol_from_external)

    df["rank_score"] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=["gene_symbol", "rank_score"])
    df = df[df["gene_symbol"].str.len() > 0]

    # Keep strongest absolute statistic for duplicated symbols.
    df["abs_rank"] = df["rank_score"].abs()
    df = (
        df.sort_values("abs_rank", ascending=False)
          .drop_duplicates("gene_symbol", keep="first")
    )

    rnk = df[["gene_symbol", "rank_score"]].sort_values(
        "rank_score", ascending=False
    )
    return rnk


def run_one(label, path):
    print(f"\n[{label}]")
    rnk = prepare_rank(label, path)
    print(f"  ranked genes: {len(rnk)}")

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
    out = OUTDIR / f"{label}_hallmark_prerank.csv"
    res.to_csv(out, index=False)

    # GSEApy may name term column "Term".
    term_col = "Term" if "Term" in res.columns else res.columns[0]
    hit = res[
        res[term_col].astype(str).str.contains("TGF", case=False, na=False)
    ].copy()

    if hit.empty:
        print("  TGF-beta Hallmark result: NOT FOUND")
        return {
            "dataset": label,
            "term": "",
            "NES": None,
            "nominal_p": None,
            "FDR_q": None,
        }

    row = hit.iloc[0]
    nes = row.get("NES")
    p = row.get("NOM p-val", row.get("NOM p-val ", None))
    fdr = row.get("FDR q-val", row.get("FDR q-val ", None))

    print(
        f"  TGF-beta Hallmark: NES={nes}, nominal_p={p}, FDR_q={fdr}"
    )

    return {
        "dataset": label,
        "term": row.get(term_col),
        "NES": nes,
        "nominal_p": p,
        "FDR_q": fdr,
    }


def main():
    print("=" * 78)
    print("ATLAS — 00L TARGETED TGF-BETA RANKED VALIDATION")
    print("=" * 78)

    rows = []
    for label, path in INPUTS.items():
        rows.append(run_one(label, path))

    summary = pd.DataFrame(rows)
    out = OUTDIR / "tgfb_ranked_validation_summary.csv"
    summary.to_csv(out, index=False)

    print("\nSummary:")
    print(summary.to_string(index=False))
    print(f"\nOutput: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
