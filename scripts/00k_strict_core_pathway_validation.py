#!/usr/bin/env python3
"""
ATLAS — 00K Strict-Core Pathway Validation

Runs direction-aware pathway over-representation analysis on the final
three-dataset strict core:
  - UP_IN_RESISTANT genes
  - DOWN_IN_RESISTANT genes

Uses Enrichr Hallmark 2020 via gseapy.

Outputs:
  results/external_validation/pathway_validation/
    strict_core_up_hallmark2020.csv
    strict_core_down_hallmark2020.csv
    strict_core_pathway_summary.csv
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INFILE = ROOT / "results/external_validation/three_dataset_strict_core_genes.csv"
OUTDIR = ROOT / "results/external_validation/pathway_validation"
OUTDIR.mkdir(parents=True, exist_ok=True)

try:
    import gseapy as gp
except ImportError:
    print("ERROR: gseapy is not installed.")
    print("Install with: python -m pip install -U gseapy")
    raise SystemExit(2)


LIBRARY = "MSigDB_Hallmark_2020"


def run_enrichr(genes, label):
    if not genes:
        return pd.DataFrame()

    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=LIBRARY,
        organism="human",
        outdir=None,
        cutoff=1.0,
    )

    res = enr.results.copy()

    if "Adjusted P-value" in res.columns:
        res = res.sort_values("Adjusted P-value", ascending=True)

    out = OUTDIR / f"strict_core_{label.lower()}_hallmark2020.csv"
    res.to_csv(out, index=False)
    return res


def main():
    df = pd.read_csv(INFILE)

    up = (
        df.loc[
            df["replication_direction"] == "UP_IN_RESISTANT",
            "gene_symbol",
        ]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    down = (
        df.loc[
            df["replication_direction"] == "DOWN_IN_RESISTANT",
            "gene_symbol",
        ]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    print("=" * 78)
    print("ATLAS — 00K STRICT-CORE PATHWAY VALIDATION")
    print("=" * 78)
    print(f"Strict-core UP genes: {len(up)}")
    print(f"Strict-core DOWN genes: {len(down)}")

    up_res = run_enrichr(up, "UP")
    down_res = run_enrichr(down, "DOWN")

    rows = []

    for direction, res in [("UP_IN_RESISTANT", up_res), ("DOWN_IN_RESISTANT", down_res)]:
        if res.empty:
            continue

        pcol = "Adjusted P-value" if "Adjusted P-value" in res.columns else None
        term_col = "Term" if "Term" in res.columns else res.columns[0]

        top = res.head(10).copy()

        print(f"\nTop Hallmark pathways — {direction}:")
        cols = [term_col]
        if pcol:
            cols.append(pcol)
        if "Combined Score" in top.columns:
            cols.append("Combined Score")
        if "Genes" in top.columns:
            cols.append("Genes")
        print(top[cols].to_string(index=False))

        for _, r in top.iterrows():
            rows.append({
                "direction": direction,
                "term": r.get(term_col),
                "adjusted_pvalue": r.get(pcol) if pcol else None,
                "combined_score": r.get("Combined Score"),
                "genes": r.get("Genes"),
            })

    summary = pd.DataFrame(rows)
    summary_out = OUTDIR / "strict_core_pathway_summary.csv"
    summary.to_csv(summary_out, index=False)

    print("\nOutputs:")
    print(OUTDIR / "strict_core_up_hallmark2020.csv")
    print(OUTDIR / "strict_core_down_hallmark2020.csv")
    print(summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
