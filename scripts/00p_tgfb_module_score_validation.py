#!/usr/bin/env python3
"""
ATLAS — 00P TGF-beta Module Score Validation

Computes sample-level expression scores for:
1) the 16-gene reproducible positive TGF-beta leading-edge module
2) the 7-gene GSE121105-specific negative leading-edge module

Datasets:
- GSE121105
- GSE237606

Method:
- library-size normalize raw counts to CPM
- log2(CPM + 1)
- z-score each gene within each dataset
- module score = mean z-score across available module genes
- compare resistant vs sensitive/parental
- for GSE237606, report overall and stratum-adjusted paired mean difference

Outputs:
  results/external_validation/pathway_validation/
    tgfb_module_sample_scores.csv
    tgfb_module_score_summary.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_rel

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/validation_expression"
PV = ROOT / "results/external_validation/pathway_validation"
PV.mkdir(parents=True, exist_ok=True)

POSITIVE = {
    "ACVR1","ARID4B","FURIN","HDAC1","LTBP2","MAP3K7","NOG","RAB31",
    "RHOA","SMAD1","SMURF2","SPTBN1","TGFB1","TGFBR1","TRIM33","XIAP"
}
NEGATIVE_121105 = {"HIPK2","ID1","ID3","SKI","SLC20A1","SMAD7","THBS1"}


def norm_symbol(x):
    s = str(x).strip()
    if "|" in s:
        s = s.split("|", 1)[0].strip()
    return s.upper()


def read_counts(path):
    if path.name.endswith(".csv.gz"):
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(path, sep="\t")

    gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene_id"})
    df["gene_symbol"] = df["gene_id"].map(norm_symbol)

    sample_cols = [c for c in df.columns if c not in {"gene_id", "gene_symbol"}]
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    # collapse duplicated symbols by summing counts
    df = df.groupby("gene_symbol", as_index=False)[sample_cols].sum()
    return df


def logcpm_z(df):
    sample_cols = [c for c in df.columns if c != "gene_symbol"]
    counts = df.set_index("gene_symbol")[sample_cols].astype(float)

    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    logcpm = np.log2(cpm + 1)

    mean = logcpm.mean(axis=1)
    sd = logcpm.std(axis=1, ddof=0).replace(0, np.nan)
    z = logcpm.sub(mean, axis=0).div(sd, axis=0).fillna(0)
    return z


def module_scores(z, genes):
    present = sorted(set(z.index) & set(genes))
    if not present:
        return pd.Series(index=z.columns, dtype=float), []
    return z.loc[present].mean(axis=0), present


def load_design(acc):
    return pd.read_csv(DATA / f"{acc}_primary_design.csv")


def map_design(acc, design):
    design = design.copy()
    design["count_column"] = design["description"].astype(str)

    if acc == "GSE121105":
        design = design[
            design["included"].astype(str).str.lower().isin(["true", "1"])
        ].copy()

    return design


def summarize(acc, scores, design, module_name, genes_present):
    d = design.set_index("count_column").loc[scores.index].copy()
    d["score"] = scores.values
    d["dataset"] = acc
    d["module"] = module_name

    r = d[d["phenotype"] == "RESISTANT"]["score"]
    s = d[d["phenotype"] == "SENSITIVE_OR_PARENTAL"]["score"]

    if len(r) and len(s):
        _, p_mwu = mannwhitneyu(r, s, alternative="two-sided")
        mean_diff = r.mean() - s.mean()
    else:
        p_mwu = np.nan
        mean_diff = np.nan

    stratum_diff = np.nan
    stratum_p = np.nan

    if acc == "GSE237606":
        d["stratum"] = d["timepoint"].astype(str) + "_" + d["stimulus"].astype(str)

        piv = (
            d.groupby(["stratum", "phenotype"])["score"]
            .mean()
            .unstack("phenotype")
            .dropna()
        )

        if {"RESISTANT", "SENSITIVE_OR_PARENTAL"}.issubset(piv.columns):
            diffs = piv["RESISTANT"] - piv["SENSITIVE_OR_PARENTAL"]
            stratum_diff = diffs.mean()
            if len(diffs) >= 2:
                _, stratum_p = ttest_rel(
                    piv["RESISTANT"],
                    piv["SENSITIVE_OR_PARENTAL"],
                )

    return d.reset_index(), {
        "dataset": acc,
        "module": module_name,
        "genes_present_n": len(genes_present),
        "genes_present": ";".join(genes_present),
        "resistant_n": len(r),
        "sensitive_parental_n": len(s),
        "mean_score_resistant": r.mean() if len(r) else np.nan,
        "mean_score_sensitive": s.mean() if len(s) else np.nan,
        "mean_difference_resistant_minus_sensitive": mean_diff,
        "mannwhitney_p": p_mwu,
        "stratum_adjusted_mean_difference": stratum_diff,
        "stratum_paired_t_p": stratum_p,
    }


def main():
    print("=" * 78)
    print("ATLAS — 00P TGF-BETA MODULE SCORE VALIDATION")
    print("=" * 78)

    configs = {
        "GSE121105": DATA / "GSE121105/GSE121105_geneCount.csv.gz",
        "GSE237606": DATA / "GSE237606/GSE237606_RawCounts.txt.gz",
    }

    sample_rows = []
    summaries = []

    for acc, path in configs.items():
        counts = read_counts(path)
        z = logcpm_z(counts)
        design = map_design(acc, load_design(acc))

        print(f"\n[{acc}]")

        for module_name, genes in [
            ("POSITIVE_TGFB_16", POSITIVE),
            ("GSE121105_NEGATIVE_TGFB_7", NEGATIVE_121105),
        ]:
            scores, present = module_scores(z, genes)

            # keep only samples in the design
            keep = [x for x in design["count_column"] if x in scores.index]
            scores = scores.loc[keep]

            sample_df, summary = summarize(
                acc, scores, design, module_name, present
            )
            sample_rows.append(sample_df)
            summaries.append(summary)

            print(
                f"  {module_name}: genes={len(present)}, "
                f"Δ(R-S)={summary['mean_difference_resistant_minus_sensitive']:.3f}, "
                f"MWU p={summary['mannwhitney_p']:.4g}"
            )

            if acc == "GSE237606":
                print(
                    f"    stratum-adjusted Δ={summary['stratum_adjusted_mean_difference']:.3f}, "
                    f"paired p={summary['stratum_paired_t_p']:.4g}"
                )

    sample_out = PV / "tgfb_module_sample_scores.csv"
    summary_out = PV / "tgfb_module_score_summary.csv"

    pd.concat(sample_rows, ignore_index=True).to_csv(sample_out, index=False)
    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(summary_out, index=False)

    print("\nSummary:")
    print(summary_df.to_string(index=False))

    print("\nOutputs:")
    print(sample_out)
    print(summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
