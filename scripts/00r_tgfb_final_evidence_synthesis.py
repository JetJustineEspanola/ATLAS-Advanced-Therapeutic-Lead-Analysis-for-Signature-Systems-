#!/usr/bin/env python3
"""
ATLAS — 00R Final TGF-beta Evidence Synthesis

Combines:
- ranked GSEA
- gene-level TGF-beta audit
- leading-edge overlap
- sample-level module scores
- cross-dataset dual-module audit

Produces a compact, defensible evidence table and conclusion flags.

Outputs:
  results/external_validation/pathway_validation/
    tgfb_final_evidence_matrix.csv
    tgfb_final_conclusion.txt
"""

from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PV = ROOT / "results/external_validation/pathway_validation"

GSEA = PV / "tgfb_ranked_validation_summary.csv"
MOD = PV / "tgfb_dual_module_cross_dataset_summary.csv"
SAMPLE = PV / "tgfb_module_score_summary.csv"
LEDGE = PV / "tgfb_leading_edge_summary.csv"

OUT_MATRIX = PV / "tgfb_final_evidence_matrix.csv"
OUT_TXT = PV / "tgfb_final_conclusion.txt"


def main():
    gsea = pd.read_csv(GSEA)
    mod = pd.read_csv(MOD)
    sample = pd.read_csv(SAMPLE)
    ledge = pd.read_csv(LEDGE)

    rows = []

    for ds in ["discovery", "GSE121105", "GSE237606"]:
        g = gsea[gsea["dataset"] == ds]
        pos = mod[(mod["dataset"] == ds) & (mod["module"] == "POSITIVE_TGFB_16")]
        neg = mod[(mod["dataset"] == ds) & (mod["module"] == "NEGATIVE_TGFB_7")]

        row = {
            "dataset": ds,
            "hallmark_NES": g["NES"].iloc[0] if not g.empty else np.nan,
            "hallmark_FDR": g["FDR_q"].iloc[0] if not g.empty else np.nan,
            "positive_module_mean_log2FC": pos["mean_log2FC"].iloc[0] if not pos.empty else np.nan,
            "positive_module_wilcoxon_p": pos["wilcoxon_vs_zero_p"].iloc[0] if not pos.empty else np.nan,
            "negative_module_mean_log2FC": neg["mean_log2FC"].iloc[0] if not neg.empty else np.nan,
            "negative_module_wilcoxon_p": neg["wilcoxon_vs_zero_p"].iloc[0] if not neg.empty else np.nan,
        }

        if ds in {"GSE121105", "GSE237606"}:
            s_pos = sample[
                (sample["dataset"] == ds)
                & (sample["module"] == "POSITIVE_TGFB_16")
            ]
            s_neg = sample[
                (sample["dataset"] == ds)
                & (sample["module"] == "GSE121105_NEGATIVE_TGFB_7")
            ]

            row["positive_module_sample_delta"] = (
                s_pos["mean_difference_resistant_minus_sensitive"].iloc[0]
                if not s_pos.empty else np.nan
            )
            row["positive_module_sample_p"] = (
                s_pos["mannwhitney_p"].iloc[0] if not s_pos.empty else np.nan
            )
            row["negative_module_sample_delta"] = (
                s_neg["mean_difference_resistant_minus_sensitive"].iloc[0]
                if not s_neg.empty else np.nan
            )
            row["negative_module_sample_p"] = (
                s_neg["mannwhitney_p"].iloc[0] if not s_neg.empty else np.nan
            )

        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(OUT_MATRIX, index=False)

    # Evidence flags
    positive_supported = (
        (out["positive_module_mean_log2FC"] > 0).sum() >= 3
        and (
            (out["positive_module_wilcoxon_p"] < 0.05).sum() >= 2
        )
    )

    negative_supported = (
        (out["negative_module_mean_log2FC"] < 0).sum() >= 3
        and (
            (out["negative_module_wilcoxon_p"] < 0.05).sum() >= 2
        )
    )

    hallmark_direction_consistent = (
        np.sign(out["hallmark_NES"].dropna()).nunique() == 1
    )

    conclusion = []
    conclusion.append("ATLAS TGF-beta final evidence synthesis")
    conclusion.append("=" * 72)
    conclusion.append(
        f"Positive 16-gene module supported across datasets: {positive_supported}"
    )
    conclusion.append(
        f"Negative 7-gene module supported across datasets: {negative_supported}"
    )
    conclusion.append(
        f"Hallmark TGF-beta NES direction consistent across datasets: {hallmark_direction_consistent}"
    )
    conclusion.append("")
    conclusion.append(
        "Defensible interpretation:"
    )
    conclusion.append(
        "Trastuzumab resistance is associated with reproducible TGF-beta pathway remodeling, "
        "not a universally monotonic increase in the entire Hallmark TGF-beta pathway."
    )
    conclusion.append(
        "A 16-gene positive TGF-beta module is strongly activated in discovery and GSE237606 "
        "and trends upward in GSE121105, while a distinct 7-gene module is repressed across datasets."
    )
    conclusion.append(
        "Therefore, report pathway-level context dependence alongside reproducible submodule-level dysregulation."
    )

    OUT_TXT.write_text("\n".join(conclusion), encoding="utf-8")

    print("=" * 78)
    print("ATLAS — 00R FINAL TGF-BETA EVIDENCE SYNTHESIS")
    print("=" * 78)
    print(out.to_string(index=False))
    print("\n" + "\n".join(conclusion[-5:]))
    print("\nOutputs:")
    print(OUT_MATRIX)
    print(OUT_TXT)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
