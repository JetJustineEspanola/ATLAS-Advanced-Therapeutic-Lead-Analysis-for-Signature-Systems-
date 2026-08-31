#!/usr/bin/env python3
"""
ATLAS — 00S Validated Resistance Evidence Export

Purpose
-------
Convert the external-validation results into a clean downstream evidence layer
for ATLAS network/drug prioritization.

Inputs
------
results/external_validation/three_dataset_strict_core_genes.csv
results/external_validation/pathway_validation/tgfb_final_evidence_matrix.csv
results/external_validation/pathway_validation/tgfb_reproducible_module_gene_table.csv
results/external_validation/pathway_validation/tgfb_dual_module_gene_effects.csv

Outputs
-------
results/external_validation/downstream/
  validated_resistance_gene_evidence.csv
  validated_tgfb_module_evidence.csv
  validated_evidence_summary.csv

Important
---------
This stage does NOT alter 04Q/04R/04U rankings yet. It creates a validated,
traceable evidence layer that can replace the older discovery-DEG fallback.
"""

from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

STRICT = ROOT / "results/external_validation/three_dataset_strict_core_genes.csv"
TGFB_FINAL = ROOT / "results/external_validation/pathway_validation/tgfb_final_evidence_matrix.csv"
TGFB_MODULE = ROOT / "results/external_validation/pathway_validation/tgfb_reproducible_module_gene_table.csv"
TGFB_DUAL = ROOT / "results/external_validation/pathway_validation/tgfb_dual_module_gene_effects.csv"

OUTDIR = ROOT / "results/external_validation/downstream"
OUTDIR.mkdir(parents=True, exist_ok=True)


def main():
    print("=" * 78)
    print("ATLAS — 00S VALIDATED RESISTANCE EVIDENCE EXPORT")
    print("=" * 78)

    strict = pd.read_csv(STRICT)

    # Standardized gene-level downstream evidence.
    out = pd.DataFrame({
        "gene_symbol": strict["gene_symbol"],
        "evidence_class": "THREE_DATASET_STRICT_CORE",
        "resistance_direction": strict["replication_direction"],
        "discovery_log2FC": strict["discovery_log2FC"],
        "discovery_padj": strict["discovery_padj"],
        "GSE121105_log2FC": strict["GSE121105_log2FC"],
        "GSE121105_padj": strict["GSE121105_padj"],
        "GSE237606_log2FC": strict["GSE237606_log2FC"],
        "GSE237606_padj": strict["GSE237606_padj"],
        "external_consensus_rank": strict["consensus_rank"],
    })

    out["validated_in_three_datasets"] = True
    out["minimum_abs_log2FC_across_three"] = out[
        ["discovery_log2FC", "GSE121105_log2FC", "GSE237606_log2FC"]
    ].abs().min(axis=1)

    # Conservative evidence strength label.
    out["validated_evidence_strength"] = np.where(
        out["minimum_abs_log2FC_across_three"] >= 2.0,
        "VERY_STRONG",
        "STRONG",
    )

    # Add TGF-beta module membership.
    tgfb_pos = pd.read_csv(TGFB_MODULE)
    positive_set = set(tgfb_pos["gene_symbol"].astype(str))

    dual = pd.read_csv(TGFB_DUAL)
    negative_set = set(
        dual.loc[
            dual["module"] == "NEGATIVE_TGFB_7",
            "gene_symbol"
        ].astype(str)
    )

    out["tgfb_positive_module_16"] = out["gene_symbol"].isin(positive_set)
    out["tgfb_negative_module_7"] = out["gene_symbol"].isin(negative_set)

    def tgfb_role(row):
        if row["tgfb_positive_module_16"]:
            return "TGF_BETA_POSITIVE_MODULE"
        if row["tgfb_negative_module_7"]:
            return "TGF_BETA_NEGATIVE_MODULE"
        return ""

    out["tgfb_module_role"] = out.apply(tgfb_role, axis=1)

    gene_out = OUTDIR / "validated_resistance_gene_evidence.csv"
    out.to_csv(gene_out, index=False)

    # TGF-beta module-specific evidence table.
    pos = pd.DataFrame({
        "gene_symbol": sorted(positive_set),
        "module": "POSITIVE_TGFB_16",
    })
    neg = pd.DataFrame({
        "gene_symbol": sorted(negative_set),
        "module": "NEGATIVE_TGFB_7",
    })
    module = pd.concat([pos, neg], ignore_index=True)

    module = module.merge(
        out[
            [
                "gene_symbol",
                "validated_in_three_datasets",
                "validated_evidence_strength",
                "resistance_direction",
            ]
        ],
        on="gene_symbol",
        how="left",
    )

    module["validated_in_three_datasets"] = (
        module["validated_in_three_datasets"].fillna(False)
    )

    module_out = OUTDIR / "validated_tgfb_module_evidence.csv"
    module.to_csv(module_out, index=False)

    tgfb_final = pd.read_csv(TGFB_FINAL)

    summary = pd.DataFrame([{
        "strict_core_gene_n": len(out),
        "strict_core_up_n": int(
            (out["resistance_direction"] == "UP_IN_RESISTANT").sum()
        ),
        "strict_core_down_n": int(
            (out["resistance_direction"] == "DOWN_IN_RESISTANT").sum()
        ),
        "very_strong_gene_n": int(
            (out["validated_evidence_strength"] == "VERY_STRONG").sum()
        ),
        "tgfb_positive_module_gene_n": len(positive_set),
        "tgfb_negative_module_gene_n": len(negative_set),
        "tgfb_positive_module_strict_core_overlap": int(
            out["tgfb_positive_module_16"].sum()
        ),
        "tgfb_negative_module_strict_core_overlap": int(
            out["tgfb_negative_module_7"].sum()
        ),
        "tgfb_hallmark_direction_consistent": False,
        "tgfb_interpretation": (
            "Reproducible TGF-beta pathway remodeling with a positive 16-gene "
            "module and a repressed 7-gene module; whole-pathway direction is "
            "context-dependent."
        ),
    }])

    summary_out = OUTDIR / "validated_evidence_summary.csv"
    summary.to_csv(summary_out, index=False)

    print(f"Strict-core genes exported: {len(out)}")
    print(
        f"  Up: {(out['resistance_direction'] == 'UP_IN_RESISTANT').sum()} | "
        f"Down: {(out['resistance_direction'] == 'DOWN_IN_RESISTANT').sum()}"
    )
    print(
        f"TGF-beta positive module overlap with strict core: "
        f"{out['tgfb_positive_module_16'].sum()}/{len(positive_set)}"
    )
    print(
        f"TGF-beta negative module overlap with strict core: "
        f"{out['tgfb_negative_module_7'].sum()}/{len(negative_set)}"
    )

    print("\nOutputs:")
    print(gene_out)
    print(module_out)
    print(summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
