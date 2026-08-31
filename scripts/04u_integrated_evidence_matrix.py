#!/usr/bin/env python3
"""
ATLAS — Stage 04U
Integrated Evidence Matrix / Final Computational Prioritization

Purpose
-------
Combine the major computational evidence layers into one transparent matrix
for final interpretation and experimental follow-up planning.

Inputs
------
04R final candidate prioritization
04T ADMET / structural assessment
04S docking
04Q network integration
04P target annotations
04O safety screening
04N regulatory / clinical evidence

Outputs
-------
results/cmap/integrated_evidence/
    ATLAS_integrated_evidence_matrix.csv
    ATLAS_experimental_validation_shortlist.csv
    ATLAS_integrated_evidence_summary.csv
    ATLAS_integrated_evidence_metadata.json

Principles
----------
- No probability of efficacy is produced.
- Safety can veto a candidate.
- Docking is supportive, not decisive.
- Regulatory absence is not penalized.
- CMap opposition remains transcriptomic evidence, not proof of resistance reversal.
- Final labels are transparent evidence categories for experimental validation.
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_04R = (
    PROJECT_ROOT / "results" / "cmap" / "final_prioritization"
    / "ATLAS_final_candidate_prioritization.csv"
)

INPUT_04T = (
    PROJECT_ROOT / "results" / "cmap" / "admet_structural"
    / "ATLAS_ADMET_structural_assessment.csv"
)

INPUT_04S = (
    PROJECT_ROOT / "results" / "cmap" / "docking"
    / "ATLAS_docking_results.csv"
)

INPUT_04Q = (
    PROJECT_ROOT / "results" / "cmap" / "network_integration"
    / "ATLAS_drug_network_prioritized.csv"
)

INPUT_04P = (
    PROJECT_ROOT / "results" / "cmap" / "drug_targets"
    / "ATLAS_CMap_drug_target_annotations.csv"
)

INPUT_04O = (
    PROJECT_ROOT / "results" / "cmap" / "safety_screening"
    / "ATLAS_CMap_safety_screening.csv"
)

INPUT_04N = (
    PROJECT_ROOT / "results" / "cmap" / "regulatory_status"
    / "ATLAS_CMap_regulatory_annotations.csv"
)

OUT_DIR = (
    PROJECT_ROOT / "results" / "cmap" / "integrated_evidence"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_MATRIX = OUT_DIR / "ATLAS_integrated_evidence_matrix.csv"
OUT_SHORTLIST = OUT_DIR / "ATLAS_experimental_validation_shortlist.csv"
OUT_SUMMARY = OUT_DIR / "ATLAS_integrated_evidence_summary.csv"
OUT_META = OUT_DIR / "ATLAS_integrated_evidence_metadata.json"


def header(text: str) -> None:
    print("\n" + "=" * 78, flush=True)
    print(text, flush=True)
    print("=" * 78, flush=True)


def clean_text(x: Any) -> str:
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def normalize_keys(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df

    out = df.copy()
    for col in ["pert_id", "pert_iname"]:
        if col in out.columns:
            out[col] = out[col].astype("string").str.strip()
    return out


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    os.replace(tmp, path)


def merge_layer(
    base: pd.DataFrame,
    layer: pd.DataFrame,
    keys: list[str] = ["pert_id", "pert_iname"],
) -> pd.DataFrame:
    if layer.empty:
        return base

    if not all(k in base.columns and k in layer.columns for k in keys):
        return base

    layer = layer.drop_duplicates(keys, keep="first").copy()

    overlap = [
        c for c in layer.columns
        if c in base.columns and c not in keys
    ]

    if overlap:
        layer = layer.drop(columns=overlap)

    return base.merge(
        layer,
        on=keys,
        how="left",
        validate="one_to_one",
    )


def yesno(cond: bool) -> str:
    return "YES" if bool(cond) else "NO"


def evidence_profile(row: pd.Series) -> dict[str, Any]:
    """
    Build interpretable evidence badges and counts.
    """
    cmap = clean_text(row.get("cmap_label"))
    safety = clean_text(row.get("safety_label"))
    target = clean_text(row.get("target_label"))
    network = clean_text(row.get("network_label"))
    regulatory = clean_text(row.get("regulatory_label"))
    developability = clean_text(
        row.get("structural_developability_category")
    )
    docking_status = clean_text(row.get("docking_status"))
    protocol = clean_text(row.get("protocol_validation"))

    strong_layers = 0
    supportive_layers = 0
    caution_layers = 0
    missing_layers = 0

    # CMap
    if cmap == "STRONG":
        strong_layers += 1
    elif cmap == "MODERATE":
        supportive_layers += 1
    else:
        missing_layers += 1

    # Safety
    if safety.startswith("PASS"):
        strong_layers += 1
    elif safety == "CAUTION":
        caution_layers += 1
    elif safety in {"INSUFFICIENT", "UNKNOWN"}:
        missing_layers += 1
    elif safety == "HIGH_RISK":
        caution_layers += 2

    # Target
    if target == "STRONG":
        strong_layers += 1
    elif target == "MODERATE":
        supportive_layers += 1
    elif target == "WEAK":
        caution_layers += 1
    else:
        missing_layers += 1

    # Network
    if network == "STRONG":
        strong_layers += 1
    elif network == "MODERATE":
        supportive_layers += 1
    elif network == "WEAK":
        caution_layers += 1
    else:
        missing_layers += 1

    # Developability
    if developability == "FAVORABLE_DEVELOPABILITY":
        strong_layers += 1
    elif developability == "ACCEPTABLE_WITH_CAUTION":
        supportive_layers += 1
    elif developability == "STRUCTURAL_REVIEW_REQUIRED":
        caution_layers += 1
    elif developability == "UNFAVORABLE_DEVELOPABILITY":
        caution_layers += 2
    else:
        missing_layers += 1

    # Docking
    if docking_status == "COMPLETED":
        if protocol == "PASS_RMSD_LE_2A":
            strong_layers += 1
        elif protocol in {
            "COMPLETED_RMSD_UNAVAILABLE",
            "CAUTION_RMSD_GT_2A",
        }:
            supportive_layers += 1
            caution_layers += 1
        else:
            supportive_layers += 1
    elif docking_status == "SKIPPED_NO_VETTED_EXPERIMENTAL_STRUCTURE":
        missing_layers += 1

    # Regulatory is context only
    if regulatory in {"HIGH", "MODERATE", "CLINICAL_CONTEXT"}:
        supportive_layers += 1

    return {
        "strong_evidence_layers_n": strong_layers,
        "supportive_evidence_layers_n": supportive_layers,
        "caution_layers_n": caution_layers,
        "missing_evidence_layers_n": missing_layers,
    }


def final_category(row: pd.Series) -> tuple[str, str]:
    """
    Final computational interpretation for experimental validation.
    """
    safety = clean_text(
        row.get("safety_screening_recommendation")
    )
    hard_flag = row.get("hard_safety_flag", False)
    hard_flag = (
        bool(hard_flag)
        if not pd.isna(hard_flag)
        else False
    )

    t_decision = clean_text(row.get("04t_decision"))
    r_decision = clean_text(row.get("final_decision"))

    strong_n = int(
        pd.to_numeric(
            pd.Series([row.get("strong_evidence_layers_n")]),
            errors="coerce",
        ).fillna(0).iloc[0]
    )
    support_n = int(
        pd.to_numeric(
            pd.Series([row.get("supportive_evidence_layers_n")]),
            errors="coerce",
        ).fillna(0).iloc[0]
    )
    caution_n = int(
        pd.to_numeric(
            pd.Series([row.get("caution_layers_n")]),
            errors="coerce",
        ).fillna(0).iloc[0]
    )

    if safety == "HIGH_RISK_DEPRIORITIZE" or hard_flag:
        return (
            "DEPRIORITIZE",
            "High-risk safety evidence overrides other computational support.",
        )

    if t_decision == "DEPRIORITIZE":
        return (
            "DEPRIORITIZE",
            "04T developability/safety assessment does not support advancement.",
        )

    # Highest-confidence gate:
    # - candidate passed the strict 04R docking gate,
    # - 04T advanced it with validated docking,
    # - no meaningful safety caution,
    # - and several independent layers converge.
    #
    # Count both strong and supportive layers because the evidence system
    # intentionally distinguishes strength without treating MODERATE support
    # as absence. Requiring 4 STRONG layers alone was too restrictive and
    # incorrectly demoted candidates such as validated sitagliptin->DPP4.
    if (
        r_decision == "DOCK_NOW"
        and t_decision == "ADVANCE_TO_04U"
        and strong_n >= 3
        and (strong_n + support_n) >= 6
        and caution_n == 0
    ):
        return (
            "PRIORITY_EXPERIMENTAL_VALIDATION",
            "Validated docking and multiple convergent strong/supportive evidence layers justify top experimental priority.",
        )

    if (
        r_decision == "DOCK_NOW"
        and t_decision in {
            "ADVANCE_TO_04U_WITH_DOCKING_CAUTION",
            "ADVANCE_TO_04U_NO_EXPERIMENTAL_DOCKING",
        }
    ):
        return (
            "EXPERIMENTAL_VALIDATION_WITH_CAUTION",
            "Candidate remains well supported, but structural evidence is incomplete or cautionary.",
        )

    if (
        r_decision == "DOCK_WITH_CAUTION"
        or caution_n >= 2
    ):
        return (
            "MANUAL_REVIEW_BEFORE_EXPERIMENT",
            "Potentially interesting candidate with meaningful caution flags.",
        )

    if strong_n + support_n >= 4:
        return (
            "SECONDARY_EXPERIMENTAL_CANDIDATE",
            "Several evidence layers support follow-up, but evidence is not as convergent as top candidates.",
        )

    return (
        "HOLD_INSUFFICIENT_CONVERGENCE",
        "Current evidence does not sufficiently converge for immediate experimental prioritization.",
    )


def experimental_priority_score(row: pd.Series) -> float:
    """
    Transparent ranking score for ordering candidates within final categories.
    Not a probability.
    """
    score = 0.0

    score += float(row.get("strong_evidence_layers_n", 0)) * 2.0
    score += float(row.get("supportive_evidence_layers_n", 0)) * 1.0
    score -= float(row.get("caution_layers_n", 0)) * 1.5
    score -= float(row.get("missing_evidence_layers_n", 0)) * 0.5

    integrated = pd.to_numeric(
        pd.Series([row.get("integrated_prioritization_score")]),
        errors="coerce",
    ).iloc[0]

    if pd.notna(integrated):
        score += min(max(integrated, -10), 20) / 4.0

    return round(score, 3)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument(
        "--max-candidates",
        type=int,
        default=25,
        help="Maximum candidates to integrate from 04T. Default: 25",
    )

    p.add_argument(
        "--experimental-top",
        type=int,
        default=5,
        help="Maximum experimental shortlist size. Default: 5",
    )

    return p.parse_args()


def main() -> int:
    args = parse_args()

    header("ATLAS — Stage 04U Integrated Evidence Matrix")

    if not INPUT_04T.exists():
        print(
            f"ERROR: missing 04T input:\n{INPUT_04T}",
            flush=True,
        )
        return 1

    base = normalize_keys(
        read_csv(INPUT_04T)
    ).head(args.max_candidates)

    r = normalize_keys(read_csv(INPUT_04R))
    q = normalize_keys(read_csv(INPUT_04Q))
    p = normalize_keys(read_csv(INPUT_04P))
    o = normalize_keys(read_csv(INPUT_04O))
    n = normalize_keys(read_csv(INPUT_04N))
    s = read_csv(INPUT_04S)

    merged = base.copy()

    # Most 04R/04Q/04P/04O/04N columns are already inherited into 04T,
    # but merging here makes 04U robust if a column was not carried forward.
    for layer in [r, q, p, o, n]:
        merged = merge_layer(merged, layer)

    if not s.empty and "pert_iname" in s.columns:
        s = s.drop_duplicates("pert_iname", keep="first").copy()

        overlap = [
            c for c in s.columns
            if c in merged.columns and c != "pert_iname"
        ]

        if overlap:
            s = s.drop(columns=overlap)

        merged = merged.merge(
            s,
            on="pert_iname",
            how="left",
        )

    print(f"\nCandidates integrated: {len(merged)}", flush=True)

    # -----------------------------------------------------------------
    # Evidence matrix
    # -----------------------------------------------------------------
    profiles = merged.apply(evidence_profile, axis=1)
    profile_df = pd.DataFrame(
        list(profiles),
        index=merged.index,
    )

    for c in profile_df.columns:
        merged[c] = profile_df[c]

    final_parts = merged.apply(final_category, axis=1)
    merged["final_evidence_category"] = [
        x[0] for x in final_parts
    ]
    merged["final_evidence_reason"] = [
        x[1] for x in final_parts
    ]

    merged["experimental_priority_score"] = merged.apply(
        experimental_priority_score,
        axis=1,
    )

    # Human-readable badges
    merged["cmap_evidence_badge"] = merged.get(
        "cmap_label",
        pd.Series("", index=merged.index),
    )

    merged["safety_evidence_badge"] = merged.get(
        "safety_label",
        pd.Series("", index=merged.index),
    )

    merged["target_evidence_badge"] = merged.get(
        "target_label",
        pd.Series("", index=merged.index),
    )

    merged["network_evidence_badge"] = merged.get(
        "network_label",
        pd.Series("", index=merged.index),
    )

    merged["developability_evidence_badge"] = merged.get(
        "structural_developability_category",
        pd.Series("", index=merged.index),
    )

    merged["docking_evidence_badge"] = np.where(
        merged.get(
            "protocol_validation",
            pd.Series("", index=merged.index),
        )
        .fillna("")
        .astype(str)
        .eq("PASS_RMSD_LE_2A"),
        "VALIDATED_DOCKING",
        np.where(
            merged.get(
                "docking_status",
                pd.Series("", index=merged.index),
            )
            .fillna("")
            .astype(str)
            .eq("COMPLETED"),
            "DOCKING_WITH_CAUTION",
            "NO_VALIDATED_DOCKING",
        ),
    )

    # -----------------------------------------------------------------
    # Ranking
    # -----------------------------------------------------------------
    category_order = {
        "PRIORITY_EXPERIMENTAL_VALIDATION": 0,
        "EXPERIMENTAL_VALIDATION_WITH_CAUTION": 1,
        "SECONDARY_EXPERIMENTAL_CANDIDATE": 2,
        "MANUAL_REVIEW_BEFORE_EXPERIMENT": 3,
        "HOLD_INSUFFICIENT_CONVERGENCE": 4,
        "DEPRIORITIZE": 5,
    }

    merged["_category_order"] = (
        merged["final_evidence_category"]
        .map(category_order)
        .fillna(99)
    )

    merged = merged.sort_values(
        [
            "_category_order",
            "experimental_priority_score",
            "integrated_prioritization_score",
        ],
        ascending=[True, False, False],
        na_position="last",
    ).reset_index(drop=True)

    merged["04u_rank"] = np.arange(
        1,
        len(merged) + 1,
    )

    merged = merged.drop(
        columns=["_category_order"]
    )

    shortlist = merged[
        merged["final_evidence_category"].isin([
            "PRIORITY_EXPERIMENTAL_VALIDATION",
            "EXPERIMENTAL_VALIDATION_WITH_CAUTION",
            "SECONDARY_EXPERIMENTAL_CANDIDATE",
        ])
    ].head(args.experimental_top).copy()

    shortlist["experimental_shortlist_rank"] = np.arange(
        1,
        len(shortlist) + 1,
    )

    # -----------------------------------------------------------------
    # Outputs
    # -----------------------------------------------------------------
    atomic_csv(merged, OUT_MATRIX)
    atomic_csv(shortlist, OUT_SHORTLIST)

    summary = (
        merged.groupby(
            "final_evidence_category",
            dropna=False,
        )
        .agg(
            compound_count=("pert_iname", "size"),
            median_experimental_priority_score=(
                "experimental_priority_score",
                "median",
            ),
            tier1_count=(
                "priority_rank",
                lambda x: int(
                    pd.to_numeric(
                        x,
                        errors="coerce",
                    ).isin([1, 2, 3, 4, 5, 6, 7]).sum()
                ),
            ),
        )
        .reset_index()
    )

    atomic_csv(summary, OUT_SUMMARY)

    metadata = {
        "stage": "04U",
        "implementation": "v2_convergent_priority_gate",
        "created_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "inputs": {
            "04R": str(INPUT_04R),
            "04T": str(INPUT_04T),
            "04S": str(INPUT_04S),
            "04Q": str(INPUT_04Q),
            "04P": str(INPUT_04P),
            "04O": str(INPUT_04O),
            "04N": str(INPUT_04N),
        },
        "candidate_count": int(len(merged)),
        "experimental_shortlist_count": int(len(shortlist)),
        "guardrails": [
            "No efficacy probability is generated.",
            "Negative CMap tau is transcriptomic opposition, not proof of resistance reversal.",
            "Safety high-risk evidence can veto a candidate.",
            "Docking is supportive evidence, not proof of binding or efficacy.",
            "Regulatory absence is not penalized.",
            "Final shortlist is intended for experimental validation planning.",
        ],
        "recommended_wet_lab_design": {
            "models": "trastuzumab-resistant HER2-positive breast cancer models",
            "arms": [
                "vehicle/control",
                "trastuzumab alone",
                "candidate alone",
                "candidate + trastuzumab",
            ],
            "goal": (
                "test whether the candidate improves response to trastuzumab "
                "in resistant HER2-positive models"
            ),
        },
        "next_stage": "experimental_validation",
    }

    atomic_json(metadata, OUT_META)

    # -----------------------------------------------------------------
    # Console
    # -----------------------------------------------------------------
    header("STAGE 04U SUMMARY")
    print(
        summary.to_string(index=False),
        flush=True,
    )

    header("TOP INTEGRATED EVIDENCE")

    show = [
        c for c in [
            "04u_rank",
            "priority_rank",
            "pert_iname",
            "cmap_evidence_badge",
            "safety_evidence_badge",
            "target_evidence_badge",
            "network_evidence_badge",
            "developability_evidence_badge",
            "docking_evidence_badge",
            "strong_evidence_layers_n",
            "supportive_evidence_layers_n",
            "caution_layers_n",
            "experimental_priority_score",
            "final_evidence_category",
        ]
        if c in merged.columns
    ]

    print(
        merged[show]
        .head(15)
        .to_string(index=False),
        flush=True,
    )

    header("EXPERIMENTAL VALIDATION SHORTLIST")

    if shortlist.empty:
        print(
            "No compounds met the current experimental-validation gate.",
            flush=True,
        )
    else:
        short_cols = [
            c for c in [
                "experimental_shortlist_rank",
                "pert_iname",
                "best_network_target",
                "validated_uniprot_accession",
                "best_affinity_kcal_mol",
                "protocol_validation",
                "experimental_priority_score",
                "final_evidence_category",
                "final_evidence_reason",
            ]
            if c in shortlist.columns
        ]

        print(
            shortlist[short_cols]
            .to_string(index=False),
            flush=True,
        )

    header("STAGE 04U COMPLETE")
    print("\nOutputs:", flush=True)
    print(f"  {OUT_MATRIX}", flush=True)
    print(f"  {OUT_SHORTLIST}", flush=True)
    print(f"  {OUT_SUMMARY}", flush=True)
    print(f"  {OUT_META}", flush=True)
    print(
        "\nNext: experimental validation in trastuzumab-resistant HER2+ models",
        flush=True,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
