#!/usr/bin/env python3
"""
ATLAS — Stage 04R
Final Candidate Prioritization / Docking Gatekeeper

Purpose
-------
Integrate the major evidence layers produced so far and produce a transparent
shortlist of compounds and drug-target pairs suitable for Stage 04S docking.

Evidence layers
---------------
04L / 04M  CMap priority and compound identity
04N        Regulatory / clinical evidence (positive context only; absence is not penalized)
04O        Safety / cytotoxicity / promiscuity / PAINS
04P        Drug-target annotation
04Q        Resistance-network support

Primary outputs
---------------
results/cmap/final_prioritization/
    ATLAS_final_candidate_prioritization.csv
    ATLAS_docking_shortlist.csv
    ATLAS_deprioritized_candidates.csv
    ATLAS_final_prioritization_summary.csv
    ATLAS_final_prioritization_metadata.json

Design principles
-----------------
- Safety can veto a candidate.
- Missing regulatory evidence does NOT penalize a candidate.
- Regulatory evidence is a small translational bonus only.
- Strong CMap evidence remains important, but cannot override severe safety risk.
- Broad/promiscuous target profiles are penalized.
- Docking is only recommended when a biologically supported target is available.
- The integrated score is a transparent ranking score, NOT a probability.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_04M = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "drug_filter"
    / "ATLAS_CMap_drug_candidates.csv"
)

INPUT_04N = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "regulatory_status"
    / "ATLAS_CMap_regulatory_annotations.csv"
)

INPUT_04O = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "safety_screening"
    / "ATLAS_CMap_safety_screening.csv"
)

INPUT_04P = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "drug_targets"
    / "ATLAS_CMap_drug_target_annotations.csv"
)

INPUT_04P_PAIRS = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "drug_targets"
    / "ATLAS_CMap_drug_target_pairs.csv"
)

INPUT_04Q = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "network_integration"
    / "ATLAS_drug_network_prioritized.csv"
)

OUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "final_prioritization"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUT_DIR / "ATLAS_final_candidate_prioritization.csv"
OUT_DOCK = OUT_DIR / "ATLAS_docking_shortlist.csv"
OUT_DROP = OUT_DIR / "ATLAS_deprioritized_candidates.csv"
OUT_SUMMARY = OUT_DIR / "ATLAS_final_prioritization_summary.csv"
OUT_META = OUT_DIR / "ATLAS_final_prioritization_metadata.json"


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

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


# ---------------------------------------------------------------------
# Evidence scoring
# ---------------------------------------------------------------------

def score_cmap(row: pd.Series) -> tuple[float, str]:
    """
    0–6 points.
    Tier is primary; consistency/mean tau lightly refines the score.
    """
    try:
        tier = int(float(row.get("priority_tier_number")))
    except Exception:
        tier = None

    if tier == 1:
        score = 6.0
        label = "STRONG"
    elif tier == 2:
        score = 4.0
        label = "MODERATE"
    elif tier == 3:
        score = 2.0
        label = "WEAK"
    elif tier == 4:
        score = 0.5
        label = "VERY_WEAK"
    else:
        score = 0.0
        label = "UNKNOWN"

    n_strong = pd.to_numeric(
        pd.Series([row.get("n_strong_negative")]),
        errors="coerce",
    ).iloc[0]

    if pd.notna(n_strong):
        if n_strong >= 3:
            score += 0.5
        elif n_strong >= 2:
            score += 0.25

    return min(score, 6.5), label


def score_safety(row: pd.Series) -> tuple[float, str, bool]:
    """
    Safety is primarily a gate, not a bonus race.

    Returns:
        score contribution
        label
        hard veto
    """
    status = clean_text(
        row.get("safety_screening_recommendation")
    )

    hard_flag_raw = row.get("hard_safety_flag", False)
    hard_flag = (
        bool(hard_flag_raw)
        if not pd.isna(hard_flag_raw)
        else False
    )

    if status == "HIGH_RISK_DEPRIORITIZE" or hard_flag:
        return -12.0, "HIGH_RISK", True

    if status == "CAUTION_MANUAL_REVIEW":
        return -2.5, "CAUTION", False

    if status == "PASS_PRELIMINARY_SCREEN":
        completeness = clean_text(
            row.get("safety_data_completeness")
        )
        if completeness == "GOOD":
            return 2.0, "PASS_GOOD_DATA", False
        return 1.0, "PASS_LIMITED_DATA", False

    if status == "INSUFFICIENT_SAFETY_DATA":
        return -1.0, "INSUFFICIENT", False

    return -1.5, "UNKNOWN", False


def score_target(row: pd.Series) -> tuple[float, str]:
    """
    0–5 points.
    """
    category = clean_text(
        row.get("target_support_category")
    )

    mapping = {
        "STRONG_TARGET_SUPPORT": (5.0, "STRONG"),
        "MODERATE_TARGET_SUPPORT": (3.0, "MODERATE"),
        "WEAK_TARGET_SUPPORT": (1.5, "WEAK"),
        "NO_TARGET_SUPPORT_FOUND": (0.0, "NONE"),
        "QUERY_ERROR": (0.0, "ERROR"),
    }

    return mapping.get(category, (0.0, "UNKNOWN"))


def score_network(row: pd.Series) -> tuple[float, str]:
    """
    0–6 points from 04Q network support.
    """
    category = clean_text(
        row.get("best_network_support_category")
    )

    raw_score = pd.to_numeric(
        pd.Series([row.get("best_network_support_score")]),
        errors="coerce",
    ).iloc[0]

    if pd.isna(raw_score):
        raw_score = 0.0

    if category == "STRONG_NETWORK_SUPPORT":
        return min(6.0, 4.0 + raw_score / 4.0), "STRONG"
    if category == "MODERATE_NETWORK_SUPPORT":
        return min(4.0, 2.0 + raw_score / 4.0), "MODERATE"
    if category == "WEAK_NETWORK_SUPPORT":
        return min(2.0, 1.0 + raw_score / 5.0), "WEAK"

    return 0.0, "NONE"


def score_regulatory(row: pd.Series) -> tuple[float, str]:
    """
    Positive-only translational context. No evidence = 0, not negative.
    Maximum +1.5.
    """
    category = clean_text(
        row.get("regulatory_evidence_category")
    )

    mapping = {
        "FDA_APPLICATION_AND_LABEL_EVIDENCE": (1.5, "HIGH"),
        "FDA_APPLICATION_RECORD_FOUND": (1.0, "MODERATE"),
        "FDA_LABEL_EVIDENCE_FOUND": (1.0, "MODERATE"),
        "CLINICAL_TRIAL_EVIDENCE_ONLY": (0.5, "CLINICAL_CONTEXT"),
        "NO_US_REGULATORY_OR_TRIAL_EVIDENCE_FOUND": (0.0, "NONE_FOUND"),
    }

    return mapping.get(category, (0.0, "UNKNOWN"))


def promiscuity_penalty(row: pd.Series) -> tuple[float, str]:
    """
    Penalize broad/non-specific pharmacology.

    Uses both PubChem assay promiscuity and ChEMBL target count.
    """
    penalty = 0.0
    reasons = []

    frac = pd.to_numeric(
        pd.Series([row.get("pubchem_active_fraction")]),
        errors="coerce",
    ).iloc[0]

    if pd.notna(frac):
        if frac >= 0.50:
            penalty -= 3.0
            reasons.append("very high assay promiscuity")
        elif frac >= 0.25:
            penalty -= 2.0
            reasons.append("high assay promiscuity")
        elif frac >= 0.15:
            penalty -= 1.0
            reasons.append("moderate assay promiscuity")

    target_n = pd.to_numeric(
        pd.Series([row.get("chembl_target_n")]),
        errors="coerce",
    ).iloc[0]

    if pd.notna(target_n):
        if target_n >= 30:
            penalty -= 3.0
            reasons.append("very broad ChEMBL target profile")
        elif target_n >= 15:
            penalty -= 2.0
            reasons.append("broad ChEMBL target profile")
        elif target_n >= 10:
            penalty -= 0.75
            reasons.append("multi-target ChEMBL profile")

    return penalty, " | ".join(reasons)


def docking_gate(row: pd.Series) -> tuple[str, str]:
    """
    Final rule-based gate for Stage 04S.
    """
    if bool(row.get("hard_safety_veto", False)):
        return (
            "DEPRIORITIZE_SAFETY",
            "Severe safety evidence overrides downstream computational support.",
        )

    safety = clean_text(row.get("safety_label"))
    target = clean_text(row.get("target_label"))
    network = clean_text(row.get("network_label"))
    cmap = clean_text(row.get("cmap_label"))
    best_target = clean_text(row.get("best_network_target"))

    if not best_target:
        return (
            "DEPRIORITIZE_WEAK_TARGET",
            "No specific network-supported protein target is available for docking.",
        )

    if target in {"NONE", "ERROR", "UNKNOWN"}:
        return (
            "DEPRIORITIZE_WEAK_TARGET",
            "Drug-target support is insufficient for a defensible docking pair.",
        )

    if network == "NONE":
        return (
            "HOLD_NETWORK_UNSUPPORTED",
            "Target evidence exists but the target lacks resistance-network support.",
        )

    if safety == "INSUFFICIENT":
        return (
            "MANUAL_REVIEW",
            "Biological evidence may be useful, but safety evidence is insufficient.",
        )

    if safety == "CAUTION":
        # Allow only strong multi-layer evidence through as a review candidate.
        if cmap == "STRONG" and target in {"STRONG", "MODERATE"} and network in {"STRONG", "MODERATE"}:
            return (
                "DOCK_WITH_CAUTION",
                "Strong multi-layer evidence justifies docking, but safety/promiscuity requires explicit caution.",
            )
        return (
            "MANUAL_REVIEW",
            "Safety caution is not offset by sufficiently strong multi-layer evidence.",
        )

    if (
        safety.startswith("PASS")
        and cmap in {"STRONG", "MODERATE"}
        and target in {"STRONG", "MODERATE"}
        and network in {"STRONG", "MODERATE"}
    ):
        return (
            "DOCK_NOW",
            "CMap, preliminary safety, target evidence, and resistance-network support converge.",
        )

    if (
        safety.startswith("PASS")
        and target in {"STRONG", "MODERATE", "WEAK"}
        and network in {"STRONG", "MODERATE", "WEAK"}
    ):
        return (
            "KEEP_FOR_REVIEW",
            "Candidate remains biologically plausible but does not meet the strongest docking gate.",
        )

    return (
        "HOLD_INSUFFICIENT_EVIDENCE",
        "Current evidence is not strong enough for docking.",
    )


# ---------------------------------------------------------------------
# Docking pair extraction
# ---------------------------------------------------------------------

def resolve_uniprot_gene(symbol: str) -> dict[str, Any]:
    """
    Resolve an exact human gene symbol to a reviewed UniProt accession.

    This is used to validate the 04Q best-network target before docking.
    We intentionally do NOT fall back to an unrelated first ChEMBL target row.
    """
    symbol = clean_text(symbol)

    empty = {
        "validated_target_symbol": symbol,
        "validated_uniprot_accession": "",
        "validated_uniprot_entry": "",
        "validated_protein_name": "",
        "target_mapping_status": "UNRESOLVED",
    }

    if not symbol:
        return empty

    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"(gene_exact:{symbol}) AND (organism_id:9606) AND (reviewed:true)",
        "format": "json",
        "fields": "accession,id,gene_names,protein_name",
        "size": 5,
    }

    try:
        r = requests.get(
            url,
            params=params,
            headers={"User-Agent": "ATLAS-04R/2.0"},
            timeout=(4, 12),
        )
        if r.status_code != 200:
            return empty
        payload = r.json()
    except Exception:
        return empty

    results = payload.get("results", []) or []
    if not results:
        return empty

    # Require the requested symbol to be present among returned gene names.
    for rec in results:
        genes = rec.get("genes", []) or []
        gene_names = set()

        for g in genes:
            primary = g.get("geneName", {}) or {}
            if primary.get("value"):
                gene_names.add(str(primary["value"]).upper())

            for syn in g.get("synonyms", []) or []:
                if syn.get("value"):
                    gene_names.add(str(syn["value"]).upper())

        if symbol.upper() not in gene_names:
            continue

        protein_desc = rec.get("proteinDescription", {}) or {}
        rec_name = protein_desc.get("recommendedName", {}) or {}
        full_name = rec_name.get("fullName", {}) or {}

        return {
            "validated_target_symbol": symbol,
            "validated_uniprot_accession": clean_text(
                rec.get("primaryAccession")
            ),
            "validated_uniprot_entry": clean_text(
                rec.get("uniProtkbId")
            ),
            "validated_protein_name": clean_text(
                full_name.get("value")
            ),
            "target_mapping_status": "VALIDATED_UNIPROT_HUMAN_REVIEWED",
        }

    return empty


def attach_best_target_details(
    ranked: pd.DataFrame,
    pairs: pd.DataFrame,
) -> pd.DataFrame:
    """
    Validate the best-network target for docking.

    Important correction:
    The previous implementation fell back to the first ChEMBL target row when
    the network target could not be matched. That could attach the wrong
    accession to the selected network target. This version never does that.

    It resolves the actual best-network gene symbol against reviewed human
    UniProt instead.
    """
    out = ranked.copy()

    cache: dict[str, dict[str, Any]] = {}
    records = []

    for _, row in out.iterrows():
        symbol = clean_text(row.get("best_network_target"))

        if symbol not in cache:
            cache[symbol] = resolve_uniprot_gene(symbol)

        records.append(cache[symbol])

    resolved = pd.DataFrame(records, index=out.index)

    for c in resolved.columns:
        out[c] = resolved[c]

    # Preserve traceable 04P ChEMBL context only when an exact preferred-name
    # match exists. No arbitrary fallback is allowed.
    exact_chembl_ids = []
    exact_chembl_accessions = []

    for _, row in out.iterrows():
        pert_id = clean_text(row.get("pert_id"))
        pert_iname = clean_text(row.get("pert_iname"))
        symbol = clean_text(row.get("best_network_target"))

        if pairs.empty or not symbol:
            exact_chembl_ids.append("")
            exact_chembl_accessions.append("")
            continue

        g = pairs[
            (pairs["pert_id"].astype(str) == pert_id)
            & (pairs["pert_iname"].astype(str) == pert_iname)
        ].copy()

        if g.empty:
            exact_chembl_ids.append("")
            exact_chembl_accessions.append("")
            continue

        mask = pd.Series(False, index=g.index)

        if "target_pref_name" in g.columns:
            mask |= (
                g["target_pref_name"]
                .fillna("")
                .astype(str)
                .str.upper()
                .eq(symbol.upper())
            )

        if "target_component_descriptions" in g.columns:
            mask |= (
                g["target_component_descriptions"]
                .fillna("")
                .astype(str)
                .str.upper()
                .str.contains(symbol.upper(), regex=False)
            )

        match = g[mask]

        if match.empty:
            exact_chembl_ids.append("")
            exact_chembl_accessions.append("")
        else:
            chosen = match.iloc[0]
            exact_chembl_ids.append(
                clean_text(chosen.get("target_chembl_id"))
            )
            exact_chembl_accessions.append(
                clean_text(chosen.get("target_accessions"))
            )

    out["exact_target_chembl_id"] = exact_chembl_ids
    out["exact_04p_target_accessions"] = exact_chembl_accessions

    return out


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument(
        "--top-docking",
        type=int,
        default=5,
        help="Maximum number of docking pairs to shortlist. Default: 5",
    )

    p.add_argument(
        "--allow-caution",
        action="store_true",
        help=(
            "Allow DOCK_WITH_CAUTION candidates into the final docking shortlist. "
            "Without this flag, only DOCK_NOW candidates are shortlisted."
        ),
    )

    return p.parse_args()


def main() -> int:
    args = parse_args()

    header("ATLAS — Stage 04R Final Candidate Prioritization")

    required = {
        "04M": INPUT_04M,
        "04O": INPUT_04O,
        "04P": INPUT_04P,
        "04Q": INPUT_04Q,
    }

    missing = [
        f"{stage}: {path}"
        for stage, path in required.items()
        if not path.exists()
    ]

    if missing:
        print("\nERROR: required upstream outputs are missing:", flush=True)
        for item in missing:
            print(f"  {item}", flush=True)
        return 1

    base = normalize_keys(read_csv(INPUT_04M))
    reg = normalize_keys(read_csv(INPUT_04N))
    safety = normalize_keys(read_csv(INPUT_04O))
    targets = normalize_keys(read_csv(INPUT_04P))
    network = normalize_keys(read_csv(INPUT_04Q))
    pairs = normalize_keys(read_csv(INPUT_04P_PAIRS))

    # Restrict to candidates that actually reached 04Q, since docking requires
    # a network-supported target context.
    if not network.empty:
        network_keys = set(
            zip(
                network["pert_id"].astype(str),
                network["pert_iname"].astype(str),
            )
        )
        base = base[
            [
                (str(r["pert_id"]), str(r["pert_iname"])) in network_keys
                for _, r in base.iterrows()
            ]
        ].copy()

    merged = base.copy()
    merged = merge_layer(merged, reg)
    merged = merge_layer(merged, safety)
    merged = merge_layer(merged, targets)
    merged = merge_layer(merged, network)

    if merged.empty:
        print("ERROR: no integrated candidates remained after merging.", flush=True)
        return 1

    # ---------------------------------------------------------------
    # Layer scores
    # ---------------------------------------------------------------
    cmap_parts = merged.apply(score_cmap, axis=1)
    merged["cmap_score"] = [x[0] for x in cmap_parts]
    merged["cmap_label"] = [x[1] for x in cmap_parts]

    safety_parts = merged.apply(score_safety, axis=1)
    merged["safety_score"] = [x[0] for x in safety_parts]
    merged["safety_label"] = [x[1] for x in safety_parts]
    merged["hard_safety_veto"] = [x[2] for x in safety_parts]

    target_parts = merged.apply(score_target, axis=1)
    merged["target_score"] = [x[0] for x in target_parts]
    merged["target_label"] = [x[1] for x in target_parts]

    network_parts = merged.apply(score_network, axis=1)
    merged["network_score"] = [x[0] for x in network_parts]
    merged["network_label"] = [x[1] for x in network_parts]

    regulatory_parts = merged.apply(score_regulatory, axis=1)
    merged["regulatory_score"] = [x[0] for x in regulatory_parts]
    merged["regulatory_label"] = [x[1] for x in regulatory_parts]

    promis_parts = merged.apply(promiscuity_penalty, axis=1)
    merged["promiscuity_penalty"] = [x[0] for x in promis_parts]
    merged["promiscuity_penalty_reason"] = [x[1] for x in promis_parts]

    # ---------------------------------------------------------------
    # Integrated score
    # ---------------------------------------------------------------
    merged["integrated_prioritization_score"] = (
        merged["cmap_score"]
        + merged["safety_score"]
        + merged["target_score"]
        + merged["network_score"]
        + merged["regulatory_score"]
        + merged["promiscuity_penalty"]
    )

    gate_parts = merged.apply(docking_gate, axis=1)
    merged["final_decision"] = [x[0] for x in gate_parts]
    merged["final_decision_reason"] = [x[1] for x in gate_parts]

    # Decision ordering before score.
    decision_order = {
        "DOCK_NOW": 0,
        "DOCK_WITH_CAUTION": 1,
        "KEEP_FOR_REVIEW": 2,
        "MANUAL_REVIEW": 3,
        "HOLD_NETWORK_UNSUPPORTED": 4,
        "HOLD_INSUFFICIENT_EVIDENCE": 5,
        "DEPRIORITIZE_WEAK_TARGET": 6,
        "DEPRIORITIZE_SAFETY": 7,
    }

    merged["_decision_order"] = (
        merged["final_decision"]
        .map(decision_order)
        .fillna(99)
    )

    merged["priority_rank"] = pd.to_numeric(
        merged["priority_rank"],
        errors="coerce",
    )

    merged = merged.sort_values(
        [
            "_decision_order",
            "integrated_prioritization_score",
            "priority_rank",
        ],
        ascending=[True, False, True],
        na_position="last",
    ).reset_index(drop=True)

    merged["final_rank"] = np.arange(1, len(merged) + 1)

    merged = attach_best_target_details(merged, pairs)

    # ---------------------------------------------------------------
    # Docking shortlist
    # ---------------------------------------------------------------
    allowed = ["DOCK_NOW"]
    if args.allow_caution:
        allowed.append("DOCK_WITH_CAUTION")

    docking = merged[
        merged["final_decision"].isin(allowed)
    ].copy()

    # Require a concrete target name.
    docking = docking[
        docking["best_network_target"]
        .fillna("")
        .astype(str)
        .str.strip()
        .ne("")
    ].copy()

    # Docking requires a validated human protein accession.
    docking = docking[
        docking["target_mapping_status"]
        .eq("VALIDATED_UNIPROT_HUMAN_REVIEWED")
        & docking["validated_uniprot_accession"]
        .fillna("")
        .astype(str)
        .str.strip()
        .ne("")
    ].copy()

    docking = docking.head(max(1, args.top_docking)).copy()
    docking["docking_rank"] = np.arange(1, len(docking) + 1)

    deprioritized = merged[
        merged["final_decision"].str.startswith("DEPRIORITIZE")
    ].copy()

    # ---------------------------------------------------------------
    # Outputs
    # ---------------------------------------------------------------
    merged = merged.drop(columns=["_decision_order"])

    atomic_csv(merged, OUT_ALL)
    atomic_csv(docking, OUT_DOCK)
    atomic_csv(deprioritized, OUT_DROP)

    summary = (
        merged.groupby("final_decision", dropna=False)
        .agg(
            compound_count=("pert_iname", "size"),
            median_integrated_score=(
                "integrated_prioritization_score",
                "median",
            ),
            tier1_count=(
                "priority_tier_number",
                lambda x: int(
                    (pd.to_numeric(x, errors="coerce") == 1).sum()
                ),
            ),
        )
        .reset_index()
        .sort_values(
            "median_integrated_score",
            ascending=False,
        )
    )

    atomic_csv(summary, OUT_SUMMARY)

    metadata = {
        "stage": "04R",
        "implementation": "v2_validated_docking_target_mapping",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "04M": str(INPUT_04M),
            "04N": str(INPUT_04N),
            "04O": str(INPUT_04O),
            "04P": str(INPUT_04P),
            "04P_pairs": str(INPUT_04P_PAIRS),
            "04Q": str(INPUT_04Q),
        },
        "candidate_count_integrated": int(len(merged)),
        "docking_shortlist_count": int(len(docking)),
        "top_docking_limit": int(args.top_docking),
        "allow_caution_in_shortlist": bool(args.allow_caution),
        "scoring": {
            "CMap": "0–6.5",
            "Safety": "+2 to -12, with hard veto",
            "Target": "0–5",
            "Network": "0–6",
            "Regulatory": "0 to +1.5 only; no negative penalty for absence",
            "Promiscuity": "0 to approximately -6",
        },
        "guardrails": [
            "Integrated score is a transparent prioritization score, not a probability.",
            "High-risk safety evidence can veto a candidate.",
            "No regulatory database hit is not treated as negative evidence.",
            "Docking is only recommended for target-supported, network-supported pairs.",
            "Docking score will not be interpreted as proof of binding or efficacy.",
        ],
        "next_stage": "04S_target_supported_molecular_docking",
    }

    atomic_json(metadata, OUT_META)

    # ---------------------------------------------------------------
    # Console
    # ---------------------------------------------------------------
    header("STAGE 04R SUMMARY")
    print(summary.to_string(index=False), flush=True)

    header("TOP INTEGRATED CANDIDATES")

    display_cols = [
        c for c in [
            "final_rank",
            "priority_rank",
            "pert_iname",
            "cmap_label",
            "safety_label",
            "target_label",
            "network_label",
            "regulatory_label",
            "promiscuity_penalty",
            "integrated_prioritization_score",
            "best_network_target",
            "final_decision",
        ]
        if c in merged.columns
    ]

    print(
        merged[display_cols]
        .head(15)
        .to_string(index=False),
        flush=True,
    )

    header("04S DOCKING SHORTLIST")

    if docking.empty:
        print(
            "\nNo candidates passed the strict DOCK_NOW gate.",
            flush=True,
        )
        print(
            "Review KEEP_FOR_REVIEW / DOCK_WITH_CAUTION candidates before relaxing the gate.",
            flush=True,
        )
    else:
        dock_cols = [
            c for c in [
                "docking_rank",
                "pert_iname",
                "best_network_target",
                "validated_target_symbol",
                "validated_uniprot_accession",
                "validated_uniprot_entry",
                "exact_target_chembl_id",
                "cmap_label",
                "safety_label",
                "target_label",
                "network_label",
                "integrated_prioritization_score",
                "final_decision",
            ]
            if c in docking.columns
        ]
        print(
            docking[dock_cols].to_string(index=False),
            flush=True,
        )

    header("STAGE 04R COMPLETE")
    print("\nOutputs:", flush=True)
    print(f"  {OUT_ALL}", flush=True)
    print(f"  {OUT_DOCK}", flush=True)
    print(f"  {OUT_DROP}", flush=True)
    print(f"  {OUT_SUMMARY}", flush=True)
    print(f"  {OUT_META}", flush=True)
    print(
        "\nNext: 04S — Target-supported molecular docking",
        flush=True,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
