#!/usr/bin/env python3
"""
ATLAS — Stage 04T
ADMET / Structural Developability Assessment

Purpose
-------
Evaluate the Stage 04R/04S candidates using transparent, local, reproducible
molecular descriptors and structural alerts.

This stage is deliberately framed as a *developability screen*, not as a
clinical ADMET predictor.

Inputs
------
04R final prioritization
04S docking results
04O safety screening
04N regulatory annotations

Outputs
-------
results/cmap/admet_structural/
    ATLAS_ADMET_structural_assessment.csv
    ATLAS_ADMET_structural_prioritized.csv
    ATLAS_ADMET_structural_summary.csv
    ATLAS_ADMET_structural_metadata.json

Key outputs
-----------
- Lipinski / Veber-style property checks
- QED
- ESOL-like aqueous solubility estimate
- structural alert counts (PAINS and Brenk)
- simple permeability / oral-likeness proxies
- docking-aware structural interpretation
- explicit decision categories for 04U integration

Guardrails
----------
- These are computational proxies, not measured PK/PD or toxicology.
- ESOL is a rough solubility estimate, not an experimental value.
- Lipinski/Veber rules are heuristics, not hard efficacy rules.
- Docking score is not proof of binding or efficacy.
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

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_04R = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "final_prioritization"
    / "ATLAS_final_candidate_prioritization.csv"
)

INPUT_04S = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "docking"
    / "ATLAS_docking_results.csv"
)

INPUT_04O = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "safety_screening"
    / "ATLAS_CMap_safety_screening.csv"
)

INPUT_04N = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "regulatory_status"
    / "ATLAS_CMap_regulatory_annotations.csv"
)

OUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "admet_structural"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUT_DIR / "ATLAS_ADMET_structural_assessment.csv"
OUT_PRIORITIZED = OUT_DIR / "ATLAS_ADMET_structural_prioritized.csv"
OUT_SUMMARY = OUT_DIR / "ATLAS_ADMET_structural_summary.csv"
OUT_META = OUT_DIR / "ATLAS_ADMET_structural_metadata.json"


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


def build_catalog(which: str) -> FilterCatalog:
    params = FilterCatalogParams()

    if which == "PAINS":
        params.AddCatalog(
            FilterCatalogParams.FilterCatalogs.PAINS
        )
    elif which == "BRENK":
        params.AddCatalog(
            FilterCatalogParams.FilterCatalogs.BRENK
        )
    else:
        raise ValueError(which)

    return FilterCatalog(params)


PAINS_CATALOG = build_catalog("PAINS")
BRENK_CATALOG = build_catalog("BRENK")


def smiles_from_row(row: pd.Series) -> str:
    candidates = [
        "canonical_smiles",
        "isomeric_smiles",
        "pubchem_canonical_smiles",
        "pubchem_isomeric_smiles",
        "smiles",
    ]

    for col in candidates:
        if col in row.index:
            value = clean_text(row.get(col))
            if value:
                return value

    return ""


def esol_logS(mol: Chem.Mol) -> float:
    """
    Delaney-style ESOL estimate.

    Delaney-style ESOL approximation:
        logS = 0.16 - 1.5 logP - 0.0062(MW)
               + 0.066 RB + 0.066 AP

    where AP = aromatic_atoms / heavy_atoms.

    Important:
    The previous version accidentally used -0.01 * MW, which made
    solubility estimates systematically too negative.

    where AP = aromatic_atoms / heavy_atoms.
    """
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    rb = Lipinski.NumRotatableBonds(mol)

    heavy = max(1, mol.GetNumHeavyAtoms())
    aromatic_atoms = sum(
        1 for atom in mol.GetAtoms()
        if atom.GetIsAromatic()
    )
    ap = aromatic_atoms / heavy

    return (
        0.16
        - 1.5 * logp
        - 0.0062 * mw
        + 0.066 * rb
        + 0.066 * ap
    )


def filter_matches(
    mol: Chem.Mol,
    catalog: FilterCatalog,
) -> tuple[int, str]:
    matches = catalog.GetMatches(mol)
    names = []

    for match in matches:
        try:
            names.append(
                clean_text(match.GetDescription())
            )
        except Exception:
            names.append("STRUCTURAL_ALERT")

    return len(matches), " | ".join(sorted(set(names)))


def compute_descriptors(
    smiles: str,
) -> dict[str, Any]:
    if not smiles:
        return {
            "structure_parse_status": "NO_SMILES",
        }

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return {
            "structure_parse_status": "INVALID_SMILES",
        }

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    fraction_csp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    formal_charge = Chem.GetFormalCharge(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    qed = QED.qed(mol)
    logS = esol_logS(mol)

    pains_n, pains_matches = filter_matches(
        mol,
        PAINS_CATALOG,
    )
    brenk_n, brenk_matches = filter_matches(
        mol,
        BRENK_CATALOG,
    )

    lipinski_violations = int(mw > 500)
    lipinski_violations += int(logp > 5)
    lipinski_violations += int(hbd > 5)
    lipinski_violations += int(hba > 10)

    veber_violations = int(rotb > 10)
    veber_violations += int(tpsa > 140)

    oral_property_pass = (
        lipinski_violations <= 1
        and veber_violations == 0
    )

    # Transparent heuristic proxies only.
    if tpsa <= 90 and logp <= 5 and mw <= 500:
        permeability_proxy = "FAVORABLE"
    elif tpsa <= 140 and mw <= 600:
        permeability_proxy = "INTERMEDIATE"
    else:
        permeability_proxy = "UNFAVORABLE"

    if logS >= -3:
        solubility_proxy = "FAVORABLE"
    elif logS >= -5:
        solubility_proxy = "INTERMEDIATE"
    else:
        solubility_proxy = "LOW"

    return {
        "structure_parse_status": "OK",
        "molecular_weight": mw,
        "clogp": logp,
        "tpsa": tpsa,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotb,
        "ring_count": rings,
        "aromatic_ring_count": aromatic_rings,
        "fraction_csp3": fraction_csp3,
        "formal_charge": formal_charge,
        "heavy_atom_count": heavy_atoms,
        "qed": qed,
        "esol_logS": logS,
        "lipinski_violations": lipinski_violations,
        "veber_violations": veber_violations,
        "oral_property_pass": oral_property_pass,
        "permeability_proxy": permeability_proxy,
        "solubility_proxy": solubility_proxy,
        "pains_alert_n_04t": pains_n,
        "pains_alerts_04t": pains_matches,
        "brenk_alert_n": brenk_n,
        "brenk_alerts": brenk_matches,
    }


def structural_score(row: pd.Series) -> tuple[float, str, str]:
    """
    Transparent developability score, not probability.

    Approximate range:
      +7 favorable
      negative = increasingly problematic
    """
    if clean_text(row.get("structure_parse_status")) != "OK":
        return -5.0, "INSUFFICIENT_STRUCTURE_DATA", "No valid molecular structure."

    score = 0.0
    reasons = []

    lip = pd.to_numeric(
        pd.Series([row.get("lipinski_violations")]),
        errors="coerce",
    ).iloc[0]
    veb = pd.to_numeric(
        pd.Series([row.get("veber_violations")]),
        errors="coerce",
    ).iloc[0]
    qed = pd.to_numeric(
        pd.Series([row.get("qed")]),
        errors="coerce",
    ).iloc[0]
    brenk = pd.to_numeric(
        pd.Series([row.get("brenk_alert_n")]),
        errors="coerce",
    ).iloc[0]
    pains = pd.to_numeric(
        pd.Series([row.get("pains_alert_n_04t")]),
        errors="coerce",
    ).iloc[0]

    perm = clean_text(row.get("permeability_proxy"))
    sol = clean_text(row.get("solubility_proxy"))

    if pd.notna(lip):
        if lip == 0:
            score += 2.0
        elif lip == 1:
            score += 1.0
            reasons.append("one Lipinski violation")
        else:
            score -= 2.0
            reasons.append("multiple Lipinski violations")

    if pd.notna(veb):
        if veb == 0:
            score += 1.0
        else:
            score -= 1.0
            reasons.append("Veber-style property concern")

    if pd.notna(qed):
        if qed >= 0.70:
            score += 2.0
        elif qed >= 0.45:
            score += 1.0
        else:
            score -= 1.0
            reasons.append("low QED")

    if perm == "FAVORABLE":
        score += 1.0
    elif perm == "UNFAVORABLE":
        score -= 1.0
        reasons.append("unfavorable permeability proxy")

    if sol == "FAVORABLE":
        score += 1.0
    elif sol == "LOW":
        score -= 1.5
        reasons.append("low ESOL solubility proxy")

    if pd.notna(pains) and pains > 0:
        score -= 1.0
        reasons.append("PAINS structural alert")

    if pd.notna(brenk):
        if brenk >= 3:
            score -= 2.0
            reasons.append("multiple Brenk structural alerts")
        elif brenk >= 1:
            score -= 0.5
            reasons.append("Brenk structural alert")

    if score >= 5:
        category = "FAVORABLE_DEVELOPABILITY"
    elif score >= 2:
        category = "ACCEPTABLE_WITH_CAUTION"
    elif score >= 0:
        category = "STRUCTURAL_REVIEW_REQUIRED"
    else:
        category = "UNFAVORABLE_DEVELOPABILITY"

    return score, category, " | ".join(reasons)


def final_04t_decision(row: pd.Series) -> tuple[str, str]:
    safety = clean_text(
        row.get("safety_screening_recommendation")
    )
    dev = clean_text(
        row.get("structural_developability_category")
    )
    docking_status = clean_text(
        row.get("docking_status")
    )
    protocol = clean_text(
        row.get("protocol_validation")
    )
    final_04r = clean_text(
        row.get("final_decision")
    )

    if safety == "HIGH_RISK_DEPRIORITIZE":
        return (
            "DEPRIORITIZE",
            "High-risk safety evidence overrides structural/docking support.",
        )

    if dev == "UNFAVORABLE_DEVELOPABILITY":
        return (
            "DEPRIORITIZE",
            "Structural developability profile is unfavorable.",
        )

    if docking_status == "COMPLETED":
        if protocol == "PASS_RMSD_LE_2A":
            if dev == "FAVORABLE_DEVELOPABILITY":
                return (
                    "ADVANCE_TO_04U_HIGH_CONFIDENCE",
                    "Validated docking plus favorable structural developability.",
                )
            return (
                "ADVANCE_TO_04U",
                "Validated docking with acceptable structural developability.",
            )

        if protocol in {
            "COMPLETED_RMSD_UNAVAILABLE",
            "CAUTION_RMSD_GT_2A",
        }:
            return (
                "ADVANCE_TO_04U_WITH_DOCKING_CAUTION",
                "Docking completed, but protocol validation is incomplete or cautionary.",
            )

    if docking_status == "SKIPPED_NO_VETTED_EXPERIMENTAL_STRUCTURE":
        if final_04r == "DOCK_NOW":
            return (
                "ADVANCE_TO_04U_NO_EXPERIMENTAL_DOCKING",
                "Candidate remains prioritized but lacks vetted experimental-structure docking.",
            )

    if dev in {
        "FAVORABLE_DEVELOPABILITY",
        "ACCEPTABLE_WITH_CAUTION",
    }:
        return (
            "KEEP_FOR_04U_REVIEW",
            "Developability is acceptable, but structural docking evidence is incomplete.",
        )

    return (
        "HOLD",
        "Evidence remains incomplete after structural assessment.",
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument(
        "--max-candidates",
        type=int,
        default=25,
        help=(
            "Assess the top N candidates from 04R. "
            "Default: 25"
        ),
    )

    return p.parse_args()


def main() -> int:
    args = parse_args()

    header("ATLAS — Stage 04T ADMET / Structural Developability Assessment")

    if not INPUT_04R.exists():
        print(
            f"ERROR: missing 04R input:\n{INPUT_04R}",
            flush=True,
        )
        return 1

    base = normalize_keys(
        read_csv(INPUT_04R)
    ).head(args.max_candidates)

    safety = normalize_keys(read_csv(INPUT_04O))
    reg = normalize_keys(read_csv(INPUT_04N))
    docking = read_csv(INPUT_04S)

    merged = base.copy()
    merged = merge_layer(merged, safety)
    merged = merge_layer(merged, reg)

    # 04S may not have pert_id, so merge by compound name.
    if not docking.empty and "pert_iname" in docking.columns:
        docking = docking.drop_duplicates(
            "pert_iname",
            keep="first",
        ).copy()

        overlap = [
            c for c in docking.columns
            if c in merged.columns and c != "pert_iname"
        ]

        if overlap:
            docking = docking.drop(columns=overlap)

        merged = merged.merge(
            docking,
            on="pert_iname",
            how="left",
        )

    print(f"\nCandidates selected: {len(merged)}", flush=True)

    # ---------------------------------------------------------------
    # Descriptor calculation
    # ---------------------------------------------------------------
    descriptor_records = []

    for _, row in merged.iterrows():
        smiles = smiles_from_row(row)
        descriptor_records.append(
            {
                "04t_smiles": smiles,
                **compute_descriptors(smiles),
            }
        )

    desc = pd.DataFrame(
        descriptor_records,
        index=merged.index,
    )

    for col in desc.columns:
        merged[col] = desc[col]

    # ---------------------------------------------------------------
    # Developability scoring
    # ---------------------------------------------------------------
    score_parts = merged.apply(
        structural_score,
        axis=1,
    )

    merged["structural_developability_score"] = [
        x[0] for x in score_parts
    ]
    merged["structural_developability_category"] = [
        x[1] for x in score_parts
    ]
    merged["structural_developability_reasons"] = [
        x[2] for x in score_parts
    ]

    decision_parts = merged.apply(
        final_04t_decision,
        axis=1,
    )

    merged["04t_decision"] = [
        x[0] for x in decision_parts
    ]
    merged["04t_decision_reason"] = [
        x[1] for x in decision_parts
    ]

    # ---------------------------------------------------------------
    # Docking-aware notes
    # ---------------------------------------------------------------
    merged["cyp3a4_interaction_context_flag"] = (
        merged.get(
            "target_symbol",
            pd.Series("", index=merged.index),
        )
        .fillna("")
        .astype(str)
        .str.upper()
        .eq("CYP3A4")
    )

    merged["cyp3a4_interaction_context_note"] = np.where(
        merged["cyp3a4_interaction_context_flag"],
        (
            "Compound was docked to CYP3A4; treat as structural interaction "
            "context only, not as proof of inhibition, induction, or metabolic liability."
        ),
        "",
    )

    # ---------------------------------------------------------------
    # Ranking
    # ---------------------------------------------------------------
    decision_order = {
        "ADVANCE_TO_04U_HIGH_CONFIDENCE": 0,
        "ADVANCE_TO_04U": 1,
        "ADVANCE_TO_04U_WITH_DOCKING_CAUTION": 2,
        "ADVANCE_TO_04U_NO_EXPERIMENTAL_DOCKING": 3,
        "KEEP_FOR_04U_REVIEW": 4,
        "HOLD": 5,
        "DEPRIORITIZE": 6,
    }

    merged["_decision_order"] = (
        merged["04t_decision"]
        .map(decision_order)
        .fillna(99)
    )

    merged = merged.sort_values(
        [
            "_decision_order",
            "structural_developability_score",
            "integrated_prioritization_score",
        ],
        ascending=[True, False, False],
        na_position="last",
    ).reset_index(drop=True)

    merged["04t_rank"] = np.arange(
        1,
        len(merged) + 1,
    )

    merged = merged.drop(
        columns=["_decision_order"]
    )

    prioritized = merged[
        merged["04t_decision"].isin([
            "ADVANCE_TO_04U_HIGH_CONFIDENCE",
            "ADVANCE_TO_04U",
            "ADVANCE_TO_04U_WITH_DOCKING_CAUTION",
            "ADVANCE_TO_04U_NO_EXPERIMENTAL_DOCKING",
            "KEEP_FOR_04U_REVIEW",
        ])
    ].copy()

    # ---------------------------------------------------------------
    # Outputs
    # ---------------------------------------------------------------
    atomic_csv(merged, OUT_ALL)
    atomic_csv(prioritized, OUT_PRIORITIZED)

    summary = (
        merged.groupby(
            "04t_decision",
            dropna=False,
        )
        .agg(
            compound_count=("pert_iname", "size"),
            median_developability_score=(
                "structural_developability_score",
                "median",
            ),
            median_qed=("qed", "median"),
            median_esol_logS=("esol_logS", "median"),
        )
        .reset_index()
    )

    atomic_csv(summary, OUT_SUMMARY)

    metadata = {
        "stage": "04T",
        "implementation": "v2_corrected_esol_coefficient",
        "created_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "inputs": {
            "04R": str(INPUT_04R),
            "04S": str(INPUT_04S),
            "04O": str(INPUT_04O),
            "04N": str(INPUT_04N),
        },
        "candidate_count": int(len(merged)),
        "methods": {
            "molecular_descriptors": "RDKit",
            "drug_likeness": "Lipinski/Veber-style heuristic checks",
            "qed": "RDKit QED",
            "solubility": "ESOL-like logS estimate",
            "structural_alerts": "RDKit PAINS + Brenk catalogs",
        },
        "guardrails": [
            "This stage is a computational developability screen, not measured ADMET.",
            "ESOL logS is an estimate, not experimental solubility.",
            "Lipinski and Veber rules are heuristics, not hard efficacy criteria.",
            "PAINS/Brenk alerts are screening flags, not proof of toxicity.",
            "CYP3A4 docking is interaction context, not proof of metabolic inhibition or induction.",
            "Docking score is not proof of binding or efficacy.",
        ],
        "next_stage": "04U_integrated_evidence_matrix",
    }

    atomic_json(metadata, OUT_META)

    # ---------------------------------------------------------------
    # Console
    # ---------------------------------------------------------------
    header("STAGE 04T SUMMARY")
    print(
        summary.to_string(index=False),
        flush=True,
    )

    header("TOP 04T CANDIDATES")

    show = [
        c for c in [
            "04t_rank",
            "priority_rank",
            "pert_iname",
            "structural_developability_category",
            "structural_developability_score",
            "molecular_weight",
            "clogp",
            "tpsa",
            "qed",
            "esol_logS",
            "lipinski_violations",
            "veber_violations",
            "pains_alert_n_04t",
            "brenk_alert_n",
            "docking_status",
            "best_affinity_kcal_mol",
            "protocol_validation",
            "04t_decision",
        ]
        if c in merged.columns
    ]

    print(
        merged[show]
        .head(15)
        .to_string(index=False),
        flush=True,
    )

    header("STAGE 04T COMPLETE")
    print("\nOutputs:", flush=True)
    print(f"  {OUT_ALL}", flush=True)
    print(f"  {OUT_PRIORITIZED}", flush=True)
    print(f"  {OUT_SUMMARY}", flush=True)
    print(f"  {OUT_META}", flush=True)
    print(
        "\nNext: 04U — Integrated evidence matrix",
        flush=True,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
