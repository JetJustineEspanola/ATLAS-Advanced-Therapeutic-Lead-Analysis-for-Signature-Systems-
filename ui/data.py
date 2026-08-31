from pathlib import Path
import json
import math
import pandas as pd
import streamlit as st

PROJECT_ROOT = Path(__file__).resolve().parents[1]

def _csv(path):
    try:
        if path.exists():
            return pd.read_csv(path)
    except Exception:
        pass
    return pd.DataFrame()

def _json(path):
    try:
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        pass
    return {}

def clean(value):
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value).strip()

def numeric(value, default=0.0):
    try:
        if value is None or pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default

def safe_col(df, name, default=""):
    if name in df.columns:
        return df[name]
    return pd.Series([default] * len(df), index=df.index)

@st.cache_data(show_spinner=False)
def load_all_data():
    integrated_dir = PROJECT_ROOT / "results" / "cmap" / "integrated_evidence"
    docking_dir = PROJECT_ROOT / "results" / "cmap" / "docking"
    admet_dir = PROJECT_ROOT / "results" / "cmap" / "admet_structural"
    final_dir = PROJECT_ROOT / "results" / "cmap" / "final_prioritization"
    reg_dir = PROJECT_ROOT / "results" / "cmap" / "regulatory_status"
    safety_dir = PROJECT_ROOT / "results" / "cmap" / "safety_screening"
    target_dir = PROJECT_ROOT / "results" / "cmap" / "drug_targets"
    network_dir = PROJECT_ROOT / "results" / "cmap" / "network_integration"

    matrix = _csv(integrated_dir / "ATLAS_integrated_evidence_matrix.csv")
    shortlist = _csv(integrated_dir / "ATLAS_experimental_validation_shortlist.csv")
    summary = _csv(integrated_dir / "ATLAS_integrated_evidence_summary.csv")
    metadata = _json(integrated_dir / "ATLAS_integrated_evidence_metadata.json")
    docking = _csv(docking_dir / "ATLAS_docking_results.csv")

    biology = [
        (
            "Resistance-associated pathway",
            "TGF-β signaling",
            "Strongest Hallmark signal in the discovery analysis; association, not causation.",
        ),
        (
            "Discovery differential expression",
            "1,405 up / 1,215 down",
            "Strict resistance-associated DEG groups used for downstream interpretation.",
        ),
        (
            "Key pathway genes",
            "SMAD2, TGFBR1, ACVR2A/2B, SMAD4, TGFB1, SMURF2, ACVR1; SMAD7 down",
            "Core TGF-β pathway components observed in the discovery analysis.",
        ),
        (
            "Additional signal",
            "PTPN11 increased",
            "Significant increase in the discovery differential-expression analysis.",
        ),
    ]

    return {
        "matrix": matrix,
        "shortlist": shortlist,
        "summary": summary,
        "metadata": metadata,
        "docking": docking,
        "biology": biology,
        "regulatory_exists": (reg_dir / "ATLAS_CMap_regulatory_annotations.csv").exists(),
        "safety_exists": (safety_dir / "ATLAS_CMap_safety_screening.csv").exists(),
        "targets_exists": (target_dir / "ATLAS_CMap_drug_target_annotations.csv").exists(),
        "network_exists": (network_dir / "ATLAS_drug_network_prioritized.csv").exists(),
        "final_exists": (final_dir / "ATLAS_final_candidate_prioritization.csv").exists(),
        "docking_exists": (docking_dir / "ATLAS_docking_results.csv").exists(),
        "admet_exists": (admet_dir / "ATLAS_ADMET_structural_assessment.csv").exists(),
        "integrated_exists": (integrated_dir / "ATLAS_integrated_evidence_matrix.csv").exists(),
    }
