
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ATLASPaths:
    root: Path

    @property
    def results(self) -> Path:
        return self.root / "results"

    @property
    def deg(self) -> Path:
        return self.results / "differential_expression" / "DEGs_resistant_vs_sensitive_annotated.csv"

    @property
    def external_summary(self) -> Path:
        return self.results / "external_validation" / "ATLAS_external_validation_summary.csv"

    @property
    def external_meta(self) -> Path:
        return self.results / "external_validation" / "ATLAS_external_validation_metadata.json"

    @property
    def cmap_priority(self) -> Path:
        return self.results / "cmap" / "prioritization" / "ATLAS_CMap_prioritized.csv"

    @property
    def drug_candidates(self) -> Path:
        return self.results / "cmap" / "drug_filter" / "ATLAS_CMap_drug_candidates.csv"

    @property
    def regulatory(self) -> Path:
        return self.results / "cmap" / "regulatory_status" / "ATLAS_CMap_regulatory_annotations.csv"

    @property
    def safety(self) -> Path:
        return self.results / "cmap" / "safety_screening" / "ATLAS_CMap_safety_screening.csv"

    @property
    def safety_priority(self) -> Path:
        return self.results / "cmap" / "safety_screening" / "ATLAS_CMap_safety_prioritized.csv"

    @property
    def targets(self) -> Path:
        return self.results / "cmap" / "drug_targets" / "ATLAS_CMap_drug_target_annotations.csv"

    @property
    def target_pairs(self) -> Path:
        return self.results / "cmap" / "drug_targets" / "ATLAS_CMap_drug_target_pairs.csv"

    @property
    def network(self) -> Path:
        return self.results / "cmap" / "network_integration" / "ATLAS_drug_network_prioritized.csv"

    @property
    def network_summary(self) -> Path:
        return self.results / "cmap" / "network_integration" / "ATLAS_network_summary.csv"

    def available_output_tables(self) -> dict[str, Path]:
        known = {
            "Discovery DEGs": self.deg,
            "External validation summary": self.external_summary,
            "Drug candidates (04M)": self.drug_candidates,
            "Regulatory annotations (04N)": self.regulatory,
            "Safety screen (04O)": self.safety,
            "Drug-target annotations (04P)": self.targets,
            "Drug-target pairs (04P)": self.target_pairs,
            "Network-prioritized drugs (04Q)": self.network,
        }
        return {k: v for k, v in known.items() if v.exists()}


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame()


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _normalize_keys(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    for c in ["pert_id", "pert_iname"]:
        if c in out.columns:
            out[c] = out[c].astype("string").str.strip()
    return out


def _merge_candidate_layer(base: pd.DataFrame, layer: pd.DataFrame) -> pd.DataFrame:
    if base.empty:
        return layer.copy()
    if layer.empty:
        return base.copy()

    common = [c for c in ["pert_id", "pert_iname"] if c in base.columns and c in layer.columns]
    if len(common) < 2:
        return base

    layer = layer.drop_duplicates(common).copy()
    overlap_nonkeys = [c for c in layer.columns if c in base.columns and c not in common]
    if overlap_nonkeys:
        layer = layer.drop(columns=overlap_nonkeys)

    return base.merge(layer, on=common, how="left")


def _coalesce_series(df: pd.DataFrame, names: list[str], default=None) -> pd.Series:
    result = pd.Series([default] * len(df), index=df.index, dtype="object")
    for name in names:
        if name not in df.columns:
            continue
        mask = result.isna() | result.eq("") | result.eq(default)
        result.loc[mask] = df.loc[mask, name]
    return result


def _classify_cmap(row: pd.Series) -> str:
    tier = row.get("priority_tier_number")
    try:
        tier = int(float(tier))
    except Exception:
        tier = None

    if tier == 1:
        return "STRONG"
    if tier == 2:
        return "MODERATE"
    if tier == 3:
        return "WEAK"
    if tier == 4:
        return "WEAK"
    return "PENDING"


def _classify_safety(row: pd.Series) -> str:
    value = str(row.get("safety_screening_recommendation", "") or "")
    if value == "PASS_PRELIMINARY_SCREEN":
        return "PASS"
    if value == "CAUTION_MANUAL_REVIEW":
        return "CAUTION"
    if value == "HIGH_RISK_DEPRIORITIZE":
        return "HIGH RISK"
    if value == "INSUFFICIENT_SAFETY_DATA":
        return "INSUFFICIENT"
    return "PENDING"


def _classify_target(row: pd.Series) -> str:
    value = str(row.get("target_support_category", "") or "")
    mapping = {
        "STRONG_TARGET_SUPPORT": "STRONG",
        "MODERATE_TARGET_SUPPORT": "MODERATE",
        "WEAK_TARGET_SUPPORT": "WEAK",
        "NO_TARGET_SUPPORT_FOUND": "INSUFFICIENT",
        "QUERY_ERROR": "INSUFFICIENT",
    }
    return mapping.get(value, "PENDING")


def _classify_network(row: pd.Series) -> str:
    value = str(row.get("best_network_support_category", "") or "")
    mapping = {
        "STRONG_NETWORK_SUPPORT": "STRONG",
        "MODERATE_NETWORK_SUPPORT": "MODERATE",
        "WEAK_NETWORK_SUPPORT": "WEAK",
        "NO_NETWORK_SUPPORT": "INSUFFICIENT",
    }
    return mapping.get(value, "PENDING")


def _classify_regulatory(row: pd.Series) -> str:
    value = str(row.get("regulatory_evidence_category", "") or "")
    if value == "FDA_APPLICATION_AND_LABEL_EVIDENCE":
        return "STRONG"
    if value in {"FDA_APPLICATION_RECORD_FOUND", "FDA_LABEL_EVIDENCE_FOUND"}:
        return "MODERATE"
    if value == "CLINICAL_TRIAL_EVIDENCE_ONLY":
        return "WEAK"
    if value == "NO_US_REGULATORY_OR_TRIAL_EVIDENCE_FOUND":
        return "INSUFFICIENT"
    return "PENDING"


def _overall_evidence(row: pd.Series) -> str:
    vals = [
        row.get("cmap_evidence"),
        row.get("safety_evidence"),
        row.get("target_evidence"),
        row.get("network_evidence"),
    ]

    if row.get("safety_evidence") == "HIGH RISK":
        return "CAUTION"

    strong = sum(v == "STRONG" for v in vals)
    moderate = sum(v in {"STRONG", "MODERATE", "PASS"} for v in vals)

    if strong >= 2 and moderate >= 3:
        return "STRONG"
    if moderate >= 2:
        return "MODERATE"
    if any(v in {"WEAK", "PASS"} for v in vals):
        return "WEAK"
    return "INSUFFICIENT"


def _decision(row: pd.Series) -> tuple[str, str]:
    safety = row.get("safety_evidence")
    cmap = row.get("cmap_evidence")
    target = row.get("target_evidence")
    network = row.get("network_evidence")

    if safety == "HIGH RISK":
        return "DEPRIORITIZE", "High-risk preliminary safety evidence outweighs downstream promise."

    if cmap == "STRONG" and target in {"STRONG", "MODERATE"} and network in {"STRONG", "MODERATE"} and safety in {"PASS", "CAUTION"}:
        return "PROCEED TO DOCKING", "Strong CMap reversal is supported by target and resistance-network evidence."

    if cmap in {"STRONG", "MODERATE"} and safety != "HIGH RISK" and target in {"STRONG", "MODERATE", "WEAK"}:
        return "PRIORITIZE", "The candidate retains multi-layer evidence and should remain in the shortlist."

    if safety == "CAUTION" or target == "INSUFFICIENT" or network == "INSUFFICIENT":
        return "MANUAL REVIEW", "The candidate has useful evidence, but a key layer is uncertain or cautionary."

    if cmap == "WEAK":
        return "HOLD — INSUFFICIENT EVIDENCE", "Current evidence does not justify escalation."

    return "HOLD — INSUFFICIENT EVIDENCE", "ATLAS is still waiting for enough downstream evidence to make a strong recommendation."


def build_candidate_table(paths: ATLASPaths) -> pd.DataFrame:
    # Start from the broadest usable candidate table.
    base = _read_csv(paths.drug_candidates)

    if base.empty:
        base = _read_csv(paths.safety)
    if base.empty:
        base = _read_csv(paths.targets)
    if base.empty:
        base = _read_csv(paths.network)

    if base.empty:
        return pd.DataFrame()

    base = _normalize_keys(base)

    for path in [paths.regulatory, paths.safety, paths.targets, paths.network]:
        layer = _normalize_keys(_read_csv(path))
        base = _merge_candidate_layer(base, layer)

    if "priority_rank" not in base.columns:
        base["priority_rank"] = np.arange(1, len(base) + 1)

    if "priority_tier" not in base.columns:
        if "priority_tier_number" in base.columns:
            base["priority_tier"] = base["priority_tier_number"].map(
                lambda x: f"Tier {int(x)}" if pd.notna(x) else "Unknown"
            )
        else:
            base["priority_tier"] = "Unknown"

    base["cmap_evidence"] = base.apply(_classify_cmap, axis=1)
    base["external_evidence"] = "GLOBAL"  # 03B validates the resistance signature, not each drug.
    base["safety_evidence"] = base.apply(_classify_safety, axis=1)
    base["target_evidence"] = base.apply(_classify_target, axis=1)
    base["network_evidence"] = base.apply(_classify_network, axis=1)
    base["regulatory_evidence"] = base.apply(_classify_regulatory, axis=1)
    base["overall_evidence"] = base.apply(_overall_evidence, axis=1)

    decisions = base.apply(_decision, axis=1)
    base["recommended_action"] = [x[0] for x in decisions]
    base["decision_reason"] = [x[1] for x in decisions]

    action_priority = {
        "PROCEED TO DOCKING": 0,
        "PRIORITIZE": 1,
        "MANUAL REVIEW": 2,
        "HOLD — INSUFFICIENT EVIDENCE": 3,
        "DEPRIORITIZE": 4,
    }
    evidence_priority = {
        "STRONG": 0,
        "MODERATE": 1,
        "WEAK": 2,
        "CAUTION": 3,
        "INSUFFICIENT": 4,
    }

    base["_action_sort"] = base["recommended_action"].map(action_priority).fillna(9)
    base["_evidence_sort"] = base["overall_evidence"].map(evidence_priority).fillna(9)
    base["priority_rank"] = pd.to_numeric(base["priority_rank"], errors="coerce")

    base = base.sort_values(
        ["_action_sort", "_evidence_sort", "priority_rank"],
        ascending=[True, True, True],
        na_position="last",
    ).reset_index(drop=True)

    base["ui_rank"] = np.arange(1, len(base) + 1)

    return base.drop(columns=["_action_sort", "_evidence_sort"])


def load_external_validation(paths: ATLASPaths) -> dict[str, Any]:
    summary = _read_csv(paths.external_summary)
    meta = _read_json(paths.external_meta)

    if summary.empty and not meta:
        return {}

    requested = meta.get("datasets_requested", [])
    successful = meta.get("datasets_successfully_analyzed", [])

    return {
        "summary": summary,
        "dataset_count": len(requested) if requested else len(summary),
        "successful_count": len(successful) if successful else len(summary),
        "metadata": meta,
    }


def _find_pathway_files(results_root: Path) -> list[Path]:
    candidates = []
    for folder_name in ["pathway_analysis", "pathways", "enrichment"]:
        folder = results_root / folder_name
        if folder.exists():
            candidates.extend(folder.rglob("*.csv"))
    return candidates


def load_resistance_summary(paths: ATLASPaths) -> dict[str, Any]:
    out: dict[str, Any] = {}

    deg = _read_csv(paths.deg)
    if not deg.empty:
        symbol_col = next(
            (c for c in ["Gene name", "gene_symbol", "Gene Symbol", "symbol"] if c in deg.columns),
            None,
        )
        if "log2FoldChange" in deg.columns and "padj" in deg.columns:
            lfc = pd.to_numeric(deg["log2FoldChange"], errors="coerce")
            padj = pd.to_numeric(deg["padj"], errors="coerce")
            sig = deg[(padj < 0.05) & (lfc.abs() >= 1)].copy()
            out["significant_deg_n"] = int(len(sig))

            if symbol_col:
                sig["_lfc"] = pd.to_numeric(sig["log2FoldChange"], errors="coerce")
                sig["_padj"] = pd.to_numeric(sig["padj"], errors="coerce")
                sig = sig.sort_values(["_padj", "_lfc"], ascending=[True, False])
                genes = []
                for _, r in sig.head(12).iterrows():
                    genes.append(
                        {
                            "Gene": str(r[symbol_col]),
                            "Direction": "UP" if r["_lfc"] > 0 else "DOWN",
                            "log2FC": round(float(r["_lfc"]), 3),
                            "FDR": float(r["_padj"]) if pd.notna(r["_padj"]) else np.nan,
                        }
                    )
                out["key_genes"] = genes

    # Try to discover a pathway table rather than hard-code a filename.
    pathway_files = _find_pathway_files(paths.results)
    best = None
    for p in pathway_files:
        try:
            df = pd.read_csv(p)
        except Exception:
            continue

        name_col = next(
            (c for c in ["Term", "pathway", "Pathway", "Name", "term"] if c in df.columns),
            None,
        )
        nes_col = next(
            (c for c in ["NES", "nes", "Normalized Enrichment Score"] if c in df.columns),
            None,
        )
        fdr_col = next(
            (c for c in ["FDR q-val", "FDR", "fdr", "padj", "Adjusted P-value"] if c in df.columns),
            None,
        )

        if name_col and (nes_col or fdr_col) and not df.empty:
            local = df.copy()
            if fdr_col:
                local["_fdr"] = pd.to_numeric(local[fdr_col], errors="coerce")
            else:
                local["_fdr"] = np.nan
            if nes_col:
                local["_nes"] = pd.to_numeric(local[nes_col], errors="coerce")
            else:
                local["_nes"] = np.nan

            if local["_fdr"].notna().any():
                row = local.sort_values(["_fdr", "_nes"], ascending=[True, False]).iloc[0]
            else:
                row = local.sort_values("_nes", ascending=False).iloc[0]

            best = {
                "name": str(row[name_col]),
                "fdr": float(row["_fdr"]) if pd.notna(row["_fdr"]) else np.nan,
                "nes": float(row["_nes"]) if pd.notna(row["_nes"]) else np.nan,
                "source": p,
            }
            break

    if best:
        out["top_pathway_name"] = best["name"]
        out["top_pathway_fdr"] = best["fdr"]
        out["top_pathway_nes"] = best["nes"]
        out["interpretation"] = (
            f"{best['name']} is one of the strongest pathway-level features associated with the resistant transcriptional state. "
            "This supports biological relevance but does not establish causality."
        )

    return out


def build_pipeline_status(paths: ATLASPaths) -> pd.DataFrame:
    rows = []

    def add(stage: str, label: str, exists: bool, detail: str, running_hint: bool = False):
        status = "Complete" if exists else ("Pending" if not running_hint else "In progress")
        rows.append({"stage": stage, "label": label, "status": status, "detail": detail})

    add("02", "Differential expression", paths.deg.exists(), "Resistant vs sensitive transcriptomic discovery.")
    add("03", "Pathway analysis", bool(_find_pathway_files(paths.results)), "Resistance-associated pathway enrichment.")
    add("03B", "External validation", paths.external_summary.exists(), "Independent GEO validation datasets.")
    add("03C", "Consensus signature", any(p.exists() for p in [
        paths.results / "consensus_signature" / "ATLAS_consensus_resistance_signature.csv",
        paths.results / "consensus_signature" / "consensus_resistance_signature.csv",
    ]), "Externally supported resistance signature.")
    add("04M", "Compound identity", paths.drug_candidates.exists(), "CMap compounds resolved to usable identities.")
    add("04N", "Regulatory / clinical evidence", paths.regulatory.exists(), "FDA/openFDA and ClinicalTrials.gov evidence.")
    add("04O", "Safety screening", paths.safety.exists(), "Cytotoxicity, promiscuity, PAINS, hazard evidence.")
    add("04P", "Drug-target annotation", paths.targets.exists(), "PubChem/ChEMBL target evidence.")
    add("04Q", "Network integration", paths.network.exists(), "STRING links between targets and resistance biology.")
    add("04R", "Final prioritization", (paths.results / "cmap" / "final_prioritization").exists(), "Integrated candidate ranking.")
    add("04S", "Docking", (paths.results / "docking").exists(), "Target-supported docking for top pairs.")
    add("04T", "ADMET", (paths.results / "admet").exists(), "Structural and ADMET assessment.")
    add("04U", "Integrated evidence", (paths.results / "integrated_evidence").exists(), "Final evidence matrix.")

    return pd.DataFrame(rows)
