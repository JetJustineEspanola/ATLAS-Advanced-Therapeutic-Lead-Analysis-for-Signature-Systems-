from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
import math
import os

import pandas as pd


def clean_value(v: Any):
    if v is None:
        return None
    if isinstance(v, float) and math.isnan(v):
        return None
    if hasattr(v, "item"):
        try:
            return v.item()
        except Exception:
            pass
    return v


def frame_records(df: pd.DataFrame, limit: int = 1000):
    records = []
    for row in df.head(limit).to_dict(orient="records"):
        records.append({str(k): clean_value(v) for k, v in row.items()})
    return records


def first_existing(paths: Iterable[Path]) -> Path | None:
    for p in paths:
        if p.exists() and p.is_file():
            return p
    return None


@dataclass
class AtlasReader:
    root: Path

    @classmethod
    def from_env(cls):
        return cls(Path(os.getenv("ATLAS_ROOT", "/home/regulus/Documents/ATLAS")).expanduser())

    def read_csv(self, path: Path) -> pd.DataFrame:
        try:
            return pd.read_csv(path, low_memory=False)
        except Exception:
            return pd.DataFrame()

    def datasets_path(self):
        return first_existing([
            self.root / "data/enriched/dataset_candidates_independence_scored.csv",
            self.root / "data/discovery/dataset_candidates_scored.csv",
            self.root / "data/discovery/dataset_candidates.csv",
        ])

    def signature_path(self):
        return first_existing([
            self.root / "results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv",
            self.root / "results/differential_expression/DEGs_resistant_vs_sensitive.csv",
        ])

    def integrated_path(self):
        return first_existing([
            self.root / "results/cmap/integrated_evidence/ATLAS_integrated_evidence_matrix.csv",
        ])

    def find_likely_csv(self, keywords: list[str]) -> Path | None:
        candidates: list[tuple[int, float, Path]] = []
        for base in [self.root / "results", self.root / "data/enriched"]:
            if not base.exists():
                continue
            for p in base.rglob("*.csv"):
                name = p.name.lower()
                score = sum(1 for k in keywords if k.lower() in name)
                if score:
                    candidates.append((score, p.stat().st_mtime, p))
        if not candidates:
            return None
        candidates.sort(reverse=True)
        return candidates[0][2]

    def dataset_rows(self):
        p = self.datasets_path()
        if not p:
            return {"rows": [], "columns": [], "source_file": None}
        df = self.read_csv(p)

        def pick(*names):
            for n in names:
                if n in df.columns:
                    return n
            return None

        mapping = {
            "dataset_id": pick("dataset_id"),
            "source": pick("source"),
            "title": pick("title"),
            "category": pick("independence_aware_category", "eligibility_category"),
            "score": pick("independence_aware_score", "eligibility_score"),
            "modality": pick("modality_class"),
            "phenotype_confidence": pick("phenotype_confidence"),
            "relationship_role": pick("relationship_role"),
            "sample_count": pick("sample_count", "sample_n"),
        }
        out = pd.DataFrame(index=df.index)
        for dest, src in mapping.items():
            out[dest] = df[src] if src else None
        return {"rows": frame_records(out, 5000), "columns": list(out.columns), "source_file": str(p)}

    def generic_rows(self, p: Path | None, limit: int = 1000):
        if not p:
            return {"rows": [], "columns": [], "source_file": None}
        df = self.read_csv(p)
        return {"rows": frame_records(df, limit), "columns": [str(c) for c in df.columns], "source_file": str(p)}

    def candidates(self):
        return self.generic_rows(self.integrated_path())

    def signature(self):
        return self.generic_rows(self.signature_path())

    def cmap(self):
        return self.generic_rows(self.integrated_path() or self.find_likely_csv(["cmap", "tau", "priorit"]))

    def docking(self):
        return self.generic_rows(self.find_likely_csv(["docking", "dock", "vina", "structural"]))

    def dashboard(self):
        warnings: list[str] = []
        drows = self.dataset_rows()["rows"]
        ddf = pd.DataFrame(drows)
        primary = []
        if not ddf.empty and "category" in ddf.columns:
            primary = ddf.loc[ddf["category"].astype(str).eq("PRIMARY_VALIDATION"), "dataset_id"].dropna().astype(str).tolist()

        sig_path = self.signature_path()
        sig_df = self.read_csv(sig_path) if sig_path else pd.DataFrame()
        integrated_path = self.integrated_path()
        integrated = self.read_csv(integrated_path) if integrated_path else pd.DataFrame()

        def find_col(df, needles):
            for c in df.columns:
                lc = str(c).lower()
                if any(n in lc for n in needles):
                    return c
            return None

        docking_col = find_col(integrated, ["docking_score", "docking score", "vina", "binding_affinity"])
        docking_count = int(integrated[docking_col].notna().sum()) if docking_col else len(self.docking()["rows"])

        metrics = [
            {"key": "signature", "label": "Signature Discovery", "value": len(sig_df), "suffix": "DEG rows", "status": "Complete" if len(sig_df) else "Waiting", "tone": "success" if len(sig_df) else "muted"},
            {"key": "cmap", "label": "Integrated Candidates", "value": len(integrated), "suffix": "rows", "status": "Available" if len(integrated) else "Waiting", "tone": "success" if len(integrated) else "muted"},
            {"key": "docking", "label": "Molecular Docking", "value": docking_count, "suffix": "evidence rows", "status": "Available" if docking_count else "Not started", "tone": "warning" if docking_count else "muted"},
            {"key": "validation", "label": "Primary Validation", "value": len(primary), "suffix": "datasets", "status": "Gated", "tone": "success" if primary else "muted"},
        ]

        total = len(ddf)
        exploratory = int((ddf.get("category", pd.Series(dtype=str)).astype(str) == "EXPLORATORY").sum()) if total else 0
        funnel = [
            {"label": "Discovered datasets", "value": total, "ratio": 1.0, "tone": "primary"},
            {"label": "Exploratory", "value": exploratory, "ratio": exploratory / total if total else 0.7, "tone": "secondary"},
            {"label": "Primary validation", "value": len(primary), "ratio": max(0.34, len(primary) / total if total else 0.34), "tone": "warning"},
        ]

        top_candidates = []
        if not integrated.empty:
            name_col = find_col(integrated, ["compound", "drug", "pert_iname", "candidate"])
            cmap_col = find_col(integrated, ["tau", "connectivity"])
            target_col = find_col(integrated, ["target"])
            score_col = find_col(integrated, ["integrated", "final_score", "priority_score"])
            status_col = find_col(integrated, ["priority", "status", "recommendation"])
            df = integrated.copy()
            if score_col:
                numeric = pd.to_numeric(df[score_col], errors="coerce")
                df = df.assign(_sort=numeric).sort_values("_sort", ascending=False, na_position="last")
            for _, row in df.head(8).iterrows():
                top_candidates.append({
                    "name": str(clean_value(row.get(name_col)) or "Unnamed") if name_col else "Unnamed",
                    "connectivity_score": clean_value(row.get(cmap_col)) if cmap_col else None,
                    "docking_score": clean_value(row.get(docking_col)) if docking_col else None,
                    "target": clean_value(row.get(target_col)) if target_col else None,
                    "final_score": clean_value(row.get(score_col)) if score_col else None,
                    "status": str(clean_value(row.get(status_col)) or "Available") if status_col else "Available",
                })

        activity_files = [p for p in [self.datasets_path(), sig_path, integrated_path, self.root / "results/pipeline_state/full_pipeline_state.json", self.root / "data/enriched/ebi_dataset_enrichment_summary.csv"] if p and p.exists()]
        activity = []
        for p in sorted(activity_files, key=lambda x: x.stat().st_mtime, reverse=True)[:6]:
            activity.append({"message": f"Updated {p.relative_to(self.root)}", "timestamp": datetime.fromtimestamp(p.stat().st_mtime, tz=timezone.utc).astimezone().isoformat(timespec="seconds"), "tone": "success"})

        if not sig_path:
            warnings.append("Differential-expression output was not found.")
        if not integrated_path:
            warnings.append("Integrated evidence matrix was not found.")

        newest = max((p.stat().st_mtime for p in activity_files), default=None)
        last_updated = datetime.fromtimestamp(newest, tz=timezone.utc).astimezone().isoformat(timespec="seconds") if newest else None
        return {"project": {"name": "ATLAS", "atlas_root": str(self.root), "last_updated": last_updated}, "metrics": metrics, "funnel": funnel, "top_candidates": top_candidates, "activity": activity, "primary_validation": primary, "warnings": warnings}
