from __future__ import annotations

from pathlib import Path
from typing import Any
import math

import pandas as pd

from atlas_reader import AtlasReader, clean_value


def _col(df: pd.DataFrame, names: list[str], contains: bool = True):
    lowered = {str(c).lower(): c for c in df.columns}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    if contains:
        for c in df.columns:
            lc = str(c).lower()
            if any(n.lower() in lc for n in names):
                return c
    return None


def _finite(v: Any):
    try:
        f = float(v)
        return f if math.isfinite(f) else None
    except Exception:
        return None


def _latest_matching(root: Path, keywords: list[str], required_cols: list[list[str]] | None = None) -> Path | None:
    if not root.exists():
        return None
    candidates: list[Path] = []
    for p in root.rglob("*.csv"):
        name = p.name.lower()
        if not any(k.lower() in name for k in keywords):
            continue
        if required_cols:
            try:
                head = pd.read_csv(p, nrows=3)
            except Exception:
                continue
            ok = all(_col(head, group) is not None for group in required_cols)
            if not ok:
                continue
        candidates.append(p)
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def build_research_statistics(reader: AtlasReader, point_limit: int = 1800):
    payload: dict[str, Any] = {
        "summary": {},
        "deg": {"available": False, "points": [], "stats": {}},
        "dataset_categories": {"available": False, "rows": []},
        "pathways": {"available": False, "rows": [], "source_file": None},
        "candidates": {"available": False, "rows": [], "source_file": None},
        "tgfb": {"available": False, "rows": [], "source_file": None},
        "notes": [
            "Negative CMap connectivity/tau is transcriptional opposition, not proof of reversal or clinical efficacy.",
            "Docking scores are structural-computational evidence, not proof of binding, efficacy, or resistance reversal.",
            "Network and target evidence support plausibility; they do not establish causal mechanism.",
            "Research charts are generated from the files currently present under ATLAS_ROOT. Missing outputs remain explicitly unavailable rather than being fabricated.",
        ],
    }

    # DEG statistics and volcano data
    sig_path = reader.signature_path()
    if sig_path:
        df = reader.read_csv(sig_path)
        lfc = _col(df, ["log2FoldChange", "log2fc", "logfc", "fold_change"])
        padj = _col(df, ["padj", "fdr", "adj_p", "adjusted_p", "qvalue"])
        gene = _col(df, ["gene_symbol", "gene name", "symbol", "gene_name", "gene", "ensembl"])
        if lfc and padj:
            work = df[[c for c in [gene, lfc, padj] if c is not None]].copy()
            work[lfc] = pd.to_numeric(work[lfc], errors="coerce")
            work[padj] = pd.to_numeric(work[padj], errors="coerce")
            work = work.dropna(subset=[lfc, padj])
            work = work[work[padj] > 0]
            sig = work[work[padj] < 0.05]
            strict = sig[sig[lfc].abs() >= 1]
            up = strict[strict[lfc] >= 1]
            down = strict[strict[lfc] <= -1]
            display = work.nsmallest(point_limit, padj).copy()
            display["neglog10_padj"] = -display[padj].clip(lower=1e-300).map(math.log10)
            points = []
            for _, r in display.iterrows():
                label = str(clean_value(r.get(gene)) or "") if gene else ""
                fc = _finite(r[lfc])
                pv = _finite(r[padj])
                y = _finite(r["neglog10_padj"])
                if fc is None or pv is None or y is None:
                    continue
                cls = "up" if (pv < .05 and fc >= 1) else "down" if (pv < .05 and fc <= -1) else "other"
                points.append({"gene": label, "log2fc": fc, "neglog10_padj": y, "padj": pv, "class": cls})
            payload["deg"] = {
                "available": True,
                "source_file": str(sig_path),
                "points": points,
                "stats": {
                    "rows": int(len(work)),
                    "fdr_lt_0_05": int(len(sig)),
                    "strict_abs_log2fc_ge_1": int(len(strict)),
                    "strict_up": int(len(up)),
                    "strict_down": int(len(down)),
                    "displayed_points": len(points),
                },
            }

    # Dataset category composition
    datasets = reader.dataset_rows().get("rows", [])
    if datasets:
        ddf = pd.DataFrame(datasets)
        if "category" in ddf.columns:
            counts = ddf["category"].fillna("UNCLASSIFIED").astype(str).value_counts()
            payload["dataset_categories"] = {
                "available": True,
                "rows": [{"category": str(k), "count": int(v)} for k, v in counts.items()],
            }
            payload["summary"]["datasets"] = int(len(ddf))
            payload["summary"]["primary_validation_datasets"] = int((ddf["category"].astype(str) == "PRIMARY_VALIDATION").sum())

    # Pathway/GSEA results
    pathway_path = _latest_matching(
        reader.root / "results",
        ["pathway", "gsea", "hallmark"],
        required_cols=[["nes", "normalized_enrichment"], ["pathway", "term", "name", "geneset"]],
    )
    if pathway_path:
        pdf = reader.read_csv(pathway_path)
        nes = _col(pdf, ["nes", "normalized_enrichment"])
        name = _col(pdf, ["pathway", "term", "name", "geneset"])
        q = _col(pdf, ["fdr", "padj", "qvalue", "adj_p"])
        if nes and name:
            work = pdf.copy()
            work[nes] = pd.to_numeric(work[nes], errors="coerce")
            work = work.dropna(subset=[nes])
            work["__abs"] = work[nes].abs()
            work = work.nlargest(16, "__abs")
            rows = []
            for _, r in work.iterrows():
                rows.append({
                    "pathway": str(r[name]),
                    "nes": clean_value(r[nes]),
                    "fdr": clean_value(r[q]) if q else None,
                })
            payload["pathways"] = {"available": True, "rows": rows, "source_file": str(pathway_path)}

    # Candidate integrated evidence
    integrated = reader.integrated_path()
    if integrated:
        cdf = reader.read_csv(integrated)
        name = _col(cdf, ["compound", "drug", "pert_iname", "candidate"])
        score = _col(cdf, ["integrated_score", "final_score", "priority_score", "integrated evidence"])
        tau = _col(cdf, ["tau", "connectivity"])
        docking = _col(cdf, ["docking_score", "vina", "binding_affinity"])
        target = _col(cdf, ["target"])
        if name:
            work = cdf.copy()
            if score:
                work[score] = pd.to_numeric(work[score], errors="coerce")
                work = work.sort_values(score, ascending=False, na_position="last")
            rows = []
            for _, r in work.head(20).iterrows():
                rows.append({
                    "candidate": str(clean_value(r.get(name)) or "Unnamed"),
                    "score": clean_value(r.get(score)) if score else None,
                    "tau": clean_value(r.get(tau)) if tau else None,
                    "docking": clean_value(r.get(docking)) if docking else None,
                    "target": clean_value(r.get(target)) if target else None,
                })
            payload["candidates"] = {"available": True, "rows": rows, "source_file": str(integrated)}
            payload["summary"]["integrated_candidates"] = int(len(cdf))

    # TGF-beta ranked validation, if a compatible output exists.
    tgfb_path = _latest_matching(
        reader.root / "results",
        ["tgfb", "tgf_beta", "tgf-beta"],
        required_cols=[["nes", "normalized_enrichment"]],
    )
    if tgfb_path:
        tdf = reader.read_csv(tgfb_path)
        nes = _col(tdf, ["nes", "normalized_enrichment"])
        dataset = _col(tdf, ["dataset", "study", "cohort", "accession"])
        fdr = _col(tdf, ["fdr", "padj", "qvalue"])
        if nes:
            rows = []
            for i, r in tdf.head(30).iterrows():
                value = _finite(r.get(nes))
                if value is None:
                    continue
                rows.append({
                    "dataset": str(clean_value(r.get(dataset)) or f"Row {i+1}") if dataset else f"Row {i+1}",
                    "nes": value,
                    "fdr": clean_value(r.get(fdr)) if fdr else None,
                })
            if rows:
                payload["tgfb"] = {"available": True, "rows": rows, "source_file": str(tgfb_path)}

    payload["summary"]["deg_rows"] = payload["deg"].get("stats", {}).get("rows", 0)
    payload["summary"]["strict_degs"] = payload["deg"].get("stats", {}).get("strict_abs_log2fc_ge_1", 0)
    return payload
