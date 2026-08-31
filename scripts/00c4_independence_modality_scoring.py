#!/usr/bin/env python3
"""
ATLAS — Stage 00C4 v2
Independence- and Modality-Aware Validation Scoring

Fixes:
- checks epigenomic/specialized modalities BEFORE generic RNA/transcriptome terms
- correctly routes CUT&Tag, ATAC-seq, MeRIP-seq, ChIP-seq away from transcriptomic validation
- recognizes scRNA-seq before bulk RNA-seq
- preserves umbrella-series caps and SRR run caps
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import duckdb
import pandas as pd

CATALOG = PROJECT_ROOT / "data/catalog/atlas_catalog.duckdb"
CONFIG = PROJECT_ROOT / "config/dataset_queries.json"
REL = PROJECT_ROOT / "data/enriched/dataset_relationship_audit.csv"
OUTDIR = PROJECT_ROOT / "data/enriched"
OUTDIR.mkdir(parents=True, exist_ok=True)


def b(v):
    if v is None:
        return False
    try:
        if pd.isna(v):
            return False
    except Exception:
        pass
    return bool(v)


def clean(v):
    if v is None:
        return ""
    try:
        if pd.isna(v):
            return ""
    except Exception:
        pass
    return str(v).strip()


def infer_modality(row):
    title = clean(row.get("title")).lower()
    text = " ".join(
        clean(row.get(c)).lower()
        for c in ["title", "assay_type", "platform", "summary"]
    )

    # Explicit title tags take precedence over mixed umbrella/summary wording.
    if "bulk" in title and any(k in title for k in ["rna-seq", "rnaseq", "rna seq"]):
        return "TRANSCRIPTOMIC_RNA_SEQ"

    if any(k in title for k in ["single-cell", "single cell", "scrna-seq", "scrna seq", "scrna"]):
        return "TRANSCRIPTOMIC_SINGLE_CELL"

    # Specific/non-transcriptomic modalities first.
    if any(k in text for k in [
        "cut & tag", "cut&tag", "cut-and-tag", "cut and tag"
    ]):
        return "EPIGENOMIC_CUT_TAG"

    if any(k in text for k in [
        "atac-seq", "atac seq", "chromatin accessibility"
    ]):
        return "EPIGENOMIC_ATAC"

    if any(k in text for k in [
        "merip-seq", "merip seq", "merip", "m6a-seq", "m6a seq"
    ]):
        return "EPITRANSCRIPTOMIC_MERIP"

    if any(k in text for k in [
        "chip-seq", "chip seq"
    ]):
        return "EPIGENOMIC_CHIP"

    # Single-cell before generic RNA-seq.
    if any(k in text for k in [
        "single-cell", "single cell", "scrna-seq", "scrna seq", "scrna"
    ]):
        return "TRANSCRIPTOMIC_SINGLE_CELL"

    if any(k in text for k in [
        "microarray", "expression array", "affymetrix"
    ]):
        return "TRANSCRIPTOMIC_MICROARRAY"

    if any(k in text for k in [
        "rna-seq", "rnaseq", "rna seq", "bulk rna",
        "transcriptome", "expression profiling by high throughput sequencing"
    ]):
        return "TRANSCRIPTOMIC_RNA_SEQ"

    return "UNKNOWN"


def base_score(row, weights):
    checks = {
        "direct_trastuzumab_resistance": b(row.get("direct_trastuzumab_resistance")),
        "her2_positive_confirmed": b(row.get("her2_positive_confirmed")),
        "resistant_sensitive_groups_defined": b(row.get("resistant_sensitive_groups_defined")),
        "biological_replication": b(row.get("biological_replication")),
        "raw_data_available": b(row.get("raw_data_available")),
        "complete_sample_metadata": b(row.get("complete_sample_metadata")),
        "independent_model_or_cohort": b(row.get("independent_model_or_cohort")),
        "pd1_pdl1_or_tgfb_relevance": b(row.get("pd1_pdl1_or_tgfb_relevance")),
    }

    score = 0
    reasons = []

    for key, flag in checks.items():
        if flag:
            score += int(weights[key])
            reasons.append(f"+{weights[key]} {key}")

    phenotype_conf = clean(row.get("phenotype_confidence")).upper()
    if phenotype_conf in {"LOW", "UNRESOLVED", ""}:
        score -= 5
        reasons.append("-5 low_or_unresolved_phenotype_confidence")

    sample_count = pd.to_numeric(
        pd.Series([row.get("sample_count")]),
        errors="coerce",
    ).iloc[0]

    if pd.notna(sample_count) and sample_count < 2:
        score -= 10
        reasons.append("-10 singleton_dataset_unit")

    return max(0, min(100, int(score))), reasons


def main():
    if not CATALOG.exists():
        print(f"ERROR: missing catalog {CATALOG}")
        return 1

    cfg = json.loads(CONFIG.read_text())
    weights = cfg["eligibility_weights"]

    con = duckdb.connect(str(CATALOG))
    df = con.execute("SELECT * FROM datasets").fetchdf()
    con.close()

    if REL.exists():
        rel = pd.read_csv(REL)[[
            "dataset_id",
            "relationship_role",
            "likely_independent_validation_unit",
            "relationship_reason",
        ]]
        df = df.merge(rel, on="dataset_id", how="left")
    else:
        df["relationship_role"] = "NOT_AUDITED"
        df["likely_independent_validation_unit"] = True
        df["relationship_reason"] = ""

    rows = []

    for _, r in df.iterrows():
        score, reasons = base_score(r, weights)
        modality = infer_modality(r)

        groups_defined = b(r.get("resistant_sensitive_groups_defined"))
        replication = b(r.get("biological_replication"))
        her2 = b(r.get("her2_positive_confirmed"))
        direct = b(r.get("direct_trastuzumab_resistance"))

        relationship_role = clean(r.get("relationship_role")) or "NO_OVERLAP_DETECTED"
        source = clean(r.get("source")).upper()
        accession = clean(r.get("source_accession")).upper()

        is_transcriptomic = modality.startswith("TRANSCRIPTOMIC_")
        is_multiomic = modality in {
            "EPIGENOMIC_ATAC",
            "EPIGENOMIC_CUT_TAG",
            "EPIGENOMIC_CHIP",
            "EPITRANSCRIPTOMIC_MERIP",
        }

        is_umbrella = relationship_role in {
            "SUPERSERIES_OR_UMBRELLA_DATASET",
            "DUPLICATE_SAMPLE_SET",
            "HIGH_OVERLAP_DATASET",
        }

        if relationship_role == "SUBSERIES_OR_NESTED_DATASET":
            reasons.append("NESTED_SUBSERIES_ALLOWED_IF_MODALITY_SPECIFIC")

        if is_umbrella:
            reasons.append("CATEGORY_CAP: umbrella/duplicate/high-overlap dataset")

        if source == "SRA" and accession.startswith("SRR"):
            reasons.append("CATEGORY_CAP: SRR run is not a study-level validation unit")

        if is_transcriptomic:
            validation_role = "TRANSCRIPTOMIC_VALIDATION"
        elif is_multiomic:
            validation_role = "SUPPORTING_MULTIOMIC_CONTEXT"
        else:
            validation_role = "UNRESOLVED_MODALITY"

        primary_gate = (
            is_transcriptomic
            and not is_umbrella
            and not (source == "SRA" and accession.startswith("SRR"))
            and groups_defined
            and replication
            and her2
            and direct
        )

        supporting_gate = (
            is_transcriptomic
            and not is_umbrella
            and not (source == "SRA" and accession.startswith("SRR"))
            and groups_defined
            and replication
            and her2
        )

        if source == "SRA" and accession.startswith("SRR"):
            category = "EXPLORATORY"
        elif is_umbrella:
            category = "EXPLORATORY"
        elif validation_role == "SUPPORTING_MULTIOMIC_CONTEXT":
            category = (
                "SUPPORTING_MULTIOMIC_CONTEXT"
                if groups_defined
                else "EXPLORATORY_MULTIOMIC"
            )
        elif primary_gate and score >= 70:
            category = "PRIMARY_VALIDATION"
        elif supporting_gate and score >= 60:
            category = "SUPPORTING_VALIDATION"
        elif is_transcriptomic and (her2 or direct) and score >= 35:
            category = "EXPLORATORY"
        else:
            category = "EXCLUDE_OR_MANUAL_REVIEW"

        rows.append({
            **r.to_dict(),
            "modality_class": modality,
            "validation_role": validation_role,
            "independence_aware_score": score,
            "independence_aware_category": category,
            "independence_aware_reasons": " | ".join(reasons),
        })

    out = pd.DataFrame(rows).sort_values(
        [
            "independence_aware_category",
            "independence_aware_score",
            "source",
            "source_accession",
        ],
        ascending=[True, False, True, True],
    )

    all_path = OUTDIR / "dataset_candidates_independence_scored.csv"
    trans_path = OUTDIR / "transcriptomic_validation_candidates.csv"
    multi_path = OUTDIR / "multiomic_context_candidates.csv"

    out.to_csv(all_path, index=False)

    trans = out[out["validation_role"] == "TRANSCRIPTOMIC_VALIDATION"].sort_values(
        ["independence_aware_score"],
        ascending=False,
    )
    trans.to_csv(trans_path, index=False)

    multi = out[out["validation_role"] == "SUPPORTING_MULTIOMIC_CONTEXT"].sort_values(
        ["independence_aware_score"],
        ascending=False,
    )
    multi.to_csv(multi_path, index=False)

    print("=" * 78)
    print("ATLAS — 00C4 v3 INDEPENDENCE/MODALITY-AWARE VALIDATION SCORING")
    print("=" * 78)

    print("\nCategory counts:")
    print(out["independence_aware_category"].value_counts().to_string())

    print("\nTop transcriptomic validation candidates:")
    cols = [
        "independence_aware_score",
        "independence_aware_category",
        "dataset_id",
        "modality_class",
        "phenotype_confidence",
        "resistant_sensitive_groups_defined",
        "biological_replication",
        "relationship_role",
        "title",
    ]
    print(trans[cols].head(20).to_string(index=False))

    print("\nSupporting multi-omic datasets:")
    if multi.empty:
        print("None")
    else:
        mcols = [
            "independence_aware_score",
            "independence_aware_category",
            "dataset_id",
            "modality_class",
            "relationship_role",
            "title",
        ]
        print(multi[mcols].head(20).to_string(index=False))

    print("\nOutputs:")
    print(all_path)
    print(trans_path)
    print(multi_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
