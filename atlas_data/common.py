from __future__ import annotations
from pathlib import Path
from typing import Any
import hashlib
import json
import re
import time
import requests
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "config" / "dataset_queries.json"
CATALOG_PATH = PROJECT_ROOT / "data" / "catalog" / "atlas_catalog.duckdb"
DISCOVERY_DIR = PROJECT_ROOT / "data" / "discovery"
MANIFEST_DIR = PROJECT_ROOT / "data" / "manifests"

for p in [CATALOG_PATH.parent, DISCOVERY_DIR, MANIFEST_DIR]:
    p.mkdir(parents=True, exist_ok=True)

def load_config():
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

def clean(x: Any) -> str:
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()

def http_get(url: str, *, params=None, timeout=(5, 30), retries=2):
    last = None
    for attempt in range(retries + 1):
        try:
            r = requests.get(
                url,
                params=params,
                timeout=timeout,
                headers={"User-Agent": "ATLAS-data-acquisition/0.1"},
            )
            r.raise_for_status()
            return r
        except Exception as e:
            last = e
            if attempt < retries:
                time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET failed: {url}: {last}")

def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()

def query_terms(text: str) -> dict[str, bool]:
    t = clean(text).lower()
    return {
        "mentions_trastuzumab": any(k in t for k in ["trastuzumab", "herceptin"]),
        "mentions_her2": any(k in t for k in ["her2", "her-2", "erbb2"]),
        "mentions_resistance": any(k in t for k in ["resistan", "refractor", "non-responder", "nonresponder"]),
        "mentions_tgfb": any(k in t for k in ["tgf-beta", "tgf beta", "tgfb", "smad", "tgfbr"]),
        "mentions_checkpoint": any(k in t for k in ["pd-l1", "pdl1", "cd274", "pd-1", "pd1", "pdcd1"]),
        "mentions_transcriptomics": any(k in t for k in ["rna-seq", "rnaseq", "expression", "transcript", "microarray", "l1000"]),
        "mentions_breast": "breast" in t,
    }

def infer_platform(text: str) -> str:
    t = clean(text).lower()
    if "single cell" in t or "single-cell" in t or "scrna" in t:
        return "SINGLE_CELL_RNA_SEQ"
    if "rna-seq" in t or "rnaseq" in t or "transcriptome sequencing" in t:
        return "BULK_RNA_SEQ"
    if "microarray" in t or "array" in t or "affymetrix" in t:
        return "MICROARRAY"
    if "l1000" in t:
        return "L1000"
    return "UNKNOWN"

def canonical_row(**kwargs):
    base = {
        "dataset_id": "",
        "source": "",
        "source_accession": "",
        "title": "",
        "summary": "",
        "organism": "",
        "assay_type": "",
        "platform": "",
        "sample_count": None,
        "publication_date": "",
        "source_url": "",
        "raw_data_available": None,
        "processed_data_available": None,
        "query_family": "",
        "query_text": "",
        "discovered_utc": pd.Timestamp.utcnow().isoformat(),
        "metadata_complete_flag": False,
        "direct_trastuzumab_resistance": False,
        "her2_positive_confirmed": False,
        "resistant_sensitive_groups_defined": False,
        "biological_replication": False,
        "complete_sample_metadata": False,
        "independent_model_or_cohort": True,
        "pd1_pdl1_or_tgfb_relevance": False,
        "phenotype_confidence": "UNRESOLVED",
        "eligibility_score": 0,
        "eligibility_category": "UNSCORED",
        "eligibility_reasons": "",
        "manual_review_required": True,
    }
    base.update(kwargs)
    text = " ".join([clean(base.get("title")), clean(base.get("summary"))])
    flags = query_terms(text)
    base["platform"] = base["platform"] or infer_platform(text)
    base["direct_trastuzumab_resistance"] = (
        flags["mentions_trastuzumab"] and flags["mentions_resistance"]
    )
    base["her2_positive_confirmed"] = flags["mentions_her2"]
    base["pd1_pdl1_or_tgfb_relevance"] = (
        flags["mentions_tgfb"] or flags["mentions_checkpoint"]
    )
    base["metadata_complete_flag"] = bool(
        base["source_accession"] and base["title"] and base["source"]
    )
    return base
