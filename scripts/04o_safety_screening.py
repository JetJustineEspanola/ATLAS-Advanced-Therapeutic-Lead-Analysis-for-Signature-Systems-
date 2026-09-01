#!/usr/bin/env python3
"""
ATLAS — Stage 04O
Cytotoxicity / Promiscuity / PAINS / Safety Screening

Primary input
-------------
results/cmap/drug_filter/ATLAS_CMap_drug_candidates.csv

Optional input
--------------
results/cmap/regulatory_status/ATLAS_CMap_regulatory_annotations.csv

The stage does NOT require 04N to be finished. If 04N exists, selected
regulatory columns are merged. If it does not exist, 04O proceeds from 04M.

Outputs
-------
results/cmap/safety_screening/
    ATLAS_CMap_safety_screening.csv
    ATLAS_CMap_safety_prioritized.csv
    ATLAS_CMap_safety_manual_review.csv
    ATLAS_CMap_safety_summary.csv
    ATLAS_CMap_safety_metadata.json
    cache/

Evidence sources
----------------
- PubChem PUG REST assay summaries (batched by CID)
- PubChem PUG REST structural properties (batched)
- PubChem PUG-View "Safety and Hazards" annotations
- Optional RDKit PAINS filters if RDKit is installed

Important guardrails
--------------------
- An active PubChem assay is not proof of clinically relevant toxicity.
- A viability/cytotoxicity assay keyword is only a screening signal.
- PAINS is an assay-interference flag, not a toxicity verdict.
- No PubChem hazard annotation does NOT mean a compound is safe.
- This is computational prioritization, not experimental validation.
"""

from __future__ import annotations

import argparse
import io
import json
import math
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from requests.adapters import HTTPAdapter


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

OPTIONAL_04N = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "regulatory_status"
    / "ATLAS_CMap_regulatory_annotations.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "safety_screening"
)
CACHE_DIR = OUTPUT_DIR / "cache"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUTPUT_DIR / "ATLAS_CMap_safety_screening.csv"
OUT_PRIORITIZED = OUTPUT_DIR / "ATLAS_CMap_safety_prioritized.csv"
OUT_MANUAL = OUTPUT_DIR / "ATLAS_CMap_safety_manual_review.csv"
OUT_SUMMARY = OUTPUT_DIR / "ATLAS_CMap_safety_summary.csv"
OUT_META = OUTPUT_DIR / "ATLAS_CMap_safety_metadata.json"

ASSAY_CACHE = CACHE_DIR / "pubchem_assay_summary_cache.csv"
STRUCTURE_CACHE = CACHE_DIR / "pubchem_structure_cache.csv"
HAZARD_CACHE = CACHE_DIR / "pubchem_hazard_cache.jsonl"


# ---------------------------------------------------------------------
# PubChem
# ---------------------------------------------------------------------

PUG_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUG_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"

USER_AGENT = "ATLAS-Safety-Screening/1.0 (research workflow)"

_thread_local = threading.local()


class RateLimiter:
    def __init__(self, rps: float):
        self.interval = 1.0 / max(float(rps), 0.01)
        self.lock = threading.Lock()
        self.next_allowed = 0.0

    def wait(self) -> None:
        with self.lock:
            now = time.monotonic()
            delay = max(0.0, self.next_allowed - now)
            if delay:
                time.sleep(delay)
                now = time.monotonic()
            self.next_allowed = max(now, self.next_allowed) + self.interval


PUG_LIMITER = RateLimiter(4.0)
PUG_VIEW_LIMITER = RateLimiter(3.0)


def get_session() -> requests.Session:
    if not hasattr(_thread_local, "session"):
        s = requests.Session()
        s.headers.update({
            "User-Agent": USER_AGENT,
            "Accept": "*/*",
        })
        adapter = HTTPAdapter(
            pool_connections=8,
            pool_maxsize=8,
            max_retries=0,
        )
        s.mount("https://", adapter)
        s.mount("http://", adapter)
        _thread_local.session = s
    return _thread_local.session


def safe_get(
    url: str,
    *,
    params: dict[str, Any] | None = None,
    limiter: RateLimiter,
    connect_timeout: float = 4.0,
    read_timeout: float = 12.0,
    retries: int = 1,
) -> requests.Response | None:
    session = get_session()

    for attempt in range(retries + 1):
        limiter.wait()
        try:
            r = session.get(
                url,
                params=params,
                timeout=(connect_timeout, read_timeout),
            )
        except requests.RequestException:
            r = None

        if r is not None and r.status_code == 200:
            return r

        if r is not None and r.status_code in {400, 404}:
            return r

        if attempt < retries:
            time.sleep(min(0.75 * (2 ** attempt), 4.0))

    return r


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def header(text: str) -> None:
    print("\n" + "=" * 76, flush=True)
    print(text, flush=True)
    print("=" * 76, flush=True)


def clean_text(x: Any) -> str:
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return re.sub(r"\s+", " ", str(x).strip())


def normalize_cid(x: Any) -> str:
    text = clean_text(x)
    if not text:
        return ""
    try:
        return str(int(float(text)))
    except Exception:
        return ""


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(tmp, path)


def chunks(items: list[str], n: int):
    for i in range(0, len(items), n):
        yield items[i:i+n]


def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def flatten_json_strings(obj: Any) -> list[str]:
    out: list[str] = []

    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(k, str):
                out.append(k)
            out.extend(flatten_json_strings(v))
    elif isinstance(obj, list):
        for v in obj:
            out.extend(flatten_json_strings(v))
    elif isinstance(obj, (str, int, float)):
        out.append(str(obj))

    return out


# ---------------------------------------------------------------------
# Candidate selection / optional 04N merge
# ---------------------------------------------------------------------

def select_candidates(df: pd.DataFrame, mode: str, tier2_top: int) -> pd.DataFrame:
    if mode == "full":
        return df.copy()

    tier = pd.to_numeric(df["priority_tier_number"], errors="coerce")
    rank = pd.to_numeric(df["priority_rank"], errors="coerce")

    t1 = df[tier == 1].copy()
    t2 = (
        df[tier == 2]
        .assign(_rank=rank[tier == 2])
        .sort_values("_rank", na_position="last")
        .head(tier2_top)
        .drop(columns="_rank")
    )

    return (
        pd.concat([t1, t2], ignore_index=True)
        .drop_duplicates(["pert_id", "pert_iname"])
        .sort_values("priority_rank", na_position="last")
        .reset_index(drop=True)
    )


def merge_optional_04n(df: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    if not OPTIONAL_04N.exists():
        return df, False

    try:
        reg = pd.read_csv(OPTIONAL_04N)
    except Exception:
        return df, False

    keys = [c for c in ["pert_id", "pert_iname"] if c in reg.columns and c in df.columns]
    if len(keys) < 2:
        return df, False

    useful = [
        "pert_id",
        "pert_iname",
        "fda_application_record_found",
        "fda_label_record_found",
        "clinical_trial_record_found",
        "regulatory_evidence_category",
        "clinical_evidence_level",
        "manual_regulatory_review_required",
    ]
    useful = [c for c in useful if c in reg.columns]

    reg = reg[useful].drop_duplicates(["pert_id", "pert_iname"])

    return (
        df.merge(
            reg,
            on=["pert_id", "pert_iname"],
            how="left",
            suffixes=("", "_04n"),
        ),
        True,
    )


# ---------------------------------------------------------------------
# Batched PubChem structure retrieval
# ---------------------------------------------------------------------

def load_structure_cache() -> pd.DataFrame:
    if STRUCTURE_CACHE.exists():
        try:
            x = pd.read_csv(STRUCTURE_CACHE, dtype={"pubchem_cid": "string"})
            if not x.empty:
                return x
        except Exception:
            pass
    return pd.DataFrame()


def fetch_structure_batch(cids: list[str]) -> pd.DataFrame:
    if not cids:
        return pd.DataFrame()

    cid_text = ",".join(cids)
    url = (
        f"{PUG_REST}/compound/cid/{cid_text}/property/"
        "CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula/"
        "CSV"
    )

    r = safe_get(
        url,
        limiter=PUG_LIMITER,
        retries=1,
    )

    if r is None or r.status_code != 200:
        return pd.DataFrame()

    try:
        x = pd.read_csv(io.StringIO(r.text))
    except Exception:
        return pd.DataFrame()

    cid_col = find_column(x, ["CID"])
    if cid_col is None:
        return pd.DataFrame()

    x["pubchem_cid"] = x[cid_col].map(normalize_cid)

    rename = {}
    for wanted, variants in {
        "canonical_smiles": ["CanonicalSMILES", "ConnectivitySMILES"],
        "isomeric_smiles": ["IsomericSMILES", "SMILES"],
        "inchikey": ["InChIKey"],
        "molecular_formula": ["MolecularFormula"],
    }.items():
        c = find_column(x, variants)
        if c:
            rename[c] = wanted

    x = x.rename(columns=rename)

    keep = [
        c for c in [
            "pubchem_cid",
            "canonical_smiles",
            "isomeric_smiles",
            "inchikey",
            "molecular_formula",
        ]
        if c in x.columns
    ]

    return x[keep].drop_duplicates("pubchem_cid")


def ensure_structures(df: pd.DataFrame, batch_size: int) -> pd.DataFrame:
    cache = load_structure_cache()
    cached_cids = set(
        cache.get("pubchem_cid", pd.Series(dtype="string"))
        .astype("string")
        .dropna()
        .tolist()
    )

    cids = sorted(set(df["pubchem_cid"].dropna().astype(str)) - cached_cids)
    new_frames = []

    if cids:
        batches = list(chunks(cids, batch_size))
        print(f"Structure enrichment: {len(cids):,} CIDs in {len(batches):,} batches", flush=True)

        for i, batch in enumerate(batches, start=1):
            x = fetch_structure_batch(batch)
            if not x.empty:
                new_frames.append(x)
            if i == 1 or i % 10 == 0 or i == len(batches):
                print(f"  structure batches {i}/{len(batches)}", flush=True)

    if new_frames:
        new = pd.concat(new_frames, ignore_index=True)
        cache = pd.concat([cache, new], ignore_index=True)
        cache = cache.drop_duplicates("pubchem_cid", keep="last")
        atomic_csv(cache, STRUCTURE_CACHE)

    if cache.empty:
        return df

    cols_to_add = [
        c for c in [
            "pubchem_cid",
            "canonical_smiles",
            "isomeric_smiles",
            "inchikey",
            "molecular_formula",
        ]
        if c in cache.columns
    ]

    return df.merge(
        cache[cols_to_add],
        on="pubchem_cid",
        how="left",
        suffixes=("", "_pubchem"),
    )


# ---------------------------------------------------------------------
# PubChem assay summaries
# ---------------------------------------------------------------------

CYTOTOXICITY_TERMS = re.compile(
    r"\b(?:"
    r"cytotox|cell\s*viability|cell\s*death|cellular\s*viability|"
    r"growth\s*inhibition|proliferation\s*inhibition|antiprolifer|"
    r"apoptos|necros|lethal|toxicity"
    r")\b",
    re.I,
)


def load_assay_cache() -> pd.DataFrame:
    if ASSAY_CACHE.exists():
        try:
            x = pd.read_csv(ASSAY_CACHE, dtype={"pubchem_cid": "string"})
            if not x.empty:
                return x
        except Exception:
            pass
    return pd.DataFrame()


def summarize_assay_table(x: pd.DataFrame) -> pd.DataFrame:
    if x.empty:
        return pd.DataFrame()

    cid_col = find_column(x, ["CID"])
    outcome_col = find_column(x, ["Activity Outcome", "ActivityOutcome"])
    assay_name_col = find_column(x, ["Assay Name", "AssayName"])
    gene_col = find_column(x, ["Target GeneID", "Target Gene ID", "GeneID"])
    accession_col = find_column(x, ["Target Accession", "TargetAccession"])
    activity_name_col = find_column(x, ["Activity Name", "ActivityName"])

    if cid_col is None:
        return pd.DataFrame()

    x = x.copy()
    x["pubchem_cid"] = x[cid_col].map(normalize_cid)

    if outcome_col:
        outcome = x[outcome_col].astype(str).str.strip().str.lower()
    else:
        outcome = pd.Series([""] * len(x), index=x.index)

    x["_active"] = outcome.eq("active")
    x["_inactive"] = outcome.eq("inactive")
    x["_inconclusive"] = outcome.str.contains("inconclusive", na=False)

    if assay_name_col:
        names = x[assay_name_col].fillna("").astype(str)
    else:
        names = pd.Series([""] * len(x), index=x.index)

    x["_cytotox_keyword"] = names.str.contains(CYTOTOXICITY_TERMS, na=False)

    if activity_name_col:
        act_names = x[activity_name_col].fillna("").astype(str)
        x["_cytotox_keyword"] |= act_names.str.contains(CYTOTOXICITY_TERMS, na=False)

    rows = []

    for cid, g in x.groupby("pubchem_cid", dropna=False):
        if not cid:
            continue

        active_n = int(g["_active"].sum())
        inactive_n = int(g["_inactive"].sum())
        interpretable = active_n + inactive_n

        cyto_rows = g[g["_cytotox_keyword"]]
        cyto_active = int(cyto_rows["_active"].sum())

        genes = []
        if gene_col:
            genes = sorted({
                clean_text(v)
                for v in g.loc[g["_active"], gene_col].tolist()
                if clean_text(v)
            })

        accessions = []
        if accession_col:
            accessions = sorted({
                clean_text(v)
                for v in g.loc[g["_active"], accession_col].tolist()
                if clean_text(v)
            })

        rows.append({
            "pubchem_cid": cid,
            "pubchem_assay_rows_n": int(len(g)),
            "pubchem_active_assay_n": active_n,
            "pubchem_inactive_assay_n": inactive_n,
            "pubchem_inconclusive_assay_n": int(g["_inconclusive"].sum()),
            "pubchem_interpretable_assay_n": interpretable,
            "pubchem_active_fraction": (
                active_n / interpretable if interpretable else math.nan
            ),
            "cytotoxicity_keyword_assay_n": int(len(cyto_rows)),
            "cytotoxicity_keyword_active_n": cyto_active,
            "active_target_gene_ids": " | ".join(genes[:100]),
            "active_target_accessions": " | ".join(accessions[:100]),
        })

    return pd.DataFrame(rows)


def fetch_assay_batch(cids: list[str]) -> pd.DataFrame:
    if not cids:
        return pd.DataFrame()

    url = (
        f"{PUG_REST}/compound/cid/{','.join(cids)}/"
        "assaysummary/CSV"
    )

    r = safe_get(
        url,
        limiter=PUG_LIMITER,
        read_timeout=20.0,
        retries=1,
    )

    if r is None or r.status_code != 200:
        return pd.DataFrame()

    try:
        raw = pd.read_csv(io.StringIO(r.text), low_memory=False)
    except Exception:
        return pd.DataFrame()

    return summarize_assay_table(raw)


def ensure_assays(df: pd.DataFrame, batch_size: int) -> pd.DataFrame:
    cache = load_assay_cache()

    cached_cids = set(
        cache.get("pubchem_cid", pd.Series(dtype="string"))
        .astype("string")
        .dropna()
        .tolist()
    )

    cids = sorted(set(df["pubchem_cid"].dropna().astype(str)) - cached_cids)
    new_frames = []

    if cids:
        batches = list(chunks(cids, batch_size))
        print(f"Assay enrichment: {len(cids):,} CIDs in {len(batches):,} batches", flush=True)

        for i, batch in enumerate(batches, start=1):
            x = fetch_assay_batch(batch)
            if not x.empty:
                new_frames.append(x)
            if i == 1 or i % 5 == 0 or i == len(batches):
                print(f"  assay batches {i}/{len(batches)}", flush=True)

    if new_frames:
        new = pd.concat(new_frames, ignore_index=True)
        cache = pd.concat([cache, new], ignore_index=True)
        cache = cache.drop_duplicates("pubchem_cid", keep="last")
        atomic_csv(cache, ASSAY_CACHE)

    if cache.empty:
        return df

    return df.merge(
        cache,
        on="pubchem_cid",
        how="left",
        suffixes=("", "_assay"),
    )


# ---------------------------------------------------------------------
# PUG-View safety/hazard annotations
# ---------------------------------------------------------------------

HAZARD_CODE_RE = re.compile(r"\bH(?:3\d{2}|4\d{2})\b", re.I)

SEVERE_HAZARD_CODES = {
    "H300", "H301", "H310", "H311", "H330", "H331",
    "H340", "H350", "H360",
    "H370", "H372",
}

HIGH_CONCERN_TEXT = re.compile(
    r"\b("
    r"carcinogen|carcinogenic|mutagen|mutagenic|"
    r"reproductive toxicity|fatal if|toxic if|"
    r"organ damage|suspected of causing cancer"
    r")\b",
    re.I,
)


def load_hazard_cache() -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}

    if not HAZARD_CACHE.exists():
        return out

    with HAZARD_CACHE.open("r", encoding="utf-8") as f:
        for line in f:
            try:
                obj = json.loads(line)
            except Exception:
                continue
            cid = normalize_cid(obj.get("pubchem_cid"))
            if cid:
                out[cid] = obj

    return out


def fetch_hazard(cid: str) -> dict[str, Any]:
    url = f"{PUG_VIEW}/data/compound/{cid}/JSON"

    r = safe_get(
        url,
        params={"heading": "Safety and Hazards"},
        limiter=PUG_VIEW_LIMITER,
        read_timeout=12.0,
        retries=1,
    )

    if r is None or r.status_code != 200:
        return {
            "pubchem_cid": cid,
            "pubchem_safety_annotation_found": False,
            "hazard_statement_codes": "",
            "severe_hazard_codes": "",
            "high_concern_hazard_text_flag": False,
            "pubchem_safety_annotation_text_chars": 0,
        }

    try:
        payload = r.json()
    except Exception:
        payload = {}

    text = " ".join(flatten_json_strings(payload))
    codes = sorted({m.upper() for m in HAZARD_CODE_RE.findall(text)})
    severe = sorted(set(codes) & SEVERE_HAZARD_CODES)

    return {
        "pubchem_cid": cid,
        "pubchem_safety_annotation_found": bool(text.strip()),
        "hazard_statement_codes": " | ".join(codes),
        "severe_hazard_codes": " | ".join(severe),
        "high_concern_hazard_text_flag": bool(HIGH_CONCERN_TEXT.search(text)),
        "pubchem_safety_annotation_text_chars": len(text),
    }


def append_hazard_cache(record: dict[str, Any]) -> None:
    with HAZARD_CACHE.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")
        f.flush()


def ensure_hazards(df: pd.DataFrame, workers: int) -> pd.DataFrame:
    cache = load_hazard_cache()
    cids = sorted(set(df["pubchem_cid"].dropna().astype(str)))
    missing = [cid for cid in cids if cid not in cache]

    if missing:
        print(f"Hazard enrichment: {len(missing):,} uncached CIDs", flush=True)

        with ThreadPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(fetch_hazard, cid): cid for cid in missing}
            completed = 0

            for fut in as_completed(futures):
                rec = fut.result()
                cid = rec["pubchem_cid"]
                cache[cid] = rec
                append_hazard_cache(rec)

                completed += 1
                if completed == 1 or completed % 25 == 0 or completed == len(missing):
                    print(f"  hazard records {completed}/{len(missing)}", flush=True)

    if not cache:
        return df

    hazard_df = pd.DataFrame(list(cache.values())).drop_duplicates(
        "pubchem_cid", keep="last"
    )

    return df.merge(
        hazard_df,
        on="pubchem_cid",
        how="left",
        suffixes=("", "_hazard"),
    )


# ---------------------------------------------------------------------
# Optional RDKit PAINS
# ---------------------------------------------------------------------

def run_pains(df: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    try:
        from rdkit import Chem
        from rdkit.Chem.FilterCatalog import (
            FilterCatalog,
            FilterCatalogParams,
        )
    except Exception:
        df["pains_evaluated"] = False
        df["pains_flag"] = pd.NA
        df["pains_filter_matches"] = "RDKIT_NOT_AVAILABLE"
        return df, False

    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    smiles_col = None
    for c in ["isomeric_smiles", "canonical_smiles"]:
        if c in df.columns:
            smiles_col = c
            break

    if smiles_col is None:
        df["pains_evaluated"] = False
        df["pains_flag"] = pd.NA
        df["pains_filter_matches"] = "NO_SMILES"
        return df, True

    evaluated = []
    flags = []
    matches_text = []

    for value in df[smiles_col].tolist():
        smi = clean_text(value)
        if not smi:
            evaluated.append(False)
            flags.append(pd.NA)
            matches_text.append("NO_SMILES")
            continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            evaluated.append(False)
            flags.append(pd.NA)
            matches_text.append("SMILES_PARSE_FAILED")
            continue

        matches = catalog.GetMatches(mol)
        labels = sorted({
            clean_text(m.GetDescription())
            for m in matches
            if clean_text(m.GetDescription())
        })

        evaluated.append(True)
        flags.append(bool(labels))
        matches_text.append(" | ".join(labels))

    df["pains_evaluated"] = evaluated
    df["pains_flag"] = flags
    df["pains_filter_matches"] = matches_text

    return df, True


# ---------------------------------------------------------------------
# Safety scoring
# ---------------------------------------------------------------------

def score_row(row: pd.Series) -> pd.Series:
    """
    Conservative screening interpretation.

    Key design change from v1:
    - Severe GHS hazard annotations can no longer receive PASS.
    - PASS requires zero identified risk signals AND at least two evaluable
      evidence layers.
    - Missing evidence is treated as uncertainty, not safety.
    - 04M research/tool classification is NOT counted as a toxicity signal.
    """
    risk = 0
    reasons: list[str] = []
    evidence_layers = 0
    force_high_risk = False
    hard_safety_flag = False

    # ---------------------------------------------------------------
    # 1) Hazard evidence
    # ---------------------------------------------------------------
    severe_codes = clean_text(row.get("severe_hazard_codes"))
    high_text_raw = row.get("high_concern_hazard_text_flag", False)
    high_text = (
        bool(high_text_raw)
        if not pd.isna(high_text_raw)
        else False
    )

    severe_set = (
        {x.strip().upper() for x in severe_codes.split(" | ") if x.strip()}
        if severe_codes
        else set()
    )

    fatal_acute = severe_set & {"H300", "H310", "H330"}
    toxic_acute = severe_set & {"H301", "H311", "H331"}
    genotox_carc_repro = severe_set & {"H340", "H350", "H360"}
    organ_toxicity = severe_set & {"H370", "H372"}

    if severe_set:
        evidence_layers += 1

        if fatal_acute:
            risk += 6
            hard_safety_flag = True
            force_high_risk = True
            reasons.append(
                "fatal-acute-toxicity GHS annotation "
                + ",".join(sorted(fatal_acute))
            )

        if genotox_carc_repro:
            risk += 6
            hard_safety_flag = True
            force_high_risk = True
            reasons.append(
                "genotoxic/carcinogenic/reproductive-hazard annotation "
                + ",".join(sorted(genotox_carc_repro))
            )

        if toxic_acute:
            risk += 3
            hard_safety_flag = True
            reasons.append(
                "acute-toxicity GHS annotation "
                + ",".join(sorted(toxic_acute))
            )

        if organ_toxicity:
            risk += 4
            hard_safety_flag = True
            reasons.append(
                "specific-organ-toxicity GHS annotation "
                + ",".join(sorted(organ_toxicity))
            )

    elif high_text:
        risk += 3
        evidence_layers += 1
        hard_safety_flag = True
        reasons.append("high-concern PubChem safety annotation text")

    elif (
        bool(row.get("pubchem_safety_annotation_found", False))
        if not pd.isna(row.get("pubchem_safety_annotation_found", False))
        else False
    ):
        # This only means the layer was evaluable; it is not evidence of safety.
        evidence_layers += 1

    # ---------------------------------------------------------------
    # 2) Cytotoxicity / viability assay evidence
    # ---------------------------------------------------------------
    cyto_active = pd.to_numeric(
        pd.Series([row.get("cytotoxicity_keyword_active_n")]),
        errors="coerce",
    ).iloc[0]

    if pd.notna(cyto_active):
        evidence_layers += 1

        if cyto_active >= 5:
            risk += 4
            reasons.append("multiple active cytotoxicity/viability assays")
        elif cyto_active >= 2:
            risk += 3
            reasons.append("repeated active cytotoxicity/viability assays")
        elif cyto_active >= 1:
            risk += 1
            reasons.append("active cytotoxicity/viability assay signal")

    # ---------------------------------------------------------------
    # 3) Assay promiscuity
    # ---------------------------------------------------------------
    frac = pd.to_numeric(
        pd.Series([row.get("pubchem_active_fraction")]),
        errors="coerce",
    ).iloc[0]
    interp = pd.to_numeric(
        pd.Series([row.get("pubchem_interpretable_assay_n")]),
        errors="coerce",
    ).iloc[0]

    if pd.notna(frac) and pd.notna(interp) and interp >= 20:
        evidence_layers += 1

        if frac >= 0.25:
            risk += 3
            reasons.append("high PubChem assay promiscuity")
        elif frac >= 0.15:
            risk += 2
            reasons.append("moderate PubChem assay promiscuity")
        elif frac >= 0.08:
            risk += 1
            reasons.append("elevated PubChem assay promiscuity")

    # ---------------------------------------------------------------
    # 4) PAINS — assay interference, NOT toxicity
    # ---------------------------------------------------------------
    pains_eval = row.get("pains_evaluated", False)
    pains_flag = row.get("pains_flag", pd.NA)

    if bool(pains_eval):
        evidence_layers += 1

        if not pd.isna(pains_flag) and bool(pains_flag):
            risk += 1
            reasons.append("PAINS structural alert (assay-interference risk)")

    # ---------------------------------------------------------------
    # 5) Translational classification (reported separately)
    # ---------------------------------------------------------------
    classification = clean_text(
        row.get("compound_classification")
    ).lower()

    nonclinical_classification = any(
        k in classification
        for k in ["research", "tool", "probe", "unresolved"]
    )

    # Do not add this to toxicity risk. It is a translational-development flag.
    if nonclinical_classification:
        reasons.append(
            "research/tool/unresolved classification "
            "(translational caution; not a toxicity signal)"
        )

    # ---------------------------------------------------------------
    # Recommendation
    # ---------------------------------------------------------------
    #
    # PASS is deliberately strict: zero detected risk AND at least two
    # independent/evaluable safety layers. Missing data cannot create a PASS.
    if force_high_risk:
        recommendation = "HIGH_RISK_DEPRIORITIZE"
    elif risk >= 6:
        recommendation = "HIGH_RISK_DEPRIORITIZE"
    elif risk > 0:
        recommendation = "CAUTION_MANUAL_REVIEW"
    elif evidence_layers >= 2:
        recommendation = "PASS_PRELIMINARY_SCREEN"
    else:
        recommendation = "INSUFFICIENT_SAFETY_DATA"

    completeness = (
        "GOOD"
        if evidence_layers >= 3
        else "LIMITED"
        if evidence_layers >= 1
        else "NONE"
    )

    return pd.Series({
        "safety_risk_score": int(risk),
        "safety_evidence_layers_n": int(evidence_layers),
        "safety_data_completeness": completeness,
        "hard_safety_flag": bool(hard_safety_flag),
        "translational_nonclinical_flag": bool(
            nonclinical_classification
        ),
        "safety_screening_recommendation": recommendation,
        "safety_screening_reasons": " | ".join(reasons),
    })


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument(
        "--mode",
        choices=["full", "priority"],
        default="priority",
        help="priority = Tier 1 + top-N Tier 2; full = all 04M candidates",
    )
    p.add_argument(
        "--tier2-top",
        type=int,
        default=100,
        help="Tier 2 candidates included in priority mode",
    )
    p.add_argument(
        "--hazard-workers",
        type=int,
        default=6,
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=20,
        help="PubChem CIDs per PUG REST batch request",
    )
    p.add_argument(
        "--skip-hazards",
        action="store_true",
        help="skip PUG-View safety/hazard enrichment",
    )
    p.add_argument(
        "--skip-assays",
        action="store_true",
        help="skip PubChem assay-summary enrichment",
    )

    args = p.parse_args()
    args.tier2_top = max(0, args.tier2_top)
    args.hazard_workers = max(1, min(args.hazard_workers, 12))
    args.batch_size = max(1, min(args.batch_size, 50))

    return args


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    header("ATLAS — Stage 04O Safety / Promiscuity / PAINS Screening")

    if not INPUT_04M.exists():
        print(f"ERROR: 04M input not found:\n{INPUT_04M}")
        return 1

    df = pd.read_csv(INPUT_04M)

    required = [
        "priority_rank",
        "priority_tier_number",
        "pert_id",
        "pert_iname",
        "pubchem_cid",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print("ERROR: missing columns:")
        for c in missing:
            print(f"  - {c}")
        return 1

    df["pubchem_cid"] = df["pubchem_cid"].map(normalize_cid)
    df.loc[df["pubchem_cid"] == "", "pubchem_cid"] = pd.NA

    df = (
        df.drop_duplicates(["pert_id", "pert_iname"])
        .sort_values("priority_rank", na_position="last")
        .reset_index(drop=True)
    )

    source_n = len(df)
    df = select_candidates(df, args.mode, args.tier2_top)

    print(f"\n04M candidates available: {source_n:,}", flush=True)
    print(f"04O mode: {args.mode}", flush=True)
    print(f"Candidates selected: {len(df):,}", flush=True)
    print(
        f"Candidates with PubChem CID: {df['pubchem_cid'].notna().sum():,}",
        flush=True,
    )

    df, merged_04n = merge_optional_04n(df)
    print(
        f"04N regulatory evidence merged: {'YES' if merged_04n else 'NO — proceeding independently'}",
        flush=True,
    )

    # PubChem structure data is useful for PAINS and downstream stages.
    df = ensure_structures(df, args.batch_size)

    if not args.skip_assays:
        df = ensure_assays(df, args.batch_size)
    else:
        print("PubChem assay enrichment skipped.", flush=True)

    if not args.skip_hazards:
        df = ensure_hazards(df, args.hazard_workers)
    else:
        print("PubChem hazard enrichment skipped.", flush=True)

    df, rdkit_available = run_pains(df)

    print(
        f"RDKit PAINS screening: {'AVAILABLE' if rdkit_available else 'NOT AVAILABLE'}",
        flush=True,
    )

    scoring = df.apply(score_row, axis=1)
    df = pd.concat([df, scoring], axis=1)

    # Strongest candidates first; within same rank, safer first.
    df = df.sort_values(
        ["safety_risk_score", "priority_rank"],
        ascending=[True, True],
        na_position="last",
    ).reset_index(drop=True)

    prioritized = df[
        df["safety_screening_recommendation"].isin([
            "PASS_PRELIMINARY_SCREEN",
            "CAUTION_MANUAL_REVIEW",
        ])
    ].copy()

    manual = df[
        df["safety_screening_recommendation"].isin([
            "CAUTION_MANUAL_REVIEW",
            "HIGH_RISK_DEPRIORITIZE",
            "INSUFFICIENT_SAFETY_DATA",
        ])
    ].copy()

    summary = (
        df.groupby(
            "safety_screening_recommendation",
            dropna=False,
        )
        .agg(
            compound_count=("pert_iname", "size"),
            tier1_count=(
                "priority_tier_number",
                lambda x: int(
                    (pd.to_numeric(x, errors="coerce") == 1).sum()
                ),
            ),
            tier2_count=(
                "priority_tier_number",
                lambda x: int(
                    (pd.to_numeric(x, errors="coerce") == 2).sum()
                ),
            ),
            median_safety_risk_score=("safety_risk_score", "median"),
        )
        .reset_index()
        .sort_values("compound_count", ascending=False)
    )

    atomic_csv(df, OUT_ALL)
    atomic_csv(prioritized, OUT_PRIORITIZED)
    atomic_csv(manual, OUT_MANUAL)
    atomic_csv(summary, OUT_SUMMARY)

    metadata = {
        "stage": "04O",
        "implementation": "conservative_safety_interpretation_v2",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_04m": str(INPUT_04M),
        "optional_04n": str(OPTIONAL_04N),
        "04n_merged": merged_04n,
        "mode": args.mode,
        "tier2_top": args.tier2_top,
        "candidate_count": int(len(df)),
        "pubchem_cid_count": int(df["pubchem_cid"].notna().sum()),
        "rdkit_pains_available": rdkit_available,
        "sources": [
            "PubChem PUG REST assay summary",
            "PubChem PUG REST compound properties",
            "PubChem PUG-View Safety and Hazards",
            "RDKit PAINS FilterCatalog when available",
        ],
        "guardrails": [
            "Active PubChem assays are screening evidence, not clinical toxicity proof.",
            "PAINS flags assay-interference liability, not toxicity.",
            "Absence of PubChem hazard annotations does not establish safety.",
        ],
        "next_stage": "04P_drug_target_annotation",
    }
    atomic_json(metadata, OUT_META)

    header("STAGE 04O SUMMARY")

    print(summary.to_string(index=False), flush=True)

    tier1 = df[pd.to_numeric(
        df["priority_tier_number"], errors="coerce"
    ) == 1].copy()

    if not tier1.empty:
        header("TIER 1 SAFETY SCREEN")
        cols = [
            c for c in [
                "priority_rank",
                "pert_iname",
                "pubchem_cid",
                "pubchem_active_fraction",
                "cytotoxicity_keyword_active_n",
                "severe_hazard_codes",
                "pains_flag",
                "safety_risk_score",
                "safety_evidence_layers_n",
                "safety_data_completeness",
                "hard_safety_flag",
                "safety_screening_recommendation",
                "safety_screening_reasons",
            ]
            if c in tier1.columns
        ]
        print(tier1[cols].to_string(index=False), flush=True)

    header("STAGE 04O COMPLETE")
    print("\nOutputs:", flush=True)
    print(f"  {OUT_ALL}", flush=True)
    print(f"  {OUT_PRIORITIZED}", flush=True)
    print(f"  {OUT_MANUAL}", flush=True)
    print(f"  {OUT_SUMMARY}", flush=True)
    print(f"  {OUT_META}", flush=True)
    print("\nNext: 04P — Drug-target annotation", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
