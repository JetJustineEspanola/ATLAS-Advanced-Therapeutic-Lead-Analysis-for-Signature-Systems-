#!/usr/bin/env python3
"""
ATLAS — Stage 04P
Drug–Target Annotation

Purpose
-------
Attach target evidence to CMap candidates that survived/entered Stage 04O.

Input preference
----------------
1. results/cmap/safety_screening/ATLAS_CMap_safety_prioritized.csv
2. results/cmap/safety_screening/ATLAS_CMap_safety_screening.csv
3. results/cmap/drug_filter/ATLAS_CMap_drug_candidates.csv

Evidence layers
---------------
A. PubChem BioAssay target accessions / gene IDs already produced by 04O.
B. ChEMBL molecule matching + activity-derived target annotations.

Outputs
-------
results/cmap/drug_targets/
    ATLAS_CMap_drug_target_annotations.csv
    ATLAS_CMap_drug_target_pairs.csv
    ATLAS_CMap_drug_target_prioritized.csv
    ATLAS_CMap_drug_target_summary.csv
    ATLAS_CMap_drug_target_metadata.json
    cache/

Important guardrails
--------------------
- A ChEMBL activity record is biochemical/pharmacological evidence, not proof
  that the target mediates trastuzumab-resistance reversal.
- PubChem BioAssay target hits can include screening artifacts or context-
  dependent activities.
- Target evidence is used to support prioritization and network analysis,
  not to claim mechanism-of-action causality.
"""

from __future__ import annotations

import argparse
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

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_PRIORITY = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "safety_screening"
    / "ATLAS_CMap_safety_prioritized.csv"
)

INPUT_ALL_04O = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "safety_screening"
    / "ATLAS_CMap_safety_screening.csv"
)

INPUT_04M = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "drug_filter"
    / "ATLAS_CMap_drug_candidates.csv"
)

OUT_DIR = PROJECT_ROOT / "results" / "cmap" / "drug_targets"
CACHE_DIR = OUT_DIR / "cache"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

OUT_ANNOT = OUT_DIR / "ATLAS_CMap_drug_target_annotations.csv"
OUT_PAIRS = OUT_DIR / "ATLAS_CMap_drug_target_pairs.csv"
OUT_PRIORITY = OUT_DIR / "ATLAS_CMap_drug_target_prioritized.csv"
OUT_SUMMARY = OUT_DIR / "ATLAS_CMap_drug_target_summary.csv"
OUT_META = OUT_DIR / "ATLAS_CMap_drug_target_metadata.json"

MOLECULE_CACHE = CACHE_DIR / "chembl_molecule_cache.jsonl"
ACTIVITY_CACHE = CACHE_DIR / "chembl_activity_cache.jsonl"
TARGET_CACHE = CACHE_DIR / "chembl_target_cache.jsonl"


# ---------------------------------------------------------------------
# ChEMBL
# ---------------------------------------------------------------------

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
USER_AGENT = "ATLAS-DrugTarget-Annotation/1.0"

_thread_local = threading.local()


class RateLimiter:
    def __init__(self, rps: float):
        self.interval = 1.0 / max(rps, 0.01)
        self.lock = threading.Lock()
        self.next_allowed = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            delay = max(0.0, self.next_allowed - now)
            if delay > 0:
                time.sleep(delay)
                now = time.monotonic()
            self.next_allowed = max(self.next_allowed, now) + self.interval


CHEMBL_LIMITER = RateLimiter(4.0)


def get_session() -> requests.Session:
    if not hasattr(_thread_local, "session"):
        s = requests.Session()
        s.headers.update({
            "User-Agent": USER_AGENT,
            "Accept": "application/json",
        })
        adapter = HTTPAdapter(pool_connections=8, pool_maxsize=8, max_retries=0)
        s.mount("https://", adapter)
        s.mount("http://", adapter)
        _thread_local.session = s
    return _thread_local.session


def safe_get(
    url: str,
    params: dict[str, Any] | None = None,
    retries: int = 1,
    connect_timeout: float = 4.0,
    read_timeout: float = 12.0,
) -> tuple[int | None, dict[str, Any] | None]:
    session = get_session()

    for attempt in range(retries + 1):
        CHEMBL_LIMITER.wait()
        try:
            r = session.get(
                url,
                params=params or {},
                timeout=(connect_timeout, read_timeout),
            )
        except requests.RequestException:
            r = None

        if r is not None and r.status_code == 200:
            try:
                return 200, r.json()
            except ValueError:
                return 200, None

        if r is not None and r.status_code in {400, 404}:
            return r.status_code, None

        if attempt < retries:
            time.sleep(min(0.75 * (2 ** attempt), 4.0))

    return (None if r is None else r.status_code), None


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


def norm_name(x: Any) -> str:
    return clean_text(x).casefold()


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(tmp, path)


def append_jsonl(path: Path, obj: dict[str, Any]) -> None:
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(obj, ensure_ascii=False) + "\n")
        f.flush()


def load_jsonl_by_key(path: Path, key_field: str) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return out

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            try:
                obj = json.loads(line)
            except Exception:
                continue
            key = clean_text(obj.get(key_field))
            if key:
                out[key] = obj
    return out


def unique_pipe(values: list[Any], limit: int = 100) -> str:
    cleaned = sorted({
        clean_text(v)
        for v in values
        if clean_text(v)
    })
    return " | ".join(cleaned[:limit])


def choose_input() -> tuple[Path, str]:
    if INPUT_PRIORITY.exists():
        return INPUT_PRIORITY, "04O_PRIORITIZED"
    if INPUT_ALL_04O.exists():
        return INPUT_ALL_04O, "04O_ALL"
    if INPUT_04M.exists():
        return INPUT_04M, "04M_FALLBACK"
    raise FileNotFoundError(
        "No valid 04P input found. Expected 04O safety output or 04M drug candidates."
    )


# ---------------------------------------------------------------------
# PubChem target evidence from 04O
# ---------------------------------------------------------------------

def parse_pipe_set(value: Any) -> set[str]:
    text = clean_text(value)
    if not text:
        return set()
    return {x.strip() for x in text.split("|") if x.strip()}


def pubchem_target_strength(row: pd.Series) -> tuple[int, str]:
    accessions = parse_pipe_set(row.get("active_target_accessions"))
    genes = parse_pipe_set(row.get("active_target_gene_ids"))

    n = max(len(accessions), len(genes))
    if n >= 5:
        return 2, "MULTIPLE_PUBCHEM_BIOASSAY_TARGETS"
    if n >= 1:
        return 1, "PUBCHEM_BIOASSAY_TARGET_PRESENT"
    return 0, ""


# ---------------------------------------------------------------------
# ChEMBL molecule matching
# ---------------------------------------------------------------------

def candidate_names(row: pd.Series | dict[str, Any]) -> list[str]:
    names = [
        clean_text(row.get("pert_iname")),
        clean_text(row.get("pubchem_title")),
    ]
    out = []
    seen = set()
    for name in names:
        key = norm_name(name)
        if name and key not in seen:
            seen.add(key)
            out.append(name)
    return out[:2]


def score_molecule_match(query_name: str, mol: dict[str, Any]) -> float:
    q = norm_name(query_name)
    pref = norm_name(mol.get("pref_name"))

    score = 0.0
    if pref == q and q:
        score += 100.0
    elif q and pref and (q in pref or pref in q):
        score += 60.0

    syns = mol.get("molecule_synonyms", []) or []
    for syn in syns:
        s = norm_name(syn.get("molecule_synonym"))
        if s == q and q:
            score += 90.0
            break

    if mol.get("molecule_chembl_id"):
        score += 5.0

    return score


def query_chembl_molecule(name: str) -> dict[str, Any]:
    status, payload = safe_get(
        f"{CHEMBL_BASE}/molecule/search.json",
        params={"q": name, "limit": 10},
        retries=1,
    )

    result = {
        "query_name": name,
        "chembl_http_status": status,
        "chembl_molecule_found": False,
        "molecule_chembl_id": "",
        "chembl_pref_name": "",
        "chembl_molecule_type": "",
        "chembl_max_phase": np.nan,
        "chembl_match_score": 0.0,
    }

    if not payload:
        return result

    molecules = payload.get("molecules", []) or []
    if not molecules:
        return result

    ranked = sorted(
        [(score_molecule_match(name, m), m) for m in molecules],
        key=lambda x: x[0],
        reverse=True,
    )

    score, best = ranked[0]
    if score < 5:
        return result

    result.update({
        "chembl_molecule_found": True,
        "molecule_chembl_id": clean_text(best.get("molecule_chembl_id")),
        "chembl_pref_name": clean_text(best.get("pref_name")),
        "chembl_molecule_type": clean_text(best.get("molecule_type")),
        "chembl_max_phase": best.get("max_phase"),
        "chembl_match_score": float(score),
    })

    return result


# ---------------------------------------------------------------------
# ChEMBL activities + targets
# ---------------------------------------------------------------------

def query_chembl_activities(
    molecule_chembl_id: str,
    min_pchembl: float,
    max_records: int,
) -> dict[str, Any]:
    status, payload = safe_get(
        f"{CHEMBL_BASE}/activity.json",
        params={
            "molecule_chembl_id": molecule_chembl_id,
            "pchembl_value__gte": min_pchembl,
            "limit": max_records,
        },
        retries=1,
        read_timeout=20.0,
    )

    result = {
        "molecule_chembl_id": molecule_chembl_id,
        "activity_http_status": status,
        "activity_record_n": 0,
        "target_chembl_ids": [],
        "assay_types": [],
        "standard_types": [],
        "pchembl_values": [],
    }

    if not payload:
        return result

    acts = payload.get("activities", []) or []

    target_ids = []
    assay_types = []
    standard_types = []
    pchembl_values = []

    for a in acts:
        tid = clean_text(a.get("target_chembl_id"))
        if tid:
            target_ids.append(tid)

        at = clean_text(a.get("assay_type"))
        if at:
            assay_types.append(at)

        st = clean_text(a.get("standard_type"))
        if st:
            standard_types.append(st)

        try:
            pv = float(a.get("pchembl_value"))
            if math.isfinite(pv):
                pchembl_values.append(pv)
        except Exception:
            pass

    result.update({
        "activity_record_n": len(acts),
        "target_chembl_ids": sorted(set(target_ids)),
        "assay_types": sorted(set(assay_types)),
        "standard_types": sorted(set(standard_types)),
        "pchembl_values": pchembl_values,
    })

    return result


def query_target(target_chembl_id: str) -> dict[str, Any]:
    status, payload = safe_get(
        f"{CHEMBL_BASE}/target/{target_chembl_id}.json",
        retries=1,
    )

    result = {
        "target_chembl_id": target_chembl_id,
        "target_http_status": status,
        "target_pref_name": "",
        "target_type": "",
        "target_organism": "",
        "target_accessions": [],
        "target_component_descriptions": [],
    }

    if not payload:
        return result

    components = payload.get("target_components", []) or []
    accessions = []
    descriptions = []

    for c in components:
        acc = clean_text(c.get("accession"))
        if acc:
            accessions.append(acc)

        desc = clean_text(c.get("component_description"))
        if desc:
            descriptions.append(desc)

    result.update({
        "target_pref_name": clean_text(payload.get("pref_name")),
        "target_type": clean_text(payload.get("target_type")),
        "target_organism": clean_text(payload.get("organism")),
        "target_accessions": sorted(set(accessions)),
        "target_component_descriptions": sorted(set(descriptions)),
    })

    return result


# ---------------------------------------------------------------------
# Annotation workflow
# ---------------------------------------------------------------------

def annotate_one_candidate(
    row_dict: dict[str, Any],
    mol_cache: dict[str, dict[str, Any]],
    act_cache: dict[str, dict[str, Any]],
    target_cache: dict[str, dict[str, Any]],
    cache_lock: threading.Lock,
    min_pchembl: float,
    max_activity_records: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:

    row = pd.Series(row_dict)

    # 1) molecule match
    molecule_match = None

    for name in candidate_names(row):
        cache_key = norm_name(name)

        with cache_lock:
            cached = mol_cache.get(cache_key)

        if cached is None:
            cached = query_chembl_molecule(name)
            with cache_lock:
                mol_cache[cache_key] = cached
                append_jsonl(
                    MOLECULE_CACHE,
                    {"cache_key": cache_key, **cached},
                )

        if cached.get("chembl_molecule_found"):
            molecule_match = cached
            break

    if molecule_match is None:
        molecule_match = {
            "chembl_molecule_found": False,
            "molecule_chembl_id": "",
            "chembl_pref_name": "",
            "chembl_molecule_type": "",
            "chembl_max_phase": np.nan,
            "chembl_match_score": 0.0,
            "query_name": "",
        }

    molecule_chembl_id = clean_text(
        molecule_match.get("molecule_chembl_id")
    )

    # 2) activity records
    if molecule_chembl_id:
        with cache_lock:
            activities = act_cache.get(molecule_chembl_id)

        if activities is None:
            activities = query_chembl_activities(
                molecule_chembl_id,
                min_pchembl=min_pchembl,
                max_records=max_activity_records,
            )
            with cache_lock:
                act_cache[molecule_chembl_id] = activities
                append_jsonl(ACTIVITY_CACHE, activities)
    else:
        activities = {
            "activity_record_n": 0,
            "target_chembl_ids": [],
            "assay_types": [],
            "standard_types": [],
            "pchembl_values": [],
        }

    # 3) target details
    target_rows: list[dict[str, Any]] = []

    for target_chembl_id in activities.get("target_chembl_ids", [])[:50]:
        with cache_lock:
            target_info = target_cache.get(target_chembl_id)

        if target_info is None:
            target_info = query_target(target_chembl_id)
            with cache_lock:
                target_cache[target_chembl_id] = target_info
                append_jsonl(TARGET_CACHE, target_info)

        pair = {
            "pert_id": clean_text(row.get("pert_id")),
            "pert_iname": clean_text(row.get("pert_iname")),
            "priority_rank": row.get("priority_rank"),
            "priority_tier_number": row.get("priority_tier_number"),
            "molecule_chembl_id": molecule_chembl_id,
            "target_chembl_id": target_chembl_id,
            "target_pref_name": clean_text(target_info.get("target_pref_name")),
            "target_type": clean_text(target_info.get("target_type")),
            "target_organism": clean_text(target_info.get("target_organism")),
            "target_accessions": unique_pipe(
                target_info.get("target_accessions", [])
            ),
            "target_component_descriptions": unique_pipe(
                target_info.get("target_component_descriptions", [])
            ),
            "target_evidence_source": "ChEMBL_activity",
            "minimum_pchembl_filter": min_pchembl,
        }
        target_rows.append(pair)

    pvals = activities.get("pchembl_values", []) or []
    median_pchembl = float(np.median(pvals)) if pvals else np.nan
    max_pchembl = float(np.max(pvals)) if pvals else np.nan

    pubchem_strength, pubchem_reason = pubchem_target_strength(row)

    chembl_target_n = len(set(
        activities.get("target_chembl_ids", []) or []
    ))

    chembl_strength = 0
    if chembl_target_n >= 1:
        chembl_strength = 1
    if chembl_target_n >= 1 and pd.notna(max_pchembl) and max_pchembl >= 7:
        chembl_strength = 2

    total_strength = pubchem_strength + chembl_strength

    if total_strength >= 3:
        support = "STRONG_TARGET_SUPPORT"
    elif total_strength >= 2:
        support = "MODERATE_TARGET_SUPPORT"
    elif total_strength >= 1:
        support = "WEAK_TARGET_SUPPORT"
    else:
        support = "NO_TARGET_SUPPORT_FOUND"

    summary = {
        "pert_id": clean_text(row.get("pert_id")),
        "pert_iname": clean_text(row.get("pert_iname")),
        "chembl_molecule_found": bool(
            molecule_match.get("chembl_molecule_found", False)
        ),
        "molecule_chembl_id": molecule_chembl_id,
        "chembl_pref_name": clean_text(
            molecule_match.get("chembl_pref_name")
        ),
        "chembl_molecule_type": clean_text(
            molecule_match.get("chembl_molecule_type")
        ),
        "chembl_max_phase": molecule_match.get("chembl_max_phase"),
        "chembl_match_score": molecule_match.get("chembl_match_score"),
        "chembl_activity_record_n": int(
            activities.get("activity_record_n", 0)
        ),
        "chembl_target_n": chembl_target_n,
        "chembl_target_ids": unique_pipe(
            activities.get("target_chembl_ids", [])
        ),
        "chembl_assay_types": unique_pipe(
            activities.get("assay_types", [])
        ),
        "chembl_standard_types": unique_pipe(
            activities.get("standard_types", [])
        ),
        "chembl_median_pchembl": median_pchembl,
        "chembl_max_pchembl": max_pchembl,
        "pubchem_target_accession_n": len(
            parse_pipe_set(row.get("active_target_accessions"))
        ),
        "pubchem_target_gene_id_n": len(
            parse_pipe_set(row.get("active_target_gene_ids"))
        ),
        "pubchem_target_support_note": pubchem_reason,
        "target_evidence_strength_score": total_strength,
        "target_support_category": support,
    }

    return summary, target_rows


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument(
        "--max-candidates",
        type=int,
        default=50,
        help=(
            "Analyze at most N top-ranked candidates from the selected input. "
            "Use 0 for all. Default: 50"
        ),
    )
    p.add_argument(
        "--workers",
        type=int,
        default=6,
    )
    p.add_argument(
        "--min-pchembl",
        type=float,
        default=5.0,
        help="Minimum ChEMBL pChEMBL value. Default: 5.0",
    )
    p.add_argument(
        "--max-activity-records",
        type=int,
        default=500,
        help="Maximum ChEMBL activity rows fetched per drug.",
    )

    args = p.parse_args()

    args.workers = max(1, min(args.workers, 12))
    args.max_candidates = max(0, args.max_candidates)
    args.max_activity_records = max(50, min(args.max_activity_records, 1000))

    return args


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    header("ATLAS — Stage 04P Drug–Target Annotation")

    try:
        input_file, input_source = choose_input()
    except Exception as exc:
        print(f"ERROR: {exc}", flush=True)
        return 1

    print(f"\nInput source: {input_source}", flush=True)
    print(f"Input file:\n{input_file}", flush=True)

    df = pd.read_csv(input_file)

    required = [
        "pert_id",
        "pert_iname",
        "priority_rank",
        "priority_tier_number",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print("\nERROR: required columns missing:", flush=True)
        for c in missing:
            print(f"  - {c}", flush=True)
        return 1

    df["priority_rank"] = pd.to_numeric(
        df["priority_rank"],
        errors="coerce",
    )

    df = (
        df.drop_duplicates(["pert_id", "pert_iname"])
        .sort_values("priority_rank", na_position="last")
        .reset_index(drop=True)
    )

    # If safety status exists, exclude only the clearest high-risk candidates.
    if "safety_screening_recommendation" in df.columns:
        before = len(df)
        df = df[
            df["safety_screening_recommendation"]
            != "HIGH_RISK_DEPRIORITIZE"
        ].copy()
        print(
            f"\nSafety-filtered: {before:,} -> {len(df):,} candidates",
            flush=True,
        )

    if args.max_candidates > 0:
        df = df.head(args.max_candidates).copy()

    print(f"Candidates selected for target annotation: {len(df):,}", flush=True)
    print(f"Minimum pChEMBL: {args.min_pchembl:g}", flush=True)
    print(f"Workers: {args.workers}", flush=True)

    mol_cache_raw = load_jsonl_by_key(MOLECULE_CACHE, "cache_key")
    act_cache = load_jsonl_by_key(ACTIVITY_CACHE, "molecule_chembl_id")
    target_cache = load_jsonl_by_key(TARGET_CACHE, "target_chembl_id")

    mol_cache = {
        k: {kk: vv for kk, vv in v.items() if kk != "cache_key"}
        for k, v in mol_cache_raw.items()
    }

    if mol_cache or act_cache or target_cache:
        print(
            f"\nCache: molecules={len(mol_cache):,}, "
            f"activities={len(act_cache):,}, "
            f"targets={len(target_cache):,}",
            flush=True,
        )

    cache_lock = threading.Lock()

    summaries: list[dict[str, Any]] = []
    pairs: list[dict[str, Any]] = []

    start = time.monotonic()

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_map = {
            executor.submit(
                annotate_one_candidate,
                row.to_dict(),
                mol_cache,
                act_cache,
                target_cache,
                cache_lock,
                args.min_pchembl,
                args.max_activity_records,
            ): row["pert_iname"]
            for _, row in df.iterrows()
        }

        done = 0
        total = len(future_map)

        for future in as_completed(future_map):
            name = future_map[future]

            try:
                summary, target_rows = future.result()
            except Exception as exc:
                print(
                    f"  WARNING: {name}: {exc!r}",
                    flush=True,
                )
                summary = {
                    "pert_id": "",
                    "pert_iname": name,
                    "target_support_category": "QUERY_ERROR",
                    "target_evidence_strength_score": 0,
                    "chembl_target_n": 0,
                }
                target_rows = []

            summaries.append(summary)
            pairs.extend(target_rows)

            done += 1

            if done == 1 or done % 10 == 0 or done == total:
                elapsed = max(time.monotonic() - start, 0.001)
                rate = done / elapsed
                remaining = total - done
                eta_min = remaining / rate / 60 if rate > 0 else math.nan

                print(
                    f"  {done}/{total} candidates annotated "
                    f"({rate:.2f}/s, ETA ~{eta_min:.1f} min)",
                    flush=True,
                )

    summary_df = pd.DataFrame(summaries)

    annotated = df.merge(
        summary_df,
        on=["pert_id", "pert_iname"],
        how="left",
    )

    # Rank within 04P:
    # stronger target support first, then original CMap priority.
    annotated = annotated.sort_values(
        ["target_evidence_strength_score", "priority_rank"],
        ascending=[False, True],
        na_position="last",
    ).reset_index(drop=True)

    # Candidates eligible for network integration.
    prioritized = annotated[
        annotated["target_support_category"].isin([
            "STRONG_TARGET_SUPPORT",
            "MODERATE_TARGET_SUPPORT",
            "WEAK_TARGET_SUPPORT",
        ])
    ].copy()

    pair_df = pd.DataFrame(pairs)

    if pair_df.empty:
        pair_df = pd.DataFrame(columns=[
            "pert_id",
            "pert_iname",
            "priority_rank",
            "priority_tier_number",
            "molecule_chembl_id",
            "target_chembl_id",
            "target_pref_name",
            "target_type",
            "target_organism",
            "target_accessions",
            "target_component_descriptions",
            "target_evidence_source",
            "minimum_pchembl_filter",
        ])

    support_summary = (
        annotated.groupby("target_support_category", dropna=False)
        .agg(
            compound_count=("pert_iname", "size"),
            tier1_count=(
                "priority_tier_number",
                lambda x: int(
                    (pd.to_numeric(x, errors="coerce") == 1).sum()
                ),
            ),
            median_chembl_target_n=("chembl_target_n", "median"),
            median_target_strength_score=(
                "target_evidence_strength_score",
                "median",
            ),
        )
        .reset_index()
        .sort_values("compound_count", ascending=False)
    )

    atomic_csv(annotated, OUT_ANNOT)
    atomic_csv(pair_df, OUT_PAIRS)
    atomic_csv(prioritized, OUT_PRIORITY)
    atomic_csv(support_summary, OUT_SUMMARY)

    metadata = {
        "stage": "04P",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_source": input_source,
        "input_file": str(input_file),
        "candidate_count": int(len(df)),
        "minimum_pchembl": float(args.min_pchembl),
        "max_activity_records_per_drug": int(args.max_activity_records),
        "workers": int(args.workers),
        "sources": [
            "PubChem BioAssay target fields propagated from 04O when available",
            "ChEMBL molecule search",
            "ChEMBL activity records",
            "ChEMBL target annotations",
        ],
        "guardrails": [
            "Target annotations are evidence, not proof of resistance mechanism.",
            "ChEMBL activity does not establish efficacy in HER2-positive disease.",
            "Only biologically supported drug-target pairs should proceed to docking.",
        ],
        "next_stage": "04Q_PPI_and_pathway_network_integration",
    }

    atomic_json(metadata, OUT_META)

    header("STAGE 04P SUMMARY")
    print(support_summary.to_string(index=False), flush=True)

    tier1 = annotated[
        pd.to_numeric(
            annotated["priority_tier_number"],
            errors="coerce",
        ) == 1
    ].copy()

    if not tier1.empty:
        header("TIER 1 TARGET EVIDENCE")

        cols = [
            c for c in [
                "priority_rank",
                "pert_iname",
                "target_support_category",
                "target_evidence_strength_score",
                "molecule_chembl_id",
                "chembl_target_n",
                "chembl_target_ids",
                "chembl_max_pchembl",
                "pubchem_target_accession_n",
                "pubchem_target_gene_id_n",
            ]
            if c in tier1.columns
        ]

        print(tier1[cols].to_string(index=False), flush=True)

    header("STAGE 04P COMPLETE")
    print("\nOutputs:", flush=True)
    print(f"  {OUT_ANNOT}", flush=True)
    print(f"  {OUT_PAIRS}", flush=True)
    print(f"  {OUT_PRIORITY}", flush=True)
    print(f"  {OUT_SUMMARY}", flush=True)
    print(f"  {OUT_META}", flush=True)
    print("\nNext: 04Q — PPI + pathway/network integration", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
