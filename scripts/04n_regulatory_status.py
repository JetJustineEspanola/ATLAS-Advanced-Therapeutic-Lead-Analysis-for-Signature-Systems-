#!/usr/bin/env python3
"""
ATLAS — Stage 04N: Optimized Regulatory / Clinical Status Verification

Input
-----
results/cmap/drug_filter/ATLAS_CMap_drug_candidates.csv

Final outputs
-------------
results/cmap/regulatory_status/
    ATLAS_CMap_regulatory_annotations.csv
    ATLAS_CMap_regulatory_review_ready.csv
    ATLAS_CMap_regulatory_manual_review.csv
    ATLAS_CMap_regulatory_summary.csv
    ATLAS_CMap_regulatory_metadata.json

Crash-safe / resume files
-------------------------
    ATLAS_CMap_regulatory_checkpoint.jsonl
    ATLAS_CMap_regulatory_partial.csv
    ATLAS_CMap_regulatory_http_cache.jsonl

Why this version is faster
--------------------------
1. Concurrent I/O with a bounded ThreadPoolExecutor.
2. Shared rate limiting so concurrency does not hammer openFDA.
3. One combined OR query per openFDA endpoint instead of 3 sequential
   field-by-field requests in the common case.
4. limit=1 for openFDA and smaller ClinicalTrials.gov payloads.
5. Short connect/read timeouts and bounded retries.
6. Persistent endpoint-result cache.
7. Append-only candidate checkpoint after every completed candidate.
8. Resume support after Ctrl+C, terminal closure, VS Code crash, etc.
9. Optional priority mode for Tier 1 + top-N Tier 2 candidates.

Scientific guardrail
--------------------
API records are evidence only. The script does NOT automatically claim that a
compound is FDA-approved for HER2-positive breast cancer, trastuzumab
resistance, or any particular indication.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import pandas as pd
import requests
from requests.adapters import HTTPAdapter


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "drug_filter"
    / "ATLAS_CMap_drug_candidates.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "regulatory_status"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

ANNOTATED_FILE = OUTPUT_DIR / "ATLAS_CMap_regulatory_annotations.csv"
REVIEW_READY_FILE = OUTPUT_DIR / "ATLAS_CMap_regulatory_review_ready.csv"
MANUAL_REVIEW_FILE = OUTPUT_DIR / "ATLAS_CMap_regulatory_manual_review.csv"
SUMMARY_FILE = OUTPUT_DIR / "ATLAS_CMap_regulatory_summary.csv"
METADATA_FILE = OUTPUT_DIR / "ATLAS_CMap_regulatory_metadata.json"

CHECKPOINT_FILE = OUTPUT_DIR / "ATLAS_CMap_regulatory_checkpoint.jsonl"
PARTIAL_FILE = OUTPUT_DIR / "ATLAS_CMap_regulatory_partial.csv"
HTTP_CACHE_FILE = OUTPUT_DIR / "ATLAS_CMap_regulatory_http_cache.jsonl"


# ---------------------------------------------------------------------------
# APIs
# ---------------------------------------------------------------------------

FDA_DRUGSFDA_URL = "https://api.fda.gov/drug/drugsfda.json"
FDA_LABEL_URL = "https://api.fda.gov/drug/label.json"
CTGOV_URL = "https://clinicaltrials.gov/api/v2/studies"

USER_AGENT = "ATLAS-Regulatory-Status/2.0 (research workflow)"

REQUIRED_COLUMNS = [
    "priority_rank",
    "priority_tier_number",
    "priority_tier",
    "pert_id",
    "pert_iname",
    "pubchem_cid",
    "pubchem_title",
    "n_negative",
    "n_strong_negative",
    "mean_tau",
    "median_tau",
    "minimum_tau",
    "compound_classification",
    "downstream_status",
]


# ---------------------------------------------------------------------------
# Global thread-safe state
# ---------------------------------------------------------------------------

_thread_local = threading.local()
_stop_event = threading.Event()

_checkpoint_lock = threading.Lock()
_cache_lock = threading.Lock()

_http_cache: dict[str, dict[str, Any]] = {}


class RateLimiter:
    """Simple global spacing limiter shared by all worker threads."""

    def __init__(self, requests_per_second: float):
        self.interval = 1.0 / max(requests_per_second, 0.01)
        self.lock = threading.Lock()
        self.next_allowed = 0.0

    def wait(self) -> None:
        with self.lock:
            now = time.monotonic()
            delay = max(0.0, self.next_allowed - now)
            if delay > 0:
                time.sleep(delay)
                now = time.monotonic()
            self.next_allowed = max(now, self.next_allowed) + self.interval


# Conservative no-key pacing. This is deliberately below an aggressive burst.
OPENFDA_LIMITER = RateLimiter(3.5)
CTGOV_LIMITER = RateLimiter(5.0)


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------

def header(text: str) -> None:
    print("\n" + "=" * 72, flush=True)
    print(text, flush=True)
    print("=" * 72, flush=True)


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return re.sub(r"\s+", " ", str(value).strip())


def normalize_name(name: str) -> str:
    return clean_text(name).casefold()


def escape_openfda_phrase(name: str) -> str:
    # Escape the two characters most likely to break a quoted query.
    return clean_text(name).replace("\\", "\\\\").replace('"', '\\"')


def get_session() -> requests.Session:
    """One Session per worker thread; requests.Session is not shared."""
    if not hasattr(_thread_local, "session"):
        session = requests.Session()
        session.headers.update({
            "User-Agent": USER_AGENT,
            "Accept": "application/json",
        })
        adapter = HTTPAdapter(
            pool_connections=4,
            pool_maxsize=4,
            max_retries=0,  # retries are handled explicitly below
        )
        session.mount("https://", adapter)
        session.mount("http://", adapter)
        _thread_local.session = session
    return _thread_local.session


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, index=False, float_format="%.6f")
    os.replace(tmp, path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    os.replace(tmp, path)


def candidate_key(pert_id: Any, pert_iname: Any) -> str:
    return f"{clean_text(pert_id)}\t{clean_text(pert_iname)}"


def candidate_search_names(row: pd.Series | dict[str, Any]) -> list[str]:
    names = [
        clean_text(row.get("pert_iname")),
        clean_text(row.get("pubchem_title")),
    ]
    seen: set[str] = set()
    out: list[str] = []

    for name in names:
        key = normalize_name(name)
        if name and key not in seen:
            seen.add(key)
            out.append(name)

    return out[:2]


# ---------------------------------------------------------------------------
# Retry + request layer
# ---------------------------------------------------------------------------

def safe_get(
    url: str,
    params: dict[str, Any],
    *,
    limiter: RateLimiter,
    timeout_connect: float,
    timeout_read: float,
    retries: int,
) -> tuple[int | None, dict[str, Any] | None]:
    """
    Fast bounded retry policy.

    - 404 / 400: do not retry
    - 429 / 5xx: short retry with Retry-After support
    - connection/timeouts: retry
    - Ctrl+C / global stop: abort quickly
    """
    session = get_session()
    attempts = max(1, retries + 1)

    for attempt in range(attempts):
        if _stop_event.is_set():
            return None, None

        limiter.wait()

        try:
            response = session.get(
                url,
                params=params,
                timeout=(timeout_connect, timeout_read),
            )
        except requests.RequestException:
            if attempt + 1 >= attempts:
                return None, None
            if _stop_event.wait(min(2.0, 0.5 * (2 ** attempt))):
                return None, None
            continue

        status = response.status_code

        if status == 200:
            try:
                return 200, response.json()
            except ValueError:
                return 200, None

        # Expected "no match" / malformed search: no expensive retry.
        if status in {400, 404}:
            return status, None

        if status in {429, 500, 502, 503, 504}:
            if attempt + 1 >= attempts:
                return status, None

            retry_after = response.headers.get("Retry-After", "")
            try:
                delay = float(retry_after)
            except (TypeError, ValueError):
                delay = 0.75 * (2 ** attempt)

            # Bound the wait so one bad API response cannot stall the run.
            delay = min(max(delay, 0.5), 8.0)
            if _stop_event.wait(delay):
                return status, None
            continue

        return status, None

    return None, None


# ---------------------------------------------------------------------------
# Persistent HTTP result cache
# ---------------------------------------------------------------------------

def load_http_cache() -> None:
    global _http_cache

    if not HTTP_CACHE_FILE.exists():
        return

    loaded = 0
    try:
        with HTTP_CACHE_FILE.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                    key = clean_text(obj.get("cache_key"))
                    result = obj.get("result")
                    if key and isinstance(result, dict):
                        _http_cache[key] = result
                        loaded += 1
                except Exception:
                    continue
    except Exception:
        return

    if loaded:
        print(f"HTTP cache loaded: {loaded:,} entries", flush=True)


def append_http_cache(cache_key: str, result: dict[str, Any]) -> None:
    obj = {
        "cache_key": cache_key,
        "saved_utc": datetime.now(timezone.utc).isoformat(),
        "result": result,
    }
    line = json.dumps(obj, ensure_ascii=False)

    with _cache_lock:
        _http_cache[cache_key] = copy.deepcopy(result)
        with HTTP_CACHE_FILE.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")
            handle.flush()


def cached_query(
    endpoint_name: str,
    name: str,
    query_fn: Callable[[str], dict[str, Any]],
) -> dict[str, Any]:
    cache_key = f"{endpoint_name}|{normalize_name(name)}"

    with _cache_lock:
        cached = _http_cache.get(cache_key)
        if cached is not None:
            return copy.deepcopy(cached)

    result = query_fn(name)
    append_http_cache(cache_key, result)
    return result


# ---------------------------------------------------------------------------
# API-specific result templates
# ---------------------------------------------------------------------------

def empty_drugsfda() -> dict[str, Any]:
    return {
        "fda_application_record_found": False,
        "fda_application_number": "",
        "fda_sponsor_name": "",
        "fda_product_number": "",
        "fda_marketing_status": "",
        "fda_dosage_form": "",
        "fda_route": "",
        "fda_active_ingredient": "",
        "fda_drugsfda_match_name": "",
        "fda_drugsfda_http_status": "",
    }


def empty_fda_label() -> dict[str, Any]:
    return {
        "fda_label_record_found": False,
        "fda_label_generic_name": "",
        "fda_label_brand_name": "",
        "fda_label_manufacturer": "",
        "fda_label_product_type": "",
        "fda_label_route": "",
        "fda_label_indications_excerpt": "",
        "fda_label_match_name": "",
        "fda_label_http_status": "",
    }


def empty_trials() -> dict[str, Any]:
    return {
        "clinical_trial_record_found": False,
        "clinical_trial_count_returned": 0,
        "clinical_trial_nct_ids": "",
        "clinical_trial_statuses": "",
        "clinical_trial_phases": "",
        "clinical_trial_conditions": "",
        "clinical_trial_match_name": "",
        "clinicaltrials_http_status": "",
    }


# ---------------------------------------------------------------------------
# Optimized endpoint queries
# ---------------------------------------------------------------------------

def query_drugsfda(
    name: str,
    timeout_connect: float,
    timeout_read: float,
    retries: int,
) -> dict[str, Any]:
    
    out = empty_drugsfda()
    if not name:
        return out

    q = escape_openfda_phrase(name)

    # One OR query instead of three serial searches.
    expression = (
        f'products.active_ingredients.name:"{q}" OR '
        f'openfda.generic_name:"{q}" OR '
        f'openfda.substance_name:"{q}"'
    )

    params: dict[str, Any] = {
        "search": expression,
        "limit": 1,
    }

    api_key = os.getenv("OPENFDA_API_KEY", "").strip()
    if api_key:
        params["api_key"] = api_key

    status, payload = safe_get(
        FDA_DRUGSFDA_URL,
        params,
        limiter=OPENFDA_LIMITER,
        timeout_connect=timeout_connect,
        timeout_read=timeout_read,
        retries=retries,
    )
    out["fda_drugsfda_http_status"] = (
        "" if status is None else str(status)
    )

    # Some openFDA parsers can reject a broad OR expression for a name.
    # A single exact active-ingredient fallback is still much cheaper than
    # the original 3-field serial loop.
    if status == 400:
        fallback = {
            "search": f'products.active_ingredients.name:"{q}"',
            "limit": 1,
        }
        if api_key:
            fallback["api_key"] = api_key

        status, payload = safe_get(
            FDA_DRUGSFDA_URL,
            fallback,
            limiter=OPENFDA_LIMITER,
            timeout_connect=timeout_connect,
            timeout_read=timeout_read,
            retries=retries,
        )
        out["fda_drugsfda_http_status"] = (
            "" if status is None else str(status)
        )

    if not payload:
        return out

    records = payload.get("results", []) or []
    if not records:
        return out

    rec = records[0]
    products = rec.get("products", []) or []
    product = products[0] if products else {}
    ingredients = product.get("active_ingredients", []) or []

    ingredient_names = [
        clean_text(x.get("name"))
        for x in ingredients
        if isinstance(x, dict)
    ]

    out.update({
        "fda_application_record_found": True,
        "fda_application_number": clean_text(
            rec.get("application_number")
        ),
        "fda_sponsor_name": clean_text(rec.get("sponsor_name")),
        "fda_product_number": clean_text(
            product.get("product_number")
        ),
        "fda_marketing_status": clean_text(
            product.get("marketing_status")
        ),
        "fda_dosage_form": clean_text(product.get("dosage_form")),
        "fda_route": clean_text(product.get("route")),
        "fda_active_ingredient": " | ".join(
            x for x in ingredient_names if x
        ),
        "fda_drugsfda_match_name": name,
    })

    return out


def query_fda_label(
    name: str,
    timeout_connect: float,
    timeout_read: float,
    retries: int,
) -> dict[str, Any]:
    out = empty_fda_label()
    if not name:
        return out

    q = escape_openfda_phrase(name)

    # One OR query instead of three serial searches.
    expression = (
        f'openfda.generic_name:"{q}" OR '
        f'openfda.substance_name:"{q}" OR '
        f'openfda.brand_name:"{q}"'
    )

    params: dict[str, Any] = {
        "search": expression,
        "limit": 1,
    }

    api_key = os.getenv("OPENFDA_API_KEY", "").strip()
    if api_key:
        params["api_key"] = api_key

    status, payload = safe_get(
        FDA_LABEL_URL,
        params,
        limiter=OPENFDA_LIMITER,
        timeout_connect=timeout_connect,
        timeout_read=timeout_read,
        retries=retries,
    )
    out["fda_label_http_status"] = (
        "" if status is None else str(status)
    )

    if status == 400:
        fallback = {
            "search": f'openfda.generic_name:"{q}"',
            "limit": 1,
        }
        if api_key:
            fallback["api_key"] = api_key

        status, payload = safe_get(
            FDA_LABEL_URL,
            fallback,
            limiter=OPENFDA_LIMITER,
            timeout_connect=timeout_connect,
            timeout_read=timeout_read,
            retries=retries,
        )
        out["fda_label_http_status"] = (
            "" if status is None else str(status)
        )

    if not payload:
        return out

    records = payload.get("results", []) or []
    if not records:
        return out

    rec = records[0]
    openfda = rec.get("openfda", {}) or {}

    def first(field: str) -> str:
        value = openfda.get(field, [])
        if isinstance(value, list):
            return clean_text(value[0]) if value else ""
        return clean_text(value)

    indication = rec.get("indications_and_usage", [])
    if isinstance(indication, list):
        indication = " ".join(clean_text(x) for x in indication)
    indication = clean_text(indication)
    if len(indication) > 700:
        indication = indication[:697] + "..."

    out.update({
        "fda_label_record_found": True,
        "fda_label_generic_name": first("generic_name"),
        "fda_label_brand_name": first("brand_name"),
        "fda_label_manufacturer": first("manufacturer_name"),
        "fda_label_product_type": first("product_type"),
        "fda_label_route": first("route"),
        "fda_label_indications_excerpt": indication,
        "fda_label_match_name": name,
    })

    return out


def query_clinical_trials(
    name: str,
    timeout_connect: float,
    timeout_read: float,
    retries: int,
) -> dict[str, Any]:
    out = empty_trials()
    if not name:
        return out

    status, payload = safe_get(
        CTGOV_URL,
        {
            "query.intr": name,
            "pageSize": 10,
            "format": "json",
        },
        limiter=CTGOV_LIMITER,
        timeout_connect=timeout_connect,
        timeout_read=timeout_read,
        retries=retries,
    )

    out["clinicaltrials_http_status"] = (
        "" if status is None else str(status)
    )

    if not payload:
        return out

    studies = payload.get("studies", []) or []
    if not studies:
        return out

    nct_ids: list[str] = []
    statuses: list[str] = []
    phases: list[str] = []
    conditions: list[str] = []

    for study in studies:
        protocol = study.get("protocolSection", {}) or {}
        ident = protocol.get("identificationModule", {}) or {}
        status_mod = protocol.get("statusModule", {}) or {}
        design = protocol.get("designModule", {}) or {}
        cond_mod = protocol.get("conditionsModule", {}) or {}

        nct = clean_text(ident.get("nctId"))
        overall = clean_text(status_mod.get("overallStatus"))

        phase_val = design.get("phases", []) or []
        if isinstance(phase_val, list):
            phase_val = ", ".join(
                clean_text(x)
                for x in phase_val
                if clean_text(x)
            )
        else:
            phase_val = clean_text(phase_val)

        cond_val = cond_mod.get("conditions", []) or []
        if isinstance(cond_val, list):
            cond_val = ", ".join(
                clean_text(x)
                for x in cond_val[:5]
                if clean_text(x)
            )
        else:
            cond_val = clean_text(cond_val)

        if nct:
            nct_ids.append(nct)
        if overall:
            statuses.append(overall)
        if phase_val:
            phases.append(phase_val)
        if cond_val:
            conditions.append(cond_val)

    out.update({
        "clinical_trial_record_found": True,
        "clinical_trial_count_returned": len(studies),
        "clinical_trial_nct_ids": " | ".join(
            sorted(set(nct_ids))
        ),
        "clinical_trial_statuses": " | ".join(
            sorted(set(statuses))
        ),
        "clinical_trial_phases": " | ".join(
            sorted(set(phases))
        ),
        "clinical_trial_conditions": " | ".join(
            sorted(set(conditions))[:20]
        ),
        "clinical_trial_match_name": name,
    })

    return out


# ---------------------------------------------------------------------------
# Evidence classification
# ---------------------------------------------------------------------------

def evidence_category(row: dict[str, Any] | pd.Series) -> str:
    application = bool(
        row.get("fda_application_record_found", False)
    )
    label = bool(row.get("fda_label_record_found", False))
    trial = bool(row.get("clinical_trial_record_found", False))

    if application and label:
        return "FDA_APPLICATION_AND_LABEL_EVIDENCE"
    if application:
        return "FDA_APPLICATION_RECORD_FOUND"
    if label:
        return "FDA_LABEL_EVIDENCE_FOUND"
    if trial:
        return "CLINICAL_TRIAL_EVIDENCE_ONLY"
    return "NO_US_REGULATORY_OR_TRIAL_EVIDENCE_FOUND"


def clinical_evidence_level(
    row: dict[str, Any] | pd.Series,
) -> str:
    application = bool(
        row.get("fda_application_record_found", False)
    )
    label = bool(row.get("fda_label_record_found", False))
    trial = bool(row.get("clinical_trial_record_found", False))

    if application and label:
        return "HIGH_REGULATORY_EVIDENCE"
    if application or label:
        return "MODERATE_REGULATORY_EVIDENCE"
    if trial:
        return "CLINICAL_INVESTIGATION_EVIDENCE"
    return "UNRESOLVED"


# ---------------------------------------------------------------------------
# Candidate processing
# ---------------------------------------------------------------------------

def process_candidate(
    row_dict: dict[str, Any],
    *,
    verification_date: str,
    timeout_connect: float,
    timeout_read: float,
    retries: int,
) -> dict[str, Any]:
    names = candidate_search_names(row_dict)

    drugsfda = empty_drugsfda()
    label = empty_fda_label()
    trials = empty_trials()

    for name in names:
        if _stop_event.is_set():
            break

        if not drugsfda["fda_application_record_found"]:
            candidate = cached_query(
                "drugsfda",
                name,
                lambda n: query_drugsfda(
                    n,
                    timeout_connect,
                    timeout_read,
                    retries,
                ),
            )
            if candidate["fda_application_record_found"]:
                drugsfda = candidate
            else:
                # Preserve last HTTP status even when unmatched.
                drugsfda["fda_drugsfda_http_status"] = candidate.get(
                    "fda_drugsfda_http_status", ""
                )

        if not label["fda_label_record_found"]:
            candidate = cached_query(
                "fda_label",
                name,
                lambda n: query_fda_label(
                    n,
                    timeout_connect,
                    timeout_read,
                    retries,
                ),
            )
            if candidate["fda_label_record_found"]:
                label = candidate
            else:
                label["fda_label_http_status"] = candidate.get(
                    "fda_label_http_status", ""
                )

        if not trials["clinical_trial_record_found"]:
            candidate = cached_query(
                "clinical_trials",
                name,
                lambda n: query_clinical_trials(
                    n,
                    timeout_connect,
                    timeout_read,
                    retries,
                ),
            )
            if candidate["clinical_trial_record_found"]:
                trials = candidate
            else:
                trials["clinicaltrials_http_status"] = candidate.get(
                    "clinicaltrials_http_status", ""
                )

        if (
            drugsfda["fda_application_record_found"]
            and label["fda_label_record_found"]
            and trials["clinical_trial_record_found"]
        ):
            break

    rec: dict[str, Any] = {
        "pert_id": clean_text(row_dict.get("pert_id")),
        "pert_iname": clean_text(row_dict.get("pert_iname")),
        "regulatory_verification_date_utc": verification_date,
        "regulatory_jurisdiction": "United States",
        "regulatory_primary_source": "FDA/openFDA",
        "clinical_trial_source": "ClinicalTrials.gov",
        **drugsfda,
        **label,
        **trials,
    }

    rec["regulatory_evidence_category"] = evidence_category(rec)
    rec["clinical_evidence_level"] = clinical_evidence_level(rec)
    rec["manual_regulatory_review_required"] = not (
        bool(rec["fda_application_record_found"])
        and bool(rec["fda_label_record_found"])
    )

    # Never auto-claim an indication.
    rec["her2_breast_cancer_indication_verified"] = False
    rec["trastuzumab_resistance_indication_verified"] = False
    rec["final_fda_approval_claim"] = (
        "NOT_ASSIGNED_AUTOMATICALLY"
    )

    return rec


# ---------------------------------------------------------------------------
# Checkpoint / resume
# ---------------------------------------------------------------------------

def load_checkpoint() -> dict[str, dict[str, Any]]:
    completed: dict[str, dict[str, Any]] = {}

    if not CHECKPOINT_FILE.exists():
        return completed

    try:
        with CHECKPOINT_FILE.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    # Ignore a truncated final line after a hard crash.
                    continue

                key = candidate_key(
                    rec.get("pert_id"),
                    rec.get("pert_iname"),
                )
                completed[key] = rec
    except Exception as exc:
        print(
            f"WARNING: checkpoint could not be read: {exc}",
            flush=True,
        )

    return completed


def append_checkpoint(rec: dict[str, Any]) -> None:
    line = json.dumps(rec, ensure_ascii=False)
    with _checkpoint_lock:
        with CHECKPOINT_FILE.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")
            handle.flush()


# ---------------------------------------------------------------------------
# Output builders
# ---------------------------------------------------------------------------

def build_annotated(
    df: pd.DataFrame,
    records: list[dict[str, Any]],
    *,
    include_pending: bool,
) -> pd.DataFrame:
    if records:
        reg = pd.DataFrame(records)
        reg = reg.drop_duplicates(
            ["pert_id", "pert_iname"],
            keep="last",
        )
    else:
        reg = pd.DataFrame(
            columns=["pert_id", "pert_iname"]
        )

    annotated = df.merge(
        reg,
        on=["pert_id", "pert_iname"],
        how="left",
        validate="one_to_one",
    )
    annotated = annotated.sort_values(
        "priority_rank",
        na_position="last",
    ).reset_index(drop=True)

    processed_keys = {
        candidate_key(r.get("pert_id"), r.get("pert_iname"))
        for r in records
    }
    annotated["regulatory_processing_status"] = [
        (
            "PROCESSED"
            if candidate_key(r["pert_id"], r["pert_iname"])
            in processed_keys
            else "PENDING"
        )
        for _, r in annotated.iterrows()
    ]

    if not include_pending:
        annotated = annotated[
            annotated["regulatory_processing_status"] == "PROCESSED"
        ].copy()

    return annotated


def build_summary(annotated: pd.DataFrame) -> pd.DataFrame:
    processed = annotated[
        annotated["regulatory_processing_status"] == "PROCESSED"
    ].copy()

    if processed.empty:
        return pd.DataFrame(
            columns=[
                "regulatory_evidence_category",
                "compound_count",
                "tier1_count",
                "tier2_count",
                "mean_cmap_tau",
                "median_cmap_tau",
            ]
        )

    rows = []
    for category, group in processed.groupby(
        "regulatory_evidence_category",
        dropna=False,
    ):
        rows.append({
            "regulatory_evidence_category": category,
            "compound_count": len(group),
            "tier1_count": int(
                (group["priority_tier_number"] == 1).sum()
            ),
            "tier2_count": int(
                (group["priority_tier_number"] == 2).sum()
            ),
            "mean_cmap_tau": group["mean_tau"].mean(),
            "median_cmap_tau": group["median_tau"].median(),
        })

    return (
        pd.DataFrame(rows)
        .sort_values("compound_count", ascending=False)
        .reset_index(drop=True)
    )


def write_partial(
    df: pd.DataFrame,
    records: list[dict[str, Any]],
) -> None:
    partial = build_annotated(
        df,
        records,
        include_pending=True,
    )
    atomic_csv(partial, PARTIAL_FILE)


def write_final_outputs(
    df: pd.DataFrame,
    records: list[dict[str, Any]],
    *,
    verification_date: str,
    args: argparse.Namespace,
    run_complete: bool,
    elapsed_seconds: float,
) -> None:
    annotated = build_annotated(
        df,
        records,
        include_pending=True,
    )

    processed = annotated[
        annotated["regulatory_processing_status"] == "PROCESSED"
    ].copy()

    review_ready = processed[
        processed["regulatory_evidence_category"]
        != "NO_US_REGULATORY_OR_TRIAL_EVIDENCE_FOUND"
    ].copy()

    manual_review = processed[
        processed[
            "manual_regulatory_review_required"
        ].fillna(True)
    ].copy()

    summary = build_summary(annotated)

    atomic_csv(annotated, ANNOTATED_FILE)
    atomic_csv(review_ready, REVIEW_READY_FILE)
    atomic_csv(manual_review, MANUAL_REVIEW_FILE)
    atomic_csv(summary, SUMMARY_FILE)

    def true_count(column: str) -> int:
        if column not in processed.columns:
            return 0
        return int(processed[column].fillna(False).astype(bool).sum())

    metadata = {
        "stage": "04N",
        "implementation": "optimized_concurrent_resumable_v2",
        "run_complete": bool(run_complete),
        "input_file": str(INPUT_FILE),
        "selected_candidates": int(len(df)),
        "candidates_processed": int(len(processed)),
        "verification_date_utc": verification_date,
        "mode": args.mode,
        "priority_tier2_top_n": int(args.tier2_top),
        "workers": int(args.workers),
        "timeout_connect_seconds": float(args.connect_timeout),
        "timeout_read_seconds": float(args.read_timeout),
        "retries": int(args.retries),
        "elapsed_seconds": round(float(elapsed_seconds), 3),
        "sources": [
            "FDA Drugs@FDA via openFDA",
            "FDA drug labeling via openFDA",
            "ClinicalTrials.gov API v2",
        ],
        "fda_application_records_found": true_count(
            "fda_application_record_found"
        ),
        "fda_label_records_found": true_count(
            "fda_label_record_found"
        ),
        "clinical_trial_records_found": true_count(
            "clinical_trial_record_found"
        ),
        "manual_regulatory_review_required": true_count(
            "manual_regulatory_review_required"
        ),
        "checkpoint_file": str(CHECKPOINT_FILE),
        "partial_file": str(PARTIAL_FILE),
        "http_cache_file": str(HTTP_CACHE_FILE),
        "important_note": (
            "Automated records are evidence only. No compound is "
            "automatically declared FDA-approved for HER2-positive "
            "breast cancer or trastuzumab resistance."
        ),
        "next_stage": (
            "04O_cytotoxicity_promiscuity_pains_safety_screening"
        ),
    }

    atomic_json(metadata, METADATA_FILE)


# ---------------------------------------------------------------------------
# CLI / selection
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "ATLAS 04N optimized regulatory/clinical evidence verification"
        )
    )

    parser.add_argument(
        "--mode",
        choices=["full", "priority"],
        default="full",
        help=(
            "full = all 04M drug candidates; "
            "priority = Tier 1 + top-N Tier 2"
        ),
    )
    parser.add_argument(
        "--tier2-top",
        type=int,
        default=100,
        help=(
            "number of Tier 2 candidates in priority mode "
            "(default: 100)"
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="concurrent candidate workers (default: 8)",
    )
    parser.add_argument(
        "--connect-timeout",
        type=float,
        default=3.5,
        help="HTTP connect timeout in seconds (default: 3.5)",
    )
    parser.add_argument(
        "--read-timeout",
        type=float,
        default=8.0,
        help="HTTP read timeout in seconds (default: 8)",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=1,
        help="retries after the first HTTP attempt (default: 1)",
    )
    parser.add_argument(
        "--snapshot-every",
        type=int,
        default=25,
        help=(
            "rewrite partial CSV every N completed candidates "
            "(default: 25)"
        ),
    )
    parser.add_argument(
        "--reset",
        action="store_true",
        help=(
            "delete checkpoint/cache/partial files before starting; "
            "final result CSVs are left untouched until replaced"
        ),
    )

    args = parser.parse_args()

    args.workers = max(1, min(int(args.workers), 16))
    args.tier2_top = max(0, int(args.tier2_top))
    args.retries = max(0, min(int(args.retries), 3))
    args.snapshot_every = max(1, int(args.snapshot_every))
    args.connect_timeout = max(1.0, float(args.connect_timeout))
    args.read_timeout = max(2.0, float(args.read_timeout))

    return args


def reset_runtime_files() -> None:
    for path in [
        CHECKPOINT_FILE,
        PARTIAL_FILE,
        HTTP_CACHE_FILE,
    ]:
        if path.exists():
            path.unlink()
            print(f"Removed: {path}", flush=True)


def select_candidates(
    df: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    if args.mode == "full":
        return df.copy()

    tier1 = df[
        df["priority_tier_number"] == 1
    ].copy()

    tier2 = (
        df[df["priority_tier_number"] == 2]
        .sort_values(
            "priority_rank",
            na_position="last",
        )
        .head(args.tier2_top)
        .copy()
    )

    selected = pd.concat(
        [tier1, tier2],
        ignore_index=True,
    )

    selected = (
        selected.drop_duplicates(
            ["pert_id", "pert_iname"]
        )
        .sort_values(
            "priority_rank",
            na_position="last",
        )
        .reset_index(drop=True)
    )

    return selected


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    header(
        "ATLAS — Stage 04N Optimized Regulatory / Clinical Verification"
    )

    print(f"\nInput:\n{INPUT_FILE}", flush=True)

    if not INPUT_FILE.exists():
        print("\nERROR: Stage 04M output not found.", flush=True)
        print(
            "Run: python scripts/04m_cmap_drug_filter.py",
            flush=True,
        )
        return 1

    if args.reset:
        reset_runtime_files()

    df = pd.read_csv(INPUT_FILE)

    print(f"\nRows loaded: {len(df):,}", flush=True)
    print(f"Columns loaded: {len(df.columns)}", flush=True)

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        print(
            "\nERROR: Required Stage 04M columns are missing:",
            flush=True,
        )
        for c in missing:
            print(f"  - {c}", flush=True)
        return 1

    df["pert_id"] = (
        df["pert_id"].astype("string").str.strip()
    )
    df["pert_iname"] = (
        df["pert_iname"].astype("string").str.strip()
    )
    df["priority_rank"] = pd.to_numeric(
        df["priority_rank"],
        errors="coerce",
    )
    df["priority_tier_number"] = pd.to_numeric(
        df["priority_tier_number"],
        errors="coerce",
    )

    df = (
        df.sort_values(
            "priority_rank",
            na_position="last",
        )
        .drop_duplicates(
            ["pert_id", "pert_iname"],
            keep="first",
        )
        .reset_index(drop=True)
    )

    df = select_candidates(df, args)

    verification_date = (
        datetime.now(timezone.utc).date().isoformat()
    )

    load_http_cache()
    completed = load_checkpoint()

    # Only keep checkpoint rows that belong to this selected run.
    selected_keys = {
        candidate_key(r["pert_id"], r["pert_iname"])
        for _, r in df.iterrows()
    }
    completed = {
        k: v
        for k, v in completed.items()
        if k in selected_keys
    }

    pending_rows: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        key = candidate_key(
            row["pert_id"],
            row["pert_iname"],
        )
        if key not in completed:
            pending_rows.append(row.to_dict())

    total = len(df)
    already_done = len(completed)
    remaining = len(pending_rows)

    print("\nRuntime configuration:", flush=True)
    print(f"  Mode: {args.mode}", flush=True)
    print(f"  Candidates selected: {total:,}", flush=True)
    print(f"  Workers: {args.workers}", flush=True)
    print(
        f"  Timeout: connect={args.connect_timeout:g}s, "
        f"read={args.read_timeout:g}s",
        flush=True,
    )
    print(f"  Retries: {args.retries}", flush=True)
    print(
        "  openFDA pacing: ~3.5 requests/sec shared",
        flush=True,
    )
    print(
        "  ClinicalTrials.gov pacing: ~5 requests/sec shared",
        flush=True,
    )

    if already_done:
        print(
            f"\nCheckpoint resumed: {already_done:,}/{total:,} "
            "already complete.",
            flush=True,
        )

    print(
        f"Remaining candidates: {remaining:,}\n",
        flush=True,
    )

    start_time = time.monotonic()

    # Write a partial snapshot immediately if resuming.
    if completed:
        write_partial(df, list(completed.values()))

    if not pending_rows:
        elapsed = time.monotonic() - start_time
        write_final_outputs(
            df,
            list(completed.values()),
            verification_date=verification_date,
            args=args,
            run_complete=True,
            elapsed_seconds=elapsed,
        )
        print("Nothing left to query; final outputs rebuilt.", flush=True)
        return 0

    executor = ThreadPoolExecutor(
        max_workers=args.workers,
        thread_name_prefix="atlas04n",
    )

    futures = {}
    for row_dict in pending_rows:
        future = executor.submit(
            process_candidate,
            row_dict,
            verification_date=verification_date,
            timeout_connect=args.connect_timeout,
            timeout_read=args.read_timeout,
            retries=args.retries,
        )
        futures[future] = row_dict

    completed_this_run = 0
    interrupted = False

    try:
        for future in as_completed(futures):
            if _stop_event.is_set():
                break

            row_dict = futures[future]

            try:
                rec = future.result()
            except Exception as exc:
                # Save a conservative error record instead of losing the row.
                rec = {
                    "pert_id": clean_text(row_dict.get("pert_id")),
                    "pert_iname": clean_text(
                        row_dict.get("pert_iname")
                    ),
                    "regulatory_verification_date_utc": (
                        verification_date
                    ),
                    "regulatory_jurisdiction": "United States",
                    "regulatory_primary_source": "FDA/openFDA",
                    "clinical_trial_source": "ClinicalTrials.gov",
                    **empty_drugsfda(),
                    **empty_fda_label(),
                    **empty_trials(),
                    "regulatory_evidence_category": (
                        "QUERY_ERROR_REQUIRES_REVIEW"
                    ),
                    "clinical_evidence_level": "UNRESOLVED",
                    "manual_regulatory_review_required": True,
                    "her2_breast_cancer_indication_verified": False,
                    "trastuzumab_resistance_indication_verified": False,
                    "final_fda_approval_claim": (
                        "NOT_ASSIGNED_AUTOMATICALLY"
                    ),
                    "regulatory_query_error": repr(exc),
                }

            key = candidate_key(
                rec.get("pert_id"),
                rec.get("pert_iname"),
            )

            completed[key] = rec
            append_checkpoint(rec)
            completed_this_run += 1

            done = len(completed)

            if (
                completed_this_run == 1
                or completed_this_run % args.snapshot_every == 0
                or done == total
            ):
                write_partial(
                    df,
                    list(completed.values()),
                )

                elapsed = max(
                    time.monotonic() - start_time,
                    0.001,
                )
                rate = completed_this_run / elapsed
                pending = total - done
                eta = (
                    pending / rate
                    if rate > 0
                    else math.nan
                )

                eta_text = (
                    f"{eta / 60:.1f} min"
                    if math.isfinite(eta)
                    else "n/a"
                )

                print(
                    f"  {done:,}/{total:,} processed "
                    f"({rate:.2f} candidates/s this run, "
                    f"ETA ~{eta_text}) — saved",
                    flush=True,
                )

    except KeyboardInterrupt:
        interrupted = True
        _stop_event.set()

        print(
            "\nCtrl+C received. Stopping new work and preserving "
            "completed results...",
            flush=True,
        )

    finally:
        if interrupted:
            for future in futures:
                future.cancel()
            executor.shutdown(
                wait=False,
                cancel_futures=True,
            )
        else:
            executor.shutdown(wait=True)

    elapsed = time.monotonic() - start_time

    # Always save what we have.
    write_partial(df, list(completed.values()))

    run_complete = len(completed) == total

    write_final_outputs(
        df,
        list(completed.values()),
        verification_date=verification_date,
        args=args,
        run_complete=run_complete,
        elapsed_seconds=elapsed,
    )

    if not run_complete:
        header("STAGE 04N PAUSED — PROGRESS SAVED")
        print(
            f"\nProcessed: {len(completed):,}/{total:,}",
            flush=True,
        )
        print(
            "Rerun the same command to resume from the checkpoint.",
            flush=True,
        )
        print(f"\nCheckpoint:\n  {CHECKPOINT_FILE}", flush=True)
        print(f"\nPartial results:\n  {PARTIAL_FILE}", flush=True)
        return 130 if interrupted else 1

    annotated = pd.read_csv(ANNOTATED_FILE)

    header("STAGE 04N REGULATORY / CLINICAL EVIDENCE SUMMARY")

    print(
        f"\nCandidates processed: {len(annotated):,}",
        flush=True,
    )

    def count_true(col: str) -> int:
        if col not in annotated:
            return 0
        return int(
            annotated[col].fillna(False).astype(bool).sum()
        )

    print(
        "FDA application records found: "
        f"{count_true('fda_application_record_found'):,}",
        flush=True,
    )
    print(
        "FDA label records found: "
        f"{count_true('fda_label_record_found'):,}",
        flush=True,
    )
    print(
        "ClinicalTrials.gov records found: "
        f"{count_true('clinical_trial_record_found'):,}",
        flush=True,
    )
    print(
        "Manual regulatory review required: "
        f"{count_true('manual_regulatory_review_required'):,}",
        flush=True,
    )

    print("\nEvidence categories:", flush=True)
    for category, count in (
        annotated["regulatory_evidence_category"]
        .value_counts(dropna=False)
        .items()
    ):
        print(f"  {category}: {count:,}", flush=True)

    tier1 = annotated[
        annotated["priority_tier_number"] == 1
    ].copy()

    header("TIER 1 REGULATORY / CLINICAL EVIDENCE")

    display_cols = [
        "priority_rank",
        "pert_iname",
        "fda_application_record_found",
        "fda_application_number",
        "fda_marketing_status",
        "fda_label_record_found",
        "clinical_trial_record_found",
        "clinical_trial_phases",
        "regulatory_evidence_category",
        "clinical_evidence_level",
        "manual_regulatory_review_required",
    ]

    if tier1.empty:
        print("\nNo Tier 1 candidates found.", flush=True)
    else:
        existing_cols = [
            c for c in display_cols if c in tier1.columns
        ]
        print(
            "\n" + tier1[existing_cols].to_string(index=False),
            flush=True,
        )

    header("STAGE 04N COMPLETE")

    print("\nOutputs:", flush=True)
    print(f"  {ANNOTATED_FILE}", flush=True)
    print(f"  {REVIEW_READY_FILE}", flush=True)
    print(f"  {MANUAL_REVIEW_FILE}", flush=True)
    print(f"  {SUMMARY_FILE}", flush=True)
    print(f"  {METADATA_FILE}", flush=True)

    print("\nCrash-safe files:", flush=True)
    print(f"  {CHECKPOINT_FILE}", flush=True)
    print(f"  {HTTP_CACHE_FILE}", flush=True)

    print(
        "\nNext: 04O — Cytotoxicity / promiscuity / "
        "PAINS / safety screening",
        flush=True,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
