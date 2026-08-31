#!/usr/bin/env python3
"""
ATLAS — Stage 04Q
PPI + Pathway / Network Integration

Purpose
-------
Connect drug targets from Stage 04P to the ATLAS resistance biology discovered
in Stages 02/03.

Primary inputs
--------------
results/cmap/drug_targets/ATLAS_CMap_drug_target_pairs.csv
results/cmap/drug_targets/ATLAS_CMap_drug_target_annotations.csv
results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv

Optional inputs
---------------
results/consensus_signature/ATLAS_consensus_resistance_signature.csv
results/pathway_analysis/... (if present)

Network source
--------------
STRING API

Outputs
-------
results/cmap/network_integration/
    ATLAS_drug_target_network_edges.csv
    ATLAS_target_resistance_gene_links.csv
    ATLAS_target_network_scores.csv
    ATLAS_drug_network_prioritized.csv
    ATLAS_network_summary.csv
    ATLAS_network_metadata.json
    cache/

Important guardrails
--------------------
- STRING edges are association evidence, not necessarily direct physical binding.
- Network proximity supports biological plausibility; it does not establish that
  a compound reverses trastuzumab resistance.
- Docking should only follow for drug-target pairs with convergent CMap, safety,
  target, and network/pathway support.
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

INPUT_PAIRS = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "drug_targets"
    / "ATLAS_CMap_drug_target_pairs.csv"
)

INPUT_ANNOT = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "drug_targets"
    / "ATLAS_CMap_drug_target_annotations.csv"
)

DISCOVERY_DEG = (
    PROJECT_ROOT
    / "results"
    / "differential_expression"
    / "DEGs_resistant_vs_sensitive_annotated.csv"
)

CONSENSUS_CANDIDATES = [
    PROJECT_ROOT / "results" / "consensus_signature" / "ATLAS_consensus_resistance_signature.csv",
    PROJECT_ROOT / "results" / "consensus_signature" / "consensus_resistance_signature.csv",
]

OUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "network_integration"
)
CACHE_DIR = OUT_DIR / "cache"

OUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

OUT_EDGES = OUT_DIR / "ATLAS_drug_target_network_edges.csv"
OUT_LINKS = OUT_DIR / "ATLAS_target_resistance_gene_links.csv"
OUT_TARGET_SCORES = OUT_DIR / "ATLAS_target_network_scores.csv"
OUT_DRUG_PRIORITY = OUT_DIR / "ATLAS_drug_network_prioritized.csv"
OUT_SUMMARY = OUT_DIR / "ATLAS_network_summary.csv"
OUT_META = OUT_DIR / "ATLAS_network_metadata.json"

STRING_ID_CACHE = CACHE_DIR / "string_identifier_cache.jsonl"
STRING_NETWORK_CACHE = CACHE_DIR / "string_network_cache.jsonl"


# ---------------------------------------------------------------------
# STRING
# ---------------------------------------------------------------------

STRING_BASE = "https://string-db.org/api"
STRING_SPECIES = 9606
USER_AGENT = "ATLAS-Network-Integration/1.0"

_thread_local = threading.local()


class RateLimiter:
    def __init__(self, rps: float):
        self.interval = 1.0 / max(float(rps), 0.01)
        self.lock = threading.Lock()
        self.next_allowed = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            delay = max(0.0, self.next_allowed - now)
            if delay > 0:
                time.sleep(delay)
                now = time.monotonic()
            self.next_allowed = max(now, self.next_allowed) + self.interval


STRING_LIMITER = RateLimiter(3.0)


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
    params: dict[str, Any] | None = None,
    *,
    retries: int = 1,
    connect_timeout: float = 4.0,
    read_timeout: float = 15.0,
) -> tuple[int | None, str | None]:
    session = get_session()

    for attempt in range(retries + 1):
        STRING_LIMITER.wait()

        try:
            r = session.get(
                url,
                params=params or {},
                timeout=(connect_timeout, read_timeout),
            )
        except requests.RequestException:
            r = None

        if r is not None and r.status_code == 200:
            return 200, r.text

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


def find_symbol_column(df: pd.DataFrame) -> str | None:
    candidates = [
        "Gene name",
        "gene_symbol",
        "Gene Symbol",
        "symbol",
        "SYMBOL",
    ]
    for c in candidates:
        if c in df.columns:
            return c
    return None


def split_pipe(value: Any) -> list[str]:
    text = clean_text(value)
    if not text:
        return []
    return [x.strip() for x in text.split("|") if x.strip()]


# ---------------------------------------------------------------------
# Resistance-gene set
# ---------------------------------------------------------------------

def load_resistance_genes(
    max_genes: int,
    padj_cutoff: float,
    abs_fc_cutoff: float,
) -> pd.DataFrame:
    # Prefer a finished 03C consensus signature if available.
    for path in CONSENSUS_CANDIDATES:
        if not path.exists():
            continue

        try:
            df = pd.read_csv(path)
        except Exception:
            continue

        symbol_col = find_symbol_column(df)
        if symbol_col is None:
            continue

        out = pd.DataFrame({
            "gene_symbol": df[symbol_col].map(clean_text),
        })

        # Try to carry consensus labels / effect sizes when present.
        for col in [
            "consensus_category",
            "atlas_log2FC",
            "log2FoldChange",
            "atlas_padj",
            "padj",
        ]:
            if col in df.columns:
                out[col] = df[col]

        out = out[out["gene_symbol"].ne("")].drop_duplicates("gene_symbol")

        # Prefer HIGH/MODERATE consensus if the column exists.
        if "consensus_category" in out.columns:
            selected = out[
                out["consensus_category"].astype(str).isin(
                    ["HIGH", "MODERATE"]
                )
            ].copy()
            if not selected.empty:
                out = selected

        return out.head(max_genes).reset_index(drop=True)

    # Fallback: discovery DEG set.
    if not DISCOVERY_DEG.exists():
        raise FileNotFoundError(
            f"Discovery DEG file not found: {DISCOVERY_DEG}"
        )

    df = pd.read_csv(DISCOVERY_DEG)
    symbol_col = find_symbol_column(df)
    if symbol_col is None:
        raise ValueError(
            "Could not identify gene-symbol column in discovery DEG file."
        )

    if "log2FoldChange" not in df.columns or "padj" not in df.columns:
        raise ValueError(
            "Discovery DEG file must contain log2FoldChange and padj."
        )

    out = pd.DataFrame({
        "gene_symbol": df[symbol_col].map(clean_text),
        "log2FoldChange": pd.to_numeric(
            df["log2FoldChange"],
            errors="coerce",
        ),
        "padj": pd.to_numeric(
            df["padj"],
            errors="coerce",
        ),
    })

    out = out[
        out["gene_symbol"].ne("")
        & out["log2FoldChange"].notna()
        & out["padj"].notna()
    ].copy()

    out = out[
        out["padj"].lt(padj_cutoff)
        & out["log2FoldChange"].abs().ge(abs_fc_cutoff)
    ].copy()

    out["abs_fc"] = out["log2FoldChange"].abs()
    out = (
        out.sort_values(
            ["padj", "abs_fc"],
            ascending=[True, False],
        )
        .drop_duplicates("gene_symbol")
        .head(max_genes)
        .drop(columns=["abs_fc"])
        .reset_index(drop=True)
    )

    return out


# ---------------------------------------------------------------------
# Target identifiers from 04P
# ---------------------------------------------------------------------

def collect_candidate_targets(
    pairs: pd.DataFrame,
    max_targets_per_drug: int,
) -> pd.DataFrame:
    required = [
        "pert_id",
        "pert_iname",
        "target_chembl_id",
        "target_pref_name",
        "target_accessions",
    ]
    missing = [c for c in required if c not in pairs.columns]
    if missing:
        raise ValueError(
            f"04P target-pair file missing columns: {missing}"
        )

    rows = []

    for (pert_id, pert_iname), group in pairs.groupby(
        ["pert_id", "pert_iname"],
        dropna=False,
    ):
        group = group.head(max_targets_per_drug)

        for _, r in group.iterrows():
            accessions = split_pipe(r.get("target_accessions"))
            pref_name = clean_text(r.get("target_pref_name"))

            # Prefer protein accessions because STRING mapping is more reliable.
            if accessions:
                for accession in accessions:
                    rows.append({
                        "pert_id": clean_text(pert_id),
                        "pert_iname": clean_text(pert_iname),
                        "target_chembl_id": clean_text(
                            r.get("target_chembl_id")
                        ),
                        "target_query_identifier": accession,
                        "target_query_type": "ACCESSION",
                        "target_pref_name": pref_name,
                    })
            elif pref_name:
                rows.append({
                    "pert_id": clean_text(pert_id),
                    "pert_iname": clean_text(pert_iname),
                    "target_chembl_id": clean_text(
                        r.get("target_chembl_id")
                    ),
                    "target_query_identifier": pref_name,
                    "target_query_type": "PREF_NAME",
                    "target_pref_name": pref_name,
                })

    return pd.DataFrame(rows).drop_duplicates(
        [
            "pert_id",
            "pert_iname",
            "target_chembl_id",
            "target_query_identifier",
        ]
    )


# ---------------------------------------------------------------------
# STRING mapping
# ---------------------------------------------------------------------

def map_identifier_to_string(identifier: str) -> dict[str, Any]:
    query_key = clean_text(identifier)

    result = {
        "query_identifier": query_key,
        "string_mapped": False,
        "string_identifier": "",
        "string_preferred_name": "",
        "string_annotation": "",
    }

    if not query_key:
        return result

    url = f"{STRING_BASE}/json/get_string_ids"
    params = {
        "identifiers": query_key,
        "species": STRING_SPECIES,
        "limit": 1,
        "echo_query": 1,
        "caller_identity": "ATLAS",
    }

    status, text = safe_get(url, params, retries=1)

    if status != 200 or not text:
        return result

    try:
        payload = json.loads(text)
    except Exception:
        return result

    if not payload:
        return result

    rec = payload[0]

    result.update({
        "string_mapped": True,
        "string_identifier": clean_text(
            rec.get("stringId")
            or rec.get("string_id")
        ),
        "string_preferred_name": clean_text(
            rec.get("preferredName")
            or rec.get("preferred_name")
        ),
        "string_annotation": clean_text(
            rec.get("annotation")
        ),
    })

    return result


def map_identifiers(
    identifiers: list[str],
    workers: int,
) -> dict[str, dict[str, Any]]:
    cache = load_jsonl_by_key(
        STRING_ID_CACHE,
        "query_identifier",
    )

    missing = [
        x for x in sorted(set(identifiers))
        if x and x not in cache
    ]

    if missing:
        print(
            f"STRING identifier mapping: {len(missing):,} uncached identifiers",
            flush=True,
        )

        lock = threading.Lock()

        with ThreadPoolExecutor(max_workers=workers) as ex:
            futures = {
                ex.submit(map_identifier_to_string, identifier): identifier
                for identifier in missing
            }

            done = 0
            for fut in as_completed(futures):
                rec = fut.result()
                key = rec["query_identifier"]

                with lock:
                    cache[key] = rec
                    append_jsonl(STRING_ID_CACHE, rec)

                done += 1
                if done == 1 or done % 25 == 0 or done == len(missing):
                    print(
                        f"  mapped {done}/{len(missing)}",
                        flush=True,
                    )

    return cache


# ---------------------------------------------------------------------
# STRING network query
# ---------------------------------------------------------------------

def fetch_string_network(
    identifiers: list[str],
    required_score: int,
) -> pd.DataFrame:
    identifiers = [x for x in identifiers if x]

    if len(identifiers) < 2:
        return pd.DataFrame()

    cache_key = "|".join(sorted(identifiers)) + f"|score={required_score}"

    cache = load_jsonl_by_key(
        STRING_NETWORK_CACHE,
        "cache_key",
    )
    if cache_key in cache:
        rows = cache[cache_key].get("rows", [])
        return pd.DataFrame(rows)

    url = f"{STRING_BASE}/tsv/network"
    params = {
        "identifiers": "%0d".join(identifiers),
        "species": STRING_SPECIES,
        "required_score": required_score,
        "network_type": "functional",
        "caller_identity": "ATLAS",
    }

    status, text = safe_get(
        url,
        params,
        retries=1,
        read_timeout=30.0,
    )

    if status != 200 or not text:
        return pd.DataFrame()

    try:
        df = pd.read_csv(
            pd.io.common.StringIO(text),
            sep="\t",
        )
    except Exception:
        return pd.DataFrame()

    append_jsonl(
        STRING_NETWORK_CACHE,
        {
            "cache_key": cache_key,
            "rows": df.to_dict("records"),
        },
    )

    return df


# ---------------------------------------------------------------------
# Network integration
# ---------------------------------------------------------------------

def classify_link(
    target_name: str,
    resistance_gene: str,
    resistance_set: set[str],
) -> str:
    if target_name == resistance_gene and target_name in resistance_set:
        return "DIRECT_TARGET_IS_RESISTANCE_GENE"
    return "STRING_NETWORK_LINK"


def build_network_links(
    network: pd.DataFrame,
    target_names: set[str],
    resistance_genes: set[str],
) -> pd.DataFrame:
    if network.empty:
        return pd.DataFrame()

    # STRING network TSV commonly contains preferredName_A/B.
    a_col = None
    b_col = None

    for candidate in [
        "preferredName_A",
        "preferred_name_A",
    ]:
        if candidate in network.columns:
            a_col = candidate
            break

    for candidate in [
        "preferredName_B",
        "preferred_name_B",
    ]:
        if candidate in network.columns:
            b_col = candidate
            break

    if a_col is None or b_col is None:
        return pd.DataFrame()

    score_col = None
    for candidate in ["score", "combined_score"]:
        if candidate in network.columns:
            score_col = candidate
            break

    rows = []

    for _, r in network.iterrows():
        a = clean_text(r[a_col])
        b = clean_text(r[b_col])

        score = (
            pd.to_numeric(pd.Series([r.get(score_col)]), errors="coerce").iloc[0]
            if score_col
            else np.nan
        )

        if a in target_names and b in resistance_genes:
            rows.append({
                "target_symbol": a,
                "resistance_gene": b,
                "string_score": score,
                "network_link_type": classify_link(
                    a, b, resistance_genes
                ),
            })

        if b in target_names and a in resistance_genes:
            rows.append({
                "target_symbol": b,
                "resistance_gene": a,
                "string_score": score,
                "network_link_type": classify_link(
                    b, a, resistance_genes
                ),
            })

    if not rows:
        return pd.DataFrame()

    return (
        pd.DataFrame(rows)
        .drop_duplicates(
            ["target_symbol", "resistance_gene"],
            keep="first",
        )
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------

def score_target_network(
    target_symbol: str,
    links: pd.DataFrame,
    resistance_set: set[str],
) -> dict[str, Any]:
    target_links = links[
        links["target_symbol"] == target_symbol
    ].copy() if not links.empty else pd.DataFrame()

    linked_n = len(target_links)
    direct = int(target_symbol in resistance_set)

    mean_score = (
        pd.to_numeric(
            target_links["string_score"],
            errors="coerce",
        ).mean()
        if linked_n
        else np.nan
    )

    max_score = (
        pd.to_numeric(
            target_links["string_score"],
            errors="coerce",
        ).max()
        if linked_n
        else np.nan
    )

    # Transparent rule-based biological proximity score.
    score = 0

    if direct:
        score += 4

    if linked_n >= 10:
        score += 4
    elif linked_n >= 5:
        score += 3
    elif linked_n >= 2:
        score += 2
    elif linked_n >= 1:
        score += 1

    if pd.notna(max_score):
        if max_score >= 0.9:
            score += 2
        elif max_score >= 0.7:
            score += 1

    if score >= 6:
        category = "STRONG_NETWORK_SUPPORT"
    elif score >= 3:
        category = "MODERATE_NETWORK_SUPPORT"
    elif score >= 1:
        category = "WEAK_NETWORK_SUPPORT"
    else:
        category = "NO_NETWORK_SUPPORT"

    return {
        "target_symbol": target_symbol,
        "target_is_resistance_gene": bool(direct),
        "linked_resistance_gene_n": int(linked_n),
        "mean_string_score_to_resistance": mean_score,
        "max_string_score_to_resistance": max_score,
        "network_support_score": int(score),
        "network_support_category": category,
        "linked_resistance_genes": (
            " | ".join(
                sorted(
                    set(
                        target_links["resistance_gene"]
                        .dropna()
                        .astype(str)
                    )
                )
            )
            if linked_n
            else ""
        ),
    }


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument(
        "--max-drugs",
        type=int,
        default=25,
        help="Top N 04P drugs to analyze. Use 0 for all. Default: 25",
    )
    p.add_argument(
        "--max-targets-per-drug",
        type=int,
        default=10,
        help="Maximum 04P target rows per drug. Default: 10",
    )
    p.add_argument(
        "--max-resistance-genes",
        type=int,
        default=200,
        help="Maximum ATLAS resistance genes in STRING network. Default: 200",
    )
    p.add_argument(
        "--deg-padj",
        type=float,
        default=0.05,
    )
    p.add_argument(
        "--deg-abs-fc",
        type=float,
        default=1.0,
    )
    p.add_argument(
        "--string-score",
        type=int,
        default=700,
        help="STRING required_score 0-1000. Default: 700 (high confidence).",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=6,
    )

    args = p.parse_args()

    args.max_drugs = max(0, args.max_drugs)
    args.max_targets_per_drug = max(1, args.max_targets_per_drug)
    args.max_resistance_genes = max(20, args.max_resistance_genes)
    args.string_score = min(max(args.string_score, 0), 1000)
    args.workers = max(1, min(args.workers, 12))

    return args


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    header("ATLAS — Stage 04Q PPI + Pathway / Network Integration")

    if not INPUT_PAIRS.exists():
        print(
            f"ERROR: Stage 04P target pairs not found:\n{INPUT_PAIRS}",
            flush=True,
        )
        return 1

    if not INPUT_ANNOT.exists():
        print(
            f"ERROR: Stage 04P annotations not found:\n{INPUT_ANNOT}",
            flush=True,
        )
        return 1

    pairs = pd.read_csv(INPUT_PAIRS)
    annot = pd.read_csv(INPUT_ANNOT)

    annot["priority_rank"] = pd.to_numeric(
        annot["priority_rank"],
        errors="coerce",
    )

    annot = annot.sort_values(
        ["target_evidence_strength_score", "priority_rank"],
        ascending=[False, True],
        na_position="last",
    ).reset_index(drop=True)

    if args.max_drugs > 0:
        selected_drugs = annot.head(args.max_drugs).copy()
    else:
        selected_drugs = annot.copy()

    selected_keys = set(
        zip(
            selected_drugs["pert_id"].astype(str),
            selected_drugs["pert_iname"].astype(str),
        )
    )

    pairs = pairs[
        [
            (str(r["pert_id"]), str(r["pert_iname"]))
            in selected_keys
            for _, r in pairs.iterrows()
        ]
    ].copy()

    print(
        f"\nSelected drugs: {len(selected_drugs):,}",
        flush=True,
    )
    print(
        f"04P target-pair rows: {len(pairs):,}",
        flush=True,
    )

    resistance = load_resistance_genes(
        max_genes=args.max_resistance_genes,
        padj_cutoff=args.deg_padj,
        abs_fc_cutoff=args.deg_abs_fc,
    )

    resistance_genes = sorted(
        set(resistance["gene_symbol"].dropna().astype(str))
    )
    resistance_set = set(resistance_genes)

    print(
        f"Resistance genes selected: {len(resistance_genes):,}",
        flush=True,
    )

    target_queries = collect_candidate_targets(
        pairs,
        max_targets_per_drug=args.max_targets_per_drug,
    )

    if target_queries.empty:
        print(
            "ERROR: No target identifiers available from 04P.",
            flush=True,
        )
        return 1

    identifiers_to_map = (
        target_queries["target_query_identifier"]
        .dropna()
        .astype(str)
        .tolist()
        + resistance_genes
    )

    id_map = map_identifiers(
        identifiers_to_map,
        workers=args.workers,
    )

    # Attach mapped target symbols.
    target_queries["string_mapped"] = target_queries[
        "target_query_identifier"
    ].map(
        lambda x: bool(
            id_map.get(clean_text(x), {}).get(
                "string_mapped", False
            )
        )
    )

    target_queries["target_symbol"] = target_queries[
        "target_query_identifier"
    ].map(
        lambda x: clean_text(
            id_map.get(clean_text(x), {}).get(
                "string_preferred_name"
            )
        )
    )

    target_queries["string_identifier"] = target_queries[
        "target_query_identifier"
    ].map(
        lambda x: clean_text(
            id_map.get(clean_text(x), {}).get(
                "string_identifier"
            )
        )
    )

    mapped_targets = target_queries[
        target_queries["string_mapped"]
        & target_queries["target_symbol"].ne("")
    ].copy()

    target_symbols = sorted(
        set(mapped_targets["target_symbol"].astype(str))
    )

    mapped_resistance = []
    for gene in resistance_genes:
        rec = id_map.get(gene, {})
        preferred = clean_text(rec.get("string_preferred_name"))
        if rec.get("string_mapped") and preferred:
            mapped_resistance.append(preferred)

    mapped_resistance = sorted(set(mapped_resistance))

    print(
        f"Mapped targets: {len(target_symbols):,}",
        flush=True,
    )
    print(
        f"Mapped resistance genes: {len(mapped_resistance):,}",
        flush=True,
    )

    network_identifiers = sorted(
        set(target_symbols) | set(mapped_resistance)
    )

    print(
        f"STRING network nodes submitted: {len(network_identifiers):,}",
        flush=True,
    )

    network = fetch_string_network(
        network_identifiers,
        required_score=args.string_score,
    )

    if network.empty:
        print(
            "WARNING: STRING returned no network edges.",
            flush=True,
        )

    atomic_csv(network, OUT_EDGES)

    links = build_network_links(
        network,
        target_names=set(target_symbols),
        resistance_genes=set(mapped_resistance),
    )

    if links.empty:
        links = pd.DataFrame(columns=[
            "target_symbol",
            "resistance_gene",
            "string_score",
            "network_link_type",
        ])

    atomic_csv(links, OUT_LINKS)

    target_score_rows = [
        score_target_network(
            target_symbol,
            links,
            set(mapped_resistance),
        )
        for target_symbol in target_symbols
    ]

    target_scores = pd.DataFrame(target_score_rows)

    if target_scores.empty:
        target_scores = pd.DataFrame(columns=[
            "target_symbol",
            "target_is_resistance_gene",
            "linked_resistance_gene_n",
            "mean_string_score_to_resistance",
            "max_string_score_to_resistance",
            "network_support_score",
            "network_support_category",
            "linked_resistance_genes",
        ])

    atomic_csv(target_scores, OUT_TARGET_SCORES)

    # Drug-level aggregation.
    tq = mapped_targets.merge(
        target_scores,
        on="target_symbol",
        how="left",
    )

    drug_rows = []

    for (pert_id, pert_iname), g in tq.groupby(
        ["pert_id", "pert_iname"],
        dropna=False,
    ):
        scores = pd.to_numeric(
            g["network_support_score"],
            errors="coerce",
        )

        best_idx = (
            scores.idxmax()
            if scores.notna().any()
            else g.index[0]
        )
        best = g.loc[best_idx]

        drug_rows.append({
            "pert_id": pert_id,
            "pert_iname": pert_iname,
            "network_target_n": int(
                g["target_symbol"].nunique()
            ),
            "best_network_target": clean_text(
                best.get("target_symbol")
            ),
            "best_network_support_score": (
                float(best.get("network_support_score"))
                if pd.notna(best.get("network_support_score"))
                else np.nan
            ),
            "best_network_support_category": clean_text(
                best.get("network_support_category")
            ),
            "best_target_linked_resistance_gene_n": int(
                best.get("linked_resistance_gene_n", 0)
                if pd.notna(
                    best.get("linked_resistance_gene_n", 0)
                )
                else 0
            ),
            "best_target_linked_resistance_genes": clean_text(
                best.get("linked_resistance_genes")
            ),
            "all_network_targets": " | ".join(
                sorted(
                    set(
                        g["target_symbol"]
                        .dropna()
                        .astype(str)
                    )
                )
            ),
        })

    drug_network = pd.DataFrame(drug_rows)

    prioritized = selected_drugs.merge(
        drug_network,
        on=["pert_id", "pert_iname"],
        how="left",
    )

    prioritized["best_network_support_score"] = pd.to_numeric(
        prioritized["best_network_support_score"],
        errors="coerce",
    ).fillna(0)

    # Integrated 04P + 04Q score.
    prioritized["target_evidence_strength_score"] = pd.to_numeric(
        prioritized["target_evidence_strength_score"],
        errors="coerce",
    ).fillna(0)

    prioritized["integrated_target_network_score"] = (
        prioritized["target_evidence_strength_score"]
        + prioritized["best_network_support_score"]
    )

    prioritized = prioritized.sort_values(
        [
            "integrated_target_network_score",
            "target_evidence_strength_score",
            "priority_rank",
        ],
        ascending=[False, False, True],
        na_position="last",
    ).reset_index(drop=True)

    atomic_csv(prioritized, OUT_DRUG_PRIORITY)

    summary = pd.DataFrame([
        {
            "selected_drugs": len(selected_drugs),
            "target_query_rows": len(target_queries),
            "mapped_target_symbols": len(target_symbols),
            "resistance_genes_requested": len(resistance_genes),
            "mapped_resistance_genes": len(mapped_resistance),
            "string_network_edges": len(network),
            "target_resistance_links": len(links),
            "targets_with_network_support": int(
                (
                    pd.to_numeric(
                        target_scores.get(
                            "network_support_score",
                            pd.Series(dtype=float),
                        ),
                        errors="coerce",
                    ) > 0
                ).sum()
            ),
        }
    ])

    atomic_csv(summary, OUT_SUMMARY)

    metadata = {
        "stage": "04Q",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_target_pairs": str(INPUT_PAIRS),
        "input_target_annotations": str(INPUT_ANNOT),
        "discovery_deg": str(DISCOVERY_DEG),
        "consensus_candidates_checked": [
            str(x) for x in CONSENSUS_CANDIDATES
        ],
        "max_drugs": args.max_drugs,
        "max_targets_per_drug": args.max_targets_per_drug,
        "max_resistance_genes": args.max_resistance_genes,
        "deg_padj_cutoff": args.deg_padj,
        "deg_abs_fc_cutoff": args.deg_abs_fc,
        "string_required_score": args.string_score,
        "species": "Homo sapiens",
        "species_taxon": STRING_SPECIES,
        "network_type": "functional",
        "sources": [
            "ATLAS Stage 04P target annotations",
            "ATLAS discovery/consensus resistance genes",
            "STRING functional association network",
        ],
        "guardrails": [
            "STRING association does not necessarily represent direct physical binding.",
            "Network proximity supports plausibility but does not prove resistance reversal.",
            "Docking candidates should require convergent target/network evidence.",
        ],
        "next_stage": "04R_final_candidate_prioritization",
    }

    atomic_json(metadata, OUT_META)

    header("STAGE 04Q SUMMARY")
    print(summary.to_string(index=False), flush=True)

    header("TOP NETWORK-SUPPORTED DRUGS")

    display_cols = [
        c for c in [
            "priority_rank",
            "pert_iname",
            "target_support_category",
            "target_evidence_strength_score",
            "best_network_target",
            "best_network_support_category",
            "best_network_support_score",
            "best_target_linked_resistance_gene_n",
            "integrated_target_network_score",
        ]
        if c in prioritized.columns
    ]

    print(
        prioritized[display_cols]
        .head(15)
        .to_string(index=False),
        flush=True,
    )

    header("STAGE 04Q COMPLETE")
    print("\nOutputs:", flush=True)
    print(f"  {OUT_EDGES}", flush=True)
    print(f"  {OUT_LINKS}", flush=True)
    print(f"  {OUT_TARGET_SCORES}", flush=True)
    print(f"  {OUT_DRUG_PRIORITY}", flush=True)
    print(f"  {OUT_SUMMARY}", flush=True)
    print(f"  {OUT_META}", flush=True)
    print(
        "\nNext: 04R — Final candidate prioritization",
        flush=True,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
