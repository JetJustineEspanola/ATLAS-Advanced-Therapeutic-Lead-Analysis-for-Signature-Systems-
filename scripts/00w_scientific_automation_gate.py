#!/usr/bin/env python3
"""
ATLAS — 00W Scientific Automation Gate v3

Key correction:
The final 00C4 independence/modality scorer is treated as the authoritative
dataset-selection layer. This gate does NOT re-apply stale earlier flags such
as manual_review_required or raw independence booleans after 00C4 has already
integrated curated phenotype, HER2, replication, relationship, and modality
evidence.

Primary rule:
    independence_aware_category == PRIMARY_VALIDATION

Secondary sanity checks:
- transcriptomic modality
- not umbrella/duplicate
- explicit resistant/sensitive groups when available
- direct trastuzumab resistance when available
- HER2 confirmed when available
- biological replication when available

The gate fails safely if fewer than --min-primary datasets remain.

Outputs
-------
results/pipeline_state/scientific_automation_gate.csv
results/pipeline_state/scientific_automation_gate.json
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "enriched" / "transcriptomic_validation_candidates.csv"

OUTDIR = ROOT / "results" / "pipeline_state"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUTDIR / "scientific_automation_gate.csv"
OUT_JSON = OUTDIR / "scientific_automation_gate.json"


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def norm_col(x):
    return str(x).strip().lower().replace(" ", "_")


def find_col(df, aliases):
    cmap = {norm_col(c): c for c in df.columns}
    for a in aliases:
        if norm_col(a) in cmap:
            return cmap[norm_col(a)]
    return None


def to_bool(x):
    if pd.isna(x):
        return False
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, float)):
        return bool(x)

    s = str(x).strip().lower()
    return s in {
        "true", "1", "yes", "y", "confirmed", "pass", "passed",
        "eligible", "high",
    }


def bool_series(s):
    return s.map(to_bool)


def contains_any(s, terms):
    t = s.fillna("").astype(str).str.lower()
    mask = pd.Series(False, index=s.index)
    for term in terms:
        mask |= t.str.contains(term.lower(), regex=False)
    return mask


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DEFAULT_INPUT))
    ap.add_argument("--min-primary", type=int, default=2)
    args = ap.parse_args()

    path = Path(args.input)
    if not path.exists():
        print(f"ERROR: candidate table not found: {path}")
        return 2

    df = pd.read_csv(path)
    if df.empty:
        print(f"ERROR: candidate table is empty: {path}")
        return 2

    accession_col = find_col(
        df,
        ["source_accession", "accession", "dataset_accession", "dataset_id"],
    )
    category_col = find_col(
        df,
        ["independence_aware_category", "validation_role", "eligibility_category"],
    )

    if accession_col is None:
        print("ERROR: no accession column found.")
        return 2

    if category_col is None:
        print("ERROR: no final validation-category column found.")
        return 2

    out = df.copy()
    out["gate_accession"] = out[accession_col].astype(str).str.strip()

    # 00C4 is the authoritative integrated selection result.
    category_text = out[category_col].fillna("").astype(str).str.upper()
    primary = category_text.eq("PRIMARY_VALIDATION")
    # tolerate harmless formatting variants
    primary |= category_text.str.replace(" ", "_", regex=False).eq("PRIMARY_VALIDATION")
    primary |= category_text.str.contains("PRIMARY_VALIDATION", regex=False)

    out["gate_primary_validation_00c4"] = primary

    # Secondary checks are applied only if the relevant column exists.
    # These are sanity checks, not a full reimplementation of 00C4.
    sanity = pd.Series(True, index=out.index)
    checks = [f"authoritative 00C4 category via {category_col}"]
    warnings = []

    modality_col = find_col(out, ["modality_class", "assay_type"])
    if modality_col:
        modality = contains_any(
            out[modality_col],
            ["rna", "transcript", "expression", "bulk"],
        )
        out["gate_transcriptomic"] = modality
        sanity &= modality
        checks.append(f"transcriptomic sanity via {modality_col}")

    relationship_col = find_col(out, ["relationship_role"])
    if relationship_col:
        bad = contains_any(
            out[relationship_col],
            ["umbrella", "duplicate", "parent_series"],
        )
        out["gate_not_umbrella_duplicate"] = ~bad
        sanity &= ~bad
        checks.append(f"relationship sanity via {relationship_col}")

    optional_bools = [
        ("gate_direct_resistance", ["direct_trastuzumab_resistance", "direct_resistance"]),
        ("gate_her2_confirmed", ["her2_positive_confirmed", "her2_confirmed"]),
        ("gate_groups_defined", ["resistant_sensitive_groups_defined", "groups_defined"]),
        ("gate_biological_replication", ["biological_replication", "has_replication"]),
    ]

    for gate_name, aliases in optional_bools:
        col = find_col(out, aliases)
        if col:
            mask = bool_series(out[col])
            out[gate_name] = mask
            sanity &= mask
            checks.append(f"{gate_name.replace('gate_', '')} sanity via {col}")
        else:
            warnings.append(
                f"Optional sanity column not found: {gate_name.replace('gate_', '')}"
            )

    out["scientific_gate_pass"] = primary & sanity

    gate_cols = [
        c for c in out.columns
        if c.startswith("gate_") and c != "gate_accession"
    ]

    def fail_reasons(row):
        failed = []
        for c in gate_cols:
            if not bool(row.get(c, False)):
                failed.append(c.replace("gate_", ""))
        return " | ".join(failed)

    out["scientific_gate_fail_reasons"] = out.apply(fail_reasons, axis=1)

    selected = out[out["scientific_gate_pass"]]
    accessions = sorted(
        selected["gate_accession"]
        .dropna()
        .astype(str)
        .loc[lambda x: x.ne("")]
        .unique()
        .tolist()
    )

    gate_pass = len(accessions) >= max(1, args.min_primary)

    out.to_csv(OUT_CSV, index=False)

    report = {
        "created_utc": utcnow(),
        "input": str(path),
        "accession_column": accession_col,
        "authoritative_category_column": category_col,
        "candidate_rows": int(len(out)),
        "minimum_primary_required": int(args.min_primary),
        "selected_primary_dataset_n": int(len(accessions)),
        "selected_primary_accessions": accessions,
        "checks_applied": checks,
        "warnings": warnings,
        "gate_pass": bool(gate_pass),
        "selection_principle": (
            "00C4 independence_aware_category is authoritative. Earlier raw/manual "
            "flags are not re-applied after integrated scoring."
        ),
    }
    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("=" * 78)
    print("ATLAS — SCIENTIFIC AUTOMATION GATE v3")
    print("=" * 78)
    print(f"Authoritative category: {category_col}")
    print(f"Candidate rows: {len(out)}")
    print(f"Usable primary datasets: {len(accessions)}")
    print("Selected:", ", ".join(accessions) if accessions else "NONE")
    print(f"Required minimum: {args.min_primary}")
    print(f"GATE: {'PASS' if gate_pass else 'FAIL'}")

    print("\nPrimary/near-primary rows:")
    interesting = (
        primary
        | category_text.str.contains("PRIMARY", regex=False)
        | category_text.str.contains("EXPLOR", regex=False)
    )
    show = out.loc[
        interesting,
        [
            "gate_accession",
            category_col,
            "scientific_gate_pass",
            "scientific_gate_fail_reasons",
        ],
    ]
    print(show.to_string(index=False) if not show.empty else "None")

    print("\nOutputs:")
    print(OUT_CSV)
    print(OUT_JSON)

    return 0 if gate_pass else 3


if __name__ == "__main__":
    raise SystemExit(main())
