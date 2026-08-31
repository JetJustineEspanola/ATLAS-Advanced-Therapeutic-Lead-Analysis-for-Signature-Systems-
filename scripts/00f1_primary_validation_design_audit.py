#!/usr/bin/env python3
"""
ATLAS — 00F1 v2 Primary Validation Design Audit

Fixes GSE237606 stimulus parsing:
  S1_Con, S1_HRG, S1_EGF
  S12_Con, S12_HRG, S12_EGF
  S24_Con, S24_HRG, S24_EGF

Outputs:
  data/validation_expression/primary_validation_design_audit.csv
  data/validation_expression/GSE121105_primary_design.csv
  data/validation_expression/GSE237606_primary_design.csv
"""

from pathlib import Path
import re
import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "data/catalog/atlas_catalog.duckdb"
OUT = ROOT / "data/validation_expression"
OUT.mkdir(parents=True, exist_ok=True)


def clean(x):
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()


def parse_237606(title: str):
    t = clean(title)

    # Titles are expected to begin with S... or R...
    phenotype = "RESISTANT" if t.upper().startswith("R") else "SENSITIVE_OR_PARENTAL"

    if "Baseline" in t:
        timepoint = "0"
        stimulus = "BASELINE"
    else:
        m = re.search(r"[SR](\d+)_([A-Za-z]+)", t)
        if m:
            timepoint = m.group(1)
            stimulus = m.group(2).upper()
        else:
            m_time = re.search(r"(\d+)\s*hr", t, re.I)
            timepoint = m_time.group(1) if m_time else "UNKNOWN"
            stimulus = "UNKNOWN"

    m_rep = re.search(r"rep\s*(\d+)", t, re.I)
    replicate = m_rep.group(1) if m_rep else "UNKNOWN"

    return phenotype, timepoint, stimulus, replicate


def main():
    con = duckdb.connect(str(DB))
    samples = con.execute(
        """
        SELECT
            dataset_id,
            sample_id,
            title,
            description,
            resistance_status,
            phenotype_confidence,
            primary_contrast_included
        FROM samples
        WHERE dataset_id IN ('GEO:GSE121105','GEO:GSE237606')
        ORDER BY dataset_id, sample_id
        """
    ).fetchdf()
    con.close()

    rows = []

    g = samples[samples["dataset_id"] == "GEO:GSE121105"].copy()
    g["included"] = g["primary_contrast_included"].fillna(False).astype(bool)
    g["phenotype"] = g["resistance_status"]
    g["timepoint"] = "NA"
    g["stimulus"] = "UNTREATED_PRIMARY_CONTRAST"
    g["replicate"] = (
        g["title"].astype(str)
        .str.extract(r"replicate(\d+)", expand=False)
        .fillna("UNKNOWN")
    )
    g.to_csv(OUT / "GSE121105_primary_design.csv", index=False)

    h = samples[samples["dataset_id"] == "GEO:GSE237606"].copy()
    parsed = h["title"].apply(parse_237606)
    h["phenotype"] = [x[0] for x in parsed]
    h["timepoint"] = [x[1] for x in parsed]
    h["stimulus"] = [x[2] for x in parsed]
    h["replicate"] = [x[3] for x in parsed]
    h["included"] = True
    h.to_csv(OUT / "GSE237606_primary_design.csv", index=False)

    print("=" * 78)
    print("ATLAS — 00F1 v2 PRIMARY VALIDATION DESIGN AUDIT")
    print("=" * 78)

    print("\n[GSE121105] curated primary contrast")
    inc = g[g["included"]]
    print(inc.groupby("phenotype").size().to_string())
    print(f"Included: {len(inc)} / {len(g)}")

    print("\n[GSE237606] phenotype x timepoint x stimulus")
    tab = (
        h.groupby(["phenotype", "timepoint", "stimulus"])
        .size()
        .reset_index(name="n")
        .sort_values(["timepoint", "stimulus", "phenotype"])
    )
    print(tab.to_string(index=False))

    pivot = tab.pivot_table(
        index=["timepoint", "stimulus"],
        columns="phenotype",
        values="n",
        fill_value=0,
        aggfunc="sum",
    ).reset_index()

    for c in ["RESISTANT", "SENSITIVE_OR_PARENTAL"]:
        if c not in pivot.columns:
            pivot[c] = 0

    pivot["balanced"] = (
        pivot["RESISTANT"] == pivot["SENSITIVE_OR_PARENTAL"]
    )
    pivot["both_present"] = (
        (pivot["RESISTANT"] > 0)
        & (pivot["SENSITIVE_OR_PARENTAL"] > 0)
    )

    print("\n[GSE237606] stratum balance")
    print(pivot.to_string(index=False))

    all_balanced = bool((pivot["balanced"] & pivot["both_present"]).all())
    print(f"\nGSE237606 fully phenotype-balanced across strata: {all_balanced}")

    rows.append({
        "dataset_id": "GEO:GSE121105",
        "included_n": len(inc),
        "resistant_n": int((inc["phenotype"] == "RESISTANT").sum()),
        "sensitive_parental_n": int((inc["phenotype"] == "SENSITIVE_OR_PARENTAL").sum()),
        "balanced_across_strata": None,
    })
    rows.append({
        "dataset_id": "GEO:GSE237606",
        "included_n": len(h),
        "resistant_n": int((h["phenotype"] == "RESISTANT").sum()),
        "sensitive_parental_n": int((h["phenotype"] == "SENSITIVE_OR_PARENTAL").sum()),
        "balanced_across_strata": all_balanced,
    })

    summary_path = OUT / "primary_validation_design_audit.csv"
    pd.DataFrame(rows).to_csv(summary_path, index=False)

    print("\nOutput:")
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
