#!/usr/bin/env python3
"""
ATLAS — Stage 03C: Consensus Trastuzumab-Resistance Signature

Input
-----
results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv
results/external_validation/ATLAS_external_gene_validation.csv

Purpose
-------
Integrate the ATLAS discovery direction with independent external datasets
without pooling raw expression matrices.

Consensus logic
---------------
HIGH confidence:
    ATLAS significant resistance DEG and same direction in >= 2 external
    datasets.

MODERATE confidence:
    ATLAS significant resistance DEG and same direction in 1 external dataset.

DISCOVERY_ONLY:
    ATLAS significant resistance DEG with no external same-direction support.

CONTRADICTORY:
    ATLAS significant resistance DEG and external evidence is predominantly
    opposite in direction.

Output
------
results/consensus_signature/
    ATLAS_consensus_resistance_genes.csv
    ATLAS_consensus_high_confidence.csv
    ATLAS_consensus_moderate_confidence.csv
    ATLAS_consensus_signature_up.gmt
    ATLAS_consensus_signature_dn.gmt
    ATLAS_consensus_summary.csv
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DISCOVERY_FILE = (
    PROJECT_ROOT
    / "results"
    / "differential_expression"
    / "DEGs_resistant_vs_sensitive_annotated.csv"
)

EXTERNAL_FILE = (
    PROJECT_ROOT
    / "results"
    / "external_validation"
    / "ATLAS_external_gene_validation.csv"
)

OUT_DIR = PROJECT_ROOT / "results" / "consensus_signature"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ALL_FILE = OUT_DIR / "ATLAS_consensus_resistance_genes.csv"
HIGH_FILE = OUT_DIR / "ATLAS_consensus_high_confidence.csv"
MOD_FILE = OUT_DIR / "ATLAS_consensus_moderate_confidence.csv"
SUMMARY_FILE = OUT_DIR / "ATLAS_consensus_summary.csv"
UP_GMT = OUT_DIR / "ATLAS_consensus_signature_up.gmt"
DN_GMT = OUT_DIR / "ATLAS_consensus_signature_dn.gmt"
META_FILE = OUT_DIR / "ATLAS_consensus_metadata.json"


def clean_symbol(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def load_discovery() -> pd.DataFrame:
    df = pd.read_csv(DISCOVERY_FILE)

    symbol_col = None
    for c in ["Gene name", "gene_symbol", "symbol", "Gene Symbol"]:
        if c in df.columns:
            symbol_col = c
            break

    if symbol_col is None:
        raise ValueError("No gene-symbol column found in discovery DEG file.")

    out = pd.DataFrame(
        {
            "gene_symbol": df[symbol_col].map(clean_symbol),
            "atlas_log2FC": pd.to_numeric(
                df["log2FoldChange"], errors="coerce"
            ),
            "atlas_padj": pd.to_numeric(
                df["padj"], errors="coerce"
            ),
        }
    )
    out = out[
        out["gene_symbol"].ne("")
        & out["atlas_log2FC"].notna()
    ]
    out["atlas_abs_fc"] = out["atlas_log2FC"].abs()
    out = (
        out.sort_values(
            ["gene_symbol", "atlas_abs_fc"],
            ascending=[True, False],
        )
        .drop_duplicates("gene_symbol")
        .drop(columns=["atlas_abs_fc"])
    )

    out["atlas_significant"] = (
        out["atlas_padj"].lt(0.05)
        & out["atlas_log2FC"].abs().ge(1)
    )
    out["atlas_direction"] = np.where(
        out["atlas_log2FC"] > 0,
        "UP",
        np.where(out["atlas_log2FC"] < 0, "DOWN", "FLAT"),
    )
    return out


def write_gmt(path: Path, name: str, genes: list[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "\t".join([name, "ATLAS 03C consensus resistance signature"] + genes)
            + "\n"
        )


def main() -> int:
    print("=" * 72)
    print("ATLAS — Stage 03C Consensus Resistance Signature")
    print("=" * 72)

    if not DISCOVERY_FILE.exists():
        print(f"ERROR: discovery file not found: {DISCOVERY_FILE}")
        return 1
    if not EXTERNAL_FILE.exists():
        print(f"ERROR: 03B file not found: {EXTERNAL_FILE}")
        return 1

    discovery = load_discovery()
    ext = pd.read_csv(EXTERNAL_FILE)

    required = {
        "validation_dataset",
        "gene_symbol",
        "external_log2FC",
    }
    missing = required - set(ext.columns)
    if missing:
        print(f"ERROR: missing 03B columns: {sorted(missing)}")
        return 1

    ext["gene_symbol"] = ext["gene_symbol"].map(clean_symbol)
    ext["external_log2FC"] = pd.to_numeric(
        ext["external_log2FC"],
        errors="coerce",
    )

    # Reduce each dataset/gene to one effect.
    ext = (
        ext.sort_values(
            "external_log2FC",
            key=lambda s: s.abs(),
            ascending=False,
        )
        .drop_duplicates(
            ["validation_dataset", "gene_symbol"]
        )
    )

    rows = []

    for _, row in discovery.iterrows():
        gene = row["gene_symbol"]
        atlas_dir = row["atlas_direction"]

        sub = ext[ext["gene_symbol"] == gene].copy()

        directions = np.where(
            sub["external_log2FC"] > 0,
            "UP",
            np.where(sub["external_log2FC"] < 0, "DOWN", "FLAT"),
        )

        same = int((directions == atlas_dir).sum())
        opposite = int(
            (
                ((atlas_dir == "UP") & (directions == "DOWN"))
                | ((atlas_dir == "DOWN") & (directions == "UP"))
            ).sum()
        )
        tested = int(len(sub))

        if not bool(row["atlas_significant"]):
            confidence = "NOT_DISCOVERY_SIGNATURE"
        elif same >= 2:
            confidence = "HIGH"
        elif same == 1:
            confidence = "MODERATE"
        elif opposite > same and opposite > 0:
            confidence = "CONTRADICTORY"
        else:
            confidence = "DISCOVERY_ONLY"

        rows.append(
            {
                "gene_symbol": gene,
                "atlas_log2FC": row["atlas_log2FC"],
                "atlas_padj": row["atlas_padj"],
                "atlas_direction": atlas_dir,
                "atlas_significant": row["atlas_significant"],
                "external_datasets_tested": tested,
                "external_same_direction": same,
                "external_opposite_direction": opposite,
                "consensus_confidence": confidence,
                "external_dataset_names": " | ".join(
                    sorted(sub["validation_dataset"].astype(str).unique())
                ),
            }
        )

    consensus = pd.DataFrame(rows)

    # Ranking keeps ATLAS effect strength but puts reproducibility first.
    confidence_order = {
        "HIGH": 0,
        "MODERATE": 1,
        "DISCOVERY_ONLY": 2,
        "CONTRADICTORY": 3,
        "NOT_DISCOVERY_SIGNATURE": 4,
    }
    consensus["confidence_order"] = consensus[
        "consensus_confidence"
    ].map(confidence_order)
    consensus["atlas_abs_log2FC"] = consensus["atlas_log2FC"].abs()

    consensus = (
        consensus.sort_values(
            [
                "confidence_order",
                "external_same_direction",
                "atlas_abs_log2FC",
                "atlas_padj",
            ],
            ascending=[True, False, False, True],
        )
        .drop(columns=["confidence_order"])
        .reset_index(drop=True)
    )

    consensus.insert(
        0,
        "consensus_rank",
        np.arange(1, len(consensus) + 1),
    )

    consensus.to_csv(ALL_FILE, index=False)

    high = consensus[
        consensus["consensus_confidence"] == "HIGH"
    ].copy()
    moderate = consensus[
        consensus["consensus_confidence"] == "MODERATE"
    ].copy()

    high.to_csv(HIGH_FILE, index=False)
    moderate.to_csv(MOD_FILE, index=False)

    signature = consensus[
        consensus["consensus_confidence"].isin(["HIGH", "MODERATE"])
    ].copy()

    up = signature[
        signature["atlas_direction"] == "UP"
    ]["gene_symbol"].drop_duplicates().tolist()

    down = signature[
        signature["atlas_direction"] == "DOWN"
    ]["gene_symbol"].drop_duplicates().tolist()

    write_gmt(
        UP_GMT,
        "ATLAS_CONSENSUS_RESISTANCE_UP",
        up,
    )
    write_gmt(
        DN_GMT,
        "ATLAS_CONSENSUS_RESISTANCE_DOWN",
        down,
    )

    summary = (
        consensus["consensus_confidence"]
        .value_counts()
        .rename_axis("consensus_confidence")
        .reset_index(name="gene_count")
    )
    summary.to_csv(SUMMARY_FILE, index=False)

    metadata = {
        "stage": "03C",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "discovery_file": str(DISCOVERY_FILE),
        "external_validation_file": str(EXTERNAL_FILE),
        "high_confidence_rule": (
            "ATLAS significant DEG and same direction in >=2 external datasets"
        ),
        "moderate_confidence_rule": (
            "ATLAS significant DEG and same direction in 1 external dataset"
        ),
        "consensus_up_genes": len(up),
        "consensus_down_genes": len(down),
        "important_note": (
            "03C integrates direction-level reproducibility across independently "
            "analyzed datasets. It does not pool raw expression matrices."
        ),
    }
    META_FILE.write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    print()
    print("Consensus counts:")
    for _, r in summary.iterrows():
        print(f"  {r['consensus_confidence']}: {int(r['gene_count']):,}")

    print()
    print(f"Consensus UP genes:   {len(up):,}")
    print(f"Consensus DOWN genes: {len(down):,}")

    print()
    print("Outputs:")
    print(f"  {ALL_FILE}")
    print(f"  {HIGH_FILE}")
    print(f"  {MOD_FILE}")
    print(f"  {UP_GMT}")
    print(f"  {DN_GMT}")
    print(f"  {SUMMARY_FILE}")

    print()
    print("=" * 72)
    print("STAGE 03C COMPLETE")
    print("=" * 72)
    print(
        "Next: compare the 03C consensus signature with the original CMap "
        "signatures and optionally run a validation CMap query."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
