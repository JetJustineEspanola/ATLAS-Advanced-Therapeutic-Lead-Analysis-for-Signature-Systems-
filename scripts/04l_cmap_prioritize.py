#!/usr/bin/env python3
"""
ATLAS — Stage 4L CMap Candidate Prioritization

Purpose
-------
Prioritize CMap compounds using cross-signature consistency.

Input
-----
results/cmap/parsed/ATLAS_CMap_cross_signature_tau.csv

The input is expected to contain:
    consensus_rank
    pert_id
    pert_iname
    ATLAS_SIG_A_TOP150
    ATLAS_SIG_B_TOP100
    ATLAS_SIG_B_TOP150
    n_signatures
    n_negative
    n_strong_negative
    mean_tau
    median_tau
    minimum_tau

Prioritization
--------------
Tier 1 — Strong consensus
    >= 2 strongly negative signatures

Tier 2 — Consistent reversal
    3 negative signatures, but < 2 strongly negative

Tier 3 — Partial reversal
    2 negative signatures

Tier 4 — Weak/single-signature
    <= 1 negative signature

Within tiers, compounds are ranked by:
    1. n_strong_negative  descending
    2. n_negative         descending
    3. median_tau         ascending
    4. mean_tau           ascending
    5. minimum_tau        ascending

Notes
-----
- Negative tau indicates reversal of the queried resistance signature.
- This stage does NOT determine FDA approval status.
- This stage does NOT perform molecular docking.
- The original CMap results are never modified.
"""

from pathlib import Path
import sys

import pandas as pd


# ============================================================
# CONFIGURATION
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    BASE_DIR
    / "results"
    / "cmap"
    / "parsed"
    / "ATLAS_CMap_cross_signature_tau.csv"
)

OUTPUT_DIR = (
    BASE_DIR
    / "results"
    / "cmap"
    / "prioritized"
)

ALL_OUTPUT = OUTPUT_DIR / "ATLAS_CMap_prioritized_all.csv"
TIER1_OUTPUT = OUTPUT_DIR / "ATLAS_CMap_tier1_strong_consensus.csv"
TIER2_OUTPUT = OUTPUT_DIR / "ATLAS_CMap_tier2_consistent.csv"
TIER3_OUTPUT = OUTPUT_DIR / "ATLAS_CMap_tier3_partial.csv"
TIER4_OUTPUT = OUTPUT_DIR / "ATLAS_CMap_tier4_weak_single.csv"
TOP_OUTPUT = OUTPUT_DIR / "ATLAS_CMap_top_candidates.csv"

# Strong negative tau threshold.
#
# This matches the existing Stage 4K concept of "strong negative"
# candidates, where tau <= -90 is treated as strongly negative.
STRONG_NEGATIVE_THRESHOLD = -90.0

# Number of candidates included in the compact top-candidates file.
TOP_N = 30


# ============================================================
# REQUIRED COLUMNS
# ============================================================

REQUIRED_COLUMNS = [
    "consensus_rank",
    "pert_id",
    "pert_iname",
    "ATLAS_SIG_A_TOP150",
    "ATLAS_SIG_B_TOP100",
    "ATLAS_SIG_B_TOP150",
    "n_signatures",
    "n_negative",
    "n_strong_negative",
    "mean_tau",
    "median_tau",
    "minimum_tau",
]

SIGNATURE_COLUMNS = [
    "ATLAS_SIG_A_TOP150",
    "ATLAS_SIG_B_TOP100",
    "ATLAS_SIG_B_TOP150",
]


# ============================================================
# HELPERS
# ============================================================

def fail(message: str) -> None:
    """Print an error and exit cleanly."""
    print()
    print("ERROR:")
    print(message)
    print()
    sys.exit(1)


def classify_tier(row: pd.Series) -> str:
    """
    Assign a CMap prioritization tier.

    Tier 1:
        >= 2 strongly negative signatures

    Tier 2:
        all 3 signatures negative, but fewer than 2 strongly negative

    Tier 3:
        exactly 2 negative signatures

    Tier 4:
        <= 1 negative signature
    """

    n_negative = int(row["n_negative"])
    n_strong_negative = int(row["n_strong_negative"])

    if n_strong_negative >= 2:
        return "Tier 1 — Strong consensus"

    if n_negative == 3:
        return "Tier 2 — Consistent reversal"

    if n_negative == 2:
        return "Tier 3 — Partial reversal"

    return "Tier 4 — Weak/single-signature"


def tier_number(tier: str) -> int:
    """Convert tier label into a numeric ordering."""
    if tier.startswith("Tier 1"):
        return 1
    if tier.startswith("Tier 2"):
        return 2
    if tier.startswith("Tier 3"):
        return 3
    return 4


# ============================================================
# MAIN
# ============================================================

def main() -> None:

    print("=" * 60)
    print("ATLAS — Stage 4L CMap Candidate Prioritization")
    print("=" * 60)

    # --------------------------------------------------------
    # Check input
    # --------------------------------------------------------

    print()
    print("Input:")
    print(INPUT_FILE)

    if not INPUT_FILE.exists():
        fail(
            f"Input file does not exist:\n{INPUT_FILE}\n\n"
            "Run Stage 4K first."
        )

    # --------------------------------------------------------
    # Load
    # --------------------------------------------------------

    print()
    print("Loading cross-signature CMap results...")

    try:
        df = pd.read_csv(INPUT_FILE)
    except Exception as exc:
        fail(f"Could not read input CSV:\n{exc}")

    print(f"Rows loaded: {len(df):,}")
    print(f"Columns loaded: {len(df.columns)}")

    # --------------------------------------------------------
    # Validate columns
    # --------------------------------------------------------

    missing = [
        col for col in REQUIRED_COLUMNS
        if col not in df.columns
    ]

    if missing:
        fail(
            "The input CSV is missing required columns:\n"
            + "\n".join(f"  - {col}" for col in missing)
        )

    # --------------------------------------------------------
    # Validate signature columns
    # --------------------------------------------------------

    print()
    print("Validating signature columns...")

    for col in SIGNATURE_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    missing_tau_values = df[SIGNATURE_COLUMNS].isna().sum()

    print("Missing tau values:")
    for col, count in missing_tau_values.items():
        print(f"  {col}: {int(count):,}")

    # --------------------------------------------------------
    # Recalculate cross-signature statistics
    # --------------------------------------------------------
    #
    # We deliberately calculate these again from the actual
    # signature columns instead of blindly trusting summary
    # columns in the input.
    #
    # This provides an audit check.
    # --------------------------------------------------------

    print()
    print("Recalculating cross-signature statistics for audit...")

    tau_matrix = df[SIGNATURE_COLUMNS]

    recalculated_n_signatures = tau_matrix.notna().sum(axis=1)
    recalculated_n_negative = (tau_matrix < 0).sum(axis=1)
    recalculated_n_strong_negative = (
        tau_matrix <= STRONG_NEGATIVE_THRESHOLD
    ).sum(axis=1)

    recalculated_mean = tau_matrix.mean(axis=1)
    recalculated_median = tau_matrix.median(axis=1)
    recalculated_minimum = tau_matrix.min(axis=1)

    # --------------------------------------------------------
    # Audit original summary values
    # --------------------------------------------------------

    audit_checks = {
        "n_signatures": recalculated_n_signatures,
        "n_negative": recalculated_n_negative,
        "n_strong_negative": recalculated_n_strong_negative,
        "mean_tau": recalculated_mean,
        "median_tau": recalculated_median,
        "minimum_tau": recalculated_minimum,
    }

    audit_failures = []

    for column, recalculated in audit_checks.items():

        original = pd.to_numeric(
            df[column],
            errors="coerce"
        )

        if column in {
            "n_signatures",
            "n_negative",
            "n_strong_negative",
        }:
            matches = original.astype("Int64").eq(
                recalculated.astype("Int64")
            )
        else:
            matches = (
                (original - recalculated).abs() <= 1e-4
            )

        mismatch_count = int((~matches.fillna(False)).sum())

        if mismatch_count > 0:
            audit_failures.append(
                (column, mismatch_count)
            )

    print()
    print("Summary-column audit:")

    if audit_failures:
        print("  WARNING: summary discrepancies detected.")

        for column, count in audit_failures:
            print(
                f"    {column}: {count:,} rows differ"
            )

        print()
        print(
            "The recalculated values will be used for "
            "prioritization."
        )

    else:
        print("  PASS — summary statistics match.")

    # --------------------------------------------------------
    # Replace summary columns with recalculated values
    # --------------------------------------------------------

    df["n_signatures"] = recalculated_n_signatures.astype(int)
    df["n_negative"] = recalculated_n_negative.astype(int)
    df["n_strong_negative"] = (
        recalculated_n_strong_negative.astype(int)
    )

    df["mean_tau"] = recalculated_mean
    df["median_tau"] = recalculated_median
    df["minimum_tau"] = recalculated_minimum

    # --------------------------------------------------------
    # Add negative-signature details
    # --------------------------------------------------------

    print()
    print("Generating signature-level classification...")

    df["strong_negative_threshold"] = STRONG_NEGATIVE_THRESHOLD

    # Number of negative signatures
    #
    # Already calculated above, but explicit Boolean fields
    # make the resulting CSV easier to audit.
    for col in SIGNATURE_COLUMNS:

        safe_name = col.replace("ATLAS_SIG_", "")

        df[f"{safe_name}_negative"] = (
            df[col] < 0
        )

        df[f"{safe_name}_strong_negative"] = (
            df[col] <= STRONG_NEGATIVE_THRESHOLD
        )

    # --------------------------------------------------------
    # Assign tiers
    # --------------------------------------------------------

    print()
    print("Assigning prioritization tiers...")

    df["priority_tier"] = df.apply(
        classify_tier,
        axis=1
    )

    df["priority_tier_number"] = df["priority_tier"].map(
        tier_number
    )

    # --------------------------------------------------------
    # Add a simple consistency score
    # --------------------------------------------------------
    #
    # This is NOT a statistical significance measure.
    #
    # It is only a transparent prioritization score:
    #
    #   +2 for each strongly negative signature
    #   +1 for each additional negative signature
    #
    # Thus:
    #   2 strong negatives = 4
    #   3 strong negatives = 6
    #   1 strong + 2 negative = 4
    #
    # The individual tau values remain the primary evidence.
    # --------------------------------------------------------

    df["cross_signature_consistency_score"] = (
        2 * df["n_strong_negative"]
        + (
            df["n_negative"]
            - df["n_strong_negative"]
        )
    )

    # --------------------------------------------------------
    # Ranking
    # --------------------------------------------------------

    print()
    print("Ranking candidates...")

    df = df.sort_values(
        by=[
            "priority_tier_number",
            "cross_signature_consistency_score",
            "n_strong_negative",
            "n_negative",
            "median_tau",
            "mean_tau",
            "minimum_tau",
        ],
        ascending=[
            True,
            False,
            False,
            False,
            True,
            True,
            True,
        ],
        kind="mergesort",
    ).reset_index(drop=True)

    # New overall priority rank
    df["priority_rank"] = range(1, len(df) + 1)

    # --------------------------------------------------------
    # Reorder important columns
    # --------------------------------------------------------

    preferred_columns = [
        "priority_rank",
        "priority_tier_number",
        "priority_tier",
        "cross_signature_consistency_score",
        "pert_id",
        "pert_iname",
        "ATLAS_SIG_A_TOP150",
        "ATLAS_SIG_B_TOP100",
        "ATLAS_SIG_B_TOP150",
        "n_signatures",
        "n_negative",
        "n_strong_negative",
        "mean_tau",
        "median_tau",
        "minimum_tau",
        "strong_negative_threshold",
        "consensus_rank",
    ]

    boolean_columns = []

    for col in SIGNATURE_COLUMNS:
        safe_name = col.replace("ATLAS_SIG_", "")
        boolean_columns.extend([
            f"{safe_name}_negative",
            f"{safe_name}_strong_negative",
        ])

    preferred_columns.extend(boolean_columns)

    remaining_columns = [
        col for col in df.columns
        if col not in preferred_columns
    ]

    df = df[
        preferred_columns + remaining_columns
    ]

    # --------------------------------------------------------
    # Create output directory
    # --------------------------------------------------------

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Save all prioritized candidates
    # --------------------------------------------------------

    print()
    print("Saving prioritized results:")

    df.to_csv(
        ALL_OUTPUT,
        index=False,
        float_format="%.6f",
    )

    print(f"  ALL:    {ALL_OUTPUT}")

    # --------------------------------------------------------
    # Split into tiers
    # --------------------------------------------------------

    tier1 = df[
        df["priority_tier_number"] == 1
    ].copy()

    tier2 = df[
        df["priority_tier_number"] == 2
    ].copy()

    tier3 = df[
        df["priority_tier_number"] == 3
    ].copy()

    tier4 = df[
        df["priority_tier_number"] == 4
    ].copy()

    tier1.to_csv(
        TIER1_OUTPUT,
        index=False,
        float_format="%.6f",
    )

    tier2.to_csv(
        TIER2_OUTPUT,
        index=False,
        float_format="%.6f",
    )

    tier3.to_csv(
        TIER3_OUTPUT,
        index=False,
        float_format="%.6f",
    )

    tier4.to_csv(
        TIER4_OUTPUT,
        index=False,
        float_format="%.6f",
    )

    print(f"  TIER 1: {TIER1_OUTPUT}")
    print(f"  TIER 2: {TIER2_OUTPUT}")
    print(f"  TIER 3: {TIER3_OUTPUT}")
    print(f"  TIER 4: {TIER4_OUTPUT}")

    # --------------------------------------------------------
    # Top candidates
    # --------------------------------------------------------
    #
    # Keep the top N from the complete prioritized ranking.
    #
    # This is a convenience file for downstream analysis,
    # NOT a claim that these are clinically suitable drugs.
    # --------------------------------------------------------

    top_candidates = df.head(TOP_N).copy()

    top_candidates.to_csv(
        TOP_OUTPUT,
        index=False,
        float_format="%.6f",
    )

    print(f"  TOP {TOP_N}: {TOP_OUTPUT}")

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 60)
    print("CMap PRIORITIZATION SUMMARY")
    print("=" * 60)

    print()
    print(f"Total compounds: {len(df):,}")
    print()

    print("Tier counts:")

    print(
        f"  Tier 1 — Strong consensus: "
        f"{len(tier1):,}"
    )

    print(
        f"  Tier 2 — Consistent reversal: "
        f"{len(tier2):,}"
    )

    print(
        f"  Tier 3 — Partial reversal: "
        f"{len(tier3):,}"
    )

    print(
        f"  Tier 4 — Weak/single-signature: "
        f"{len(tier4):,}"
    )

    # --------------------------------------------------------
    # Strong consensus candidates
    # --------------------------------------------------------

    print()
    print("Strong consensus candidates:")
    print()

    if tier1.empty:
        print("  None.")
    else:
        display_columns = [
            "priority_rank",
            "pert_id",
            "pert_iname",
            "ATLAS_SIG_A_TOP150",
            "ATLAS_SIG_B_TOP100",
            "ATLAS_SIG_B_TOP150",
            "n_negative",
            "n_strong_negative",
            "mean_tau",
            "median_tau",
            "minimum_tau",
        ]

        print(
            tier1[display_columns]
            .head(20)
            .to_string(index=False)
        )

    # --------------------------------------------------------
    # Top 30
    # --------------------------------------------------------

    print()
    print("=" * 60)
    print(f"TOP {TOP_N} PRIORITIZED CMap CANDIDATES")
    print("=" * 60)
    print()

    display_columns = [
        "priority_rank",
        "priority_tier_number",
        "pert_id",
        "pert_iname",
        "ATLAS_SIG_A_TOP150",
        "ATLAS_SIG_B_TOP100",
        "ATLAS_SIG_B_TOP150",
        "n_negative",
        "n_strong_negative",
        "mean_tau",
        "median_tau",
        "minimum_tau",
    ]

    print(
        top_candidates[display_columns].to_string(
            index=False
        )
    )

    # --------------------------------------------------------
    # Interpretation reminder
    # --------------------------------------------------------

    print()
    print("=" * 60)
    print("INTERPRETATION")
    print("=" * 60)

    print(
        "\nNegative tau indicates that a compound's CMap "
        "expression signature tends to oppose the queried "
        "resistance signature."
    )

    print(
        "\nTier 1 candidates show strong reversal in at least "
        "two signatures and therefore receive the highest "
        "cross-signature priority."
    )

    print(
        "\nThis stage does NOT establish drug efficacy, "
        "FDA approval, clinical suitability, or mechanism."
    )

    print(
        "\nFDA/drug-status filtering and molecular docking "
        "should be performed as separate downstream analyses."
    )

    print()
    print("=" * 60)
    print("STAGE 4L CMap PRIORITIZATION COMPLETE")
    print("=" * 60)

    print()
    print("Primary output:")
    print(ALL_OUTPUT)

    print()
    print("Top-candidate output:")
    print(TOP_OUTPUT)

    print()


if __name__ == "__main__":
    main()