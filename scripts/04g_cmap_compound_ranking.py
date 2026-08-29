from pathlib import Path

import pandas as pd


# ============================================================
# ATLAS — Stage 4G: CMap Compound-Level Ranking
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PARSED_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "parsed"
)

INPUT_FILE = (
    PARSED_DIR
    / "TOP100_cmap_negative_compounds.csv"
)

OUTPUT_FILE = (
    PARSED_DIR
    / "TOP100_cmap_compound_level_ranking.csv"
)


# ------------------------------------------------------------
# Load negative compound signatures
# ------------------------------------------------------------

print("=" * 60)
print("ATLAS — CMap Compound-Level Ranking")
print("=" * 60)

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Input file not found:\n{INPUT_FILE}"
    )

results = pd.read_csv(
    INPUT_FILE
)

print(
    f"\nNegative compound signatures: "
    f"{len(results):,}"
)


# ------------------------------------------------------------
# Basic validation
# ------------------------------------------------------------

required_columns = [
    "pert_id",
    "pert_iname",
    "connectivity",
    "cell_id",
]

missing = [
    column
    for column in required_columns
    if column not in results.columns
]

if missing:
    raise ValueError(
        f"Missing required columns: {missing}"
    )


# ------------------------------------------------------------
# Clean compound names
# ------------------------------------------------------------

results["pert_iname"] = (
    results["pert_iname"]
    .astype("string")
    .str.strip()
)

results["pert_id"] = (
    results["pert_id"]
    .astype("string")
    .str.strip()
)


# ------------------------------------------------------------
# Keep only valid negative scores
# ------------------------------------------------------------

results = results[
    results["connectivity"].notna()
    & (results["connectivity"] < 0)
].copy()


# ============================================================
# Compound-level aggregation
# ============================================================
#
# One compound can have multiple signatures because it may
# have been tested:
#
#   - in different cell lines
#   - at different doses
#   - at different times
#
# We therefore summarize the repeated observations.
# ============================================================

compound_ranking = (
    results
    .groupby(
        [
            "pert_id",
            "pert_iname",
        ],
        dropna=False,
    )
    .agg(
        n_signatures=(
            "connectivity",
            "count",
        ),

        median_connectivity=(
            "connectivity",
            "median",
        ),

        mean_connectivity=(
            "connectivity",
            "mean",
        ),

        best_connectivity=(
            "connectivity",
            "min",
        ),

        weakest_connectivity=(
            "connectivity",
            "max",
        ),

        sd_connectivity=(
            "connectivity",
            "std",
        ),

        n_cells=(
            "cell_id",
            "nunique",
        ),

        cell_lines=(
            "cell_id",
            lambda x: ";".join(
                sorted(
                    set(
                        str(v)
                        for v in x.dropna()
                    )
                )
            ),
        ),

        doses=(
            "pert_idose",
            lambda x: ";".join(
                sorted(
                    set(
                        str(v)
                        for v in x.dropna()
                    )
                )
            ),
        ),

        times=(
            "pert_itime",
            lambda x: ";".join(
                sorted(
                    set(
                        str(v)
                        for v in x.dropna()
                    )
                )
            ),
        ),
    )
    .reset_index()
)


# ------------------------------------------------------------
# Replace undefined SD for single observations
# ------------------------------------------------------------

compound_ranking["sd_connectivity"] = (
    compound_ranking["sd_connectivity"]
    .fillna(0)
)


# ------------------------------------------------------------
# Count how consistently each compound is negative
# ------------------------------------------------------------

negative_counts = (
    results
    .groupby(
        [
            "pert_id",
            "pert_iname",
        ],
        dropna=False,
    )["connectivity"]
    .agg(
        n_negative_signatures="count",
    )
    .reset_index()
)

compound_ranking = compound_ranking.merge(
    negative_counts,
    on=[
        "pert_id",
        "pert_iname",
    ],
    how="left",
)


# ------------------------------------------------------------
# Rank compounds primarily by median connectivity
# ------------------------------------------------------------

compound_ranking = (
    compound_ranking
    .sort_values(
        [
            "median_connectivity",
            "best_connectivity",
        ],
        ascending=[
            True,
            True,
        ],
    )
    .reset_index(drop=True)
)


compound_ranking.insert(
    0,
    "rank",
    range(
        1,
        len(compound_ranking) + 1,
    ),
)


# ============================================================
# Print top compounds
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "TOP 30 COMPOUNDS BY MEDIAN NEGATIVE CONNECTIVITY"
)

print(
    "=" * 60
)

display_columns = [
    "rank",
    "pert_id",
    "pert_iname",
    "n_signatures",
    "n_cells",
    "median_connectivity",
    "mean_connectivity",
    "best_connectivity",
    "weakest_connectivity",
]

print(
    compound_ranking[
        display_columns
    ]
    .head(30)
    .to_string(
        index=False
    )
)


# ============================================================
# Save compound-level ranking
# ============================================================

compound_ranking.to_csv(
    OUTPUT_FILE,
    index=False,
)

print(
    f"\nSaved compound-level ranking:\n"
    f"{OUTPUT_FILE}"
)


# ============================================================
# Create a high-confidence subset
# ============================================================
#
# For now, define a preliminary robust subset as:
#
#   at least 2 independent CMap signatures
#   AND
#   median connectivity <= -0.30
#
# This is NOT the final drug-selection criterion.
# It is only a reproducible screening layer.
# ============================================================

robust_candidates = compound_ranking[
    (compound_ranking["n_signatures"] >= 2)
    & (
        compound_ranking["median_connectivity"]
        <= -0.30
    )
].copy()


robust_file = (
    PARSED_DIR
    / "TOP100_cmap_preliminary_robust_candidates.csv"
)

robust_candidates.to_csv(
    robust_file,
    index=False,
)

print(
    "\nPreliminary robust candidates:",
    len(robust_candidates),
)

print(
    f"Saved:\n{robust_file}"
)


# ============================================================
# Complete
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "CMap COMPOUND-LEVEL RANKING COMPLETE"
)

print(
    "=" * 60
)