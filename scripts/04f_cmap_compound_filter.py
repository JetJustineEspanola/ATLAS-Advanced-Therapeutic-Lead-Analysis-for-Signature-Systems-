from pathlib import Path

import pandas as pd


# ============================================================
# ATLAS — Stage 4F: CMap Compound Filtering
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
    / "TOP100_cmap_connectivity_all.csv"
)

OUTPUT_FILE = (
    PARSED_DIR
    / "TOP100_cmap_compounds.csv"
)

NEGATIVE_OUTPUT_FILE = (
    PARSED_DIR
    / "TOP100_cmap_negative_compounds.csv"
)


# ------------------------------------------------------------
# Load results
# ------------------------------------------------------------

print("=" * 60)
print("ATLAS — CMap Compound Filtering")
print("=" * 60)

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Input file not found:\n{INPUT_FILE}"
    )

results = pd.read_csv(
    INPUT_FILE
)

print(
    f"\nTotal CMap signatures: "
    f"{len(results):,}"
)


# ------------------------------------------------------------
# Inspect perturbation types
# ------------------------------------------------------------

print(
    "\nPerturbation types:"
)

print(
    results["pert_type"]
    .value_counts()
)


# ------------------------------------------------------------
# Keep compounds only
# ------------------------------------------------------------

compounds = results[
    results["pert_type"] == "trt_cp"
].copy()

print(
    "\nCompound perturbation signatures:",
    f"{len(compounds):,}"
)


# ------------------------------------------------------------
# Sort by connectivity
# ------------------------------------------------------------

compounds = compounds.sort_values(
    "connectivity",
    ascending=True,
)


# ------------------------------------------------------------
# Negative compounds
# ------------------------------------------------------------

negative_compounds = compounds[
    compounds["connectivity"] < 0
].copy()

print(
    "Negative compound signatures:",
    f"{len(negative_compounds):,}"
)


# ------------------------------------------------------------
# Show top compounds
# ------------------------------------------------------------

print(
    "\n" + "=" * 60
)

print(
    "TOP 25 NEGATIVELY CONNECTED COMPOUNDS"
)

print(
    "=" * 60
)

print(
    negative_compounds[
        [
            "pert_id",
            "pert_iname",
            "connectivity",
            "cell_id",
            "pert_idose",
            "pert_itime",
        ]
    ]
    .head(25)
    .to_string(index=False)
)


# ------------------------------------------------------------
# Save compound tables
# ------------------------------------------------------------

compounds.to_csv(
    OUTPUT_FILE,
    index=False,
)

negative_compounds.to_csv(
    NEGATIVE_OUTPUT_FILE,
    index=False,
)


print(
    f"\nSaved all compound signatures:\n"
    f"{OUTPUT_FILE}"
)

print(
    f"Saved negative compound signatures:\n"
    f"{NEGATIVE_OUTPUT_FILE}"
)


# ------------------------------------------------------------
# Unique compound summary
# ------------------------------------------------------------

compound_summary = (
    negative_compounds
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
        worst_connectivity=(
            "connectivity",
            "max",
        ),
        n_negative_signatures=(
            "connectivity",
            lambda x: (x < 0).sum(),
        ),
    )
    .reset_index()
    .sort_values(
        "median_connectivity",
        ascending=True,
    )
)


summary_file = (
    PARSED_DIR
    / "TOP100_cmap_compound_summary.csv"
)

compound_summary.to_csv(
    summary_file,
    index=False,
)

print(
    f"\nSaved unique-compound summary:\n"
    f"{summary_file}"
)


# ------------------------------------------------------------
# Complete
# ------------------------------------------------------------

print(
    "\n" + "=" * 60
)

print(
    "CMap COMPOUND FILTERING COMPLETE"
)

print(
    "=" * 60
)