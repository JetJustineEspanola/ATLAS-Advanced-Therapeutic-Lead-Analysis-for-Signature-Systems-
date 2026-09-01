from pathlib import Path

import h5py
import numpy as np
import pandas as pd


# ============================================================
# ATLAS — Stage 4F: Parse CMap Tau Scores
# ============================================================
#
# ps_pert_summary.gctx structure for this query:
#
#   matrix = 1 query × 8,798 perturbagens
#   ROW    = perturbagens
#   COL    = query
#
# Therefore:
#
#   matrix[0, :]
#       corresponds to the 8,798 perturbagens.
#
# Tau is the background-adjusted CMap score in [-100, +100].
# Negative tau = opposing transcriptional relationship.
# ============================================================


PROJECT_ROOT = Path(__file__).resolve().parents[1]

JOB_ROOT = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "raw"
    / "job_6a92e7670262700013c4e5e9"
)

RESULTS_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "parsed"
)

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Locate extracted CMap job
# ============================================================

job_dirs = list(
    JOB_ROOT.glob(
        "my_analysis.sig_gutc_tool.*"
    )
)

if not job_dirs:
    raise FileNotFoundError(
        f"Could not find extracted CMap job in:\n{JOB_ROOT}"
    )

JOB_DIR = job_dirs[0]


# ============================================================
# Locate tau file
# ============================================================

TAU_FILE = (
    JOB_DIR
    / "matrices"
    / "gutc"
    / "ps_pert_summary.gctx"
)

if not TAU_FILE.exists():
    raise FileNotFoundError(
        f"Tau file not found:\n{TAU_FILE}"
    )


print("=" * 60)
print("ATLAS — CMap Tau Result Parser")
print("=" * 60)

print("\nTau file:")
print(TAU_FILE)


# ============================================================
# Helper
# ============================================================

def decode_value(value):
    """Safely decode GCTX metadata values."""

    if isinstance(value, bytes):
        return value.decode(
            "utf-8",
            errors="replace",
        )

    if isinstance(value, np.bytes_):
        return bytes(value).decode(
            "utf-8",
            errors="replace",
        )

    return str(value)


# ============================================================
# Read GCTX
# ============================================================

print("\nOpening GCTX...")

with h5py.File(
    TAU_FILE,
    "r",
) as h5:

    matrix = h5[
        "/0/DATA/0/matrix"
    ][:]

    row_ids = [
        decode_value(value)
        for value in h5[
            "/0/META/ROW/id"
        ][:]
    ]

    col_ids = [
        decode_value(value)
        for value in h5[
            "/0/META/COL/id"
        ][:]
    ]

    pert_ids = [
        decode_value(value)
        for value in h5[
            "/0/META/ROW/pert_id"
        ][:]
    ]

    pert_inames = [
        decode_value(value)
        for value in h5[
            "/0/META/ROW/pert_iname"
        ][:]
    ]

    pert_types = [
        decode_value(value)
        for value in h5[
            "/0/META/ROW/pert_type"
        ][:]
    ]


# ============================================================
# Inspect structure
# ============================================================

print(
    "\nMatrix shape:",
    matrix.shape,
)

print(
    "Number of perturbagens:",
    len(row_ids),
)

print(
    "Number of queries:",
    len(col_ids),
)


# ============================================================
# Validate expected structure
# ============================================================

if matrix.ndim != 2:
    raise RuntimeError(
        f"Expected a 2D matrix, got {matrix.ndim}D."
    )

if matrix.shape[0] != 1:
    raise RuntimeError(
        "This parser expects one query row in this CMap result."
    )

if matrix.shape[1] != len(row_ids):
    raise RuntimeError(
        "Matrix column count does not match perturbagen metadata."
    )

if len(row_ids) != len(pert_ids):
    raise RuntimeError(
        "Perturbagen ID count does not match row count."
    )

if len(row_ids) != len(pert_inames):
    raise RuntimeError(
        "Perturbagen name count does not match row count."
    )

if len(row_ids) != len(pert_types):
    raise RuntimeError(
        "Perturbagen type count does not match row count."
    )

if len(col_ids) != 1:
    raise RuntimeError(
        f"Expected one query, found {len(col_ids)}."
    )


# ============================================================
# Extract tau values
# ============================================================

query_name = col_ids[0]

tau_values = matrix[0, :]

print(
    "\nQuery:",
    query_name,
)

print(
    "Tau values extracted:",
    len(tau_values),
)


# ============================================================
# Build result table
# ============================================================

tau_results = pd.DataFrame(
    {
        "query": query_name,
        "pert_id": pert_ids,
        "pert_iname": pert_inames,
        "pert_type": pert_types,
        "tau": tau_values,
    }
)


# ============================================================
# Basic QC
# ============================================================

print(
    "\nTau statistics:"
)

print(
    tau_results["tau"].describe()
)

missing_tau = (
    tau_results["tau"]
    .isna()
    .sum()
)

print(
    "\nMissing tau values:",
    missing_tau,
)


# ============================================================
# Check tau range
# ============================================================

finite_tau = tau_results[
    tau_results["tau"].notna()
]["tau"]

if not finite_tau.empty:

    tau_min = finite_tau.min()
    tau_max = finite_tau.max()

    print(
        "\nObserved tau range:",
        tau_min,
        "to",
        tau_max,
    )

    if tau_min < -100 or tau_max > 100:
        print(
            "WARNING: observed tau values fall outside "
            "the expected [-100, +100] range."
        )


# ============================================================
# Keep compounds
# ============================================================

compound_tau = tau_results[
    tau_results["pert_type"] == "trt_cp"
].copy()

print(
    "\nCompound perturbagens:",
    f"{len(compound_tau):,}",
)


# ============================================================
# Negative tau compounds
# ============================================================

negative_tau = compound_tau[
    compound_tau["tau"] < 0
].copy()

negative_tau = negative_tau.sort_values(
    "tau",
    ascending=True,
)

print(
    "Negative tau compound signatures:",
    f"{len(negative_tau):,}",
)


# ============================================================
# Strong negative tau
# ============================================================
#
# CLUE documentation:
#   tau <= -90
#   is generally considered a strong score for
#   further investigation.
# ============================================================

strong_negative_tau = compound_tau[
    compound_tau["tau"] <= -90
].copy()

strong_negative_tau = (
    strong_negative_tau
    .sort_values(
        "tau",
        ascending=True,
    )
)

print(
    "Strong negative tau compounds (tau <= -90):",
    f"{len(strong_negative_tau):,}",
)


# ============================================================
# Print top compounds
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "TOP 30 COMPOUNDS BY MOST NEGATIVE TAU"
)

print(
    "=" * 60
)

display_columns = [
    "pert_id",
    "pert_iname",
    "pert_type",
    "tau",
]

print(
    negative_tau[
        display_columns
    ]
    .head(30)
    .to_string(
        index=False
    )
)


# ============================================================
# Print strong candidates
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "STRONG NEGATIVE TAU COMPOUNDS (tau <= -90)"
)

print(
    "=" * 60
)

if strong_negative_tau.empty:

    print(
        "No compound perturbagens reached tau <= -90."
    )

else:

    print(
        strong_negative_tau[
            display_columns
        ]
        .head(30)
        .to_string(
            index=False
        )
    )


# ============================================================
# Save all perturbagen tau results
# ============================================================

all_tau_file = (
    RESULTS_DIR
    / "TOP100_cmap_tau_all_perturbagens.csv"
)

tau_results.to_csv(
    all_tau_file,
    index=False,
)

print(
    f"\nSaved all tau results:\n"
    f"{all_tau_file}"
)


# ============================================================
# Save compound tau results
# ============================================================

compound_tau_file = (
    RESULTS_DIR
    / "TOP100_cmap_tau_compounds.csv"
)

compound_tau.to_csv(
    compound_tau_file,
    index=False,
)

print(
    f"Saved compound tau results:\n"
    f"{compound_tau_file}"
)


# ============================================================
# Save negative tau compounds
# ============================================================

negative_tau_file = (
    RESULTS_DIR
    / "TOP100_cmap_negative_tau_compounds.csv"
)

negative_tau.to_csv(
    negative_tau_file,
    index=False,
)

print(
    f"Saved negative tau compounds:\n"
    f"{negative_tau_file}"
)


# ============================================================
# Save strong candidates
# ============================================================

strong_tau_file = (
    RESULTS_DIR
    / "TOP100_cmap_strong_negative_tau.csv"
)

strong_negative_tau.to_csv(
    strong_tau_file,
    index=False,
)

print(
    f"Saved strong negative tau candidates:\n"
    f"{strong_tau_file}"
)


# ============================================================
# Summary
# ============================================================

summary = pd.DataFrame(
    {
        "query": [query_name],
        "total_perturbagens": [
            len(tau_results)
        ],
        "compound_perturbagens": [
            len(compound_tau)
        ],
        "negative_tau_compounds": [
            len(negative_tau)
        ],
        "strong_negative_tau_compounds": [
            len(strong_negative_tau)
        ],
        "minimum_tau": [
            tau_results["tau"].min()
        ],
        "maximum_tau": [
            tau_results["tau"].max()
        ],
    }
)

summary_file = (
    RESULTS_DIR
    / "TOP100_cmap_tau_summary.csv"
)

summary.to_csv(
    summary_file,
    index=False,
)

print(
    f"Saved summary:\n"
    f"{summary_file}"
)


# ============================================================
# Complete
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "CMap TAU PARSING COMPLETE"
)

print(
    "=" * 60
)