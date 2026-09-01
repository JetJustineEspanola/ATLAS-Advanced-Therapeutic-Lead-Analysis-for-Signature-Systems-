from pathlib import Path

import h5py
import numpy as np
import pandas as pd


# ============================================================
# ATLAS — Stage 4E: Parse CMap Connectivity Results
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

CMAP_ROOT = (
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
)

PARSED_DIR = (
    RESULTS_DIR
    / "parsed"
)

PARSED_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Locate CMap job directory
# ============================================================

job_directories = list(
    CMAP_ROOT.glob(
        "my_analysis.sig_gutc_tool.*"
    )
)

if not job_directories:
    raise FileNotFoundError(
        "Could not locate extracted CMap job directory."
    )

JOB_DIR = job_directories[0]

CS_FILE = (
    JOB_DIR
    / "matrices"
    / "gutc"
    / "cs_sig.gctx"
)

if not CS_FILE.exists():
    raise FileNotFoundError(
        f"Connectivity matrix not found:\n{CS_FILE}"
    )


print("=" * 60)
print("ATLAS — CMap Connectivity Result Parser")
print("=" * 60)

print("\nJob directory:")
print(JOB_DIR)

print("\nConnectivity file:")
print(CS_FILE)


# ============================================================
# Helper functions
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


def read_gctx_vector(
    h5file,
    dataset_path,
):
    """Read and decode a one-dimensional GCTX metadata dataset."""

    values = h5file[dataset_path][:]

    return [
        decode_value(value)
        for value in values
    ]


# ============================================================
# Read GCTX
# ============================================================

print("\nOpening GCTX...")

with h5py.File(
    CS_FILE,
    "r",
) as h5:

    matrix = h5[
        "/0/DATA/0/matrix"
    ][:]

    print(
        "\nMatrix shape:",
        matrix.shape,
    )

    # --------------------------------------------------------
    # Query ID
    # --------------------------------------------------------

    query_ids = read_gctx_vector(
        h5,
        "/0/META/COL/id",
    )

    # --------------------------------------------------------
    # Perturbation signature IDs
    # --------------------------------------------------------

    signature_ids = read_gctx_vector(
        h5,
        "/0/META/ROW/id",
    )

    # --------------------------------------------------------
    # Perturbagen metadata
    # --------------------------------------------------------

    pert_id = read_gctx_vector(
        h5,
        "/0/META/ROW/pert_id",
    )

    pert_iname = read_gctx_vector(
        h5,
        "/0/META/ROW/pert_iname",
    )

    pert_type = read_gctx_vector(
        h5,
        "/0/META/ROW/pert_type",
    )

    cell_id = read_gctx_vector(
        h5,
        "/0/META/ROW/cell_id",
    )

    pert_idose = read_gctx_vector(
        h5,
        "/0/META/ROW/pert_idose",
    )

    pert_itime = read_gctx_vector(
        h5,
        "/0/META/ROW/pert_itime",
    )

    is_touchstone = read_gctx_vector(
        h5,
        "/0/META/ROW/is_touchstone",
    )

    is_exemplar = read_gctx_vector(
        h5,
        "/0/META/ROW/is_exemplar",
    )


# ============================================================
# Validate matrix orientation
# ============================================================

if matrix.ndim != 2:
    raise RuntimeError(
        f"Expected a 2D matrix, got {matrix.ndim}D."
    )

if matrix.shape[0] != len(query_ids):
    raise RuntimeError(
        "Number of matrix rows does not match "
        "number of query IDs."
    )

if matrix.shape[1] != len(signature_ids):
    raise RuntimeError(
        "Number of matrix columns does not match "
        "number of perturbation signatures."
    )


# ============================================================
# This CMap result has:
#
#   rows    = queries
#   columns = perturbation signatures
#
# We have one query, so extract that row.
# ============================================================

query_id = query_ids[0]

connectivity_scores = matrix[0, :]

print(
    "\nQuery:",
    query_id,
)

print(
    "Number of perturbation signatures:",
    len(connectivity_scores),
)


# ============================================================
# Build perturbagen result table
# ============================================================

results = pd.DataFrame(
    {
        "query": query_id,
        "signature_id": signature_ids,
        "pert_id": pert_id,
        "pert_iname": pert_iname,
        "pert_type": pert_type,
        "cell_id": cell_id,
        "pert_idose": pert_idose,
        "pert_itime": pert_itime,
        "is_touchstone": is_touchstone,
        "is_exemplar": is_exemplar,
        "connectivity": connectivity_scores,
    }
)


# ============================================================
# Validate values
# ============================================================

print(
    "\nConnectivity statistics:"
)

print(
    results["connectivity"].describe()
)


missing_scores = (
    results["connectivity"]
    .isna()
    .sum()
)

print(
    "\nMissing connectivity scores:",
    missing_scores,
)


# ============================================================
# Sort by negative connectivity
# ============================================================

negative_results = (
    results[
        results["connectivity"] < 0
    ]
    .sort_values(
        "connectivity",
        ascending=True,
    )
    .copy()
)


# ============================================================
# Print top candidates
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "TOP 25 MOST NEGATIVELY CONNECTED PERTURBATIONS"
)

print(
    "=" * 60
)

display_columns = [
    "pert_id",
    "pert_iname",
    "pert_type",
    "cell_id",
    "connectivity",
]

print(
    negative_results[
        display_columns
    ]
    .head(25)
    .to_string(
        index=False
    )
)


# ============================================================
# Save complete results
# ============================================================

full_file = (
    PARSED_DIR
    / "TOP100_cmap_connectivity_all.csv"
)

results.to_csv(
    full_file,
    index=False,
)

print(
    f"\nSaved full connectivity results:\n"
    f"{full_file}"
)


# ============================================================
# Save negative-connectivity results
# ============================================================

negative_file = (
    PARSED_DIR
    / "TOP100_cmap_negative_connectivity.csv"
)

negative_results.to_csv(
    negative_file,
    index=False,
)

print(
    f"Saved negative-connectivity results:\n"
    f"{negative_file}"
)


# ============================================================
# Save top 100 candidates
# ============================================================

top100_file = (
    PARSED_DIR
    / "TOP100_cmap_top_negative_candidates.csv"
)

negative_results.head(100).to_csv(
    top100_file,
    index=False,
)

print(
    f"Saved top 100 candidates:\n"
    f"{top100_file}"
)


# ============================================================
# Summary
# ============================================================

summary = pd.DataFrame(
    {
        "query": [query_id],
        "total_perturbation_signatures": [
            len(results)
        ],
        "negative_connectivity_signatures": [
            len(negative_results)
        ],
        "most_negative_score": [
            negative_results[
                "connectivity"
            ].min()
        ],
    }
)

summary_file = (
    PARSED_DIR
    / "TOP100_cmap_summary.csv"
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
    "CMap CONNECTIVITY PARSING COMPLETE"
)

print(
    "=" * 60
)