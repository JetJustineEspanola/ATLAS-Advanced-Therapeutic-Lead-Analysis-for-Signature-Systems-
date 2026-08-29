from pathlib import Path

import h5py
import numpy as np
import pandas as pd


# ============================================================
# ATLAS — Stage 4K: Parse All CMap Tau Results
# ============================================================
#
# Combines perturbagen-level tau results from:
#
#   1. ATLAS_SIG_B_TOP100
#   2. ATLAS_SIG_B_TOP150
#   3. ATLAS_SIG_A_TOP150
#
# Each result is read from:
#
#   matrices/gutc/ps_pert_summary.gctx
#
# The output is a common compound-level table that can
# subsequently be used for cross-signature consensus ranking.
# ============================================================


PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)

RAW_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "raw"
)

PARSED_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "parsed"
)

PARSED_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Job definitions
# ============================================================

JOBS = [
    {
        "signature": "ATLAS_SIG_B_TOP100",
        "job_id": "6a92e7670262700013c4e5e9",
    },
    {
        "signature": "ATLAS_SIG_B_TOP150",
        "job_id": "6a92f7a00262700013c4e5f1",
    },
    {
        "signature": "ATLAS_SIG_A_TOP150",
        "job_id": "6a92f7a109168d001496ada7",
    },
]


# ============================================================
# Helpers
# ============================================================

def decode_value(value):
    """Safely decode GCTX byte/string metadata."""

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


def locate_job_directory(job_id):
    """Locate the extracted CMap job directory."""

    root = (
        RAW_DIR
        / f"job_{job_id}"
    )

    if not root.exists():

        raise FileNotFoundError(
            f"Extracted job directory not found:\n{root}"
        )

    candidates = list(
        root.glob(
            "my_analysis.sig_gutc_tool.*"
        )
    )

    if not candidates:

        raise FileNotFoundError(
            f"No extracted CMap analysis found in:\n{root}"
        )

    return candidates[0]


def parse_tau_file(
    signature,
    job_id,
):
    """
    Parse one ps_pert_summary.gctx file.

    Returns one row per perturbagen.
    """

    job_dir = locate_job_directory(
        job_id
    )

    tau_file = (
        job_dir
        / "matrices"
        / "gutc"
        / "ps_pert_summary.gctx"
    )

    if not tau_file.exists():

        raise FileNotFoundError(
            f"Tau file not found:\n{tau_file}"
        )

    print(
        f"\n{'=' * 60}"
    )

    print(
        f"Parsing: {signature}"
    )

    print(
        f"Job ID: {job_id}"
    )

    print(
        f"File: {tau_file}"
    )

    with h5py.File(
        tau_file,
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


    print(
        "Matrix shape:",
        matrix.shape,
    )

    print(
        "Perturbagens:",
        len(row_ids),
    )

    print(
        "Queries:",
        len(col_ids),
    )


    # --------------------------------------------------------
    # Validate
    # --------------------------------------------------------

    if matrix.ndim != 2:

        raise RuntimeError(
            f"Unexpected matrix dimensions: "
            f"{matrix.shape}"
        )

    if matrix.shape[0] != 1:

        raise RuntimeError(
            "Expected one query dimension."
        )

    if matrix.shape[1] != len(row_ids):

        raise RuntimeError(
            "Matrix columns do not match "
            "perturbagen metadata."
        )

    if len(col_ids) != 1:

        raise RuntimeError(
            "Expected exactly one ATLAS query."
        )


    query_name = col_ids[0]

    tau_values = matrix[0, :]


    # --------------------------------------------------------
    # Build dataframe
    # --------------------------------------------------------

    results = pd.DataFrame(
        {
            "signature": signature,
            "job_id": job_id,
            "query": query_name,
            "pert_id": pert_ids,
            "pert_iname": pert_inames,
            "pert_type": pert_types,
            "tau": tau_values,
        }
    )


    # --------------------------------------------------------
    # Convert tau to numeric
    # --------------------------------------------------------

    results["tau"] = pd.to_numeric(
        results["tau"],
        errors="coerce",
    )


    # --------------------------------------------------------
    # Basic QC
    # --------------------------------------------------------

    print(
        "Missing tau:",
        results["tau"].isna().sum(),
    )

    finite_tau = results[
        results["tau"].notna()
    ]["tau"]

    if not finite_tau.empty:

        print(
            "Tau range:",
            finite_tau.min(),
            "to",
            finite_tau.max(),
        )


    # --------------------------------------------------------
    # Keep compounds
    # --------------------------------------------------------

    compounds = results[
        results["pert_type"] == "trt_cp"
    ].copy()

    print(
        "Compound perturbagens:",
        f"{len(compounds):,}",
    )


    return results, compounds


# ============================================================
# Parse all signatures
# ============================================================

all_results = []
all_compounds = []


for job in JOBS:

    results, compounds = parse_tau_file(
        job["signature"],
        job["job_id"],
    )

    all_results.append(
        results
    )

    all_compounds.append(
        compounds
    )


# ============================================================
# Combine full perturbagen results
# ============================================================

all_tau = pd.concat(
    all_results,
    ignore_index=True,
)

all_compounds_df = pd.concat(
    all_compounds,
    ignore_index=True,
)


print(
    "\n" + "=" * 60
)

print(
    "COMBINED CMap RESULTS"
)

print(
    "=" * 60
)

print(
    "Total perturbagen rows:",
    f"{len(all_tau):,}",
)

print(
    "Total compound rows:",
    f"{len(all_compounds_df):,}",
)


# ============================================================
# Save individual long-format results
# ============================================================

all_tau_file = (
    PARSED_DIR
    / "ATLAS_CMap_all_tau_long.csv"
)

all_tau.to_csv(
    all_tau_file,
    index=False,
)

print(
    f"\nSaved all tau results:\n"
    f"{all_tau_file}"
)


compound_long_file = (
    PARSED_DIR
    / "ATLAS_CMap_compound_tau_long.csv"
)

all_compounds_df.to_csv(
    compound_long_file,
    index=False,
)

print(
    f"Saved compound tau results:\n"
    f"{compound_long_file}"
)


# ============================================================
# Create cross-signature matrix
# ============================================================
#
# Rows:
#   compounds
#
# Columns:
#   one tau per ATLAS signature
#
# If a compound does not occur in a given signature,
# its value remains NaN.
# ============================================================

cross_signature = (
    all_compounds_df[
        [
            "pert_id",
            "pert_iname",
            "signature",
            "tau",
        ]
    ]
    .drop_duplicates(
        subset=[
            "pert_id",
            "signature",
        ],
        keep="first",
    )
    .pivot(
        index=[
            "pert_id",
            "pert_iname",
        ],
        columns="signature",
        values="tau",
    )
    .reset_index()
)


# Remove pandas column-axis name
cross_signature.columns.name = None


# ============================================================
# Ensure expected columns exist
# ============================================================

expected_signatures = [
    "ATLAS_SIG_B_TOP100",
    "ATLAS_SIG_B_TOP150",
    "ATLAS_SIG_A_TOP150",
]

for signature in expected_signatures:

    if signature not in cross_signature.columns:

        cross_signature[
            signature
        ] = np.nan


# ============================================================
# Calculate cross-signature metrics
# ============================================================

tau_columns = expected_signatures


cross_signature[
    "n_signatures"
] = (
    cross_signature[
        tau_columns
    ]
    .notna()
    .sum(axis=1)
)


cross_signature[
    "n_negative"
] = (
    cross_signature[
        tau_columns
    ]
    .lt(0)
    .sum(axis=1)
)


cross_signature[
    "n_strong_negative"
] = (
    cross_signature[
        tau_columns
    ]
    .le(-90)
    .sum(axis=1)
)


cross_signature[
    "mean_tau"
] = (
    cross_signature[
        tau_columns
    ]
    .mean(
        axis=1,
        skipna=True,
    )
)


cross_signature[
    "median_tau"
] = (
    cross_signature[
        tau_columns
    ].median(
        axis=1,
        skipna=True,
    )
)


cross_signature[
    "minimum_tau"
] = (
    cross_signature[
        tau_columns
    ].min(
        axis=1,
        skipna=True,
    )
)


# ============================================================
# Rank by cross-signature consistency
# ============================================================
#
# Priority:
#
#   1. Number of strong negative signatures
#   2. Number of negative signatures
#   3. Median tau
#   4. Mean tau
#
# This is a screening rank, not a clinical score.
# ============================================================

cross_signature = (
    cross_signature
    .sort_values(
        [
            "n_strong_negative",
            "n_negative",
            "median_tau",
            "mean_tau",
        ],
        ascending=[
            False,
            False,
            True,
            True,
        ],
    )
    .reset_index(
        drop=True
    )
)


cross_signature.insert(
    0,
    "consensus_rank",
    range(
        1,
        len(cross_signature) + 1,
    ),
)


# ============================================================
# Save cross-signature table
# ============================================================

cross_signature_file = (
    PARSED_DIR
    / "ATLAS_CMap_cross_signature_tau.csv"
)

cross_signature.to_csv(
    cross_signature_file,
    index=False,
)

print(
    f"\nSaved cross-signature tau matrix:\n"
    f"{cross_signature_file}"
)


# ============================================================
# Strong consensus candidates
# ============================================================
#
# Preliminary criterion:
#
#   present in >= 2 signatures
#   AND
#   negative tau in >= 2 signatures
#
# This is intentionally a screening criterion.
# ============================================================

consensus_candidates = cross_signature[
    (
        cross_signature[
            "n_signatures"
        ] >= 2
    )
    &
    (
        cross_signature[
            "n_negative"
        ] >= 2
    )
].copy()


consensus_file = (
    PARSED_DIR
    / "ATLAS_CMap_consensus_candidates.csv"
)

consensus_candidates.to_csv(
    consensus_file,
    index=False,
)

print(
    "\nPreliminary consensus candidates:",
    f"{len(consensus_candidates):,}",
)

print(
    f"Saved:\n{consensus_file}"
)


# ============================================================
# Strongest consensus candidates
# ============================================================

strong_consensus = cross_signature[
    (
        cross_signature[
            "n_strong_negative"
        ] >= 2
    )
].copy()


strong_consensus_file = (
    PARSED_DIR
    / "ATLAS_CMap_strong_consensus.csv"
)

strong_consensus.to_csv(
    strong_consensus_file,
    index=False,
)

print(
    "\nCandidates with strong negative tau "
    "in >=2 signatures:",
    f"{len(strong_consensus):,}",
)

print(
    f"Saved:\n{strong_consensus_file}"
)


# ============================================================
# Display top consensus candidates
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "TOP 30 CROSS-SIGNATURE CMap CANDIDATES"
)

print(
    "=" * 60
)

display_columns = [
    "consensus_rank",
    "pert_id",
    "pert_iname",
    "ATLAS_SIG_B_TOP100",
    "ATLAS_SIG_B_TOP150",
    "ATLAS_SIG_A_TOP150",
    "n_signatures",
    "n_negative",
    "n_strong_negative",
    "mean_tau",
    "median_tau",
    "minimum_tau",
]

print(
    cross_signature[
        display_columns
    ]
    .head(30)
    .to_string(
        index=False
    )
)


# ============================================================
# Summary
# ============================================================

summary = pd.DataFrame(
    {
        "signature": expected_signatures,
        "compound_rows": [
            len(
                all_compounds_df[
                    all_compounds_df[
                        "signature"
                    ] == signature
                ]
            )
            for signature in expected_signatures
        ],
    }
)


summary_file = (
    PARSED_DIR
    / "ATLAS_CMap_cross_signature_summary.csv"
)

summary.to_csv(
    summary_file,
    index=False,
)


# ============================================================
# Complete
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "STAGE 4K CMap CROSS-SIGNATURE PARSING COMPLETE"
)

print(
    "=" * 60
)

print(
    "\nPrimary output:"
)

print(
    cross_signature_file
)

print(
    "\nConsensus output:"
)

print(
    consensus_file
)