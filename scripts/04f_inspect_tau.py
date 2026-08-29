from pathlib import Path

import h5py
import numpy as np


# ============================================================
# Inspect CMap ps_pert_summary.gctx
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

JOB_ROOT = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "raw"
    / "job_6a92e7670262700013c4e5e9"
)

job_dirs = list(
    JOB_ROOT.glob(
        "my_analysis.sig_gutc_tool.*"
    )
)

if not job_dirs:
    raise FileNotFoundError(
        "CMap job directory not found."
    )

JOB_DIR = job_dirs[0]

TAU_FILE = (
    JOB_DIR
    / "matrices"
    / "gutc"
    / "ps_pert_summary.gctx"
)

print("=" * 60)
print("ATLAS — Inspect ps_pert_summary.gctx")
print("=" * 60)

print("\nFile:")
print(TAU_FILE)


def decode(value):
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


with h5py.File(TAU_FILE, "r") as h5:

    print("\nDatasets:")

    def visitor(name, obj):

        if isinstance(
            obj,
            h5py.Dataset,
        ):
            print(
                f"  {name}"
                f" | shape={obj.shape}"
                f" | dtype={obj.dtype}"
            )

    h5.visititems(visitor)

    matrix = h5[
        "/0/DATA/0/matrix"
    ][:]

    print(
        "\nMatrix shape:",
        matrix.shape,
    )

    # ROW metadata
    row_id = [
        decode(x)
        for x in h5[
            "/0/META/ROW/id"
        ][:]
    ]

    row_pert_id = [
        decode(x)
        for x in h5[
            "/0/META/ROW/pert_id"
        ][:]
    ]

    row_pert_iname = [
        decode(x)
        for x in h5[
            "/0/META/ROW/pert_iname"
        ][:]
    ]

    # COL metadata
    col_id = [
        decode(x)
        for x in h5[
            "/0/META/COL/id"
        ][:]
    ]


print("\nROW count:", len(row_id))
print("COL count:", len(col_id))

print("\nFirst 10 ROW IDs:")
print(row_id[:10])

print("\nFirst 10 ROW perturbagen IDs:")
print(row_pert_id[:10])

print("\nFirst 10 ROW perturbagen names:")
print(row_pert_iname[:10])

print("\nCOL IDs:")
print(col_id)

print(
    "\nFirst matrix values:",
    matrix.flatten()[:20],
)

print(
    "\nMatrix minimum:",
    matrix.min(),
)

print(
    "Matrix maximum:",
    matrix.max(),
)