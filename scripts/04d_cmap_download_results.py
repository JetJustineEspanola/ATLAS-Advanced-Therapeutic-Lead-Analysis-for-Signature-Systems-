from pathlib import Path
import os
import json
import tarfile

import requests
from dotenv import load_dotenv


# ============================================================
# ATLAS — CMap Result Download
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

load_dotenv(
    dotenv_path=PROJECT_ROOT / ".env"
)

API_KEY = os.getenv("CLUE_API_KEY")

if not API_KEY:
    raise RuntimeError(
        "CLUE_API_KEY was not found in ATLAS/.env"
    )


# ------------------------------------------------------------
# Job information
# ------------------------------------------------------------

JOB_ID = "6a92e7670262700013c4e5e9"

CMAP_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
)

RAW_DIR = (
    CMAP_DIR
    / "raw"
)

RAW_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ------------------------------------------------------------
# Ask CLUE for the completed job
# ------------------------------------------------------------

status_url = (
    "https://api.clue.io/api/jobs/findByJobId/"
    + JOB_ID
)

response = requests.get(
    status_url,
    headers={
        "user_key": API_KEY,
        "Accept": "application/json",
    },
    timeout=60,
)

print("=" * 60)
print("ATLAS — CMap Result Download")
print("=" * 60)

print("\nHTTP status:", response.status_code)

if not response.ok:
    print(response.text[:3000])
    raise RuntimeError(
        "Could not retrieve CMap job information."
    )

job_data = response.json()

status = job_data.get("status")
download_status = job_data.get(
    "download_status"
)
download_url = job_data.get(
    "download_url"
)

print("Job status:", status)
print("Download status:", download_status)


# ------------------------------------------------------------
# Validate job completion
# ------------------------------------------------------------

if status != "completed":
    raise RuntimeError(
        f"CMap job is not completed. "
        f"Current status: {status}"
    )

if not download_url:
    raise RuntimeError(
        "CMap reports completion but did not provide "
        "a download URL."
    )


# CLUE may return a URL beginning with //.
if download_url.startswith("//"):
    download_url = "https:" + download_url

print("\nDownload URL received.")


# ------------------------------------------------------------
# Download archive
# ------------------------------------------------------------

archive_path = (
    RAW_DIR
    / f"cmap_{JOB_ID}.tar.gz"
)

print(
    f"\nDownloading results to:\n"
    f"{archive_path}"
)

download_response = requests.get(
    download_url,
    timeout=300,
    stream=True,
)

print(
    "Download HTTP status:",
    download_response.status_code
)

if not download_response.ok:
    print(
        download_response.text[:2000]
    )

    raise RuntimeError(
        "CMap result download failed."
    )


with open(
    archive_path,
    "wb",
) as handle:

    for chunk in download_response.iter_content(
        chunk_size=1024 * 1024
    ):

        if chunk:
            handle.write(chunk)


print(
    f"Downloaded archive: "
    f"{archive_path.stat().st_size:,} bytes"
)


# ------------------------------------------------------------
# Extract archive
# ------------------------------------------------------------

extract_dir = (
    RAW_DIR
    / f"job_{JOB_ID}"
)

extract_dir.mkdir(
    parents=True,
    exist_ok=True,
)

print(
    f"\nExtracting to:\n"
    f"{extract_dir}"
)


with tarfile.open(
    archive_path,
    "r:gz",
) as archive:

    # Security check against path traversal.
    for member in archive.getmembers():

        target = (
            extract_dir
            / member.name
        ).resolve()

        if not str(target).startswith(
            str(extract_dir.resolve())
        ):
            raise RuntimeError(
                "Unsafe path detected in CMap archive."
            )

    archive.extractall(
        extract_dir
    )


print("Extraction complete.")


# ------------------------------------------------------------
# Save safe job metadata
# ------------------------------------------------------------

safe_metadata = {
    "job_id": JOB_ID,
    "status": status,
    "download_status": download_status,
    "downloaded_archive": str(
        archive_path
    ),
    "extract_directory": str(
        extract_dir
    ),
}

metadata_file = (
    CMAP_DIR
    / "TOP100_completed_job.json"
)

with open(
    metadata_file,
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        safe_metadata,
        handle,
        indent=2,
    )


# ------------------------------------------------------------
# Show important output files
# ------------------------------------------------------------

print("\nCMap output files:")

for path in sorted(
    extract_dir.rglob("*")
):

    if path.is_file():
        print(
            f"  {path.relative_to(extract_dir)}"
        )


print("\n" + "=" * 60)
print("CMap RESULT DOWNLOAD COMPLETE")
print("=" * 60)

print(
    f"\nResults directory:\n"
    f"{extract_dir}"
)