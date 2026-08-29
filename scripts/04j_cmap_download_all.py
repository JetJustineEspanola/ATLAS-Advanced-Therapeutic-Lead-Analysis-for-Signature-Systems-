from pathlib import Path
import json
import os
import tarfile
import time

import requests
from dotenv import load_dotenv


# ============================================================
# ATLAS — Stage 4J: Download All Completed CMap Jobs
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

load_dotenv(
    PROJECT_ROOT / ".env"
)

API_KEY = os.getenv("CLUE_API_KEY")

if not API_KEY:
    raise RuntimeError(
        "CLUE_API_KEY was not found in ATLAS/.env"
    )


CMAP_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
)

JOB_DIR = (
    CMAP_DIR
    / "jobs"
)

RAW_DIR = (
    CMAP_DIR
    / "raw"
)

RAW_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

MANIFEST_FILE = (
    JOB_DIR
    / "cmap_job_manifest.json"
)

if not MANIFEST_FILE.exists():
    raise FileNotFoundError(
        f"Job manifest not found:\n{MANIFEST_FILE}"
    )


# ============================================================
# Load jobs
# ============================================================

with open(
    MANIFEST_FILE,
    "r",
    encoding="utf-8",
) as handle:
    manifest = json.load(handle)


jobs = []

jobs.extend(
    manifest.get(
        "existing_completed_queries",
        [],
    )
)

jobs.extend(
    manifest.get(
        "new_submissions",
        [],
    )
)


# Remove duplicate job IDs
unique_jobs = {}

for job in jobs:

    job_id = job.get("job_id")

    if job_id:
        unique_jobs[job_id] = job

jobs = list(
    unique_jobs.values()
)


# ============================================================
# Helpers
# ============================================================

def already_extracted(job_id):
    """
    Determine whether this job was already extracted.
    """

    job_root = (
        RAW_DIR
        / f"job_{job_id}"
    )

    if not job_root.exists():
        return False

    return any(
        job_root.glob(
            "my_analysis.sig_gutc_tool.*"
        )
    )


def archive_is_valid(archive_path):
    """
    Check whether the archive is a readable tar.gz file.
    """

    if not archive_path.exists():
        return False

    if archive_path.stat().st_size == 0:
        return False

    try:

        with tarfile.open(
            archive_path,
            "r:gz",
        ) as archive:

            archive.getmembers()

        return True

    except (
        tarfile.TarError,
        EOFError,
        OSError,
    ):

        return False


# ============================================================
# Main
# ============================================================

print("=" * 60)
print("ATLAS — Download All Completed CMap Jobs")
print("=" * 60)


download_records = []


for job in jobs:

    signature = job.get(
        "signature",
        "UNKNOWN",
    )

    job_id = job.get(
        "job_id"
    )

    if not job_id:
        continue

    print(
        f"\n{signature}"
    )

    print(
        f"  Job ID: {job_id}"
    )


    # --------------------------------------------------------
    # Already extracted
    # --------------------------------------------------------

    if already_extracted(
        job_id
    ):

        print(
            "  Already extracted — skipping."
        )

        download_records.append(
            {
                "signature": signature,
                "job_id": job_id,
                "status": "already_extracted",
            }
        )

        continue


    archive_path = (
        RAW_DIR
        / f"cmap_{job_id}.tar.gz"
    )


    # --------------------------------------------------------
    # If an earlier download exists, test it first
    # --------------------------------------------------------

    if archive_path.exists():

        size_mb = (
            archive_path.stat().st_size
            / 1024
            / 1024
        )

        print(
            f"  Existing archive: "
            f"{size_mb:.2f} MB"
        )

        if archive_is_valid(
            archive_path
        ):

            print(
                "  Archive is valid."
            )

        else:

            print(
                "  Existing archive is incomplete."
            )

            print(
                "  Removing incomplete archive."
            )

            archive_path.unlink()


    # --------------------------------------------------------
    # Get completed job information
    # --------------------------------------------------------

    status_url = (
        "https://api.clue.io/api/jobs/findByJobId/"
        + job_id
    )

    response = requests.get(
        status_url,
        headers={
            "user_key": API_KEY,
            "Accept": "application/json",
        },
        timeout=60,
    )

    print(
        f"  Status HTTP: "
        f"{response.status_code}"
    )

    if not response.ok:

        print(
            "  Could not retrieve job status."
        )

        download_records.append(
            {
                "signature": signature,
                "job_id": job_id,
                "status": "status_request_failed",
            }
        )

        continue


    data = response.json()

    status = data.get(
        "status"
    )

    download_url = data.get(
        "download_url"
    )


    print(
        f"  Job status: {status}"
    )


    if status != "completed":

        print(
            "  Job is not completed."
        )

        download_records.append(
            {
                "signature": signature,
                "job_id": job_id,
                "status": status,
            }
        )

        continue


    if not download_url:

        print(
            "  No download URL returned."
        )

        download_records.append(
            {
                "signature": signature,
                "job_id": job_id,
                "status": "no_download_url",
            }
        )

        continue


    if download_url.startswith("//"):

        download_url = (
            "https:"
            + download_url
        )


    # --------------------------------------------------------
    # Download
    # --------------------------------------------------------

    print(
        "\n  Starting download..."
    )

    download_response = requests.get(
        download_url,
        stream=True,
        timeout=300,
    )

    print(
        "  Download HTTP:",
        download_response.status_code,
    )

    if not download_response.ok:

        print(
            "  Download failed."
        )

        download_records.append(
            {
                "signature": signature,
                "job_id": job_id,
                "status": "download_failed",
            }
        )

        continue


    total_bytes = (
        int(
            download_response.headers.get(
                "Content-Length",
                0,
            )
        )
        or None
    )

    downloaded_bytes = 0

    start_time = time.time()

    temp_path = (
        RAW_DIR
        / f"cmap_{job_id}.tar.gz.part"
    )

    print(
        f"  Saving to:\n"
        f"  {archive_path}"
    )


    with open(
        temp_path,
        "wb",
    ) as handle:

        for chunk in download_response.iter_content(
            chunk_size=1024 * 1024
        ):

            if not chunk:
                continue

            handle.write(
                chunk
            )

            downloaded_bytes += len(
                chunk
            )

            elapsed = (
                time.time()
                - start_time
            )

            speed_mb = (
                downloaded_bytes
                / 1024
                / 1024
                / max(elapsed, 0.001)
            )

            downloaded_mb = (
                downloaded_bytes
                / 1024
                / 1024
            )

            if total_bytes:

                total_mb = (
                    total_bytes
                    / 1024
                    / 1024
                )

                percent = (
                    downloaded_bytes
                    / total_bytes
                    * 100
                )

                print(
                    f"\r  Progress: "
                    f"{downloaded_mb:8.2f}/"
                    f"{total_mb:8.2f} MB "
                    f"({percent:6.2f}%) "
                    f"{speed_mb:6.2f} MB/s",
                    end="",
                    flush=True,
                )

            else:

                print(
                    f"\r  Downloaded: "
                    f"{downloaded_mb:.2f} MB "
                    f"({speed_mb:.2f} MB/s)",
                    end="",
                    flush=True,
                )


    print()


    # --------------------------------------------------------
    # Replace temporary file
    # --------------------------------------------------------

    temp_path.replace(
        archive_path
    )

    print(
        f"  Download complete: "
        f"{archive_path.stat().st_size:,} bytes"
    )


    # --------------------------------------------------------
    # Validate archive
    # --------------------------------------------------------

    print(
        "  Validating archive..."
    )

    if not archive_is_valid(
        archive_path
    ):

        print(
            "  ERROR: archive validation failed."
        )

        download_records.append(
            {
                "signature": signature,
                "job_id": job_id,
                "status": "invalid_archive",
            }
        )

        continue


    print(
        "  Archive validation passed."
    )


    # --------------------------------------------------------
    # Extract
    # --------------------------------------------------------

    extract_dir = (
        RAW_DIR
        / f"job_{job_id}"
    )

    extract_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        f"  Extracting to:\n"
        f"  {extract_dir}"
    )


    with tarfile.open(
        archive_path,
        "r:gz",
    ) as archive:

        root_path = (
            extract_dir
            .resolve()
        )

        for member in archive.getmembers():

            target = (
                extract_dir
                / member.name
            ).resolve()

            if not str(
                target
            ).startswith(
                str(root_path)
            ):

                raise RuntimeError(
                    "Unsafe path detected."
                )

        archive.extractall(
            extract_dir
        )


    print(
        "  Extraction complete."
    )


    download_records.append(
        {
            "signature": signature,
            "job_id": job_id,
            "status": "downloaded",
            "archive": str(
                archive_path
            ),
            "extract_dir": str(
                extract_dir
            ),
        }
    )


# ============================================================
# Save manifest
# ============================================================

output_file = (
    JOB_DIR
    / "cmap_download_manifest.json"
)

with open(
    output_file,
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        download_records,
        handle,
        indent=2,
    )


# ============================================================
# Summary
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "CMap DOWNLOAD SUMMARY"
)

print(
    "=" * 60
)

for record in download_records:

    print(
        f"{record['signature']}: "
        f"{record['status']}"
    )

print(
    f"\nSaved:\n{output_file}"
)

print(
    "\nStage 4J download complete."
)