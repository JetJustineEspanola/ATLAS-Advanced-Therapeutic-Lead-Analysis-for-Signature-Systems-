from pathlib import Path
import json
import os

import requests
from dotenv import load_dotenv


# ============================================================
# ATLAS — Stage 4H: Check All CMap Jobs
# ============================================================

PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)

load_dotenv(
    dotenv_path=PROJECT_ROOT / ".env"
)

API_KEY = os.getenv(
    "CLUE_API_KEY"
)

if not API_KEY:
    raise RuntimeError(
        "CLUE_API_KEY was not found in ATLAS/.env"
    )


MANIFEST_FILE = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "jobs"
    / "cmap_job_manifest.json"
)


if not MANIFEST_FILE.exists():
    raise FileNotFoundError(
        f"Job manifest not found:\n{MANIFEST_FILE}"
)


with open(
    MANIFEST_FILE,
    "r",
    encoding="utf-8",
) as handle:

    manifest = json.load(handle)


# ============================================================
# Collect jobs
# ============================================================

jobs = []

for job in manifest.get(
    "existing_completed_queries",
    []
):

    jobs.append(job)


for job in manifest.get(
    "new_submissions",
    []
):

    jobs.append(job)


# Remove entries without job IDs
jobs = [
    job
    for job in jobs
    if job.get("job_id")
]


# ============================================================
# Check status
# ============================================================

print("=" * 60)
print("ATLAS — CMap Job Status")
print("=" * 60)

updated_jobs = []


for job in jobs:

    signature = job.get(
        "signature",
        "UNKNOWN",
    )

    job_id = job.get(
        "job_id"
    )

    url = (
        "https://api.clue.io/api/jobs/findByJobId/"
        + job_id
    )

    response = requests.get(
        url,
        headers={
            "user_key": API_KEY,
            "Accept": "application/json",
        },
        timeout=60,
    )

    print(
        f"\n{signature}"
    )

    print(
        f"  Job ID: {job_id}"
    )

    print(
        f"  HTTP status: {response.status_code}"
    )

    try:

        data = response.json()

    except ValueError:

        print(
            "  Could not parse JSON response."
        )

        print(
            response.text[:1000]
        )

        updated_jobs.append(
            {
                **job,
                "http_status": response.status_code,
                "status": None,
            }
        )

        continue


    status = data.get(
        "status"
    )

    download_status = data.get(
        "download_status"
    )

    print(
        f"  Status: {status}"
    )

    if download_status is not None:

        print(
            f"  Download status: "
            f"{download_status}"
        )


    updated_jobs.append(
        {
            **job,
            "http_status": response.status_code,
            "status": status,
            "download_status": download_status,
        }
    )


# ============================================================
# Save updated status manifest
# ============================================================

status_manifest = {
    "jobs": updated_jobs
}


STATUS_FILE = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "jobs"
    / "cmap_job_status_all.json"
)


with open(
    STATUS_FILE,
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        status_manifest,
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
    "CMap JOB SUMMARY"
)

print(
    "=" * 60
)


for job in updated_jobs:

    print(
        f"{job.get('signature')}: "
        f"{job.get('status')}"
    )


completed = sum(
    1
    for job in updated_jobs
    if job.get("status") == "completed"
)

active = sum(
    1
    for job in updated_jobs
    if job.get("status") in {
        "pending",
        "submitted",
        "running",
        "queued",
    }
)

failed = sum(
    1
    for job in updated_jobs
    if job.get("status") in {
        "failed",
        "error",
        "cancelled",
    }
)

print(
    f"\nCompleted: {completed}"
)

print(
    f"Active:    {active}"
)

print(
    f"Failed:    {failed}"
)

print(
    f"\nSaved:\n{STATUS_FILE}"
)