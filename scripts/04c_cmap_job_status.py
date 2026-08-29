from pathlib import Path
import os
import json

import requests
from dotenv import load_dotenv


# ============================================================
# ATLAS — CMap Job Status
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
# CMap job
# ------------------------------------------------------------

JOB_ID = "6a92e7670262700013c4e5e9"

url = (
    "https://api.clue.io/api/jobs/findByJobId/"
    + JOB_ID
)


# ------------------------------------------------------------
# Request job status
# ------------------------------------------------------------

response = requests.get(
    url,
    headers={
        "user_key": API_KEY,
        "Accept": "application/json",
    },
    timeout=60,
)


print("=" * 60)
print("ATLAS — CMap Job Status")
print("=" * 60)

print("\nJob ID:")
print(JOB_ID)

print("\nHTTP status:")
print(response.status_code)


# ------------------------------------------------------------
# Parse response
# ------------------------------------------------------------

try:

    data = response.json()

except ValueError:

    print("\nNon-JSON response:")
    print(response.text[:3000])

    raise SystemExit(
        "Could not parse CMap response."
    )


# ------------------------------------------------------------
# Display safely
# ------------------------------------------------------------

# Do NOT print the entire response because CLUE may
# include account/API information.

safe_output = {
    "status": data.get("status"),
    "job_id": data.get("job_id"),
    "message": data.get("message"),
}

print("\nCMap status response:")
print(
    json.dumps(
        safe_output,
        indent=2,
    )
)


# ------------------------------------------------------------
# Save safe status record
# ------------------------------------------------------------

CMAP_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
)

CMAP_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

status_file = (
    CMAP_DIR
    / "TOP100_job_status.json"
)

with open(
    status_file,
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        safe_output,
        handle,
        indent=2,
    )

print(
    f"\nSaved safe status record:\n"
    f"{status_file}"
)