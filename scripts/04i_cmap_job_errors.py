from pathlib import Path
import json
import os

import requests
from dotenv import load_dotenv


# ============================================================
# ATLAS — Stage 4I: Inspect CMap Job Errors
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


# ============================================================
# Jobs to inspect
# ============================================================

JOBS = [
    {
        "signature": "ATLAS_SIG_B_TOP150",
        "job_id": "6a92f3310262700013c4e5ef",
    },
    {
        "signature": "ATLAS_SIG_A_TOP150",
        "job_id": "6a92f33309168d001496ada5",
    },
]


# ============================================================
# Remove sensitive fields
# ============================================================

def sanitize(data):
    """
    Recursively remove sensitive API/account fields.
    """

    if isinstance(data, dict):

        safe = {}

        for key, value in data.items():

            key_lower = str(key).lower()

            if key_lower in {
                "api_key",
                "user_key",
                "user_id",
                "authorization",
                "token",
            }:
                continue

            safe[key] = sanitize(
                value
            )

        return safe

    if isinstance(data, list):

        return [
            sanitize(item)
            for item in data
        ]

    return data


# ============================================================
# Query CLUE
# ============================================================

for job in JOBS:

    signature = job["signature"]
    job_id = job["job_id"]

    print("\n" + "=" * 60)
    print(signature)
    print("=" * 60)

    print("\nJob ID:")
    print(job_id)

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
        "\nHTTP status:",
        response.status_code,
    )

    try:

        data = response.json()

    except ValueError:

        print(
            "\nNon-JSON response:"
        )

        print(
            response.text[:5000]
        )

        continue


    safe_data = sanitize(
        data
    )

    print(
        "\nSafe CLUE response:"
    )

    print(
        json.dumps(
            safe_data,
            indent=2,
        )
    )