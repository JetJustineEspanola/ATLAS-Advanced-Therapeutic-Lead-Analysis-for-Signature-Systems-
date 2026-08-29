from pathlib import Path
import os
import json

import requests
from dotenv import load_dotenv


# ============================================================
# ATLAS — CMap API Submission Test
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Explicitly load ATLAS/.env
load_dotenv(
    dotenv_path=PROJECT_ROOT / ".env"
)

API_KEY = os.getenv("CLUE_API_KEY")

if not API_KEY:
    raise RuntimeError(
        "CLUE_API_KEY was not found in ATLAS/.env"
    )


# ------------------------------------------------------------
# CMap files
# ------------------------------------------------------------

CMAP_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
)

UP_FILE = (
    CMAP_DIR
    / "ATLAS_SIG_B_TOP100_up.gmt"
)

DOWN_FILE = (
    CMAP_DIR
    / "ATLAS_SIG_B_TOP100_dn.gmt"
)

# ------------------------------------------------------------
# Validate files
# ------------------------------------------------------------

if not UP_FILE.exists():
    raise FileNotFoundError(
        f"UP file not found:\n{UP_FILE}"
    )

if not DOWN_FILE.exists():
    raise FileNotFoundError(
        f"DOWN file not found:\n{DOWN_FILE}"
    )


# ------------------------------------------------------------
# Read GMT files
# ------------------------------------------------------------

with open(
    UP_FILE,
    "r",
    encoding="utf-8",
) as handle:
    up_gmt = handle.read().strip()


with open(
    DOWN_FILE,
    "r",
    encoding="utf-8",
) as handle:
    down_gmt = handle.read().strip()


print("=" * 60)
print("ATLAS — CMap API Submission Test")
print("=" * 60)

print("\nUP file:")
print(UP_FILE)

print("\nDOWN file:")
print(DOWN_FILE)

print("\nUP GMT:")
print(up_gmt[:200])

print("\nDOWN GMT:")
print(down_gmt[:200])


# ------------------------------------------------------------
# Construct JSON payload
# ------------------------------------------------------------

payload = {
    "tool_id": "sig_gutc_tool",
    "name": "ATLAS_BT474_TrastuzumabResistance_TOP100",
    "data_type": "L1000",
    "dataset": "Touchstone",
    "ignoreWarnings": True,
    "uptag-cmapfile": up_gmt,
    "dntag-cmapfile": down_gmt,
}


# ------------------------------------------------------------
# Submit query
# ------------------------------------------------------------

url = "https://api.clue.io/api/jobs"

response = requests.post(
    url,
    headers={
        "user_key": API_KEY,
        "Content-Type": "application/json",
        "Accept": "application/json",
    },
    json=payload,
    timeout=120,
)


# ------------------------------------------------------------
# Print response
# ------------------------------------------------------------

print("\nHTTP status:", response.status_code)

print("\nResponse:")

try:
    response_json = response.json()

    print(
        json.dumps(
            response_json,
            indent=2,
        )
    )

except ValueError:

    print(
        response.text[:5000]
    )

    response_json = None


# ------------------------------------------------------------
# Save successful response
# ------------------------------------------------------------

if response.ok:

    output_file = (
        CMAP_DIR
        / "TOP100_submission_response.json"
    )

    with open(
        output_file,
        "w",
        encoding="utf-8",
    ) as handle:

        json.dump(
            response_json,
            handle,
            indent=2,
        )

    print(
        f"\nSaved response:\n{output_file}"
    )

    print(
        "\nCMap query submitted successfully."
    )

else:

    print(
        "\nCMap query was not submitted successfully."
    )