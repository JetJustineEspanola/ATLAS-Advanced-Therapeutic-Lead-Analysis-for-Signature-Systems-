from pathlib import Path
import json
import os

import requests
from dotenv import load_dotenv


# ============================================================
# ATLAS — Stage 4G: CMap Multi-Signature Submission
# ============================================================
#
# Existing completed query:
#
#   ATLAS_SIG_B_TOP100
#   Job ID: 6a92e7670262700013c4e5e9
#
# This script submits ONLY the new Entrez-based signatures:
#
#   ATLAS_SIG_B_TOP150
#   ATLAS_SIG_A_TOP150
#
# The old failed TOP150 jobs are NOT reused.
#
# ============================================================


# ============================================================
# 1. Paths
# ============================================================

PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)

load_dotenv(
    dotenv_path=PROJECT_ROOT / ".env"
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

JOB_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 2. API key
# ============================================================

API_KEY = os.getenv(
    "CLUE_API_KEY"
)

if not API_KEY:
    raise RuntimeError(
        "CLUE_API_KEY was not found in ATLAS/.env"
    )


# ============================================================
# 3. Signatures to submit
# ============================================================

SIGNATURES = [
    {
        "name": "ATLAS_SIG_B_TOP150",
        "up_file": (
            "ATLAS_SIG_B_TOP150_up_entrez.gmt"
        ),
        "down_file": (
            "ATLAS_SIG_B_TOP150_dn_entrez.gmt"
        ),
        "role": "secondary robustness",
    },
    {
        "name": "ATLAS_SIG_A_TOP150",
        "up_file": (
            "ATLAS_SIG_A_TOP150_up_entrez.gmt"
        ),
        "down_file": (
            "ATLAS_SIG_A_TOP150_dn_entrez.gmt"
        ),
        "role": "alternative robustness",
    },
]


# ============================================================
# 4. Read GMT
# ============================================================

def read_gmt(
    filepath: Path,
):
    """
    Read a one-line GMT file.
    """

    with open(
        filepath,
        "r",
        encoding="utf-8",
    ) as handle:

        line = (
            handle
            .readline()
            .rstrip("\n")
        )

    fields = line.split("\t")

    if len(fields) < 3:
        raise ValueError(
            f"Invalid GMT file:\n{filepath}"
        )

    name = fields[0]

    description = fields[1]

    genes = [
        gene.strip()
        for gene in fields[2:]
        if gene.strip()
    ]

    return (
        name,
        description,
        genes,
    )


# ============================================================
# 5. Load and validate one signature
# ============================================================

def load_signature(
    definition,
):
    """
    Load an Entrez-based UP/DOWN GMT pair.
    """

    up_path = (
        CMAP_DIR
        / definition["up_file"]
    )

    down_path = (
        CMAP_DIR
        / definition["down_file"]
    )

    if not up_path.exists():
        raise FileNotFoundError(
            f"Missing UP file:\n{up_path}"
        )

    if not down_path.exists():
        raise FileNotFoundError(
            f"Missing DOWN file:\n{down_path}"
        )

    up_name, _, up_genes = read_gmt(
        up_path
    )

    down_name, _, down_genes = read_gmt(
        down_path
    )

    expected = definition["name"]

    if up_name != (
        f"{expected}_UP"
    ):
        raise ValueError(
            f"Unexpected UP GMT name:\n"
            f"Expected: {expected}_UP\n"
            f"Found:    {up_name}"
        )

    if down_name != (
        f"{expected}_DN"
    ):
        raise ValueError(
            f"Unexpected DOWN GMT name:\n"
            f"Expected: {expected}_DN\n"
            f"Found:    {down_name}"
        )

    # --------------------------------------------------------
    # Verify they actually look like Entrez IDs.
    # --------------------------------------------------------

    all_genes = (
        up_genes
        + down_genes
    )

    non_numeric = [
        gene
        for gene in all_genes
        if not gene.isdigit()
    ]

    if non_numeric:

        raise ValueError(
            f"{expected} contains non-Entrez "
            f"entries. Examples:\n"
            f"{non_numeric[:10]}"
        )

    return {
        "name": expected,
        "role": definition["role"],
        "up_genes": up_genes,
        "down_genes": down_genes,
    }


# ============================================================
# 6. Serialize GMT
# ============================================================

def serialize_gmt(
    name,
    genes,
):
    """
    Convert genes into one-line GMT text.
    """

    genes = list(
        dict.fromkeys(
            str(gene).strip()
            for gene in genes
            if str(gene).strip()
        )
    )

    return "\t".join(
        [
            name,
            "ATLAS",
            *genes,
        ]
    )


# ============================================================
# 7. Submit one query
# ============================================================

def submit_query(
    signature,
):
    """
    Submit one L1000/Touchstone query.
    """

    up_gmt = serialize_gmt(
        f"{signature['name']}_UP",
        signature["up_genes"],
    )

    down_gmt = serialize_gmt(
        f"{signature['name']}_DN",
        signature["down_genes"],
    )

    payload = {
        "tool_id": "sig_gutc_tool",
        "name": (
            "ATLAS_BT474_TrastuzumabResistance_"
            + signature["name"]
        ),
        "data_type": "L1000",
        "dataset": "Touchstone",
        "ignoreWarnings": True,
        "uptag-cmapfile": up_gmt,
        "dntag-cmapfile": down_gmt,
    }

    response = requests.post(
        "https://api.clue.io/api/jobs",
        headers={
            "user_key": API_KEY,
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        json=payload,
        timeout=120,
    )

    return response


# ============================================================
# 8. Extract useful response information
# ============================================================

def extract_submission_info(
    response,
):
    """
    Extract status/job ID/errors without preserving the
    full CLUE response, which may contain sensitive account
    information.
    """

    try:

        data = response.json()

    except ValueError:

        return {
            "job_id": None,
            "status": None,
            "errors": [
                response.text[:2000]
            ],
        }

    result = data.get(
        "result",
        {}
    )

    if not isinstance(
        result,
        dict,
    ):
        result = {}

    job_id = (
        result.get("job_id")
        or data.get("job_id")
    )

    status = data.get(
        "status"
    )

    errors = []

    # --------------------------------------------------------
    # Top-level system errors
    # --------------------------------------------------------

    system = data.get(
        "system",
        {}
    )

    if isinstance(
        system,
        dict,
    ):

        for error in system.get(
            "errors",
            [],
        ):

            if isinstance(
                error,
                dict,
            ):

                errors.append(
                    error.get(
                        "text",
                        str(error),
                    )
                )

            else:

                errors.append(
                    str(error)
                )

    # --------------------------------------------------------
    # UP errors
    # --------------------------------------------------------

    up = data.get(
        "up",
        {}
    )

    if isinstance(
        up,
        dict,
    ):

        for error in up.get(
            "errors",
            [],
        ):

            if isinstance(
                error,
                dict,
            ):

                errors.append(
                    "UP: "
                    + error.get(
                        "text",
                        str(error),
                    )
                )

            else:

                errors.append(
                    "UP: "
                    + str(error)
                )

    # --------------------------------------------------------
    # DOWN errors
    # --------------------------------------------------------

    down = data.get(
        "down",
        {}
    )

    if isinstance(
        down,
        dict,
    ):

        for error in down.get(
            "errors",
            [],
        ):

            if isinstance(
                error,
                dict,
            ):

                errors.append(
                    "DOWN: "
                    + error.get(
                        "text",
                        str(error),
                    )
                )

            else:

                errors.append(
                    "DOWN: "
                    + str(error)
                )

    return {
        "job_id": job_id,
        "status": status,
        "errors": errors,
    }


# ============================================================
# 9. Main
# ============================================================

print("=" * 60)
print("ATLAS — CMap Valid Multi-Signature Submission")
print("=" * 60)


submission_records = []


for definition in SIGNATURES:

    signature = load_signature(
        definition
    )

    print(
        f"\n{signature['name']}"
    )

    print(
        f"  Role: {signature['role']}"
    )

    print(
        f"  UP Entrez genes:   "
        f"{len(signature['up_genes'])}"
    )

    print(
        f"  DOWN Entrez genes: "
        f"{len(signature['down_genes'])}"
    )

    print(
        "\n  Submitting..."
    )

    response = submit_query(
        signature
    )

    print(
        "  HTTP status:",
        response.status_code,
    )

    info = extract_submission_info(
        response
    )

    job_id = info["job_id"]
    status = info["status"]
    errors = info["errors"]

    # --------------------------------------------------------
    # Success is determined by job ID.
    # --------------------------------------------------------

    submitted = (
        job_id is not None
    )

    if submitted:

        print(
            f"  Job ID: {job_id}"
        )

        print(
            f"  Status: {status}"
        )

    else:

        print(
            "  No job ID returned."
        )

        if errors:

            print(
                "  CLUE errors:"
            )

            for error in errors:

                print(
                    f"    - {error}"
                )

        else:

            print(
                "  No detailed error returned."
            )

    submission_records.append(
        {
            "signature": signature["name"],
            "role": signature["role"],
            "up_genes": len(
                signature["up_genes"]
            ),
            "down_genes": len(
                signature["down_genes"]
            ),
            "http_status": response.status_code,
            "status": status,
            "job_id": job_id,
            "submitted": submitted,
            "errors": errors,
        }
    )


# ============================================================
# 10. Save safe manifest
# ============================================================

manifest = {
    "existing_completed_queries": [
        {
            "signature": "ATLAS_SIG_B_TOP100",
            "job_id": "6a92e7670262700013c4e5e9",
            "status": "completed",
        }
    ],
    "new_submissions": submission_records,
}


manifest_file = (
    JOB_DIR
    / "cmap_job_manifest.json"
)

with open(
    manifest_file,
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        manifest,
        handle,
        indent=2,
    )


# ============================================================
# 11. Summary
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "CMap MULTI-SIGNATURE SUBMISSION SUMMARY"
)

print(
    "=" * 60
)

for record in submission_records:

    print(
        f"{record['signature']}: "
        f"submitted={record['submitted']} | "
        f"HTTP={record['http_status']} | "
        f"status={record['status']} | "
        f"job_id={record['job_id']}"
    )

print(
    f"\nSaved job manifest:\n"
    f"{manifest_file}"
)

print(
    "\nStage 4G submission complete."
)