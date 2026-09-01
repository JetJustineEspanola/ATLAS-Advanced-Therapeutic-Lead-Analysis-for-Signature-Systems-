#!/usr/bin/env python3

"""
ATLAS — Stage 4M: CMap Compound Identity and Drug Filtering

Purpose
-------
Take the prioritized CMap compounds produced by Stage 4L and prepare
them for downstream drug/regulatory verification.

This stage:

1. Loads:
       results/cmap/prioritized/ATLAS_CMap_prioritized_all.csv

2. Preserves all CMap scores and rankings from Stage 4L.

3. Resolves compound identity against PubChem where possible.

4. Retrieves:
       - PubChem CID
       - molecular formula
       - molecular weight
       - canonical SMILES
       - isomeric SMILES
       - InChI
       - InChIKey
       - PubChem title

5. Flags compounds that appear to be:
       - recognizable chemical compounds
       - internal/research compound identifiers
       - unresolved compounds

6. Generates candidate sets for the next ATLAS stage.

IMPORTANT
---------
This script DOES NOT determine FDA approval.

A PubChem match does NOT mean:
    - FDA-approved
    - clinically effective
    - suitable for breast cancer
    - safe
    - validated against trastuzumab resistance

Regulatory/drug status will be independently verified in Stage 4N.
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests


# ============================================================
# 1. Project paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

CMAP_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
)

PRIORITIZED_DIR = (
    CMAP_DIR
    / "prioritized"
)

OUTPUT_DIR = (
    CMAP_DIR
    / "drug_filter"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 2. Input
# ============================================================

INPUT_FILE = (
    PRIORITIZED_DIR
    / "ATLAS_CMap_prioritized_all.csv"
)


# ============================================================
# 3. Outputs
# ============================================================

ANNOTATED_FILE = (
    OUTPUT_DIR
    / "ATLAS_CMap_compound_annotations.csv"
)

DRUG_CANDIDATES_FILE = (
    OUTPUT_DIR
    / "ATLAS_CMap_drug_candidates.csv"
)

RESEARCH_COMPOUNDS_FILE = (
    OUTPUT_DIR
    / "ATLAS_CMap_research_compounds.csv"
)

MANUAL_REVIEW_FILE = (
    OUTPUT_DIR
    / "ATLAS_CMap_manual_review.csv"
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "ATLAS_CMap_filter_summary.csv"
)

METADATA_FILE = (
    OUTPUT_DIR
    / "ATLAS_CMap_filter_metadata.json"
)


# ============================================================
# 4. PubChem configuration
# ============================================================

PUBCHEM_BASE = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
)

REQUEST_TIMEOUT = 20

# Small pause to avoid hammering PubChem.
REQUEST_DELAY = 0.08

MAX_RETRIES = 3


# ============================================================
# 5. Expected Stage 4L columns
# ============================================================

REQUIRED_COLUMNS = [
    "priority_rank",
    "priority_tier_number",
    "priority_tier",
    "pert_id",
    "pert_iname",
    "ATLAS_SIG_A_TOP150",
    "ATLAS_SIG_B_TOP100",
    "ATLAS_SIG_B_TOP150",
    "n_signatures",
    "n_negative",
    "n_strong_negative",
    "mean_tau",
    "median_tau",
    "minimum_tau",
]


# ============================================================
# 6. HTTP session
# ============================================================

SESSION = requests.Session()

SESSION.headers.update(
    {
        "User-Agent": (
            "ATLAS-CMap-Drug-Annotation/1.0 "
            "(bioinformatics research workflow)"
        ),
        "Accept": "application/json",
    }
)


# ============================================================
# 7. Helper functions
# ============================================================

def print_header(text: str) -> None:

    print()
    print("=" * 60)
    print(text)
    print("=" * 60)


def clean_text(value: Any) -> str:

    if pd.isna(value):
        return ""

    return re.sub(
        r"\s+",
        " ",
        str(value).strip(),
    )


def normalize_name(value: Any) -> str:

    text = clean_text(value).lower()

    text = (
        text
        .replace("β", "beta")
        .replace("–", "-")
        .replace("—", "-")
    )

    return text


# ============================================================
# 8. Research/internal compound patterns
# ============================================================

RESEARCH_PATTERNS = [

    # Broad / LINCS compound identifiers
    r"^brd-[a-z0-9]+$",

    # Common medicinal chemistry / tool-compound prefixes
    r"^unc[-_ ]?\d+",
    r"^sb[-_ ]?\d+",
    r"^pd[-_ ]?\d+",
    r"^vu[-_ ]?\d+",
    r"^ql[-_ ]?",
    r"^alw[-_ ]?",
    r"^sj[-_ ]?\d+",
    r"^ro[-_ ]?\d+",
    r"^arc[-_ ]?\d+",
    r"^ci[-_ ]?\d+",
    r"^way[-_ ]?\d+",
    r"^dmp[-_ ]?\d+",
    r"^gs[-_ ]?\d+",

    # Generic anonymous IDs
    r"^compound[-_ ]?\d+",
    r"^cmpd[-_ ]?\d+",
    r"^cpd[-_ ]?\d+",
]


def looks_like_research_identifier(
    name: str,
) -> bool:

    normalized = normalize_name(name)

    if not normalized:
        return False

    for pattern in RESEARCH_PATTERNS:

        if re.search(
            pattern,
            normalized,
            flags=re.IGNORECASE,
        ):

            return True

    return False


# ============================================================
# 9. Safe PubChem request
# ============================================================

def safe_request(
    url: str,
) -> dict[str, Any] | None:

    for attempt in range(
        1,
        MAX_RETRIES + 1,
    ):

        try:

            response = SESSION.get(
                url,
                timeout=REQUEST_TIMEOUT,
            )

        except requests.RequestException:

            if attempt == MAX_RETRIES:
                return None

            time.sleep(
                attempt
            )

            continue

        if response.status_code == 200:

            try:
                return response.json()

            except ValueError:
                return None

        # 404 simply means not found.
        if response.status_code == 404:
            return None

        # Rate limiting / temporary server errors.
        if response.status_code in {
            429,
            500,
            502,
            503,
            504,
        }:

            if attempt == MAX_RETRIES:
                return None

            time.sleep(
                attempt * 2
            )

            continue

        return None

    return None


# ============================================================
# 10. PubChem name lookup
# ============================================================

def get_pubchem_properties(
    compound_name: str,
) -> dict[str, Any]:

    result = {

        "pubchem_cid": None,

        "pubchem_title": None,

        "molecular_formula": None,

        "molecular_weight": None,

        "canonical_smiles": None,

        "isomeric_smiles": None,

        "inchi": None,

        "inchikey": None,

        "pubchem_status": "not_found",
    }

    compound_name = clean_text(
        compound_name
    )

    if not compound_name:
        return result


    encoded_name = requests.utils.quote(
        compound_name,
        safe="",
    )


    properties = ",".join(
        [
            "Title",
            "MolecularFormula",
            "MolecularWeight",
            "CanonicalSMILES",
            "IsomericSMILES",
            "InChI",
            "InChIKey",
        ]
    )


    url = (
        f"{PUBCHEM_BASE}"
        f"/compound/name/"
        f"{encoded_name}"
        f"/property/"
        f"{properties}"
        f"/JSON"
    )


    data = safe_request(
        url
    )


    if data is None:
        return result


    try:

        records = (
            data[
                "PropertyTable"
            ][
                "Properties"
            ]
        )

    except (
        KeyError,
        TypeError,
    ):

        return result


    if not records:
        return result


    record = records[0]


    result[
        "pubchem_cid"
    ] = record.get(
        "CID"
    )


    result[
        "pubchem_title"
    ] = record.get(
        "Title"
    )


    result[
        "molecular_formula"
    ] = record.get(
        "MolecularFormula"
    )


    result[
        "molecular_weight"
    ] = record.get(
        "MolecularWeight"
    )


    # PubChem occasionally changes property-key naming.
    result[
        "canonical_smiles"
    ] = (
        record.get(
            "ConnectivitySMILES"
        )
        or
        record.get(
            "CanonicalSMILES"
        )
    )


    result[
        "isomeric_smiles"
    ] = (
        record.get(
            "SMILES"
        )
        or
        record.get(
            "IsomericSMILES"
        )
    )


    result[
        "inchi"
    ] = record.get(
        "InChI"
    )


    result[
        "inchikey"
    ] = record.get(
        "InChIKey"
    )


    result[
        "pubchem_status"
    ] = "resolved"


    return result


# ============================================================
# 11. Compound classification
# ============================================================

def classify_compound(
    pert_iname: str,
    pubchem_status: str,
) -> tuple[str, str]:

    """
    Conservative Stage 4M classification.

    This is an identity/technical classification,
    NOT a regulatory-status classification.
    """

    if looks_like_research_identifier(
        pert_iname
    ):

        return (

            "research_or_tool_compound",

            (
                "Compound name resembles an internal, "
                "medicinal-chemistry, or research identifier."
            ),
        )


    if pubchem_status == "resolved":

        return (

            "identified_compound",

            (
                "Compound identity and structure were "
                "resolved through PubChem."
            ),
        )


    return (

        "manual_identity_review",

        (
            "No confident PubChem identity was found "
            "from the CMap perturbagen name."
        ),
    )


# ============================================================
# 12. Downstream eligibility
# ============================================================

def assign_downstream_status(
    row: pd.Series,
) -> str:

    """
    Decide whether the compound can move to Stage 4N.

    Stage 4N will perform actual drug/regulatory verification.
    """

    classification = row[
        "compound_classification"
    ]


    if (
        classification
        == "research_or_tool_compound"
    ):

        return (
            "RESEARCH_COMPOUND_REVIEW"
        )


    if (
        classification
        == "manual_identity_review"
    ):

        return (
            "MANUAL_IDENTITY_REVIEW"
        )


    # Identified chemical with multi-signature evidence.
    if (
        row["n_negative"] >= 2
    ):

        return (
            "CANDIDATE_FOR_DRUG_STATUS_REVIEW"
        )


    return (
        "LOWER_PRIORITY_IDENTIFIED_COMPOUND"
    )


# ============================================================
# 13. Main
# ============================================================

def main() -> int:

    print_header(
        "ATLAS — Stage 4M CMap Compound Identity & Drug Filtering"
    )


    # --------------------------------------------------------
    # Input
    # --------------------------------------------------------

    print()
    print("Input:")
    print(INPUT_FILE)


    if not INPUT_FILE.exists():

        print()
        print(
            "ERROR: Stage 4L output was not found."
        )

        print(
            "Run:"
        )

        print(
            "python scripts/04l_cmap_prioritize.py"
        )

        return 1


    # --------------------------------------------------------
    # Load Stage 4L
    # --------------------------------------------------------

    print()
    print(
        "Loading prioritized CMap compounds..."
    )


    df = pd.read_csv(
        INPUT_FILE
    )


    print(
        f"Rows loaded: {len(df):,}"
    )

    print(
        f"Columns loaded: {len(df.columns)}"
    )


    # --------------------------------------------------------
    # Validate schema
    # --------------------------------------------------------

    missing_columns = [

        column

        for column
        in REQUIRED_COLUMNS

        if column
        not in df.columns
    ]


    if missing_columns:

        print()
        print(
            "ERROR: Required Stage 4L columns are missing:"
        )

        for column in missing_columns:

            print(
                f"  - {column}"
            )

        return 1


    # --------------------------------------------------------
    # Numeric conversion
    # --------------------------------------------------------

    numeric_columns = [

        "priority_rank",

        "priority_tier_number",

        "ATLAS_SIG_A_TOP150",

        "ATLAS_SIG_B_TOP100",

        "ATLAS_SIG_B_TOP150",

        "n_signatures",

        "n_negative",

        "n_strong_negative",

        "mean_tau",

        "median_tau",

        "minimum_tau",
    ]


    for column in numeric_columns:

        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )


    # --------------------------------------------------------
    # Compound identity cleanup
    # --------------------------------------------------------

    df["pert_id"] = (
        df["pert_id"]
        .astype("string")
        .str.strip()
    )


    df["pert_iname"] = (
        df["pert_iname"]
        .astype("string")
        .str.strip()
    )


    before = len(
        df
    )


    df = (
        df
        .drop_duplicates(
            subset=[
                "pert_id",
            ],
            keep="first",
        )
        .copy()
    )


    print()
    print(
        f"Unique perturbagens: {len(df):,}"
    )


    if before != len(df):

        print(
            "Duplicate perturbagens removed:",
            f"{before - len(df):,}",
        )


    # --------------------------------------------------------
    # PubChem annotation
    # --------------------------------------------------------

    print()
    print(
        "Resolving compound identities with PubChem..."
    )

    print(
        "This may take several minutes for thousands of compounds."
    )

    print()


    annotations = []


    total = len(
        df
    )


    for counter, (
        index,
        row,
    ) in enumerate(
        df.iterrows(),
        start=1,
    ):


        pert_id = clean_text(
            row["pert_id"]
        )


        pert_iname = clean_text(
            row["pert_iname"]
        )


        properties = (
            get_pubchem_properties(
                pert_iname
            )
        )


        classification, reason = (
            classify_compound(

                pert_iname,

                properties[
                    "pubchem_status"
                ],
            )
        )


        annotation = {

            "pert_id": pert_id,

            "pert_iname": pert_iname,

            **properties,

            "compound_classification":
                classification,

            "classification_reason":
                reason,

            # Stage 4N intentionally owns these.
            "regulatory_status":
                "NOT_VERIFIED",

            "regulatory_source":
                None,

            "regulatory_verification_date":
                None,
        }


        annotations.append(
            annotation
        )


        if (
            counter == 1
            or counter % 100 == 0
            or counter == total
        ):

            print(
                f"  {counter:,}/{total:,} compounds processed"
            )


        time.sleep(
            REQUEST_DELAY
        )


    annotation_df = pd.DataFrame(
        annotations
    )


    # --------------------------------------------------------
    # Merge
    # --------------------------------------------------------

    annotated = df.merge(

        annotation_df,

        on=[
            "pert_id",
            "pert_iname",
        ],

        how="left",

        validate="one_to_one",
    )


    # --------------------------------------------------------
    # Derived structure flags
    # --------------------------------------------------------

    annotated[
        "has_pubchem_identity"
    ] = (

        annotated[
            "pubchem_status"
        ]
        == "resolved"
    )


    annotated[
        "has_structure"
    ] = (

        annotated[
            "inchikey"
        ]
        .fillna("")
        .astype(str)
        .str.strip()
        .ne("")
    )


    annotated[
        "multi_signature_negative"
    ] = (

        annotated[
            "n_negative"
        ]
        >= 2
    )


    annotated[
        "strong_cmap_support"
    ] = (

        annotated[
            "n_strong_negative"
        ]
        >= 1
    )


    # --------------------------------------------------------
    # Stage 4N status
    # --------------------------------------------------------

    annotated[
        "downstream_status"
    ] = annotated.apply(

        assign_downstream_status,

        axis=1,
    )


    # --------------------------------------------------------
    # Preserve 04L ranking
    # --------------------------------------------------------

    annotated = (

        annotated

        .sort_values(
            [
                "priority_rank",
            ],
            ascending=True,
            na_position="last",
        )

        .reset_index(
            drop=True
        )
    )


    # ========================================================
    # 14. Output groups
    # ========================================================

    drug_candidates = annotated[

        annotated[
            "downstream_status"
        ]
        == "CANDIDATE_FOR_DRUG_STATUS_REVIEW"

    ].copy()


    research_compounds = annotated[

        annotated[
            "downstream_status"
        ]
        == "RESEARCH_COMPOUND_REVIEW"

    ].copy()


    manual_review = annotated[

        annotated[
            "downstream_status"
        ]
        == "MANUAL_IDENTITY_REVIEW"

    ].copy()


    # ========================================================
    # 15. Save outputs
    # ========================================================

    print()
    print(
        "Saving Stage 4M outputs..."
    )


    annotated.to_csv(

        ANNOTATED_FILE,

        index=False,

        float_format="%.6f",
    )


    drug_candidates.to_csv(

        DRUG_CANDIDATES_FILE,

        index=False,

        float_format="%.6f",
    )


    research_compounds.to_csv(

        RESEARCH_COMPOUNDS_FILE,

        index=False,

        float_format="%.6f",
    )


    manual_review.to_csv(

        MANUAL_REVIEW_FILE,

        index=False,

        float_format="%.6f",
    )


    # ========================================================
    # 16. Summary
    # ========================================================

    summary_rows = []


    for status, group in (
        annotated
        .groupby(
            "downstream_status",
            dropna=False,
        )
    ):

        summary_rows.append(
            {

                "downstream_status":
                    status,

                "compound_count":
                    len(group),

                "tier1_count":
                    int(
                        (
                            group[
                                "priority_tier_number"
                            ]
                            == 1
                        ).sum()
                    ),

                "mean_tau":
                    group[
                        "mean_tau"
                    ].mean(),

                "median_tau":
                    group[
                        "median_tau"
                    ].median(),

                "minimum_tau":
                    group[
                        "minimum_tau"
                    ].min(),
            }
        )


    summary = pd.DataFrame(
        summary_rows
    )


    if not summary.empty:

        summary = summary.sort_values(

            "compound_count",

            ascending=False,
        )


    summary.to_csv(

        SUMMARY_FILE,

        index=False,

        float_format="%.6f",
    )


    metadata = {

        "stage":
            "04M",

        "input_file":
            str(INPUT_FILE),

        "input_compounds":
            int(len(df)),

        "annotated_compounds":
            int(len(annotated)),

        "pubchem_resolved":
            int(
                annotated[
                    "has_pubchem_identity"
                ].sum()
            ),

        "candidate_for_drug_status_review":
            int(
                len(
                    drug_candidates
                )
            ),

        "research_compounds":
            int(
                len(
                    research_compounds
                )
            ),

        "manual_review":
            int(
                len(
                    manual_review
                )
            ),

        "regulatory_status_verified":
            False,

        "next_stage":
            "04N_regulatory_and_drug_status_verification",
    }


    with open(

        METADATA_FILE,

        "w",

        encoding="utf-8",

    ) as handle:

        json.dump(

            metadata,

            handle,

            indent=2,
        )


    # ========================================================
    # 17. Console summary
    # ========================================================

    print_header(
        "STAGE 4M COMPOUND FILTER SUMMARY"
    )


    print()
    print(
        f"Total prioritized compounds: "
        f"{len(annotated):,}"
    )


    print()
    print(
        "PubChem identity resolution:"
    )


    print(
        "  Resolved:",
        f"{annotated['has_pubchem_identity'].sum():,}",
    )


    print(
        "  Unresolved:",
        f"{(~annotated['has_pubchem_identity']).sum():,}",
    )


    print()
    print(
        "Compound classification:"
    )


    classification_counts = (

        annotated[
            "compound_classification"
        ]
        .value_counts(
            dropna=False
        )
    )


    for classification, count in (
        classification_counts.items()
    ):

        print(
            f"  {classification}: "
            f"{count:,}"
        )


    print()
    print(
        "Downstream groups:"
    )


    downstream_counts = (

        annotated[
            "downstream_status"
        ]
        .value_counts(
            dropna=False
        )
    )


    for status, count in (
        downstream_counts.items()
    ):

        print(
            f"  {status}: "
            f"{count:,}"
        )


    # --------------------------------------------------------
    # Tier 1 overview
    # --------------------------------------------------------

    tier1 = annotated[

        annotated[
            "priority_tier_number"
        ]
        == 1

    ].copy()


    print_header(
        "TIER 1 COMPOUND ANNOTATION"
    )


    display_columns = [

        "priority_rank",

        "pert_id",

        "pert_iname",

        "mean_tau",

        "median_tau",

        "n_negative",

        "n_strong_negative",

        "pubchem_cid",

        "pubchem_title",

        "compound_classification",

        "downstream_status",
    ]


    if tier1.empty:

        print()
        print(
            "No Tier 1 compounds found."
        )

    else:

        print()
        print(

            tier1[
                display_columns
            ]
            .to_string(
                index=False
            )
        )


    # ========================================================
    # 18. Complete
    # ========================================================

    print_header(
        "STAGE 4M COMPLETE"
    )


    print()
    print(
        "Primary annotation output:"
    )

    print(
        ANNOTATED_FILE
    )


    print()
    print(
        "Candidates for Stage 4N:"
    )

    print(
        DRUG_CANDIDATES_FILE
    )


    print()
    print(
        "Research/tool compounds:"
    )

    print(
        RESEARCH_COMPOUNDS_FILE
    )


    print()
    print(
        "Manual identity review:"
    )

    print(
        MANUAL_REVIEW_FILE
    )


    print()
    print(
        "Summary:"
    )

    print(
        SUMMARY_FILE
    )


    print()
    print(
        "IMPORTANT:"
    )

    print(
        "No FDA/regulatory status has been inferred."
    )

    print(
        "Stage 4N will independently verify "
        "drug/regulatory status."
    )


    return 0


# ============================================================
# 19. Entry point
# ============================================================

if __name__ == "__main__":

    sys.exit(
        main()
    )