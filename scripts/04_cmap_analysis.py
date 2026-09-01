from pathlib import Path

from dotenv import load_dotenv
import pandas as pd
import mygene


# ============================================================
# ATLAS — Stage 4: CMap Signature Construction
# ============================================================
#
# Purpose:
#   1. Build multiple resistance signatures from Stage 2 DEGs.
#   2. Save HUGO-symbol UP/DOWN GMT files.
#   3. Convert CMap-query signatures to Entrez Gene IDs.
#   4. Save mapping/audit reports.
#
# CMap submission is handled separately by:
#
#   04g_cmap_submit_all.py
#
# IMPORTANT:
#   CMap API submissions use Entrez IDs.
#
# File naming:
#
#   <signature>_up.gmt
#   <signature>_dn.gmt
#
# Entrez versions:
#
#   <signature>_up_entrez.gmt
#   <signature>_dn_entrez.gmt
# ============================================================


# ============================================================
# 1. Project paths
# ============================================================

PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)

load_dotenv(
    dotenv_path=PROJECT_ROOT / ".env"
)

DEG_DIR = (
    PROJECT_ROOT
    / "results"
    / "differential_expression"
)

PATHWAY_DIR = (
    PROJECT_ROOT
    / "results"
    / "pathway_analysis"
)

CMAP_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
)

CMAP_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 2. Input files
# ============================================================

DEG_FILE = (
    DEG_DIR
    / "DEGs_resistant_vs_sensitive_annotated.csv"
)

TARGETED_FILE = (
    PATHWAY_DIR
    / "targeted_gene_results.csv"
)


# ============================================================
# 3. Load DEG results
# ============================================================

print("=" * 60)
print("ATLAS — Stage 4 CMap Signature Construction")
print("=" * 60)

print("\nLoading DEG results:")
print(DEG_FILE)

if not DEG_FILE.exists():
    raise FileNotFoundError(
        f"DEG file not found:\n{DEG_FILE}"
    )

deg = pd.read_csv(
    DEG_FILE,
    index_col=0,
)

print(
    f"Genes loaded: {len(deg):,}"
)


# ============================================================
# 4. Validate DEG columns
# ============================================================

required_columns = [
    "Gene name",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
]

missing_columns = [
    column
    for column in required_columns
    if column not in deg.columns
]

if missing_columns:
    raise ValueError(
        "Missing required DEG columns:\n"
        + str(missing_columns)
    )


# ============================================================
# 5. Clean gene symbols
# ============================================================

deg["Gene name"] = (
    deg["Gene name"]
    .astype("string")
    .str.strip()
    .str.upper()
)

deg = deg.dropna(
    subset=["Gene name"]
)

deg = deg[
    deg["Gene name"] != ""
].copy()


# ============================================================
# 6. Keep one row per gene symbol
# ============================================================

deg = (
    deg
    .sort_values(
        "padj",
        na_position="last",
    )
    .drop_duplicates(
        subset=["Gene name"],
        keep="first",
    )
)

print(
    "Genes with unique usable symbols: "
    f"{len(deg):,}"
)


# ============================================================
# 7. Strict resistance signature
# ============================================================

strict_up = deg[
    (deg["padj"] < 0.05)
    & (deg["log2FoldChange"] >= 1)
].copy()

strict_down = deg[
    (deg["padj"] < 0.05)
    & (deg["log2FoldChange"] <= -1)
].copy()

print("\nStrict resistance signature:")
print(
    f"  UP:   {len(strict_up):,}"
)

print(
    f"  DOWN: {len(strict_down):,}"
)


# ============================================================
# 8. Wald-statistic signatures
# ============================================================

def make_top_signature(
    dataframe: pd.DataFrame,
    n: int,
):
    """
    Select top Resistant-up and Resistant-down genes
    using the Wald statistic.
    """

    up = (
        dataframe[
            dataframe["log2FoldChange"] > 0
        ]
        .sort_values(
            "stat",
            ascending=False,
        )
        .head(n)
        .copy()
    )

    down = (
        dataframe[
            dataframe["log2FoldChange"] < 0
        ]
        .sort_values(
            "stat",
            ascending=True,
        )
        .head(n)
        .copy()
    )

    return up, down


top100_up, top100_down = (
    make_top_signature(
        deg,
        100,
    )
)

top150_wald_up, top150_wald_down = (
    make_top_signature(
        deg,
        150,
    )
)

print("\nTop-100 Wald signature:")
print(
    f"  UP:   {len(top100_up)}"
)

print(
    f"  DOWN: {len(top100_down)}"
)

print("\nTop-150 Wald signature:")
print(
    f"  UP:   {len(top150_wald_up)}"
)

print(
    f"  DOWN: {len(top150_wald_down)}"
)


# ============================================================
# 9. Alternative TOP150 signature
# ============================================================
#
# Statistically significant genes are ranked by absolute
# effect size (log2 fold change).
# ============================================================

significant_deg = deg[
    deg["padj"] < 0.05
].copy()

alternative_up = (
    significant_deg[
        significant_deg["log2FoldChange"] > 0
    ]
    .sort_values(
        "log2FoldChange",
        ascending=False,
    )
    .head(150)
    .copy()
)

alternative_down = (
    significant_deg[
        significant_deg["log2FoldChange"] < 0
    ]
    .sort_values(
        "log2FoldChange",
        ascending=True,
    )
    .head(150)
    .copy()
)

print("\nTop-150 log2FC signature:")
print(
    f"  UP:   {len(alternative_up)}"
)

print(
    f"  DOWN: {len(alternative_down)}"
)


# ============================================================
# 10. TGF-beta-associated signature
# ============================================================

tgfb_up = pd.DataFrame(
    columns=deg.columns
)

tgfb_down = pd.DataFrame(
    columns=deg.columns
)

if TARGETED_FILE.exists():

    targeted = pd.read_csv(
        TARGETED_FILE
    )

    required_targeted_columns = [
        "pathway",
        "Gene name",
    ]

    missing_targeted = [
        column
        for column in required_targeted_columns
        if column not in targeted.columns
    ]

    if missing_targeted:
        raise ValueError(
            "targeted_gene_results.csv is missing columns:\n"
            + str(missing_targeted)
        )

    targeted["Gene name"] = (
        targeted["Gene name"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    tgfb_targets = (
        targeted[
            targeted["pathway"]
            == "TGF-beta signaling"
        ]["Gene name"]
        .dropna()
        .drop_duplicates()
        .tolist()
    )

    tgfb_signature = deg[
        deg["Gene name"].isin(
            tgfb_targets
        )
    ].copy()

    tgfb_up = tgfb_signature[
        (tgfb_signature["padj"] < 0.05)
        & (tgfb_signature["log2FoldChange"] > 0)
    ].copy()

    tgfb_down = tgfb_signature[
        (tgfb_signature["padj"] < 0.05)
        & (tgfb_signature["log2FoldChange"] < 0)
    ].copy()

else:

    print(
        "\nWARNING: targeted_gene_results.csv "
        "was not found."
    )

print("\nTGF-beta-associated signature:")
print(
    f"  UP:   {len(tgfb_up)}"
)

print(
    f"  DOWN: {len(tgfb_down)}"
)


# ============================================================
# 11. GMT writer
# ============================================================

def write_gmt(
    filepath: Path,
    name: str,
    genes,
):
    """
    Write a one-line GMT file.

    Format:
        gene_set_name <TAB> description <TAB> gene1 ...
    """

    clean_genes = []

    for gene in genes:

        if pd.isna(gene):
            continue

        gene = str(gene).strip()

        if not gene:
            continue

        clean_genes.append(
            gene
        )

    clean_genes = list(
        dict.fromkeys(
            clean_genes
        )
    )

    with open(
        filepath,
        "w",
        encoding="utf-8",
    ) as handle:

        handle.write(
            "\t".join(
                [
                    name,
                    "ATLAS",
                    *clean_genes,
                ]
            )
            + "\n"
        )

    return len(clean_genes)


# ============================================================
# 12. Write signature pair
# ============================================================

def write_signature(
    prefix: str,
    up_genes,
    down_genes,
):
    """
    Write HUGO-symbol UP/DOWN files.

    Files:
        <prefix>_up.gmt
        <prefix>_dn.gmt

    Internal GMT names:
        <prefix>_UP
        <prefix>_DN
    """

    up_file = (
        CMAP_DIR
        / f"{prefix}_up.gmt"
    )

    down_file = (
        CMAP_DIR
        / f"{prefix}_dn.gmt"
    )

    up_count = write_gmt(
        up_file,
        f"{prefix}_UP",
        up_genes,
    )

    down_count = write_gmt(
        down_file,
        f"{prefix}_DN",
        down_genes,
    )

    print(
        f"\n{prefix}"
    )

    print(
        f"  UP genes:   {up_count}"
    )

    print(
        f"  DOWN genes: {down_count}"
    )

    print(
        f"  UP file:    {up_file.name}"
    )

    print(
        f"  DOWN file:  {down_file.name}"
    )

    return {
        "signature": prefix,
        "up_count": up_count,
        "down_count": down_count,
        "up_file": str(up_file),
        "down_file": str(down_file),
    }


# ============================================================
# 13. Generate all signatures
# ============================================================

signature_records = []


signature_records.append(
    write_signature(
        "ATLAS_SIG_A_STRICT",
        strict_up["Gene name"],
        strict_down["Gene name"],
    )
)


signature_records.append(
    write_signature(
        "ATLAS_SIG_B_TOP100",
        top100_up["Gene name"],
        top100_down["Gene name"],
    )
)


signature_records.append(
    write_signature(
        "ATLAS_SIG_B_TOP150",
        top150_wald_up["Gene name"],
        top150_wald_down["Gene name"],
    )
)


signature_records.append(
    write_signature(
        "ATLAS_SIG_A_TOP150",
        alternative_up["Gene name"],
        alternative_down["Gene name"],
    )
)


signature_records.append(
    write_signature(
        "ATLAS_SIG_C_TGFB",
        tgfb_up["Gene name"],
        tgfb_down["Gene name"],
    )
)


# ============================================================
# 14. Save signature manifest
# ============================================================

manifest = pd.DataFrame(
    signature_records
)

manifest_file = (
    CMAP_DIR
    / "signature_manifest.csv"
)

manifest.to_csv(
    manifest_file,
    index=False,
)

print(
    f"\nSaved signature manifest:\n"
    f"{manifest_file}"
)


# ============================================================
# 15. Save signature definitions
# ============================================================

definitions = pd.DataFrame(
    {
        "signature": [
            "ATLAS_SIG_A_STRICT",
            "ATLAS_SIG_B_TOP100",
            "ATLAS_SIG_B_TOP150",
            "ATLAS_SIG_A_TOP150",
            "ATLAS_SIG_C_TGFB",
        ],
        "description": [
            (
                "Significant resistance DEGs "
                "with |log2FC| >= 1"
            ),
            (
                "Top 100 Resistant-up and "
                "Resistant-down genes by Wald statistic"
            ),
            (
                "Top 150 Resistant-up and "
                "Resistant-down genes by Wald statistic"
            ),
            (
                "Significant genes ranked by "
                "absolute log2 fold change"
            ),
            (
                "Significant TGF-beta-associated "
                "resistance genes"
            ),
        ],
        "cmap_role": [
            "Reference; too large for direct CMap query",
            "Primary; already completed",
            "Secondary robustness query",
            "Alternative robustness query",
            "Exploratory; currently too small",
        ],
    }
)

definitions_file = (
    CMAP_DIR
    / "signature_definitions.csv"
)

definitions.to_csv(
    definitions_file,
    index=False,
)

print(
    f"Saved signature definitions:\n"
    f"{definitions_file}"
)


# ============================================================
# 16. HUGO → Entrez conversion
# ============================================================

def symbols_to_entrez(
    gene_symbols,
):
    """
    Convert HUGO symbols to Entrez Gene IDs.

    Returns an audit table containing:
        query_symbol
        matched_symbol
        entrez_id
        mapping_status
    """

    symbols = (
        pd.Series(gene_symbols)
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
        .drop_duplicates()
        .tolist()
    )

    if not symbols:

        return pd.DataFrame(
            columns=[
                "query_symbol",
                "matched_symbol",
                "entrez_id",
                "mapping_status",
            ]
        )

    print(
        f"\nMapping {len(symbols):,} gene symbols "
        "to Entrez IDs..."
    )

    mg = mygene.MyGeneInfo()

    query_results = mg.querymany(
        symbols,
        scopes="symbol",
        fields="entrezgene,symbol",
        species="human",
        as_dataframe=False,
        returnall=False,
        verbose=False,
    )

    rows = []

    for result in query_results:

        query_symbol = str(
            result.get(
                "query",
                "",
            )
        ).strip().upper()

        matched_symbol = str(
            result.get(
                "symbol",
                "",
            )
        ).strip().upper()

        if result.get(
            "notfound",
            False,
        ):

            rows.append(
                {
                    "query_symbol": query_symbol,
                    "matched_symbol": matched_symbol,
                    "entrez_id": None,
                    "mapping_status": "not_found",
                }
            )

            continue

        entrez = result.get(
            "entrezgene"
        )

        if isinstance(
            entrez,
            list,
        ):

            entrez = (
                entrez[0]
                if entrez
                else None
            )

        if entrez is None:

            rows.append(
                {
                    "query_symbol": query_symbol,
                    "matched_symbol": matched_symbol,
                    "entrez_id": None,
                    "mapping_status": "no_entrez_id",
                }
            )

        else:

            rows.append(
                {
                    "query_symbol": query_symbol,
                    "matched_symbol": matched_symbol,
                    "entrez_id": str(entrez),
                    "mapping_status": "mapped",
                }
            )

    return pd.DataFrame(
        rows
    )


# ============================================================
# 17. Convert a complete signature to Entrez
# ============================================================

def convert_signature_to_entrez(
    prefix: str,
    up_symbols,
    down_symbols,
):
    """
    Convert one HUGO UP/DOWN signature into Entrez IDs.

    Saves:
        <prefix>_up_entrez.gmt
        <prefix>_dn_entrez.gmt
        <prefix>_entrez_mapping.csv
    """

    up_symbols = (
        pd.Series(up_symbols)
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
        .drop_duplicates()
        .tolist()
    )

    down_symbols = (
        pd.Series(down_symbols)
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
        .drop_duplicates()
        .tolist()
    )

    all_symbols = (
        up_symbols
        + down_symbols
    )

    mapping = symbols_to_entrez(
        all_symbols
    )

    mapping["direction"] = (
        mapping["query_symbol"].apply(
            lambda symbol:
            "UP"
            if symbol in up_symbols
            else "DOWN"
        )
    )

    mapped = mapping[
        mapping["mapping_status"] == "mapped"
    ].copy()

    mapped = (
        mapped
        .drop_duplicates(
            subset=[
                "query_symbol",
                "direction",
            ],
            keep="first",
        )
    )

    mapped_up = (
        mapped[
            mapped["direction"] == "UP"
        ]["entrez_id"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    mapped_down = (
        mapped[
            mapped["direction"] == "DOWN"
        ]["entrez_id"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    # --------------------------------------------------------
    # Save mapping report
    # --------------------------------------------------------

    mapping_file = (
        CMAP_DIR
        / f"{prefix}_entrez_mapping.csv"
    )

    mapping.to_csv(
        mapping_file,
        index=False,
    )

    # --------------------------------------------------------
    # Save Entrez GMT files
    # --------------------------------------------------------

    up_file = (
        CMAP_DIR
        / f"{prefix}_up_entrez.gmt"
    )

    down_file = (
        CMAP_DIR
        / f"{prefix}_dn_entrez.gmt"
    )

    write_gmt(
        up_file,
        f"{prefix}_UP",
        mapped_up,
    )

    write_gmt(
        down_file,
        f"{prefix}_DN",
        mapped_down,
    )

    print(
        f"\nEntrez conversion — {prefix}"
    )

    print(
        f"  Original UP:        {len(up_symbols)}"
    )

    print(
        f"  Original DOWN:      {len(down_symbols)}"
    )

    print(
        f"  Mapped UP Entrez:   {len(mapped_up)}"
    )

    print(
        f"  Mapped DOWN Entrez: {len(mapped_down)}"
    )

    print(
        f"  Mapping report:     {mapping_file.name}"
    )

    print(
        f"  UP file:            {up_file.name}"
    )

    print(
        f"  DOWN file:          {down_file.name}"
    )

    return {
        "signature": prefix,
        "original_up": len(up_symbols),
        "original_down": len(down_symbols),
        "mapped_up": len(mapped_up),
        "mapped_down": len(mapped_down),
        "mapping_file": str(mapping_file),
        "up_file": str(up_file),
        "down_file": str(down_file),
    }


# ============================================================
# 18. Convert CMap-queryable signatures to Entrez
# ============================================================

entrez_records = []


# TOP100
entrez_records.append(
    convert_signature_to_entrez(
        "ATLAS_SIG_B_TOP100",
        top100_up["Gene name"],
        top100_down["Gene name"],
    )
)


# TOP150 Wald
entrez_records.append(
    convert_signature_to_entrez(
        "ATLAS_SIG_B_TOP150",
        top150_wald_up["Gene name"],
        top150_wald_down["Gene name"],
    )
)


# TOP150 log2FC
entrez_records.append(
    convert_signature_to_entrez(
        "ATLAS_SIG_A_TOP150",
        alternative_up["Gene name"],
        alternative_down["Gene name"],
    )
)


# ============================================================
# 19. Save Entrez conversion summary
# ============================================================

entrez_summary = pd.DataFrame(
    entrez_records
)

entrez_summary_file = (
    CMAP_DIR
    / "cmap_entrez_conversion_summary.csv"
)

entrez_summary.to_csv(
    entrez_summary_file,
    index=False,
)

print(
    f"\nSaved Entrez conversion summary:\n"
    f"{entrez_summary_file}"
)


# ============================================================
# 20. Final CMap signature validation
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "FINAL CMap SIGNATURE VALIDATION"
)

print(
    "=" * 60
)

for record in entrez_records:

    print(
        f"\n{record['signature']}"
    )

    print(
        f"  Original UP:        "
        f"{record['original_up']}"
    )

    print(
        f"  Original DOWN:      "
        f"{record['original_down']}"
    )

    print(
        f"  Entrez UP:          "
        f"{record['mapped_up']}"
    )

    print(
        f"  Entrez DOWN:        "
        f"{record['mapped_down']}"
    )


# ============================================================
# 21. Complete
# ============================================================

print(
    "\n" + "=" * 60
)

print(
    "STAGE 4 SIGNATURE CONSTRUCTION COMPLETE"
)

print(
    "=" * 60
)

print(
    "\nReady for CMap API submission:"
)

print(
    "  ATLAS_SIG_B_TOP100"
)

print(
    "  ATLAS_SIG_B_TOP150"
)

print(
    "  ATLAS_SIG_A_TOP150"
)

print(
    "\nNot submitted by this script:"
)

print(
    "  ATLAS_SIG_A_STRICT"
)

print(
    "  ATLAS_SIG_C_TGFB"
)