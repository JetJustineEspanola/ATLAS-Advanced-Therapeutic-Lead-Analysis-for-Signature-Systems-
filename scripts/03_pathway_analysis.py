from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp


# ============================================================
# ATLAS — Stage 3: Pathway and Gene-Set Analysis
# ============================================================

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RESULTS_DIR = (
    PROJECT_ROOT
    / "results"
    / "differential_expression"
)

PATHWAY_DIR = (
    PROJECT_ROOT
    / "results"
    / "pathway_analysis"
)

PATHWAY_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------

DEG_FILE = (
    RESULTS_DIR
    / "DEGs_resistant_vs_sensitive_annotated.csv"
)

SIGNIFICANT_FILE = (
    RESULTS_DIR
    / "significant_DEGs_annotated.csv"
)


# ------------------------------------------------------------
# 3. Load DEG results
# ------------------------------------------------------------

print("=" * 60)
print("ATLAS — Stage 3 Pathway Analysis")
print("=" * 60)

print(f"\nLoading DEG results:\n{DEG_FILE}")

if not DEG_FILE.exists():
    raise FileNotFoundError(
        f"Could not find DEG file:\n{DEG_FILE}"
    )

deg = pd.read_csv(
    DEG_FILE,
    index_col=0,
)

print(f"DEG rows loaded: {len(deg):,}")

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
        "Missing required columns: "
        + str(missing_columns)
    )


# ------------------------------------------------------------
# 4. Clean gene symbols
# ------------------------------------------------------------

deg["Gene name"] = (
    deg["Gene name"]
    .astype("string")
    .str.strip()
    .str.upper()
)

# Remove missing gene symbols
deg_symbols = deg.dropna(
    subset=["Gene name"]
).copy()

# Remove blank gene symbols
deg_symbols = deg_symbols[
    deg_symbols["Gene name"] != ""
].copy()

# Remove duplicate gene symbols
deg_symbols = (
    deg_symbols
    .sort_values("padj", na_position="last")
    .drop_duplicates(
        subset=["Gene name"],
        keep="first",
    )
)

print(
    "\nGenes with usable symbols:",
    f"{len(deg_symbols):,}"
)


# ------------------------------------------------------------
# 5. Create a ranked list for GSEA
# ------------------------------------------------------------

# Use the Wald statistic as the ranking metric.
#
# Positive statistic:
#   higher expression in Resistant
#
# Negative statistic:
#   higher expression in Sensitive

ranked_genes = (
    deg_symbols[
        ["Gene name", "stat"]
    ]
    .dropna()
    .rename(
        columns={
            "Gene name": "gene",
            "stat": "score",
        }
    )
    .sort_values(
        "score",
        ascending=False,
    )
)

print(
    "\nRanked genes available for GSEA:",
    f"{len(ranked_genes):,}"
)


# ------------------------------------------------------------
# 6. Save ranked gene list
# ------------------------------------------------------------

ranked_file = PATHWAY_DIR / "ranked_genes_for_gsea.rnk"

ranked_genes.to_csv(
    ranked_file,
    sep="\t",
    index=False,
    header=False,
)

print(f"Saved ranked gene list: {ranked_file}")


# ------------------------------------------------------------
# 7. Download / retrieve relevant gene-set libraries
# ------------------------------------------------------------

print("\nChecking available Enrichr libraries...")

libraries = gp.get_library_name(
    organism="Human"
)

print(
    f"Available Human libraries: {len(libraries):,}"
)

# We will use libraries that are appropriate for
# broad pathway interpretation.
preferred_libraries = [
    "MSigDB_Hallmark_2020",
    "Reactome_2022",
]

available_libraries = [
    library
    for library in preferred_libraries
    if library in libraries
]

print("\nSelected libraries:")

for library in available_libraries:
    print(f"  - {library}")

if not available_libraries:
    raise RuntimeError(
        "None of the expected pathway libraries "
        "are available from Enrichr."
    )


# ------------------------------------------------------------
# 8. Run pre-ranked GSEA
# ------------------------------------------------------------

print("\nRunning pre-ranked GSEA...")

gsea_results = gp.prerank(
    rnk=ranked_genes,
    gene_sets=available_libraries,
    min_size=10,
    max_size=500,
    permutation_num=1000,
    seed=42,
    verbose=True,
)


# ------------------------------------------------------------
# 9. Save GSEA results
# ------------------------------------------------------------

gsea_table = gsea_results.res2d.copy()

gsea_file = (
    PATHWAY_DIR
    / "gsea_results.csv"
)

gsea_table.to_csv(
    gsea_file,
    index=False,
)

print(
    f"\nSaved GSEA results: {gsea_file}"
)


# ------------------------------------------------------------
# 10. Load significant DEGs
# ------------------------------------------------------------

if SIGNIFICANT_FILE.exists():

    significant = pd.read_csv(
        SIGNIFICANT_FILE,
        index_col=0,
    )

    significant["Gene name"] = (
        significant["Gene name"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    significant = significant.dropna(
        subset=["Gene name"]
    )

else:

    print(
        "\nWARNING: significant DEG file not found."
    )

    significant = deg_symbols[
        deg_symbols["padj"] < 0.05
    ].copy()


# ------------------------------------------------------------
# 11. Split Resistant-up and Resistant-down genes
# ------------------------------------------------------------

up_genes = (
    significant[
        (significant["padj"] < 0.05)
        & (significant["log2FoldChange"] >= 1)
    ]["Gene name"]
    .dropna()
    .drop_duplicates()
    .tolist()
)

down_genes = (
    significant[
        (significant["padj"] < 0.05)
        & (significant["log2FoldChange"] <= -1)
    ]["Gene name"]
    .dropna()
    .drop_duplicates()
    .tolist()
)

print(
    "\nSignature-sized DEG groups:"
)

print(
    "Resistant-upregulated:",
    f"{len(up_genes):,}"
)

print(
    "Resistant-downregulated:",
    f"{len(down_genes):,}"
)


# ------------------------------------------------------------
# 12. Over-representation analysis
# ------------------------------------------------------------

def run_enrichment(
    gene_list,
    label,
):
    """
    Run over-representation analysis using Enrichr.
    """

    if len(gene_list) == 0:
        print(
            f"\nSkipping {label}: empty gene list."
        )
        return None

    print(
        f"\nRunning enrichment for {label}..."
    )

    enrichment = gp.enrichr(
        gene_list=gene_list,
        gene_sets=available_libraries,
        organism="human",
        outdir=None,
    )

    table = enrichment.results.copy()

    output_file = (
        PATHWAY_DIR
        / f"enriched_{label}.csv"
    )

    table.to_csv(
        output_file,
        index=False,
    )

    print(
        f"Saved: {output_file}"
    )

    return table


up_results = run_enrichment(
    up_genes,
    "upregulated",
)

down_results = run_enrichment(
    down_genes,
    "downregulated",
)

# ------------------------------------------------------------
# 13. Targeted pathway/component analysis
# ------------------------------------------------------------

# ============================================================
# TGF-beta targeted analysis
# ============================================================

tgfb_genes = [
    "TGFB1",
    "TGFBR1",
    "TGFBR2",
    "SMAD2",
    "SMAD3",
    "SMAD4",
    "SMAD7",
    "SMURF1",
    "SMURF2",
    "ACVR1",
    "ACVR2A",
    "ACVR2B",
    "BMP2",
    "BMPR1A",
    "BMPR2",
]

# ============================================================
# PD-1 / PD-L1 targeted analysis
# ============================================================

# Core genes explicitly associated with the PD-1/PD-L1 axis.
# Reactome identifies PDCD1 as PD-1, CD274 as PD-L1,
# and PDCD1LG2 as PD-L2.
pd1_pdl1_genes = [
    "PDCD1",
    "CD274",
    "PDCD1LG2",
    "PTPN6",
    "PTPN11",
]


def summarize_target_genes(
    deg_table,
    target_genes,
    pathway_name,
):
    """
    Extract target genes from the DEG table and report
    their differential-expression statistics.
    """

    subset = deg_table[
        deg_table["Gene name"].isin(target_genes)
    ].copy()

    subset = subset[
        [
            "Gene name",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
        ]
    ].sort_values(
        "padj",
        na_position="last",
    )

    subset.insert(
        0,
        "pathway",
        pathway_name,
    )

    return subset


tgfb_gene_results = summarize_target_genes(
    deg_symbols,
    tgfb_genes,
    "TGF-beta signaling",
)

pd1_pdl1_gene_results = summarize_target_genes(
    deg_symbols,
    pd1_pdl1_genes,
    "PD-1/PD-L1 axis",
)


# ------------------------------------------------------------
# 14. Save targeted gene results
# ------------------------------------------------------------

targeted_gene_results = pd.concat(
    [
        tgfb_gene_results,
        pd1_pdl1_gene_results,
    ],
    ignore_index=True,
)

targeted_gene_file = (
    PATHWAY_DIR
    / "targeted_gene_results.csv"
)

targeted_gene_results.to_csv(
    targeted_gene_file,
    index=False,
)

print(
    f"\nSaved targeted gene results: "
    f"{targeted_gene_file}"
)


# ------------------------------------------------------------
# 15. Print targeted gene summary
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("TARGETED GENE SUMMARY")
print("=" * 60)

print("\nTGF-beta-related genes:")
print(
    tgfb_gene_results.to_string(index=False)
)

print("\nPD-1/PD-L1-related genes:")
print(
    pd1_pdl1_gene_results.to_string(index=False)
)


# ------------------------------------------------------------
# 16. Targeted pathway interpretation
# ------------------------------------------------------------

def pathway_status(
    gene_results,
    pathway_name,
    padj_threshold=0.05,
):
    """
    Summarize whether targeted genes show significant
    differential expression.
    """

    if gene_results.empty:
        return {
            "pathway": pathway_name,
            "genes_detected": 0,
            "significant_genes": 0,
            "status": "No target genes detected",
        }

    significant_genes = gene_results[
        gene_results["padj"] < padj_threshold
    ]

    return {
        "pathway": pathway_name,
        "genes_detected": len(gene_results),
        "significant_genes": len(significant_genes),
        "status": (
            "Target genes show significant differential expression"
            if len(significant_genes) > 0
            else "No significant target genes"
        ),
    }


targeted_summary = pd.DataFrame(
    [
        pathway_status(
            tgfb_gene_results,
            "TGF-beta signaling",
        ),
        pathway_status(
            pd1_pdl1_gene_results,
            "PD-1/PD-L1 axis",
        ),
    ]
)

targeted_summary_file = (
    PATHWAY_DIR
    / "targeted_pathway_summary.csv"
)

targeted_summary.to_csv(
    targeted_summary_file,
    index=False,
)

print(
    f"\nSaved targeted pathway summary: "
    f"{targeted_summary_file}"
)


# ------------------------------------------------------------
# 17. Final Stage 3 summary
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("STAGE 3 PATHWAY ANALYSIS SUMMARY")
print("=" * 60)

print("\nTGF-beta:")
print(
    f"  Target genes detected: "
    f"{len(tgfb_gene_results)}"
)

print(
    f"  Significant target genes: "
    f"{(tgfb_gene_results['padj'] < 0.05).sum()}"
)

print("\nPD-1/PD-L1:")
print(
    f"  Target genes detected: "
    f"{len(pd1_pdl1_gene_results)}"
)

print(
    f"  Significant target genes: "
    f"{(pd1_pdl1_gene_results['padj'] < 0.05).sum()}"
)

print("\nStage 3 targeted analysis complete.")