"""
02_differential_expression.py

1. Load gene_counts.csv
2. Load sample metadata
3. Filter uninformative/very-low-count genes
4. Prepare the count matrix appropriately
5. Run PyDESeq2
6. Compare Resistant vs Sensitive
7. Save DEG results
8. Generate MA plot
9. Generate volcano plot
10. Save an analysis summary
"""

from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import pandas as pd
import numpy as np

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ============================================================
# ATLAS — Stage 2: Differential Expression Analysis
# ============================================================

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR = PROJECT_ROOT / "results" / "differential_expression"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

COUNTS_FILE = PROCESSED_DIR / "gene_counts.csv"


# ------------------------------------------------------------
# 2. Load gene-level counts
# ------------------------------------------------------------

print("=" * 60)
print("ATLAS — Stage 2 Differential Expression")
print("=" * 60)

print(f"\nLoading: {COUNTS_FILE}")

if not COUNTS_FILE.exists():
    raise FileNotFoundError(
        f"Could not find gene count matrix: {COUNTS_FILE}"
    )

gene_counts = pd.read_csv(
    COUNTS_FILE,
    index_col=0,
)

print(f"Count matrix shape: {gene_counts.shape}")


# ------------------------------------------------------------
# 3. Define sample metadata
# ------------------------------------------------------------

metadata = pd.DataFrame(
    {
        "condition": [
            "Sensitive",
            "Sensitive",
            "Sensitive",
            "Resistant",
            "Resistant",
            "Resistant",
        ]
    },
    index=["TS1", "TS2", "TS3", "TR1", "TR2", "TR3"],
)

print("\nSample metadata:")
print(metadata)


# ------------------------------------------------------------
# 4. Validate samples
# ------------------------------------------------------------

if list(gene_counts.columns) != list(metadata.index):
    raise ValueError(
        "Sample columns in gene_counts.csv do not match metadata."
    )

print("\nSample validation PASSED.")


# ------------------------------------------------------------
# 5. Validate values
# ------------------------------------------------------------

if gene_counts.isna().any().any():
    raise ValueError("Count matrix contains missing values.")

if (gene_counts < 0).any().any():
    raise ValueError("Count matrix contains negative values.")

print("Missing-value check PASSED.")
print("Negative-value check PASSED.")

# ------------------------------------------------------------
# 6. Filter genes with very low counts
# ------------------------------------------------------------

min_count = 10
min_samples = 3

keep_genes = (
    (gene_counts >= min_count)
    .sum(axis=1)
    >= min_samples
)

gene_counts_filtered = gene_counts.loc[keep_genes].copy()

print("\nLow-count filtering:")
print(f"Original genes:        {gene_counts.shape[0]:,}")
print(f"Genes retained:        {gene_counts_filtered.shape[0]:,}")
print(f"Genes removed:         {(~keep_genes).sum():,}")
print(
    f"Filter: >= {min_count} counts "
    f"in >= {min_samples} samples"
)

# ------------------------------------------------------------
# 7. Prepare integer counts for PyDESeq2
# ------------------------------------------------------------

rounded_counts = np.rint(
    gene_counts_filtered
).astype(int)

print("\nCount preparation:")
print(
    "Fractional values before rounding:",
    np.sum(
        ~np.isclose(
            gene_counts_filtered.to_numpy(),
            np.round(gene_counts_filtered.to_numpy())
        )
    )
)

print(
    "Fractional values after rounding:",
    np.sum(
        ~np.isclose(
            rounded_counts.to_numpy(),
            np.round(rounded_counts.to_numpy())
        )
    )
)

print(
    "Minimum rounded count:",
    rounded_counts.min().min()
)

print(
    "Maximum rounded count:",
    rounded_counts.max().max()
)

# ------------------------------------------------------------
# 8. Save prepared counts and metadata
# ------------------------------------------------------------

prepared_counts_file = (
    PROCESSED_DIR / "deseq_counts.csv"
)

prepared_metadata_file = (
    PROCESSED_DIR / "sample_metadata.csv"
)

rounded_counts.to_csv(
    prepared_counts_file
)

metadata.to_csv(
    prepared_metadata_file
)

print("\nPrepared DESeq2 inputs saved:")
print(f"  Counts:   {prepared_counts_file}")
print(f"  Metadata: {prepared_metadata_file}")

# Verify saved inputs
check_deseq_counts = pd.read_csv(
    prepared_counts_file,
    index_col=0
)

check_metadata = pd.read_csv(
    prepared_metadata_file,
    index_col=0
)

print("\nPrepared counts shape:", check_deseq_counts.shape)
print("Metadata shape:", check_metadata.shape)
print("Counts contain only integers:",
      np.issubdtype(check_deseq_counts.dtypes.iloc[0], np.integer))

# ------------------------------------------------------------
# 9. Prepare data for PyDESeq2
# ------------------------------------------------------------

# PyDESeq2 expects:
#   rows    = samples
#   columns = genes
#   values  = integer counts

counts_for_deseq = check_deseq_counts.T.copy()

# Make sure metadata follows exactly the same sample order
metadata_for_deseq = check_metadata.loc[
    counts_for_deseq.index
].copy()

print("\nPyDESeq2 input shape:")
print("  Counts:", counts_for_deseq.shape)
print("  Metadata:", metadata_for_deseq.shape)

print("\nPyDESeq2 sample order:")
print(counts_for_deseq.index.tolist())

print("\nMetadata:")
print(metadata_for_deseq)

# ------------------------------------------------------------
# 10. Create PyDESeq2 dataset
# ------------------------------------------------------------

dds = DeseqDataSet(
    counts=counts_for_deseq,
    metadata=metadata_for_deseq,
    design="~condition",
    refit_cooks=True,
    n_cpus=1,
)

print("\nPyDESeq2 dataset created successfully.")

# ------------------------------------------------------------
# 11. Run differential-expression analysis
# ------------------------------------------------------------

print("\nRunning PyDESeq2...")

dds.deseq2()

print("PyDESeq2 analysis completed.")

# ------------------------------------------------------------
# 12. Resistant vs Sensitive contrast
# ------------------------------------------------------------

stat_res = DeseqStats(
    dds,
    contrast=("condition", "Resistant", "Sensitive"),
    n_cpus=1,
)

stat_res.summary()

results = stat_res.results_df.copy()

print("\nDEG results shape:", results.shape)

print("\nTop 10 results:")
print(results.head(10))

# ------------------------------------------------------------
# 13. Sort and inspect DEG results
# ------------------------------------------------------------

results_sorted = results.sort_values(
    by="padj",
    na_position="last"
)

print("\nTop 20 genes by adjusted p-value:")
print(results_sorted.head(20))

# ------------------------------------------------------------
# 14. DEG summary
# ------------------------------------------------------------

# Significant based on adjusted p-value
significant = results[
    results["padj"] < 0.05
].copy()

# Significant AND biologically substantial effect
significant_fc = significant[
    significant["log2FoldChange"].abs() >= 1
].copy()

# Resistant-upregulated genes
upregulated = significant[
    significant["log2FoldChange"] >= 1
].copy()

# Resistant-downregulated genes
downregulated = significant[
    significant["log2FoldChange"] <= -1
].copy()

print("\n" + "=" * 60)
print("DEG SUMMARY")
print("=" * 60)

print(f"Genes tested:                {len(results):,}")
print(f"Significant (padj < 0.05):   {len(significant):,}")
print(
    f"Significant + |log2FC| >= 1: "
    f"{len(significant_fc):,}"
)
print(f"Upregulated in Resistant:    {len(upregulated):,}")
print(f"Downregulated in Resistant:  {len(downregulated):,}")

print("\nTop 20 genes by absolute log2 fold change:")

top_effect_genes = results.loc[
    results["log2FoldChange"].abs().sort_values(
        ascending=False
    ).index
].head(20)

print(top_effect_genes)

# ------------------------------------------------------------
# 15. Add gene names
# ------------------------------------------------------------

tx2gene_full = pd.read_csv(
    PROCESSED_DIR.parent / "raw" / "tx2gene.tsv",
    sep="\t"
)

gene_annotation = (
    tx2gene_full[
        ["Gene stable ID", "Gene name"]
    ]
    .drop_duplicates()
    .dropna(subset=["Gene stable ID"])
)

gene_annotation = (
    gene_annotation
    .drop_duplicates(subset=["Gene stable ID"])
    .set_index("Gene stable ID")
)

results_annotated = results.join(
    gene_annotation,
    how="left"
)

results_annotated = results_annotated[
    [
        "Gene name",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
    ]
]

print("\nTop 20 significant genes with names:")

print(
    results_annotated
    .sort_values("padj", na_position="last")
    .head(20)
)

# ------------------------------------------------------------
# 16. Volcano plot
# ------------------------------------------------------------

volcano = results.copy()

# Avoid log10(0)
volcano["padj_plot"] = volcano["padj"].clip(lower=1e-300)

volcano["neg_log10_padj"] = -np.log10(
    volcano["padj_plot"]
)

fig, ax = plt.subplots(figsize=(9, 7))

ax.scatter(
    volcano["log2FoldChange"],
    volcano["neg_log10_padj"],
    s=8,
    alpha=0.5,
)

ax.axvline(
    1,
    linestyle="--",
    linewidth=1,
)

ax.axvline(
    -1,
    linestyle="--",
    linewidth=1,
)

ax.axhline(
    -np.log10(0.05),
    linestyle="--",
    linewidth=1,
)

ax.set_xlabel("log2 Fold Change")
ax.set_ylabel("-log10 Adjusted p-value")
ax.set_title("Resistant vs Sensitive — Differential Expression")

plt.tight_layout()

volcano_file = RESULTS_DIR / "volcano_plot.png"
plt.savefig(volcano_file, dpi=300)
plt.close()

print(f"Saved: {volcano_file}")

# ------------------------------------------------------------
# 17. MA plot
# ------------------------------------------------------------

ma = results.copy()

fig, ax = plt.subplots(figsize=(9, 7))

ax.scatter(
    np.log2(ma["baseMean"] + 1),
    ma["log2FoldChange"],
    s=8,
    alpha=0.5,
)

ax.axhline(
    0,
    linestyle="--",
    linewidth=1,
)

ax.axhline(
    1,
    linestyle="--",
    linewidth=1,
)

ax.axhline(
    -1,
    linestyle="--",
    linewidth=1,
)

ax.set_xlabel("log2 Mean Expression")
ax.set_ylabel("log2 Fold Change")
ax.set_title("MA Plot — Resistant vs Sensitive")

plt.tight_layout()

ma_file = RESULTS_DIR / "ma_plot.png"
plt.savefig(ma_file, dpi=300)
plt.close()

print(f"Saved: {ma_file}")

results_annotated.to_csv(
    RESULTS_DIR / "DEGs_resistant_vs_sensitive_annotated.csv"
)

significant_annotated = results_annotated.loc[
    significant.index
]

significant_annotated.to_csv(
    RESULTS_DIR / "significant_DEGs_annotated.csv"
)

