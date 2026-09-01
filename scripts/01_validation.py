"""
01_validation.py
│
├── Check input files
├── Read Salmon files
├── Validate transcript IDs
├── Read tx2gene
├── Validate mapping uniqueness
├── Calculate transcript match rate
├── Construct gene-level matrix
├── Check library sizes
├── Check sample correlations
├── Perform PCA
├── Final QC
└── Save gene_counts.csv
"""

from pathlib import Path
import gzip

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA


# ============================================================
# ATLAS — Stage 1: Data Validation
# ============================================================

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

# This file is:
# ATLAS/scripts/01_validation.py
#
# parents[1] therefore points to:
# ATLAS/

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DIR = PROJECT_ROOT / "data" / "raw"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
QC_DIR = RESULTS_DIR / "qc"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# 2. Define samples
# ------------------------------------------------------------

samples = {
    "TS1": RAW_DIR / "GSM9067960_TS1_quant.sf.gz",
    "TS2": RAW_DIR / "GSM9067961_TS2_quant.sf.gz",
    "TS3": RAW_DIR / "GSM9067962_TS3_quant.sf.gz",
    "TR1": RAW_DIR / "GSM9067963_TR1_quant.sf.gz",
    "TR2": RAW_DIR / "GSM9067964_TR2_quant.sf.gz",
    "TR3": RAW_DIR / "GSM9067965_TR3_quant.sf.gz",
}

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


# ------------------------------------------------------------
# 3. Utility functions
# ------------------------------------------------------------

def strip_version(ids: pd.Series) -> pd.Series:
    """
    Remove Ensembl transcript version suffixes.

    Example:
        ENST00000456328.2 -> ENST00000456328
    """
    return (
        ids.astype(str)
        .str.split(".")
        .str[0]
    )


def load_quant_file(filepath: Path) -> pd.DataFrame:
    """Load a Salmon quant.sf.gz file."""
    return pd.read_csv(
        filepath,
        sep="\t",
        compression="gzip",
    )


def load_gene_counts(
    sample_name: str,
    filepath: Path,
    tx2gene_mapping: pd.DataFrame,
    mapping_set: set,
) -> tuple[pd.Series, dict]:
    """
    Convert transcript-level Salmon counts into gene-level counts.
    """

    quant = pd.read_csv(
        filepath,
        sep="\t",
        compression="gzip",
        usecols=["Name", "NumReads"],
    )

    # Remove transcript version suffix
    quant["transcript_id"] = strip_version(quant["Name"])

    # Check transcript mapping
    matched = quant["transcript_id"].isin(mapping_set)

    total_transcripts = len(quant)
    matched_transcripts = int(matched.sum())
    unmatched_transcripts = total_transcripts - matched_transcripts

    match_rate = (
        matched_transcripts / total_transcripts * 100
        if total_transcripts > 0
        else 0
    )

    print(
        f"{sample_name}: "
        f"{matched_transcripts:,}/{total_transcripts:,} transcripts mapped "
        f"({match_rate:.2f}%)"
    )

    # Join transcript IDs to gene IDs
    quant = quant.merge(
        tx2gene_mapping,
        left_on="transcript_id",
        right_on="Transcript stable ID",
        how="left",
    )

    # Keep transcripts that have a gene mapping
    quant = quant.dropna(subset=["Gene stable ID"])

    # Aggregate transcript counts to gene counts
    gene_counts = (
        quant.groupby("Gene stable ID")["NumReads"]
        .sum()
        .rename(sample_name)
    )

    validation = {
        "sample": sample_name,
        "total_transcripts": total_transcripts,
        "matched_transcripts": matched_transcripts,
        "unmatched_transcripts": unmatched_transcripts,
        "match_rate_percent": match_rate,
        "gene_count_rows": len(gene_counts),
    }

    return gene_counts, validation


# ============================================================
# 4. Input validation
# ============================================================

print("=" * 60)
print("ATLAS — Stage 1 Data Validation")
print("=" * 60)

print(f"\nProject root: {PROJECT_ROOT}")
print(f"Raw data:     {RAW_DIR}")


expected_files = list(samples.values()) + [RAW_DIR / "tx2gene.tsv"]

print("\nChecking input files...")

missing_files = []

for filepath in expected_files:
    if filepath.exists():
        print(f"OK       {filepath.name}")
    else:
        print(f"MISSING  {filepath.name}")
        missing_files.append(filepath)

if missing_files:
    raise FileNotFoundError(
        "One or more required input files are missing."
    )

print("\nInput validation PASSED.")
 

# ============================================================
# 5. Validate metadata
# ============================================================

print("\nChecking sample metadata...")

assert list(samples.keys()) == list(metadata.index), (
    "Sample IDs do not match metadata."
)

print(metadata)
print("\nMetadata validation PASSED.")


# ============================================================
# 6. Load and validate tx2gene mapping
# ============================================================

print("\nLoading tx2gene mapping...")

tx2gene = pd.read_csv(
    RAW_DIR / "tx2gene.tsv",
    sep="\t",
)

required_columns = [
    "Gene stable ID",
    "Transcript stable ID",
    "Gene name",
]

missing_columns = [
    col for col in required_columns
    if col not in tx2gene.columns
]

if missing_columns:
    raise ValueError(
        f"tx2gene.tsv is missing columns: {missing_columns}"
    )

print(f"Mapping rows: {len(tx2gene):,}")
print("Columns:", tx2gene.columns.tolist())


# Clean transcript IDs
tx2gene_clean = tx2gene[
    ["Gene stable ID", "Transcript stable ID"]
].copy()

tx2gene_clean["Transcript stable ID"] = strip_version(
    tx2gene_clean["Transcript stable ID"]
)

# Check duplicate transcript mappings
total_mapping_rows = len(tx2gene_clean)
unique_transcripts = tx2gene_clean["Transcript stable ID"].nunique()
duplicate_rows = total_mapping_rows - unique_transcripts

print(f"Unique transcript IDs: {unique_transcripts:,}")
print(f"Duplicate transcript rows: {duplicate_rows:,}")

if duplicate_rows != 0:
    raise ValueError(
        "Duplicate transcript IDs detected in tx2gene mapping."
    )

# Check whether a transcript maps to multiple genes
transcript_gene_counts = (
    tx2gene_clean
    .groupby("Transcript stable ID")["Gene stable ID"]
    .nunique()
)

multi_gene_transcripts = transcript_gene_counts[
    transcript_gene_counts > 1
]

print(
    "Transcripts mapping to multiple genes:",
    len(multi_gene_transcripts)
)

if len(multi_gene_transcripts) != 0:
    raise ValueError(
        "At least one transcript maps to multiple genes."
    )

mapping_set = set(
    tx2gene_clean["Transcript stable ID"]
)

print("\nMapping validation PASSED.")


# ============================================================
# 7. Build gene-level count matrix
# ============================================================

print("\nBuilding gene-level count matrix...")

gene_count_series = []
validation_results = []

for sample_name, filepath in samples.items():

    gene_counts, validation = load_gene_counts(
        sample_name=sample_name,
        filepath=filepath,
        tx2gene_mapping=tx2gene_clean,
        mapping_set=mapping_set,
    )

    gene_count_series.append(gene_counts)
    validation_results.append(validation)


gene_counts = (
    pd.concat(gene_count_series, axis=1)
    .fillna(0)
)

gene_counts.index.name = "Gene stable ID"

print("\nGene count matrix created.")
print(f"Genes:   {gene_counts.shape[0]:,}")
print(f"Samples: {gene_counts.shape[1]}")


# ============================================================
# 8. Matrix validation
# ============================================================

print("\nValidating gene count matrix...")

duplicate_genes = gene_counts.index.duplicated().sum()
missing_values = gene_counts.isna().sum().sum()

expected_samples = [
    "TS1", "TS2", "TS3",
    "TR1", "TR2", "TR3"
]

assert duplicate_genes == 0, "Duplicate gene IDs detected."
assert missing_values == 0, "Missing values detected."
assert list(gene_counts.columns) == expected_samples, (
    "Unexpected sample columns."
)

print(f"Duplicate gene IDs: {duplicate_genes}")
print(f"Missing values:     {missing_values}")
print(f"Sample columns:     {gene_counts.columns.tolist()}")

print("\nGene matrix validation PASSED.")


# ============================================================
# 9. Save validation report
# ============================================================

validation_df = pd.DataFrame(validation_results)

validation_report = QC_DIR / "transcript_mapping_validation.csv"

validation_df.to_csv(
    validation_report,
    index=False,
)

print(f"\nSaved: {validation_report}")


# ============================================================
# 10. Library-size QC
# ============================================================

print("\nCalculating library sizes...")

library_sizes = (
    gene_counts.sum()
    .sort_values(ascending=False)
)

print("\nTotal estimated counts per sample:")
print(library_sizes)


library_sizes_path = QC_DIR / "library_sizes.csv"

library_sizes.rename(
    "total_estimated_counts"
).to_csv(
    library_sizes_path,
    header=True,
)

# Plot
fig, ax = plt.subplots(figsize=(8, 5))

library_sizes.plot(
    kind="bar",
    ax=ax,
)

ax.set_ylabel("Total estimated counts")
ax.set_xlabel("Sample")
ax.set_title("Library Size by Sample")
ax.tick_params(axis="x", rotation=0)

plt.tight_layout()

library_plot = QC_DIR / "library_sizes.png"
plt.savefig(library_plot, dpi=300)
plt.close()

print(f"Saved: {library_plot}")


# ============================================================
# 11. Sample correlation QC
# ============================================================

print("\nCalculating sample correlations...")

log_counts = np.log2(gene_counts + 1)

correlation = log_counts.corr()

print("\nSample correlation matrix:")
print(correlation.round(6))

correlation_path = QC_DIR / "sample_correlation.csv"

correlation.to_csv(correlation_path)


# Correlation heatmap
fig, ax = plt.subplots(figsize=(7, 6))

im = ax.imshow(
    correlation,
    vmin=0,
    vmax=1,
)

ax.set_xticks(range(len(correlation.columns)))
ax.set_yticks(range(len(correlation.index)))

ax.set_xticklabels(correlation.columns)
ax.set_yticklabels(correlation.index)

for i in range(len(correlation.index)):
    for j in range(len(correlation.columns)):
        ax.text(
            j,
            i,
            f"{correlation.iloc[i, j]:.3f}",
            ha="center",
            va="center",
        )

plt.colorbar(
    im,
    ax=ax,
    label="Pearson correlation",
)

ax.set_title("Sample-to-Sample Correlation")

plt.tight_layout()

correlation_plot = QC_DIR / "sample_correlation.png"
plt.savefig(correlation_plot, dpi=300)
plt.close()

print(f"Saved: {correlation_plot}")


# ============================================================
# 12. PCA
# ============================================================

print("\nRunning PCA...")

# Samples must be rows and genes must be columns
X = log_counts.T

pca = PCA(n_components=2)

pca_result = pca.fit_transform(X)

pc1_variance = pca.explained_variance_ratio_[0]
pc2_variance = pca.explained_variance_ratio_[1]

print(f"PC1 explained variance: {pc1_variance:.6f}")
print(f"PC2 explained variance: {pc2_variance:.6f}")
print(
    f"Total explained variance: "
    f"{pc1_variance + pc2_variance:.6f}"
)

pca_df = pd.DataFrame(
    pca_result,
    columns=["PC1", "PC2"],
    index=gene_counts.columns,
)

pca_df["Condition"] = metadata.loc[
    pca_df.index,
    "condition"
]

print("\nPCA coordinates:")
print(pca_df)


# Save PCA coordinates
pca_coordinates = QC_DIR / "pca_coordinates.csv"

pca_df.to_csv(pca_coordinates)

# PCA plot
fig, ax = plt.subplots(figsize=(8, 6))

for condition in ["Sensitive", "Resistant"]:

    subset = pca_df[
        pca_df["Condition"] == condition
    ]

    ax.scatter(
        subset["PC1"],
        subset["PC2"],
        label=condition,
        s=80,
    )

    for sample, row in subset.iterrows():
        ax.annotate(
            sample,
            (row["PC1"], row["PC2"]),
            xytext=(5, 5),
            textcoords="offset points",
        )


ax.set_xlabel(
    f"PC1 ({pc1_variance * 100:.2f}% variance)"
)

ax.set_ylabel(
    f"PC2 ({pc2_variance * 100:.2f}% variance)"
)

ax.set_title("PCA of Gene Expression")
ax.legend()

plt.tight_layout()

pca_plot = QC_DIR / "pca.png"
plt.savefig(pca_plot, dpi=300)
plt.close()

print(f"Saved: {pca_plot}")


# ============================================================
# 13. Final matrix statistics
# ============================================================

zero_count_genes = (
    gene_counts.sum(axis=1) == 0
).sum()

detected_genes = (
    gene_counts.sum(axis=1) > 0
).sum()


print("\nFinal gene-count matrix statistics:")
print(f"Genes:                    {gene_counts.shape[0]:,}")
print(f"Samples:                  {gene_counts.shape[1]}")
print(f"Genes with zero counts:   {zero_count_genes:,}")
print(f"Genes with >=1 count:     {detected_genes:,}")
print(f"Missing values:           {missing_values}")
print(f"Duplicate genes:          {duplicate_genes}")


# ============================================================
# 14. Save gene count matrix
# ============================================================

gene_counts_file = (
    PROCESSED_DIR / "gene_counts.csv"
)

gene_counts.to_csv(
    gene_counts_file
)

print(f"\nSaved gene-count matrix: {gene_counts_file}")


# ============================================================
# 15. Reload the saved matrix
# ============================================================

print("\nReloading saved gene-count matrix...")

check_counts = pd.read_csv(
    gene_counts_file,
    index_col=0,
)

reloaded_missing = check_counts.isna().sum().sum()
reloaded_duplicates = check_counts.index.duplicated().sum()

print(f"Reloaded shape:     {check_counts.shape}")
print(f"Missing values:     {reloaded_missing}")
print(f"Duplicate gene IDs: {reloaded_duplicates}")

assert check_counts.shape == gene_counts.shape
assert reloaded_missing == 0
assert reloaded_duplicates == 0

print("\nSaved-file validation PASSED.")


# ============================================================
# 16. Final report
# ============================================================

summary = pd.DataFrame(
    {
        "check": [
            "Salmon files found",
            "tx2gene mapping rows",
            "Transcript match rate",
            "Duplicate transcript mappings",
            "Transcripts mapping to multiple genes",
            "Gene matrix rows",
            "Gene matrix columns",
            "Missing gene counts",
            "Duplicate gene IDs",
            "Genes with zero counts",
            "Genes with detected counts",
            "PCA PC1 variance",
            "PCA PC2 variance",
            "PCA PC1+PC2 variance",
        ],
        "result": [
            len(samples),
            len(tx2gene),
            f"{validation_df['match_rate_percent'].min():.2f}%",
            duplicate_rows,
            len(multi_gene_transcripts),
            gene_counts.shape[0],
            gene_counts.shape[1],
            missing_values,
            duplicate_genes,
            zero_count_genes,
            detected_genes,
            f"{pc1_variance * 100:.2f}%",
            f"{pc2_variance * 100:.2f}%",
            f"{(pc1_variance + pc2_variance) * 100:.2f}%",
        ],
    }
)

summary_file = QC_DIR / "stage1_validation_summary.csv"

summary.to_csv(
    summary_file,
    index=False,
)

print(f"\nSaved: {summary_file}")


# ============================================================
# DONE
# ============================================================

print("\n" + "=" * 60)
print("STAGE 1 VALIDATION COMPLETE")
print("=" * 60)

print("\nOutputs:")
print(f"  Gene counts:      {gene_counts_file}")
print(f"  Validation report: {validation_report}")
print(f"  QC results:        {QC_DIR}")
print("\nStage 1 PASSED.")