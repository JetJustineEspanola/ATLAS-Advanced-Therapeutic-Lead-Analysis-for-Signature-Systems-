# Project Folder Structure

Below is the full project folder structure for this workspace (ATLAS):

```
.
├── README.md
├── SOP.md
├── data
│   ├── processed
│   │   ├── deseq_counts.csv
│   │   ├── gene_counts.csv
│   │   └── sample_metadata.csv
│   └── raw
│       └── tx2gene.tsv
├── notebooks
│   ├── 00_Definition_of_Terms.ipynb
│   ├── 01_stage1_validation.ipynb
│   └── 02_differential_expression.ipynb
├── results
│   ├── differential_expression
│   │   ├── DEGs_resistant_vs_sensitive_annotated.csv
│   │   ├── DEGs_resistant_vs_sensitive.csv
│   │   └── significant_DEGs_annotated.csv
│   └── qc
│       ├── library_sizes.csv
│       ├── pca_coordinates.csv
│       ├── sample_correlation.csv
│       ├── stage1_validation_summary.csv
│       └── transcript_mapping_validation.csv
├── scripts
│   ├── 01_validation.py
│   └── 02_differential_expression.py
└── src

```

Generated on: 2026-08-29
