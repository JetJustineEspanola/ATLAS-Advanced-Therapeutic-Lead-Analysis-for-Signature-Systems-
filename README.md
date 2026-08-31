# ATLAS Data Acquisition Phase 1 — 00A to 00C

This package adds a scalable dataset-discovery and validation layer before the
existing ATLAS pipeline.

## Stages

### 00A Dataset Discovery
Queries:
- NCBI GEO DataSets
- NCBI SRA
- NCI GDC / TCGA-BRCA transcriptome files
- cBioPortal studies

Outputs:
- `data/discovery/dataset_candidates.csv`
- `data/discovery/dataset_candidates.parquet` when a Parquet engine is available
- `data/manifests/00a_discovery_manifest.json`

### 00B Metadata Catalog
Creates:
- `data/catalog/atlas_catalog.duckdb`

Tables:
- datasets
- samples
- files
- processing_runs
- qc_results
- validation_results

### 00C Dataset Eligibility
Applies the ATLAS transparent eligibility rubric and produces:
- `data/discovery/dataset_candidates_scored.csv`
- `data/discovery/dataset_candidates_scored.parquet`

## Install

From the ATLAS project root:

```bash
source .venv/bin/activate
python -m pip install -U duckdb pyarrow requests pandas
```

Copy the package contents into the project root.

## Run

```bash
python -u scripts/run_data_acquisition_phase1.py \
  --sources geo,sra,gdc,cbioportal \
  --retmax 25 \
  --query-family all
```

For a faster first test:

```bash
python -u scripts/run_data_acquisition_phase1.py \
  --sources geo,sra \
  --retmax 10 \
  --query-family trastuzumab_resistance
```

## Important limitation

00C is deliberately conservative. Search-result metadata cannot reliably infer:
- true biological replicate structure,
- resistant-vs-sensitive group definitions,
- exact trastuzumab exposure,
- clone structure,
- complete sample phenotype.

Those fields remain `False` or unresolved until a later metadata-enrichment /
manual-review stage. This prevents search text from being mistaken for strong
experimental evidence.

## Next build

00D should:
- fetch full study/sample metadata for shortlisted accessions,
- populate the `samples` and `files` tables,
- determine biological vs technical replicates,
- generate resumable download manifests,
- calculate checksums,
- route datasets to raw or processed acquisition workflows.
