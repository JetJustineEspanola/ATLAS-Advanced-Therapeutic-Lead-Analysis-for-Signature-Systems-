# ATLAS 00D — Metadata Enrichment + Download Planning

Run after 00A–00C.

```bash
source .venv/bin/activate

python -u scripts/00d_metadata_enrichment.py \
  --min-score 40 \
  --sources GEO,SRA \
  --max-datasets 30 \
  --workers 6
```

Then re-score with sample-level evidence:

```bash
python -u scripts/00c2_rescore_after_enrichment.py
```

Outputs:

```text
data/enriched/
├── sample_metadata.csv
├── sample_metadata.parquet
├── download_manifest.csv
├── download_manifest.parquet
├── dataset_enrichment_summary.csv
└── dataset_candidates_rescored.csv
```

00D does **not download FASTQ files yet**. It produces a download plan first.

This is intentional: ATLAS should validate phenotype and replicate structure before
large raw-data acquisition.

Next stage:
- 00E download manager + checksums
- then platform-specific QC / harmonization
