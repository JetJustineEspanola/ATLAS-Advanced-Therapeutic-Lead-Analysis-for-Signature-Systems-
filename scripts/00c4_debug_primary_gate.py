#!/usr/bin/env python3
from pathlib import Path
import duckdb

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "data/catalog/atlas_catalog.duckdb"

TARGETS = [
    "GEO:GSE121105",
    "GEO:GSE237606",
    "GEO:GSE114575",
    "GEO:GSE123754",
]

con = duckdb.connect(str(DB))
df = con.execute(
    '''
    SELECT
        dataset_id,
        source,
        source_accession,
        title,
        direct_trastuzumab_resistance,
        her2_positive_confirmed,
        resistant_sensitive_groups_defined,
        biological_replication,
        raw_data_available,
        complete_sample_metadata,
        independent_model_or_cohort,
        phenotype_confidence,
        sample_count
    FROM datasets
    WHERE dataset_id IN (?, ?, ?, ?)
    ORDER BY dataset_id
    ''',
    TARGETS,
).fetchdf()
con.close()

print("=" * 100)
print("ATLAS — 00C4 PRIMARY GATE DEBUG")
print("=" * 100)
print(df.to_string(index=False))

print("\nPrimary gate requirements:")
print("direct_trastuzumab_resistance == True")
print("her2_positive_confirmed == True")
print("resistant_sensitive_groups_defined == True")
print("biological_replication == True")
print("plus transcriptomic / non-umbrella / non-SRR conditions in 00C4")
