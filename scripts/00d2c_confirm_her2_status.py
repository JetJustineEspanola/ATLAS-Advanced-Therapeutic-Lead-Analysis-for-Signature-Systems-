#!/usr/bin/env python3
from pathlib import Path
import duckdb

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "data/catalog/atlas_catalog.duckdb"

updates = {
    "GEO:GSE121105": True,
    "GEO:GSE123754": True,
}

con = duckdb.connect(str(DB))

for dataset_id, value in updates.items():
    con.execute(
        '''
        UPDATE datasets
        SET her2_positive_confirmed = ?
        WHERE dataset_id = ?
        ''',
        [bool(value), dataset_id],
    )

df = con.execute(
    '''
    SELECT
        dataset_id,
        her2_positive_confirmed,
        direct_trastuzumab_resistance,
        resistant_sensitive_groups_defined,
        biological_replication,
        phenotype_confidence
    FROM datasets
    WHERE dataset_id IN ('GEO:GSE121105','GEO:GSE123754')
    ORDER BY dataset_id
    '''
).fetchdf()

con.close()

print("Updated HER2 confirmation flags:")
print(df.to_string(index=False))
