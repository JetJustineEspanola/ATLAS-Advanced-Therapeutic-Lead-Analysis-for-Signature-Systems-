#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
from pathlib import Path
import pandas as pd
import duckdb

from atlas_data.common import PROJECT_ROOT, CATALOG_PATH, DISCOVERY_DIR

SCHEMA = """
CREATE TABLE IF NOT EXISTS datasets (
    dataset_id VARCHAR PRIMARY KEY,
    source VARCHAR,
    source_accession VARCHAR,
    title VARCHAR,
    summary VARCHAR,
    organism VARCHAR,
    assay_type VARCHAR,
    platform VARCHAR,
    sample_count BIGINT,
    publication_date VARCHAR,
    source_url VARCHAR,
    raw_data_available BOOLEAN,
    processed_data_available BOOLEAN,
    query_family VARCHAR,
    query_text VARCHAR,
    discovered_utc VARCHAR,
    metadata_complete_flag BOOLEAN,
    direct_trastuzumab_resistance BOOLEAN,
    her2_positive_confirmed BOOLEAN,
    resistant_sensitive_groups_defined BOOLEAN,
    biological_replication BOOLEAN,
    complete_sample_metadata BOOLEAN,
    independent_model_or_cohort BOOLEAN,
    pd1_pdl1_or_tgfb_relevance BOOLEAN,
    phenotype_confidence VARCHAR,
    eligibility_score INTEGER,
    eligibility_category VARCHAR,
    eligibility_reasons VARCHAR,
    manual_review_required BOOLEAN
);

CREATE TABLE IF NOT EXISTS samples (
    dataset_id VARCHAR,
    sample_id VARCHAR,
    biological_group VARCHAR,
    resistance_status VARCHAR,
    treatment VARCHAR,
    cell_line VARCHAR,
    patient_id VARCHAR,
    clone_id VARCHAR,
    replicate_id VARCHAR,
    replicate_type VARCHAR,
    qc_status VARCHAR,
    qc_reason VARCHAR
);

CREATE TABLE IF NOT EXISTS files (
    dataset_id VARCHAR,
    file_id VARCHAR,
    file_name VARCHAR,
    file_type VARCHAR,
    local_path VARCHAR,
    remote_url VARCHAR,
    size_bytes BIGINT,
    md5 VARCHAR,
    sha256 VARCHAR,
    download_status VARCHAR,
    source_release VARCHAR
);

CREATE TABLE IF NOT EXISTS processing_runs (
    run_id VARCHAR,
    dataset_id VARCHAR,
    stage VARCHAR,
    tool VARCHAR,
    tool_version VARCHAR,
    parameters_json VARCHAR,
    started_utc VARCHAR,
    ended_utc VARCHAR,
    status VARCHAR,
    log_path VARCHAR
);

CREATE TABLE IF NOT EXISTS qc_results (
    dataset_id VARCHAR,
    sample_id VARCHAR,
    metric VARCHAR,
    value DOUBLE,
    unit VARCHAR,
    qc_status VARCHAR,
    qc_reason VARCHAR
);

CREATE TABLE IF NOT EXISTS validation_results (
    dataset_id VARCHAR,
    validation_type VARCHAR,
    metric VARCHAR,
    value DOUBLE,
    interpretation VARCHAR
);
"""

def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input",
        default=str(DISCOVERY_DIR / "dataset_candidates.csv"),
    )
    args = p.parse_args()

    inp = Path(args.input)
    if not inp.exists():
        print(f"ERROR: input not found: {inp}")
        return 1

    df = pd.read_csv(inp)

    con = duckdb.connect(str(CATALOG_PATH))
    con.execute(SCHEMA)

    # Register dataframe and upsert all discovery rows.
    con.register("incoming_df", df)
    cols = [r[1] for r in con.execute("PRAGMA table_info('datasets')").fetchall()]
    common = [c for c in cols if c in df.columns]

    insert_cols = ", ".join(common)
    select_cols = ", ".join(f'"{c}"' for c in common)

    # DuckDB supports INSERT OR REPLACE.
    con.execute(
        f"INSERT OR REPLACE INTO datasets ({insert_cols}) "
        f"SELECT {select_cols} FROM incoming_df"
    )

    count = con.execute("SELECT COUNT(*) FROM datasets").fetchone()[0]
    by_source = con.execute(
        "SELECT source, COUNT(*) n FROM datasets GROUP BY source ORDER BY n DESC"
    ).fetchdf()

    con.close()

    print("=== 00B COMPLETE ===")
    print(f"Catalog: {CATALOG_PATH}")
    print(f"Datasets: {count}")
    print(by_source.to_string(index=False))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
