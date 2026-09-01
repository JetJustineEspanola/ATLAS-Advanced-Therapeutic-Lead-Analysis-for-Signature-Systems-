#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
OUTDIR = RESULTS / "pipeline_state"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_JSON = OUTDIR / "dataset_volume_summary.json"
OUT_SUMMARY_CSV = OUTDIR / "dataset_volume_summary.csv"
OUT_ACCESSION_CSV = OUTDIR / "dataset_volume_by_accession.csv"
OUT_FILES_CSV = OUTDIR / "dataset_volume_files.csv"
OUT_TXT = OUTDIR / "dataset_volume_summary.txt"

ACCESSION_RE = re.compile(r"\b(?:GSE|GSM|SRR|SRP|PRJNA|PRJEB)\d+\b", re.I)

CANDIDATE_TABLES = [
    DATA / "enriched" / "dataset_candidates_independence_scored.csv",
    DATA / "enriched" / "transcriptomic_validation_candidates.csv",
    DATA / "enriched" / "dataset_candidates.csv",
    DATA / "catalog" / "dataset_candidates.csv",
    DATA / "catalog" / "candidates.csv",
]

SIZE_COLUMN_HINTS = {
    "bytes": ["bytes", "size_bytes", "file_size_bytes", "download_size_bytes", "total_bytes"],
    "mb": ["size_mb", "file_size_mb", "download_size_mb", "total_mb"],
    "gb": ["size_gb", "file_size_gb", "download_size_gb", "total_gb"],
}


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def human_bytes(n):
    n = float(n or 0)
    units = ["B", "KB", "MB", "GB", "TB"]
    i = 0
    while n >= 1024 and i < len(units) - 1:
        n /= 1024
        i += 1
    return f"{n:.2f} {units[i]}"


def find_accession(text):
    m = ACCESSION_RE.search(str(text))
    return m.group(0).upper() if m else ""


def classify_path(path):
    rel = str(path.relative_to(ROOT)).replace("\\", "/")
    if rel.startswith("data/validation_expression/"):
        return "validation_expression"
    if rel.startswith("data/raw/"):
        return "raw_data"
    if rel.startswith("data/processed/"):
        return "processed_data"
    if rel.startswith("data/enriched/"):
        return "metadata_enriched"
    if rel.startswith("data/catalog/"):
        return "metadata_catalog"
    if rel.startswith("data/"):
        return "other_data"
    if rel.startswith("results/external_validation/"):
        return "external_validation_results"
    if rel.startswith("results/cmap/"):
        return "cmap_results"
    return "other"


def scan_files():
    rows = []
    roots = [DATA, RESULTS / "external_validation", RESULTS / "cmap"]

    for base in roots:
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if not path.is_file():
                continue
            try:
                size = path.stat().st_size
            except OSError:
                continue

            rel = path.relative_to(ROOT)
            rows.append({
                "path": str(rel),
                "category": classify_path(path),
                "accession": find_accession(str(rel)),
                "bytes": int(size),
                "mb": size / (1024 ** 2),
                "gb": size / (1024 ** 3),
                "extension": "".join(path.suffixes).lower(),
            })

    return pd.DataFrame(rows)


def load_candidate_metadata():
    frames = []
    for path in CANDIDATE_TABLES:
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        df = df.copy()
        df["_source_table"] = str(path.relative_to(ROOT))
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    return pd.concat(frames, ignore_index=True, sort=False)


def accession_column(df):
    for c in [
        "source_accession",
        "accession",
        "dataset_accession",
        "geo_accession",
        "series_accession",
        "dataset_id",
    ]:
        if c in df.columns:
            return c
    return None


def remote_size_from_catalog(df):
    if df.empty:
        return 0, []

    cols_lower = {str(c).lower(): c for c in df.columns}
    used = []
    total_bytes = 0.0

    for unit, hints in SIZE_COLUMN_HINTS.items():
        for hint in hints:
            if hint not in cols_lower:
                continue
            col = cols_lower[hint]
            vals = pd.to_numeric(df[col], errors="coerce").dropna()
            if vals.empty:
                continue

            if unit == "bytes":
                total_bytes += vals.sum()
            elif unit == "mb":
                total_bytes += vals.sum() * (1024 ** 2)
            elif unit == "gb":
                total_bytes += vals.sum() * (1024 ** 3)

            used.append(col)

    return int(total_bytes), used


def main():
    files = scan_files()
    candidates = load_candidate_metadata()

    if files.empty:
        files = pd.DataFrame(columns=[
            "path", "category", "accession", "bytes", "mb", "gb", "extension"
        ])

    local_total = int(files["bytes"].sum())
    data_only = files[files["path"].str.startswith("data/")]
    local_data_total = int(data_only["bytes"].sum()) if not data_only.empty else 0

    comparison_files = files[
        files["category"].isin(["validation_expression", "raw_data", "processed_data"])
    ]
    comparison_bytes = int(comparison_files["bytes"].sum()) if not comparison_files.empty else 0

    discovered_accessions = set()
    candidate_rows = 0
    source_tables = []

    if not candidates.empty:
        candidate_rows = len(candidates)
        source_tables = sorted(candidates["_source_table"].dropna().astype(str).unique())
        acc_col = accession_column(candidates)
        if acc_col:
            for value in candidates[acc_col].dropna().astype(str):
                acc = find_accession(value) or value.strip()
                if acc:
                    discovered_accessions.add(acc)

    downloaded_accessions = set(
        files.loc[files["accession"].ne(""), "accession"].astype(str)
    )

    remote_bytes, remote_size_columns = remote_size_from_catalog(candidates)

    cat = (
        files.groupby("category", dropna=False)["bytes"]
        .agg(["count", "sum"])
        .reset_index()
        .rename(columns={"count": "file_count", "sum": "bytes"})
    )

    if not cat.empty:
        cat["mb"] = cat["bytes"] / (1024 ** 2)
        cat["gb"] = cat["bytes"] / (1024 ** 3)
        cat["human_size"] = cat["bytes"].map(human_bytes)

    accession_rows = []

    local_acc = files[files["accession"].ne("")]
    if not local_acc.empty:
        grouped = local_acc.groupby("accession")["bytes"].agg(["count", "sum"]).reset_index()
        for _, r in grouped.iterrows():
            acc = str(r["accession"])
            b = int(r["sum"])
            accession_rows.append({
                "accession": acc,
                "discovered_in_metadata": acc in discovered_accessions,
                "local_file_count": int(r["count"]),
                "local_bytes": b,
                "local_mb": b / (1024 ** 2),
                "local_gb": b / (1024 ** 3),
                "local_human_size": human_bytes(b),
            })

    seen = {r["accession"] for r in accession_rows}
    for acc in sorted(discovered_accessions - seen):
        accession_rows.append({
            "accession": acc,
            "discovered_in_metadata": True,
            "local_file_count": 0,
            "local_bytes": 0,
            "local_mb": 0.0,
            "local_gb": 0.0,
            "local_human_size": "0.00 B",
        })

    by_accession = pd.DataFrame(accession_rows)

    summary_df = pd.DataFrame([
        {
            "metric": "candidate_metadata_rows_scanned",
            "value": candidate_rows,
            "unit": "rows",
            "note": "May contain the same dataset in multiple metadata tables.",
        },
        {
            "metric": "unique_dataset_accessions_discovered",
            "value": len(discovered_accessions),
            "unit": "datasets",
            "note": "Unique accessions sifted at metadata level.",
        },
        {
            "metric": "unique_accessions_with_local_files",
            "value": len(downloaded_accessions),
            "unit": "datasets",
            "note": "Accessions inferred from local file paths.",
        },
        {
            "metric": "local_data_directory_bytes",
            "value": local_data_total,
            "unit": "bytes",
            "note": human_bytes(local_data_total),
        },
        {
            "metric": "comparison_download_bytes",
            "value": comparison_bytes,
            "unit": "bytes",
            "note": human_bytes(comparison_bytes),
        },
        {
            "metric": "tracked_local_bytes_data_plus_validation_results",
            "value": local_total,
            "unit": "bytes",
            "note": human_bytes(local_total),
        },
        {
            "metric": "catalog_advertised_remote_bytes",
            "value": remote_bytes,
            "unit": "bytes",
            "note": human_bytes(remote_bytes) if remote_bytes else "Unavailable in current metadata schema",
        },
    ])

    files.to_csv(OUT_FILES_CSV, index=False)
    by_accession.to_csv(OUT_ACCESSION_CSV, index=False)
    summary_df.to_csv(OUT_SUMMARY_CSV, index=False)

    report = {
        "created_utc": utcnow(),
        "candidate_tables_scanned": source_tables,
        "candidate_metadata_rows_scanned": candidate_rows,
        "unique_dataset_accessions_discovered": len(discovered_accessions),
        "unique_accessions_with_local_files": len(downloaded_accessions),
        "local_data_directory": {
            "bytes": local_data_total,
            "mb": local_data_total / (1024 ** 2),
            "gb": local_data_total / (1024 ** 3),
            "human": human_bytes(local_data_total),
        },
        "comparison_downloads": {
            "bytes": comparison_bytes,
            "mb": comparison_bytes / (1024 ** 2),
            "gb": comparison_bytes / (1024 ** 3),
            "human": human_bytes(comparison_bytes),
        },
        "tracked_local_total": {
            "bytes": local_total,
            "mb": local_total / (1024 ** 2),
            "gb": local_total / (1024 ** 3),
            "human": human_bytes(local_total),
        },
        "catalog_advertised_remote_size": {
            "bytes": remote_bytes,
            "mb": remote_bytes / (1024 ** 2),
            "gb": remote_bytes / (1024 ** 3),
            "human": human_bytes(remote_bytes) if remote_bytes else None,
            "source_columns": remote_size_columns,
            "note": (
                "Only populated when metadata contains explicit size columns. "
                "Metadata discovery does not mean the full remote dataset was downloaded."
            ),
        },
    }

    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")

    lines = [
        "=" * 88,
        "ATLAS — DATASET VOLUME / STORAGE METADATA",
        "=" * 88,
        f"Candidate metadata rows scanned:   {candidate_rows:,}",
        f"Unique dataset accessions sifted:  {len(discovered_accessions):,}",
        f"Accessions with local files:       {len(downloaded_accessions):,}",
        "",
        f"Local data/ size:                  {human_bytes(local_data_total)}",
        f"Downloaded comparison data:        {human_bytes(comparison_bytes)}",
        f"Tracked local data + result files: {human_bytes(local_total)}",
        "",
        (
            f"Catalog-advertised remote size:     {human_bytes(remote_bytes)}"
            if remote_bytes
            else "Catalog-advertised remote size:     unavailable in current metadata"
        ),
        "",
        "Important:",
        "- 'Sifted' means metadata records ATLAS examined, not bytes downloaded.",
        "- Local byte totals are exact filesystem sizes.",
        "- Remote totals are reported only if the source metadata already includes explicit sizes.",
        "",
        "Largest local categories:",
    ]

    if not cat.empty:
        for _, r in cat.sort_values("bytes", ascending=False).head(10).iterrows():
            lines.append(
                f"- {r['category']}: {human_bytes(int(r['bytes']))} "
                f"across {int(r['file_count']):,} files"
            )
    else:
        lines.append("- none")

    lines.extend([
        "",
        "Outputs:",
        str(OUT_JSON),
        str(OUT_SUMMARY_CSV),
        str(OUT_ACCESSION_CSV),
        str(OUT_FILES_CSV),
    ])

    OUT_TXT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
