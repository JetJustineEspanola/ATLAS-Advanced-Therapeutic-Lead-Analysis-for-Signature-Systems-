#!/usr/bin/env python3
"""
ATLAS — Stage 03B: External Dataset Validation

Goal
----
Validate the discovery trastuzumab-resistance signature against independent
HER2+ trastuzumab-resistance transcriptomic datasets.

This script is deliberately split into reproducible sub-steps:
  03B.1 Dataset registry
  03B.2 Download GEO metadata / processed data
  03B.3 Build sample-group manifests
  03B.4 Analyze each external dataset independently
  03B.5 Compare external DE direction/effect with the ATLAS discovery DEGs
  03B.6 Write a cross-dataset validation summary

Important
---------
Raw expression matrices from different technologies are NOT merged together.
Each external dataset is analyzed separately and only the gene-level results
are compared across datasets.

Known initial validation datasets
---------------------------------
GSE89216  — BT-474 trastuzumab-sensitive vs resistant, Affymetrix array
GSE15043  — BT-474 parental vs Herceptin-resistant clones, Affymetrix array
GSE55005  — BT-474 vs BTR50 resistant, RNA-seq processed comparison file

Usage
-----
python scripts/03b_external_validation.py

Optional:
python scripts/03b_external_validation.py --download-only
python scripts/03b_external_validation.py --analyze-only
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import math
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests
from scipy import stats

try:
    import GEOparse
except ImportError:
    GEOparse = None

try:
    import xlrd
except ImportError:
    xlrd = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DISCOVERY_DEG_FILE = (
    PROJECT_ROOT
    / "results"
    / "differential_expression"
    / "DEGs_resistant_vs_sensitive_annotated.csv"
)

OUT_ROOT = PROJECT_ROOT / "results" / "external_validation"
RAW_DIR = OUT_ROOT / "raw"
PROCESSED_DIR = OUT_ROOT / "processed"

REGISTRY_FILE = OUT_ROOT / "dataset_registry.csv"
CROSS_GENE_FILE = OUT_ROOT / "ATLAS_external_gene_validation.csv"
SUMMARY_FILE = OUT_ROOT / "ATLAS_external_validation_summary.csv"
METADATA_FILE = OUT_ROOT / "ATLAS_external_validation_metadata.json"

for d in [OUT_ROOT, RAW_DIR, PROCESSED_DIR]:
    d.mkdir(parents=True, exist_ok=True)


SESSION = requests.Session()
SESSION.headers.update(
    {"User-Agent": "ATLAS-External-Validation/1.0"}
)

REQUEST_TIMEOUT = 60
MAX_RETRIES = 3


# ------------------------------------------------------------------
# Curated registry.
#
# This is intentionally explicit and auditable. Dataset discovery can be
# automated later, but inclusion/exclusion should remain transparent.
# ------------------------------------------------------------------

DATASETS = [
    {
        "accession": "GSE89216",
        "title": "Identification of resistance mechanisms to anti-HER2 antibody in breast cancer cell line BT-474",
        "organism": "Homo sapiens",
        "model": "BT-474",
        "platform_type": "microarray",
        "platform": "Affymetrix Human Gene 2.0 ST Array",
        "contrast": "resistant untreated vs sensitive untreated",
        "resistant_samples": "GSM2361001;GSM2361002",
        "sensitive_samples": "GSM2361005;GSM2361006",
        "series_matrix_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89216/matrix/GSE89216_series_matrix.txt.gz",
        "processed_comparison_url": "",
        "include": True,
        "priority": "primary",
        "notes": "Direct BT-474 acquired trastuzumab-resistance comparison; untreated resistant and untreated sensitive samples.",
    },
    {
        "accession": "GSE15043",
        "title": "Gene expression profiles of Herceptin-resistant breast cancer cells",
        "organism": "Homo sapiens",
        "model": "BT-474",
        "platform_type": "microarray",
        "platform": "Affymetrix Human Genome U133 Plus 2.0 Array",
        "contrast": "four resistant clones vs parental BT474",
        "resistant_samples": "GSM375721;GSM375722;GSM375723;GSM375724;GSM375725;GSM375726;GSM375727;GSM375728",
        "sensitive_samples": "GSM375719;GSM375720",
        "series_matrix_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15043/matrix/GSE15043_series_matrix.txt.gz",
        "processed_comparison_url": "",
        "include": True,
        "priority": "primary",
        "notes": "Direct parental BT474 vs four independent Herceptin-resistant clones.",
    },
    {
        "accession": "GSE55005",
        "title": "mRNA profiling reveals determinants of trastuzumab efficiency in HER2-positive breast cancer",
        "organism": "Homo sapiens",
        "model": "BT-474 / BTR50",
        "platform_type": "rna_seq_processed",
        "platform": "Illumina HiSeq 2000",
        "contrast": "BTR50 resistant vs BT474 parental",
        "resistant_samples": "GSM1327857",
        "sensitive_samples": "GSM1327853",
        "series_matrix_url": "",
        "processed_comparison_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55005/suppl/GSE55005_nexprs_btr50_bt474.xls.gz",
        "include": True,
        "priority": "primary",
        "notes": "Processed study-provided BTR50 vs BT474 comparison. Single sample per state, so use study-provided effect/statistics rather than a new replicate-level t-test.",
    },
]


def safe_get(url: str) -> bytes | None:
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = SESSION.get(url, timeout=REQUEST_TIMEOUT)
        except requests.RequestException:
            if attempt == MAX_RETRIES:
                return None
            time.sleep(attempt)
            continue

        if r.status_code == 200:
            return r.content

        if r.status_code in {429, 500, 502, 503, 504}:
            if attempt == MAX_RETRIES:
                return None
            time.sleep(attempt * 2)
            continue

        return None
    return None


def write_registry() -> pd.DataFrame:
    registry = pd.DataFrame(DATASETS)
    registry.to_csv(REGISTRY_FILE, index=False)
    return registry


def download_dataset(row: pd.Series) -> dict[str, Any]:
    acc = row["accession"]
    ds_dir = RAW_DIR / acc
    ds_dir.mkdir(parents=True, exist_ok=True)

    status = {
        "accession": acc,
        "series_matrix_downloaded": False,
        "processed_comparison_downloaded": False,
        "download_error": "",
    }

    if row.get("series_matrix_url"):
        path = ds_dir / f"{acc}_series_matrix.txt.gz"
        if not path.exists():
            content = safe_get(str(row["series_matrix_url"]))
            if content is None:
                status["download_error"] += "series_matrix_failed;"
            else:
                path.write_bytes(content)
        status["series_matrix_downloaded"] = path.exists()

    if row.get("processed_comparison_url"):
        path = ds_dir / Path(str(row["processed_comparison_url"])).name
        if not path.exists():
            content = safe_get(str(row["processed_comparison_url"]))
            if content is None:
                status["download_error"] += "processed_comparison_failed;"
            else:
                path.write_bytes(content)
        status["processed_comparison_downloaded"] = path.exists()

    return status


def read_geo_series_matrix(path: Path) -> tuple[pd.DataFrame, list[str]]:
    """
    Read the numerical table from a GEO series-matrix file.
    Returns expression matrix indexed by ID_REF and sample GSM columns.
    """
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    start = None
    end = None
    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            start = i + 1
        elif line.startswith("!series_matrix_table_end"):
            end = i
            break

    if start is None or end is None:
        raise ValueError(f"Could not find GEO matrix table markers in {path}")

    table_text = "".join(lines[start:end])
    df = pd.read_csv(io.StringIO(table_text), sep="\t", quotechar='"')
    df.columns = [str(c).strip('"') for c in df.columns]
    id_col = df.columns[0]
    df[id_col] = df[id_col].astype(str).str.strip('"')
    df = df.set_index(id_col)

    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    return df, list(df.columns)


def read_geo_platform_annotation(accession: str, platform_id: str) -> pd.DataFrame | None:
    """Retrieve and cache GEO platform annotation robustly."""
    cache_dir = RAW_DIR / "_geo_platform_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{platform_id}.annotation.tsv.gz"

    if cache_file.exists():
        try:
            x = pd.read_csv(cache_file, sep="\t", dtype=str, compression="gzip", low_memory=False)
            if not x.empty:
                print(f"  {platform_id} annotation loaded from cache ({len(x):,} rows)")
                return x
        except Exception:
            pass

    # GPL570 has an official GEO annotation file.  It contains metadata
    # preamble lines before the actual 21-column table, so locate the ID header
    # rather than using comment='#' (some metadata lines are not comments).
    if platform_id == "GPL570":
        # GEO bucket naming is historically inconsistent for older GPL IDs.
        # Try the current GPL5nnn bucket first, then the legacy GPLnnn path.
        urls = [
            "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL5nnn/GPL570/annot/GPL570.annot.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
        ]
        for url in urls:
            content = safe_get(url)
            if content is None:
                continue
            try:
                raw = gzip.decompress(content).decode("utf-8", errors="replace")
            except Exception:
                continue
            lines = raw.splitlines()
            header_i = None
            for i, line in enumerate(lines):
                fields = line.split("\t")
                if fields and fields[0].strip().strip('"') == "ID" and len(fields) > 5:
                    header_i = i
                    break
            if header_i is not None:
                x = pd.read_csv(
                    io.StringIO("\n".join(lines[header_i:])),
                    sep="\t",
                    dtype=str,
                    low_memory=False,
                )
                x.columns = [str(c).strip().strip('"') for c in x.columns]
                x.to_csv(cache_file, sep="\t", index=False, compression="gzip")
                print(
                    f"  GPL570 annotation downloaded from NCBI "
                    f"({len(x):,} rows)"
                )
                return x

    # GPL16686 does not expose the same *.annot.gz endpoint consistently.
    # Fetch the GPL family SOFT archive and extract the GPL16686 platform table.
    if platform_id == "GPL16686":
        candidates = [
            "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL16nnn/GPL16686/soft/GPL16686_family.soft.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL16nnn/GPL16686/soft/GPL16686.soft.gz",
        ]
        for url in candidates:
            content = safe_get(url)
            if content is None:
                continue
            try:
                raw = gzip.decompress(content).decode("utf-8", errors="replace")
            except Exception:
                continue
            lines = raw.splitlines()
            begin = end_i = None
            for i, line in enumerate(lines):
                if line.startswith("!platform_table_begin"):
                    begin = i + 1
                elif begin is not None and line.startswith("!platform_table_end"):
                    end_i = i
                    break
            if begin is not None and end_i is not None:
                x = pd.read_csv(
                    io.StringIO("\n".join(lines[begin:end_i])),
                    sep="\t", dtype=str, low_memory=False
                )
                x.columns = [str(c).strip().strip('"') for c in x.columns]
                x.to_csv(cache_file, sep="\t", index=False, compression="gzip")
                print(f"  GPL16686 annotation downloaded from NCBI SOFT ({len(x):,} rows)")
                return x

    raise RuntimeError(
        f"{accession}: could not obtain a usable local annotation table for {platform_id}"
    )


PLATFORM_IDS = {
    "GSE89216": "GPL16686",
    "GSE15043": "GPL570",
}


def choose_symbol_column(annotation: pd.DataFrame) -> str | None:
    """
    Find the most likely gene-symbol field in a GEO platform table.
    """
    exact = [
        "Gene Symbol",
        "Gene symbol",
        "GENE_SYMBOL",
        "Symbol",
        "SYMBOL",
        "gene_assignment",
    ]
    for c in exact:
        if c in annotation.columns:
            return c

    for c in annotation.columns:
        cl = str(c).strip().lower()
        if "gene symbol" in cl or cl == "symbol":
            return c

    for c in annotation.columns:
        cl = str(c).strip().lower()
        if "gene_assignment" in cl:
            return c

    return None


def clean_gene_symbol(value: Any) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip().strip('"')
    if not text:
        return ""

    # Standard GEO "Gene Symbol" fields.
    if "///" in text:
        text = text.split("///")[0].strip()
    if " // " in text:
        parts = [p.strip() for p in text.split(" // ")]
        # gene_assignment often begins with an accession/Entrez ID and the
        # symbol appears in the second field.
        if len(parts) >= 2 and re.fullmatch(r"[A-Za-z0-9._-]+", parts[1] or ""):
            text = parts[1]
        else:
            text = parts[0]

    return text.strip()


def microarray_deg(
    accession: str,
    matrix: pd.DataFrame,
    resistant: list[str],
    sensitive: list[str],
    annotation: pd.DataFrame,
) -> pd.DataFrame:
    missing = [x for x in resistant + sensitive if x not in matrix.columns]
    if missing:
        raise ValueError(f"{accession}: sample columns missing from series matrix: {missing}")

    ann_id = "ID" if "ID" in annotation.columns else annotation.columns[0]
    symbol_col = choose_symbol_column(annotation)
    if symbol_col is None:
        raise ValueError(
            f"{accession}: could not find a gene-symbol column in {list(annotation.columns)}"
        )

    ann = annotation[[ann_id, symbol_col]].copy()
    ann[ann_id] = ann[ann_id].astype(str)
    ann["gene_symbol"] = ann[symbol_col].map(clean_gene_symbol)
    ann = ann[ann["gene_symbol"].ne("")].drop_duplicates(subset=[ann_id])

    x = matrix.reset_index()
    x = x.rename(columns={x.columns[0]: ann_id})
    x[ann_id] = x[ann_id].astype(str)
    x = x.merge(ann[[ann_id, "gene_symbol"]], on=ann_id, how="inner")

    # GEO processed Affymetrix values are already normalized/log-scale.
    #
    # GSE15043 contains two RNA preparations for each of four resistant clones.
    # Treating all 8 resistant arrays as fully independent would pseudo-replicate
    # the clone effect. Average the replicate pair within each resistant clone,
    # then compare the four clone-level values against the two parental arrays.
    if accession == "GSE15043":
        clone_pairs = [
            ["GSM375721", "GSM375722"],  # BT/HerR1.0C
            ["GSM375723", "GSM375724"],  # BT/HerR1.0E
            ["GSM375725", "GSM375726"],  # BT/HerR0.2D
            ["GSM375727", "GSM375728"],  # BT/HerR0.2J
        ]
        clone_cols = []
        for idx, pair in enumerate(clone_pairs, start=1):
            col = f"resistant_clone_{idx}_mean"
            x[col] = x[pair].mean(axis=1)
            clone_cols.append(col)
        resistant_for_test = clone_cols
        x["resistant_mean"] = x[clone_cols].mean(axis=1)
    else:
        resistant_for_test = resistant
        x["resistant_mean"] = x[resistant].mean(axis=1)

    x["sensitive_mean"] = x[sensitive].mean(axis=1)
    x["external_log2FC"] = x["resistant_mean"] - x["sensitive_mean"]

    pvals = []
    for _, r in x.iterrows():
        rv = pd.to_numeric(
            r[resistant_for_test], errors="coerce"
        ).dropna().to_numpy(float)
        sv = pd.to_numeric(
            r[sensitive], errors="coerce"
        ).dropna().to_numpy(float)
        if len(rv) >= 2 and len(sv) >= 2:
            try:
                p = stats.ttest_ind(rv, sv, equal_var=False, nan_policy="omit").pvalue
            except Exception:
                p = np.nan
        else:
            p = np.nan
        pvals.append(p)

    x["pvalue"] = pvals

    # Benjamini-Hochberg FDR.
    p = x["pvalue"].to_numpy(float)
    valid = np.isfinite(p)
    padj = np.full(len(x), np.nan)
    if valid.any():
        pv = p[valid]
        order = np.argsort(pv)
        ranked = pv[order]
        n = len(ranked)
        q = ranked * n / np.arange(1, n + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q = np.clip(q, 0, 1)
        restored = np.empty(n)
        restored[order] = q
        padj[np.where(valid)[0]] = restored
    x["padj"] = padj

    # Collapse multiple probes to gene symbol by strongest absolute effect.
    x["abs_fc"] = x["external_log2FC"].abs()
    x = (
        x.sort_values(["gene_symbol", "abs_fc"], ascending=[True, False])
        .drop_duplicates("gene_symbol")
        .drop(columns=["abs_fc"])
    )

    x.insert(0, "dataset", accession)
    return x[
        [
            "dataset",
            "gene_symbol",
            "external_log2FC",
            "pvalue",
            "padj",
            "resistant_mean",
            "sensitive_mean",
        ]
    ]


def analyze_gse89216_or_15043(row: pd.Series) -> pd.DataFrame:
    acc = row["accession"]
    path = RAW_DIR / acc / f"{acc}_series_matrix.txt.gz"
    if not path.exists():
        raise FileNotFoundError(path)

    matrix, _ = read_geo_series_matrix(path)
    annotation = read_geo_platform_annotation(acc, PLATFORM_IDS[acc])
    if annotation is None:
        raise RuntimeError(f"{acc}: failed to retrieve platform annotation")

    resistant = str(row["resistant_samples"]).split(";")
    sensitive = str(row["sensitive_samples"]).split(";")

    return microarray_deg(
        accession=acc,
        matrix=matrix,
        resistant=resistant,
        sensitive=sensitive,
        annotation=annotation,
    )


def analyze_gse55005(row: pd.Series) -> pd.DataFrame:
    """
    Parse GSE55005_nexprs_btr50_bt474.xls.gz.

    Verified workbook schema:
        first column = gene symbol
        BTR50        = trastuzumab-resistant expression
        BT474        = parental/sensitive expression

    Values are positive normalized expression values. We compute a log2
    expression ratio with a small pseudocount:

        external_log2FC = log2((BTR50 + 1) / (BT474 + 1))

    There is only one untreated resistant and one untreated parental profile,
    so no replicate-based p-value/FDR is fabricated. GSE55005 contributes
    directional/effect-size validation only.
    """
    acc = row["accession"]
    path = RAW_DIR / acc / "GSE55005_nexprs_btr50_bt474.xls.gz"

    if not path.exists():
        raise FileNotFoundError(path)

    if xlrd is None:
        raise RuntimeError(
            "GSE55005 requires xlrd. Install with: "
            "python -m pip install xlrd"
        )

    with gzip.open(path, "rb") as handle:
        workbook_bytes = handle.read()

    book = xlrd.open_workbook(file_contents=workbook_bytes)
    sheet = book.sheet_by_index(0)

    if sheet.nrows < 2 or sheet.ncols < 3:
        raise ValueError(
            f"{acc}: unexpected workbook dimensions "
            f"{sheet.nrows} x {sheet.ncols}"
        )

    headers = [
        str(sheet.cell_value(0, c)).strip()
        for c in range(sheet.ncols)
    ]

    # The first header is blank in the source workbook.
    headers[0] = "gene_symbol"

    rows = [
        [sheet.cell_value(r, c) for c in range(sheet.ncols)]
        for r in range(1, sheet.nrows)
    ]
    df = pd.DataFrame(rows, columns=headers)

    print(f"  GSE55005 workbook dimensions: {len(df):,} genes x {len(df.columns)} columns")
    print(f"  GSE55005 workbook columns: {list(df.columns)}")

    required = {"gene_symbol", "BTR50", "BT474"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{acc}: missing expected workbook columns: {sorted(missing)}"
        )

    gene = df["gene_symbol"].map(clean_gene_symbol)
    resistant = pd.to_numeric(df["BTR50"], errors="coerce")
    sensitive = pd.to_numeric(df["BT474"], errors="coerce")

    # Positive normalized expression values -> log2 ratio.
    pseudocount = 1.0
    log2fc = np.log2(
        (resistant + pseudocount)
        / (sensitive + pseudocount)
    )

    out = pd.DataFrame(
        {
            "dataset": acc,
            "gene_symbol": gene,
            "external_log2FC": log2fc,
            "pvalue": np.nan,
            "padj": np.nan,
            "resistant_mean": resistant,
            "sensitive_mean": sensitive,
        }
    )

    out = out[
        out["gene_symbol"].ne("")
        & out["external_log2FC"].notna()
        & np.isfinite(out["external_log2FC"])
    ].copy()

    out["abs_fc"] = out["external_log2FC"].abs()
    out = (
        out.sort_values(
            ["gene_symbol", "abs_fc"],
            ascending=[True, False],
        )
        .drop_duplicates("gene_symbol")
        .drop(columns=["abs_fc"])
    )

    return out


def load_discovery() -> pd.DataFrame:
    if not DISCOVERY_DEG_FILE.exists():
        raise FileNotFoundError(
            f"Discovery DEG file not found: {DISCOVERY_DEG_FILE}"
        )

    df = pd.read_csv(DISCOVERY_DEG_FILE)

    symbol_col = None
    for c in ["Gene name", "gene_symbol", "symbol", "Gene Symbol"]:
        if c in df.columns:
            symbol_col = c
            break

    if symbol_col is None:
        raise ValueError(
            "Could not identify a gene-symbol column in the discovery DEG file."
        )

    required = ["log2FoldChange", "padj"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Discovery DEG file missing required column: {c}")

    out = pd.DataFrame(
        {
            "gene_symbol": df[symbol_col].map(clean_gene_symbol),
            "atlas_log2FC": pd.to_numeric(df["log2FoldChange"], errors="coerce"),
            "atlas_padj": pd.to_numeric(df["padj"], errors="coerce"),
        }
    )
    out = out[out["gene_symbol"].ne("") & out["atlas_log2FC"].notna()]
    out["atlas_abs_fc"] = out["atlas_log2FC"].abs()
    out = (
        out.sort_values(["gene_symbol", "atlas_abs_fc"], ascending=[True, False])
        .drop_duplicates("gene_symbol")
        .drop(columns=["atlas_abs_fc"])
    )
    return out


def validate_against_discovery(
    discovery: pd.DataFrame,
    external: pd.DataFrame,
    accession: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    merged = discovery.merge(external, on="gene_symbol", how="inner")

    merged["atlas_direction"] = np.where(
        merged["atlas_log2FC"] > 0,
        "UP",
        np.where(merged["atlas_log2FC"] < 0, "DOWN", "FLAT"),
    )
    merged["external_direction"] = np.where(
        merged["external_log2FC"] > 0,
        "UP",
        np.where(merged["external_log2FC"] < 0, "DOWN", "FLAT"),
    )
    merged["direction_agreement"] = (
        merged["atlas_direction"] == merged["external_direction"]
    )

    atlas_sig = (
        merged["atlas_padj"].lt(0.05)
        & merged["atlas_log2FC"].abs().ge(1)
    )

    # External statistical significance is only defined when the external
    # dataset actually has replicate-derived adjusted p-values. GSE55005 has
    # one untreated profile per state and therefore contributes effect/direction
    # evidence only.
    if merged["padj"].notna().any():
        ext_sig = (
            merged["external_log2FC"].abs().ge(1)
            & merged["padj"].lt(0.05)
        )
        both_sig = atlas_sig & ext_sig
    else:
        ext_sig = pd.Series(False, index=merged.index)
        both_sig = pd.Series(False, index=merged.index)

    corr_all = np.nan
    if len(merged) >= 3:
        corr_all = stats.spearmanr(
            merged["atlas_log2FC"],
            merged["external_log2FC"],
            nan_policy="omit",
        ).statistic

    corr_sig = np.nan
    sig_subset = merged[both_sig]
    if len(sig_subset) >= 3:
        corr_sig = stats.spearmanr(
            sig_subset["atlas_log2FC"],
            sig_subset["external_log2FC"],
            nan_policy="omit",
        ).statistic

    # Whole-transcriptome sign agreement can be dominated by genes with
    # near-zero external effects. Report an additional informative-effect metric.
    informative_effect = merged["external_log2FC"].abs().ge(0.5)
    atlas_sig_direction = merged.loc[atlas_sig, "direction_agreement"]

    metrics = {
        "dataset": accession,
        "overlapping_genes": int(len(merged)),
        "direction_agreement_n": int(merged["direction_agreement"].sum()),
        "direction_agreement_fraction": float(
            merged["direction_agreement"].mean()
        ) if len(merged) else np.nan,
        "informative_external_effect_n": int(informative_effect.sum()),
        "informative_effect_direction_agreement_fraction": float(
            merged.loc[informative_effect, "direction_agreement"].mean()
        ) if informative_effect.any() else np.nan,
        "atlas_significant_overlap_n": int(atlas_sig.sum()),
        "atlas_significant_direction_agreement_fraction": float(
            atlas_sig_direction.mean()
        ) if len(atlas_sig_direction) else np.nan,
        "atlas_sig_external_strong_effect_n": int(
            (atlas_sig & merged["external_log2FC"].abs().ge(1)).sum()
        ),
        "atlas_sig_external_strong_effect_agreement_fraction": float(
            merged.loc[
                atlas_sig & merged["external_log2FC"].abs().ge(1),
                "direction_agreement"
            ].mean()
        ) if (atlas_sig & merged["external_log2FC"].abs().ge(1)).any() else np.nan,
        "both_significant_n": int(both_sig.sum()),
        "both_significant_direction_agreement_n": int(
            merged.loc[both_sig, "direction_agreement"].sum()
        ),
        "both_significant_direction_agreement_fraction": float(
            merged.loc[both_sig, "direction_agreement"].mean()
        ) if both_sig.any() else np.nan,
        "spearman_all_overlap": corr_all,
        "spearman_both_significant": corr_sig,
        "validation_evidence_type": (
            "DIRECTIONAL_EFFECT_ONLY"
            if accession == "GSE55005"
            else (
                "CLONE_AWARE_REPLICATE_STATISTICAL_AND_DIRECTIONAL"
                if accession == "GSE15043"
                else "EXPLORATORY_REPLICATE_STATISTICAL_AND_DIRECTIONAL"
            )
        ),
    }

    merged.insert(0, "validation_dataset", accession)
    return merged, metrics


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args()

    print("=" * 72)
    print("ATLAS — Stage 03B External Dataset Validation")
    print("=" * 72)

    registry = write_registry()
    print(f"Dataset registry written: {REGISTRY_FILE}")

    download_status = []

    if not args.analyze_only:
        print("\nDownloading selected GEO processed data...")
        for _, row in registry[registry["include"] == True].iterrows():  # noqa: E712
            print(f"  {row['accession']} ...")
            status = download_dataset(row)
            download_status.append(status)
            print(
                f"    series_matrix={status['series_matrix_downloaded']} "
                f"processed_comparison={status['processed_comparison_downloaded']} "
                f"errors={status['download_error'] or 'none'}"
            )

        pd.DataFrame(download_status).to_csv(
            OUT_ROOT / "download_status.csv",
            index=False,
        )

    if args.download_only:
        print("\n03B download-only mode complete.")
        return 0

    discovery = load_discovery()

    all_validation_rows = []
    metrics_rows = []
    analysis_errors = []

    print("\nAnalyzing each external dataset independently...")

    for _, row in registry[registry["include"] == True].iterrows():  # noqa: E712
        acc = row["accession"]
        print(f"\n[{acc}]")

        try:
            if acc in {"GSE89216", "GSE15043"}:
                external = analyze_gse89216_or_15043(row)
            elif acc == "GSE55005":
                external = analyze_gse55005(row)
            else:
                raise ValueError(f"No analysis adapter implemented for {acc}")

            ds_dir = PROCESSED_DIR / acc
            ds_dir.mkdir(parents=True, exist_ok=True)

            external.to_csv(
                ds_dir / "external_deg_results.csv",
                index=False,
            )

            validated, metrics = validate_against_discovery(
                discovery,
                external,
                acc,
            )

            validated.to_csv(
                ds_dir / "validation_gene_comparison.csv",
                index=False,
            )

            pd.DataFrame([metrics]).to_csv(
                ds_dir / "validation_metrics.csv",
                index=False,
            )

            all_validation_rows.append(validated)
            metrics_rows.append(metrics)

            print(f"  External genes: {len(external):,}")
            print(f"  Overlap with ATLAS: {metrics['overlapping_genes']:,}")
            print(
                "  Direction agreement (all overlap): "
                f"{metrics['direction_agreement_fraction']:.3f}"
            )
            if not math.isnan(
                metrics["informative_effect_direction_agreement_fraction"]
            ):
                print(
                    "  Direction agreement (|external effect| >= 0.5): "
                    f"{metrics['informative_effect_direction_agreement_fraction']:.3f}"
                )
            if not math.isnan(
                metrics["atlas_significant_direction_agreement_fraction"]
            ):
                print(
                    "  ATLAS-significant direction agreement: "
                    f"{metrics['atlas_significant_direction_agreement_fraction']:.3f}"
                )
            print(
                "  Both-significant agreement: "
                f"{metrics['both_significant_direction_agreement_fraction']:.3f}"
                if not math.isnan(metrics["both_significant_direction_agreement_fraction"])
                else "  Both-significant agreement: n/a"
            )

        except Exception as exc:
            analysis_errors.append(
                {"dataset": acc, "error": repr(exc)}
            )
            print(f"  ERROR: {exc}")

    if all_validation_rows:
        all_validation = pd.concat(
            all_validation_rows,
            ignore_index=True,
        )
        all_validation.to_csv(CROSS_GENE_FILE, index=False)

    summary = pd.DataFrame(metrics_rows)
    summary.to_csv(SUMMARY_FILE, index=False)

    if analysis_errors:
        pd.DataFrame(analysis_errors).to_csv(
            OUT_ROOT / "analysis_errors.csv",
            index=False,
        )

    metadata = {
        "stage": "03B",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "discovery_deg_file": str(DISCOVERY_DEG_FILE),
        "datasets_requested": [
            x["accession"] for x in DATASETS if x["include"]
        ],
        "datasets_successfully_analyzed": [
            x["dataset"] for x in metrics_rows
        ],
        "analysis_errors": analysis_errors,
        "methodological_note": (
            "External datasets are analyzed independently. Raw expression "
            "values are not pooled across platforms. Validation is based on "
            "gene-level direction/effect concordance and overlap."
        ),
        "next_stage": "03C consensus resistance signature",
    }
    METADATA_FILE.write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    print("\n" + "=" * 72)
    print("STAGE 03B COMPLETE")
    print("=" * 72)
    print(f"Summary: {SUMMARY_FILE}")
    print(f"Cross-gene validation: {CROSS_GENE_FILE}")
    print(f"Metadata: {METADATA_FILE}")

    if analysis_errors:
        print(
            "\nSome datasets require schema/manual adapter review. "
            "See analysis_errors.csv."
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
