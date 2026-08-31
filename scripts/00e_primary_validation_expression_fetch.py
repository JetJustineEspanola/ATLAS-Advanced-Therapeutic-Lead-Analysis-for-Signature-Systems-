#!/usr/bin/env python3
"""
ATLAS — 00E Primary Validation Expression Fetch

Downloads processed/raw-count expression files for the two current
PRIMARY_VALIDATION cohorts:
  - GSE121105
  - GSE237606

This stage intentionally prefers GEO processed count files over re-downloading
FASTQ/SRA because cross-study DE validation can be performed from deposited
gene-level counts while avoiding unnecessary large raw-read downloads.

Outputs:
  data/validation_expression/GSE121105/
  data/validation_expression/GSE237606/
  data/validation_expression/primary_expression_manifest.csv
"""

from __future__ import annotations

from pathlib import Path
import argparse
import hashlib
import html
import re
import shutil
import sys
import tarfile
from urllib.parse import urljoin

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/validation_expression"
OUT.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "ATLAS-research-pipeline/1.0"}


def geo_stub(acc: str) -> str:
    return re.sub(r"\d{1,3}$", "nnn", acc)


def suppl_url(acc: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_stub(acc)}/{acc}/suppl/"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def list_links(url: str) -> list[str]:
    r = requests.get(url, headers=UA, timeout=60)
    r.raise_for_status()
    links = re.findall(r'href=["\']([^"\']+)["\']', r.text, flags=re.I)
    names = []
    for link in links:
        name = html.unescape(link.split("?")[0]).rstrip("/").split("/")[-1]
        if name and name not in {".", ".."}:
            names.append(name)
    return sorted(set(names))


def download(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  cached: {dest.name}")
        return dest

    print(f"  downloading: {dest.name}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, headers=UA, timeout=120, stream=True) as r:
        r.raise_for_status()
        with tmp.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    tmp.replace(dest)
    return dest


def safe_extract_tar(path: Path, dest: Path):
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "r:*") as tf:
        base = dest.resolve()
        members = []
        for m in tf.getmembers():
            target = (dest / m.name).resolve()
            if base not in target.parents and target != base:
                raise RuntimeError(f"Unsafe tar path: {m.name}")
            members.append(m)
        tf.extractall(dest, members=members)


def select_files(acc: str, names: list[str]) -> list[str]:
    # Prefer count-oriented processed data and raw archives.
    wanted = []
    for n in names:
        low = n.lower()
        if any(k in low for k in [
            "count", "raw", "expr", "expression"
        ]) and any(low.endswith(ext) for ext in [
            ".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz",
            ".tar", ".tar.gz", ".tgz"
        ]):
            wanted.append(n)

    # Known strong choice for GSE121105.
    if acc == "GSE121105":
        preferred = [n for n in names if "genecount" in n.lower()]
        if preferred:
            return preferred

    return wanted


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--accessions",
        default="GSE121105,GSE237606",
        help="Comma-separated GEO Series accessions",
    )
    args = ap.parse_args()

    accessions = [x.strip().upper() for x in args.accessions.split(",") if x.strip()]
    manifest = []

    print("=" * 78)
    print("ATLAS — 00E PRIMARY VALIDATION EXPRESSION FETCH")
    print("=" * 78)

    for acc in accessions:
        print(f"\n[{acc}]")
        d = OUT / acc
        d.mkdir(parents=True, exist_ok=True)

        url = suppl_url(acc)
        try:
            names = list_links(url)
        except Exception as e:
            print(f"  ERROR listing GEO supplementary directory: {e}")
            continue

        selected = select_files(acc, names)

        print(f"  supplementary files listed: {len(names)}")
        print(f"  selected expression/count files: {len(selected)}")

        if not selected:
            print("  WARNING: no count-oriented supplementary files detected.")
            print(f"  Inspect manually: {url}")
            continue

        for name in selected:
            try:
                path = download(urljoin(url, name), d / name)
                record = {
                    "accession": acc,
                    "filename": name,
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "extracted": False,
                }

                low = name.lower()
                if low.endswith((".tar", ".tar.gz", ".tgz")):
                    extract_dir = d / "extracted"
                    print(f"  extracting: {name}")
                    safe_extract_tar(path, extract_dir)
                    record["extracted"] = True
                    record["extract_dir"] = str(extract_dir)

                manifest.append(record)

            except Exception as e:
                print(f"  ERROR {name}: {e}")

    mdf = pd.DataFrame(manifest)
    manifest_path = OUT / "primary_expression_manifest.csv"
    mdf.to_csv(manifest_path, index=False)

    print("\n" + "=" * 78)
    print("00E COMPLETE")
    print("=" * 78)
    if not mdf.empty:
        cols = ["accession", "filename", "size_bytes", "extracted"]
        print(mdf[cols].to_string(index=False))
    else:
        print("No files downloaded.")

    print(f"\nManifest: {manifest_path}")
    print(f"Data root: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
