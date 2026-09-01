#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def run(name, args):
    cmd = [sys.executable, "-u", str(ROOT/"scripts"/name), *args]
    print("\n" + "="*78)
    print(" ".join(cmd))
    print("="*78)
    return subprocess.call(cmd, cwd=ROOT)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sources", default="geo,sra,gdc,cbioportal")
    p.add_argument("--retmax", type=int, default=25)
    p.add_argument("--query-family", default="all")
    args = p.parse_args()

    if run("00a_dataset_discovery.py", [
        "--sources", args.sources,
        "--retmax", str(args.retmax),
        "--query-family", args.query_family,
    ]) != 0:
        return 1

    if run("00b_metadata_catalog.py", []) != 0:
        return 1

    if run("00c_dataset_eligibility.py", []) != 0:
        return 1

    print("\nATLAS 00A-00C acquisition phase complete.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
