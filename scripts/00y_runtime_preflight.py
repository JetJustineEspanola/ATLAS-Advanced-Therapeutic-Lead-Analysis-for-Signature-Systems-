#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib
import json
import os
import re
import shutil
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
OUTDIR = ROOT / "results" / "pipeline_state"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_JSON = OUTDIR / "runtime_preflight.json"
OUT_CSV = OUTDIR / "runtime_preflight.csv"

REQUIRED_IMPORTS = [
    "pandas", "numpy", "scipy", "requests", "duckdb",
    "pydeseq2", "gseapy", "rdkit", "vina", "openbabel",
]

HOSTS = [
    ("NCBI", "eutils.ncbi.nlm.nih.gov", 443),
    ("STRING", "string-db.org", 443),
    ("UniProt", "rest.uniprot.org", 443),
    ("PubChem", "pubchem.ncbi.nlm.nih.gov", 443),
    ("ClinicalTrials", "clinicaltrials.gov", 443),
]

MIN_FREE_GB = 3.0


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def runner_scripts():
    runner = SCRIPTS / "run_atlas_full_auto.py"
    if not runner.exists():
        return []
    text = runner.read_text(encoding="utf-8", errors="replace")
    return re.findall(
        r'Stage\(\s*"[^"]+"\s*,\s*"([^"]+\.py)"',
        text,
        flags=re.MULTILINE,
    )


def check_import(name):
    try:
        importlib.import_module(name)
        return True, ""
    except Exception as exc:
        return False, f"{type(exc).__name__}: {exc}"


def tcp_check(host, port=443, timeout=8):
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True, ""
    except Exception as exc:
        return False, f"{type(exc).__name__}: {exc}"


def discover_cmap_env_vars():
    path = SCRIPTS / "04g_cmap_submit_all.py"
    if not path.exists():
        return [], False

    text = path.read_text(encoding="utf-8", errors="replace")

    patterns = [
        r'os\.getenv\(\s*["\']([^"\']+)["\']',
        r'os\.environ\.get\(\s*["\']([^"\']+)["\']',
        r'os\.environ\[\s*["\']([^"\']+)["\']\s*\]',
    ]

    names = set()
    for pat in patterns:
        names.update(re.findall(pat, text))

    hardcoded = bool(
        re.search(
            r'(?i)(api[_-]?key|token|authorization)\s*=\s*["\'][^"\']{12,}["\']',
            text,
        )
    )

    return sorted(names), hardcoded


def main():
    rows = []
    hard_fail = False
    warnings = []

    py_ok = sys.version_info >= (3, 10)
    rows.append({
        "check": "python_version",
        "target": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "status": "PASS" if py_ok else "FAIL",
        "detail": "requires Python >= 3.10",
    })
    hard_fail |= not py_ok

    for mod in REQUIRED_IMPORTS:
        ok, detail = check_import(mod)
        rows.append({
            "check": "python_import",
            "target": mod,
            "status": "PASS" if ok else "FAIL",
            "detail": detail,
        })
        hard_fail |= not ok

    refs = runner_scripts()
    if not refs:
        rows.append({
            "check": "runner_stage_discovery",
            "target": "run_atlas_full_auto.py",
            "status": "FAIL",
            "detail": "no runner stages detected",
        })
        hard_fail = True
    else:
        for name in refs:
            exists = (SCRIPTS / name).exists()
            rows.append({
                "check": "runner_script",
                "target": name,
                "status": "PASS" if exists else "FAIL",
                "detail": "",
            })
            hard_fail |= not exists

    for directory in [ROOT / "results", ROOT / "logs", ROOT / "data"]:
        try:
            directory.mkdir(parents=True, exist_ok=True)
            probe = directory / ".atlas_write_probe"
            probe.write_text("ok", encoding="utf-8")
            probe.unlink()
            ok, detail = True, ""
        except Exception as exc:
            ok, detail = False, f"{type(exc).__name__}: {exc}"

        rows.append({
            "check": "writable_directory",
            "target": str(directory),
            "status": "PASS" if ok else "FAIL",
            "detail": detail,
        })
        hard_fail |= not ok

    usage = shutil.disk_usage(ROOT)
    free_gb = usage.free / (1024 ** 3)
    disk_ok = free_gb >= MIN_FREE_GB
    rows.append({
        "check": "disk_space",
        "target": str(ROOT),
        "status": "PASS" if disk_ok else "FAIL",
        "detail": f"{free_gb:.2f} GiB free; minimum {MIN_FREE_GB:.1f} GiB",
    })
    hard_fail |= not disk_ok

    for label, host, port in HOSTS:
        ok, detail = tcp_check(host, port)
        rows.append({
            "check": "network_connectivity",
            "target": label,
            "status": "PASS" if ok else "FAIL",
            "detail": host if ok else f"{host}: {detail}",
        })
        hard_fail |= not ok

    env_vars, hardcoded = discover_cmap_env_vars()
    if env_vars:
        for name in env_vars:
            set_ok = bool(os.environ.get(name))
            rows.append({
                "check": "cmap_environment",
                "target": name,
                "status": "PASS" if set_ok else "FAIL",
                "detail": "set" if set_ok else "environment variable not set",
            })
            hard_fail |= not set_ok
    else:
        msg = (
            "No CMap credential environment variable was automatically detected "
            "in 04g_cmap_submit_all.py; verify its authentication method."
        )
        warnings.append(msg)
        rows.append({
            "check": "cmap_environment",
            "target": "auto-detect",
            "status": "WARN",
            "detail": msg,
        })

    if hardcoded:
        msg = (
            "Possible hard-coded credential detected in 04g_cmap_submit_all.py. "
            "Value was not read or printed; migrate credentials to environment variables."
        )
        warnings.append(msg)
        rows.append({
            "check": "credential_hygiene",
            "target": "04g_cmap_submit_all.py",
            "status": "WARN",
            "detail": msg,
        })

    fields = ["check", "target", "status", "detail"]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    report = {
        "created_utc": utcnow(),
        "project_root": str(ROOT),
        "overall_status": "FAIL" if hard_fail else "PASS",
        "warnings": warnings,
        "checks": rows,
    }
    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("=" * 88)
    print("ATLAS — RUNTIME / NETWORK / DEPENDENCY PREFLIGHT")
    print("=" * 88)

    for r in rows:
        print(f"{r['status']:4s}  {r['check']:24s}  {r['target']}")
        if r["detail"] and r["status"] != "PASS":
            print(f"      {r['detail']}")

    print()
    print(f"OVERALL: {'FAIL' if hard_fail else 'PASS'}")
    print(f"JSON: {OUT_JSON}")
    print(f"CSV:  {OUT_CSV}")

    return 2 if hard_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
