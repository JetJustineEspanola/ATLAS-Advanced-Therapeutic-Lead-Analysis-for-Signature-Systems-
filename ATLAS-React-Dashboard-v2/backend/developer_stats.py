from __future__ import annotations

from pathlib import Path
import csv
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone


def _service_state(name: str):
    try:
        p = subprocess.run(
            ["systemctl", "--user", "show", name, "-p", "ActiveState", "-p", "SubState", "-p", "MainPID"],
            capture_output=True,
            text=True,
            timeout=2,
        )
        values = {}
        for line in p.stdout.splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                values[k] = v
        return {
            "available": p.returncode == 0,
            "active_state": values.get("ActiveState", "unknown"),
            "sub_state": values.get("SubState", "unknown"),
            "main_pid": values.get("MainPID", "0"),
        }
    except Exception as exc:
        return {"available": False, "active_state": "unknown", "sub_state": str(exc), "main_pid": "0"}


def _queue_summary(root: Path):
    path = root / "results/pipeline_state/dataset_queue/dataset_queue.csv"
    summary = {"source_file": str(path), "exists": path.exists(), "rows": 0, "status_counts": {}, "recent": []}
    if not path.exists():
        return summary
    try:
        with path.open(newline="", encoding="utf-8", errors="replace") as f:
            rows = list(csv.DictReader(f))
        summary["rows"] = len(rows)
        status_key = None
        if rows:
            for k in rows[0].keys():
                if k.lower() in {"status", "queue_status", "state"}:
                    status_key = k
                    break
        counts = {}
        if status_key:
            for row in rows:
                s = row.get(status_key) or "UNKNOWN"
                counts[s] = counts.get(s, 0) + 1
        summary["status_counts"] = counts
        summary["recent"] = rows[-10:]
    except Exception as exc:
        summary["error"] = str(exc)
    return summary


def _output_registry(root: Path):
    entries = []
    watched = [
        "results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv",
        "data/enriched/dataset_candidates_independence_scored.csv",
        "results/cmap/integrated_evidence/ATLAS_integrated_evidence_matrix.csv",
        "results/pipeline_state/dataset_queue/dataset_queue.csv",
    ]
    for rel in watched:
        p = root / rel
        entries.append({
            "path": rel,
            "exists": p.exists(),
            "size_mb": round(p.stat().st_size / (1024**2), 3) if p.exists() else None,
            "modified": datetime.fromtimestamp(p.stat().st_mtime, tz=timezone.utc).astimezone().isoformat(timespec="seconds") if p.exists() else None,
        })
    return entries


def build_developer_statistics(root: Path):
    st = os.statvfs(root if root.exists() else Path.home())
    total = st.f_blocks * st.f_frsize
    free = st.f_bavail * st.f_frsize
    service = _service_state("atlas-dataset-queue.service")
    return {
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "pid": os.getpid(),
            "cwd": os.getcwd(),
            "atlas_root": str(root),
        },
        "filesystem": {
            "total_gb": round(total / (1024**3), 2),
            "free_gb": round(free / (1024**3), 2),
            "used_gb": round((total - free) / (1024**3), 2),
            "free_percent": round((free / total * 100), 1) if total else None,
        },
        "queue_service": service,
        "dataset_queue": _queue_summary(root),
        "outputs": _output_registry(root),
        "api": [
            "GET /api/health",
            "GET /api/dashboard",
            "GET /api/datasets",
            "GET /api/signature",
            "GET /api/candidates",
            "GET /api/cmap",
            "GET /api/docking",
            "GET /api/research/statistics",
            "GET /api/developer/statistics",
            "GET /api/settings",
            "PUT /api/settings",
        ],
    }
