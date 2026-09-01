#!/usr/bin/env python3
"""
ATLAS Pipeline Orchestrator

Runs ATLAS stages in order, skips completed stages by default, records logs,
and can resume from a chosen stage.

Examples
--------
Run downstream pipeline:
    python -u scripts/run_atlas_pipeline.py --from 04N --to 04U

Force rerun from safety onward:
    python -u scripts/run_atlas_pipeline.py --from 04O --to 04U --force

Run only 04R-04U:
    python -u scripts/run_atlas_pipeline.py --from 04R --to 04U

Launch Streamlit after completion:
    python -u scripts/run_atlas_pipeline.py --from 04N --to 04U --ui
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT_ROOT / "scripts"
LOG_DIR = PROJECT_ROOT / "results" / "pipeline_logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

STAGES = [
    {
        "id": "04N",
        "script": "04n_regulatory_status.py",
        "args": ["--mode", "full", "--workers", "8"],
        "output": PROJECT_ROOT / "results/cmap/regulatory_status/ATLAS_CMap_regulatory_annotations.csv",
    },
    {
        "id": "04O",
        "script": "04o_safety_screening.py",
        "args": ["--mode", "priority", "--tier2-top", "100"],
        "output": PROJECT_ROOT / "results/cmap/safety_screening/ATLAS_CMap_safety_screening.csv",
    },
    {
        "id": "04P",
        "script": "04p_drug_target_annotation.py",
        "args": ["--max-candidates", "50", "--workers", "6"],
        "output": PROJECT_ROOT / "results/cmap/drug_targets/ATLAS_CMap_drug_target_annotations.csv",
    },
    {
        "id": "04Q",
        "script": "04q_network_integration.py",
        "args": [
            "--max-drugs", "25",
            "--max-targets-per-drug", "10",
            "--max-resistance-genes", "200",
            "--string-score", "700",
            "--workers", "6",
        ],
        "output": PROJECT_ROOT / "results/cmap/network_integration/ATLAS_drug_network_prioritized.csv",
    },
    {
        "id": "04R",
        "script": "04r_final_candidate_prioritization.py",
        "args": ["--top-docking", "5"],
        "output": PROJECT_ROOT / "results/cmap/final_prioritization/ATLAS_docking_shortlist.csv",
    },
    {
        "id": "04S",
        "script": "04s_target_supported_docking.py",
        "args": ["--exhaustiveness", "16", "--num-modes", "9", "--max-candidates", "5"],
        "output": PROJECT_ROOT / "results/cmap/docking/ATLAS_docking_results.csv",
    },
    {
        "id": "04T",
        "script": "04t_admet_structural_assessment.py",
        "args": ["--max-candidates", "25"],
        "output": PROJECT_ROOT / "results/cmap/admet_structural/ATLAS_ADMET_structural_assessment.csv",
    },
    {
        "id": "04U",
        "script": "04u_integrated_evidence_matrix.py",
        "args": ["--max-candidates", "25", "--experimental-top", "5"],
        "output": PROJECT_ROOT / "results/cmap/integrated_evidence/ATLAS_integrated_evidence_matrix.csv",
    },
]

STAGE_IDS = [s["id"] for s in STAGES]


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def print_header(text: str) -> None:
    print("\n" + "=" * 78, flush=True)
    print(text, flush=True)
    print("=" * 78, flush=True)


def run_stage(stage: dict, force: bool) -> dict:
    sid = stage["id"]
    script = SCRIPTS / stage["script"]
    output = Path(stage["output"])

    if not script.exists():
        return {
            "stage": sid,
            "status": "MISSING_SCRIPT",
            "script": str(script),
            "output": str(output),
            "started_utc": now(),
            "ended_utc": now(),
        }

    if output.exists() and not force:
        print(f"[{sid}] SKIP — output already exists", flush=True)
        return {
            "stage": sid,
            "status": "SKIPPED_EXISTING_OUTPUT",
            "script": str(script),
            "output": str(output),
            "started_utc": now(),
            "ended_utc": now(),
        }

    cmd = [sys.executable, "-u", str(script), *stage["args"]]
    log_path = LOG_DIR / f"{sid}.log"

    print(f"[{sid}] RUN", flush=True)
    print(" ".join(cmd), flush=True)

    started = now()

    with log_path.open("w", encoding="utf-8") as log:
        proc = subprocess.Popen(
            cmd,
            cwd=str(PROJECT_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        assert proc.stdout is not None

        for line in proc.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()

        code = proc.wait()

    ended = now()

    status = "COMPLETED" if code == 0 and output.exists() else "FAILED"

    return {
        "stage": sid,
        "status": status,
        "returncode": code,
        "script": str(script),
        "output": str(output),
        "log": str(log_path),
        "started_utc": started,
        "ended_utc": ended,
    }


def launch_ui() -> None:
    app = PROJECT_ROOT / "app.py"

    if not app.exists():
        print("UI not launched: app.py not found.", flush=True)
        return

    print_header("Launching Streamlit UI")
    os.execvp(
        sys.executable,
        [
            sys.executable,
            "-m",
            "streamlit",
            "run",
            str(app),
        ],
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument(
        "--from",
        dest="from_stage",
        choices=STAGE_IDS,
        default="04N",
    )
    p.add_argument(
        "--to",
        dest="to_stage",
        choices=STAGE_IDS,
        default="04U",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Rerun selected stages even when outputs already exist.",
    )
    p.add_argument(
        "--ui",
        action="store_true",
        help="Launch Streamlit after successful completion.",
    )

    return p.parse_args()


def main() -> int:
    args = parse_args()

    start_idx = STAGE_IDS.index(args.from_stage)
    end_idx = STAGE_IDS.index(args.to_stage)

    if start_idx > end_idx:
        print("ERROR: --from stage must come before --to stage.", flush=True)
        return 2

    selected = STAGES[start_idx:end_idx + 1]

    print_header("ATLAS Pipeline Automation")
    print(f"Project: {PROJECT_ROOT}", flush=True)
    print(f"Stages: {args.from_stage} -> {args.to_stage}", flush=True)
    print(f"Force rerun: {args.force}", flush=True)

    manifest = {
        "started_utc": now(),
        "project_root": str(PROJECT_ROOT),
        "from_stage": args.from_stage,
        "to_stage": args.to_stage,
        "force": args.force,
        "stages": [],
    }

    for stage in selected:
        result = run_stage(stage, force=args.force)
        manifest["stages"].append(result)

        if result["status"] in {"FAILED", "MISSING_SCRIPT"}:
            print_header(f"PIPELINE STOPPED AT {stage['id']}")
            print(json.dumps(result, indent=2), flush=True)

            manifest["ended_utc"] = now()
            manifest["status"] = "FAILED"

            manifest_path = LOG_DIR / "latest_manifest.json"
            manifest_path.write_text(
                json.dumps(manifest, indent=2),
                encoding="utf-8",
            )
            return 1

    manifest["ended_utc"] = now()
    manifest["status"] = "COMPLETED"

    manifest_path = LOG_DIR / "latest_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print_header("ATLAS PIPELINE COMPLETE")

    for result in manifest["stages"]:
        print(
            f"{result['stage']:>4}  {result['status']}",
            flush=True,
        )

    print(f"\nManifest: {manifest_path}", flush=True)
    print(f"Logs: {LOG_DIR}", flush=True)

    if args.ui:
        launch_ui()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
