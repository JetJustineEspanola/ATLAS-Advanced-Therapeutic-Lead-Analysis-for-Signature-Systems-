#!/usr/bin/env python3
"""
ATLAS — Full End-to-End Automation Runner

Runs the ATLAS pipeline from online dataset discovery through final integrated
evidence, with:
- ordered stages
- checkpoint/resume support
- per-stage logs
- skip-on-success behavior
- optional full refresh
- asynchronous CMap polling
- safe stop on failed stages
- ability to resume from any stage

Designed around the current scripts/ folder layout.

Usage
-----
python -u scripts/run_atlas_full_auto.py --full
python -u scripts/run_atlas_full_auto.py --full --refresh-data
python -u scripts/run_atlas_full_auto.py --from-stage 00g --to-stage 04u
python -u scripts/run_atlas_full_auto.py --list-stages
python -u scripts/run_atlas_full_auto.py --rerun-stage 04q

Notes
-----
- This orchestrator assumes each underlying script can run with its defaults.
- It does NOT fabricate success: any nonzero exit code stops the pipeline.
- CMap submission/check/download is treated specially because jobs may be async.
- Existing successful checkpoints are skipped unless --refresh-data,
  --force, or --rerun-stage is used.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"

STATE_DIR = ROOT / "results" / "pipeline_state"
LOG_DIR = ROOT / "logs" / "full_pipeline"
STATE_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

STATE_FILE = STATE_DIR / "full_pipeline_state.json"

PYTHON = sys.executable


@dataclass
class Stage:
    key: str
    script: str
    group: str
    description: str
    checkpoint: Optional[str] = None
    args: tuple[str, ...] = ()
    special: Optional[str] = None


STAGES = [
    # ------------------------------------------------------------------
    # Online discovery / metadata / validation study selection
    # ------------------------------------------------------------------
    Stage(
        "00y", "00y_runtime_preflight.py", "preflight",
        "Runtime, network, dependency, and credential preflight",
        "results/pipeline_state/runtime_preflight.json",
    ),

    Stage(
        "00a", "00a_dataset_discovery.py", "acquisition",
        "Discover candidate datasets online",
        "data/catalog/dataset_candidates.csv",
    ),
    Stage(
        "00a2", "00a2_ebi_external_discovery.py", "acquisition",
        "Discover BioStudies/ArrayExpress candidates and merge external leads",
        args=("--merge-existing",),
    ),
    Stage(
        "00b", "00b_metadata_catalog.py", "acquisition",
        "Build/update metadata catalog",
        "data/catalog/atlas_metadata.duckdb",
    ),
    Stage(
        "00c", "00c_dataset_eligibility.py", "acquisition",
        "Initial dataset eligibility scoring",
    ),
    Stage(
        "00d", "00d_metadata_enrichment.py", "acquisition",
        "Enrich candidate dataset metadata online",
    ),
    Stage(
        "00d_ebi", "00d_ebi_metadata_enrichment.py", "acquisition",
        "Enrich BioStudies/ArrayExpress sample metadata and SDRF records",
    ),
    Stage(
        "00d1", "00d1_phenotype_audit.py", "phenotype",
        "Audit resistance phenotype metadata",
    ),
    Stage(
        "00d2", "00d2_conservative_phenotype_classifier.py", "phenotype",
        "Conservative phenotype classification",
    ),
    Stage(
        "00d2b", "00d2b_curated_study_phenotype_mapping.py", "phenotype",
        "Apply curated phenotype mappings",
    ),
    Stage(
        "00d2c", "00d2c_confirm_her2_status.py", "phenotype",
        "Confirm HER2-positive context",
    ),
    Stage(
        "00d3", "00d3_dataset_relationship_audit.py", "phenotype",
        "Audit study overlap / umbrella relationships",
    ),
    Stage(
        "00c4", "00c4_independence_modality_scoring.py", "phenotype",
        "Score modality, independence, and primary-validation eligibility",
        "data/enriched/transcriptomic_validation_candidates.csv",
    ),

    # ------------------------------------------------------------------
    # Expression acquisition + primary external validation
    # ------------------------------------------------------------------
    # ATLAS_SCIENTIFIC_GATE_STAGE_V1
    Stage(
        "00w", "00w_scientific_automation_gate.py", "phenotype",
        "Scientific automation gate before expression validation",
        "results/pipeline_state/scientific_automation_gate.json",
    ),

    Stage(
        "00e", "00e_primary_validation_expression_fetch.py", "validation",
        "Download processed primary-validation expression matrices",
        "data/validation_expression/primary_expression_manifest.csv",
    ),
    # ATLAS_DATA_VOLUME_STAGE_V1
    Stage(
        "00aa", "00aa_dataset_volume_metadata.py", "metadata",
        "Record dataset counts and local/remote storage volume metadata",
        "results/pipeline_state/dataset_volume_summary.json",
    ),

    Stage(
        "00f", "00f_primary_validation_matrix_inspection.py", "validation",
        "Inspect and map expression matrices",
    ),
    Stage(
        "00f1", "00f1_primary_validation_design_audit.py", "validation",
        "Audit primary-validation experimental designs",
    ),
    Stage(
        "00g", "00g_primary_validation_differential_expression.py", "validation",
        "Run primary-validation differential expression",
        "results/external_validation/primary_validation_DE_summary.csv",
    ),
    Stage(
        "00h", "00h_cross_study_deg_concordance.py", "validation",
        "Compute cross-study DEG concordance",
        "results/external_validation/replicated_primary_DEGs.csv",
    ),
    Stage(
        "00i", "00i_consensus_resistance_signature.py", "validation",
        "Build external consensus resistance signature",
        "results/external_validation/consensus_resistance_signature.csv",
    ),
    Stage(
        "00j", "00j_discovery_consensus_validation.py", "validation",
        "Validate external consensus against discovery",
        "results/external_validation/three_dataset_strict_core_genes.csv",
    ),

    # ------------------------------------------------------------------
    # Mechanism / pathway validation
    # ------------------------------------------------------------------
    Stage(
        "00k", "00k_strict_core_pathway_validation.py", "mechanism",
        "Strict-core pathway validation",
    ),
    Stage(
        "00l", "00l_tgfb_ranked_validation.py", "mechanism",
        "Ranked TGF-beta pathway validation",
    ),
    Stage(
        "00m", "00m_tgfb_gene_level_audit.py", "mechanism",
        "TGF-beta gene-level direction audit",
    ),
    Stage(
        "00n", "00n_tgfb_leading_edge_comparison.py", "mechanism",
        "TGF-beta leading-edge comparison",
    ),
    Stage(
        "00o", "00o_tgfb_reproducible_module_validation.py", "mechanism",
        "Validate reproducible positive TGF-beta module",
    ),
    Stage(
        "00p", "00p_tgfb_module_score_validation.py", "mechanism",
        "Sample-level TGF-beta module scoring",
    ),
    Stage(
        "00qv", "00q_tgfb_dual_module_cross_dataset_audit.py", "mechanism",
        "Cross-dataset dual-module TGF-beta audit",
    ),
    Stage(
        "00r", "00r_tgfb_final_evidence_synthesis.py", "mechanism",
        "Final TGF-beta evidence synthesis",
        "results/external_validation/pathway_validation/tgfb_final_evidence_matrix.csv",
    ),
    Stage(
        "00s", "00s_validated_resistance_evidence_export.py", "mechanism",
        "Export validated resistance evidence for downstream integration",
        "results/external_validation/downstream/validated_resistance_gene_evidence.csv",
    ),

    # ------------------------------------------------------------------
    # CMap / perturbational evidence
    # ------------------------------------------------------------------
    Stage(
        "04g", "04g_cmap_submit_all.py", "cmap",
        "Submit all CMap jobs",
        special="cmap_submit",
    ),
    Stage(
        "04h", "04h_cmap_check_all_jobs.py", "cmap",
        "Poll CMap jobs until complete",
        special="cmap_poll",
    ),
    Stage(
        "04j", "04j_cmap_download_all.py", "cmap",
        "Download completed CMap results",
    ),
    Stage(
        "04k", "04k_cmap_parse_all_tau.py", "cmap",
        "Parse all CMap tau results",
    ),
    Stage(
        "04l", "04l_cmap_prioritize.py", "cmap",
        "Pure CMap prioritization",
    ),
    Stage(
        "04m", "04m_cmap_drug_filter.py", "drug",
        "Compound identity / drug filtering",
    ),
    Stage(
        "04n", "04n_regulatory_status.py", "drug",
        "Regulatory and clinical-trial status",
    ),
    Stage(
        "04o", "04o_safety_screening.py", "drug",
        "Safety / cytotoxicity / PAINS screening",
    ),
    Stage(
        "04p", "04p_drug_target_annotation.py", "drug",
        "Drug-target annotation",
    ),
    Stage(
        "04q", "04q_network_integration.py", "network",
        "Validated resistance-network integration",
        "results/cmap/network_integration/ATLAS_drug_network_prioritized.csv",
        args=("--max-resistance-genes", "242"),
    ),
    Stage(
        "04r", "04r_final_candidate_prioritization.py", "prioritization",
        "Final multi-layer candidate prioritization",
    ),
    Stage(
        "04s", "04s_target_supported_docking.py", "docking",
        "Target-supported docking",
    ),
    Stage(
        "04t", "04t_admet_structural_assessment.py", "docking",
        "ADMET / structural assessment",
    ),
    Stage(
        "04u", "04u_integrated_evidence_matrix.py", "final",
        "Build final integrated evidence matrix",
        "results/cmap/integrated_evidence/ATLAS_integrated_evidence_matrix.csv",
    ),
]


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def load_state():
    if not STATE_FILE.exists():
        return {"created_utc": utcnow(), "stages": {}}
    try:
        return json.loads(STATE_FILE.read_text(encoding="utf-8"))
    except Exception:
        return {"created_utc": utcnow(), "stages": {}}


def save_state(state):
    tmp = STATE_FILE.with_suffix(".tmp")
    tmp.write_text(json.dumps(state, indent=2), encoding="utf-8")
    tmp.replace(STATE_FILE)


def checkpoint_exists(stage: Stage):
    if not stage.checkpoint:
        return False
    return (ROOT / stage.checkpoint).exists()


def run_command(stage: Stage, extra_args=()):
    script_path = SCRIPTS / stage.script
    if not script_path.exists():
        raise FileNotFoundError(f"Missing script: {script_path}")

    cmd = [PYTHON, "-u", str(script_path), *stage.args, *extra_args]

    log_path = LOG_DIR / f"{stage.key}_{stage.script.replace('.py', '')}.log"

    print("\n" + "=" * 88)
    print(f"[{stage.key}] {stage.description}")
    print("=" * 88)
    print("Command:", " ".join(cmd))
    print("Log:", log_path)

    with log_path.open("w", encoding="utf-8") as log:
        proc = subprocess.Popen(
            cmd,
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=os.environ.copy(),
        )

        captured = []
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="")
            log.write(line)
            captured.append(line)

        rc = proc.wait()

    return rc, "".join(captured), log_path


def cmap_is_pending(output: str):
    text = output.lower()

    pending_tokens = [
        "pending",
        "running",
        "queued",
        "submitted",
        "processing",
        "not complete",
        "not completed",
        "in progress",
    ]

    done_tokens = [
        "all jobs complete",
        "all jobs completed",
        "all completed",
        "complete: true",
        "status: completed",
        "status=completed",
    ]

    if any(x in text for x in done_tokens):
        return False

    return any(x in text for x in pending_tokens)


def run_cmap_poll(stage: Stage, poll_minutes: int, max_wait_hours: int):
    deadline = time.time() + max_wait_hours * 3600
    attempt = 0

    while True:
        attempt += 1
        print(f"\nCMap status poll attempt {attempt}")

        rc, output, log_path = run_command(stage)

        if rc != 0:
            return rc, output, log_path

        if not cmap_is_pending(output):
            return rc, output, log_path

        if time.time() >= deadline:
            print(
                f"ERROR: CMap jobs still appear pending after "
                f"{max_wait_hours} hours."
            )
            return 124, output, log_path

        print(f"CMap jobs still pending; sleeping {poll_minutes} minutes.")
        time.sleep(poll_minutes * 60)



# ATLAS_DEPENDENCY_AWARE_RESUME_V2
def stage_script_changed_since_success(stage: Stage, previous: dict) -> bool:
    """True if this stage script changed after its last successful run."""
    ended = previous.get("ended_utc")
    if not ended:
        return False

    script_path = SCRIPTS / stage.script
    if not script_path.exists():
        return True

    try:
        ended_dt = datetime.fromisoformat(ended.replace("Z", "+00:00"))
        script_dt = datetime.fromtimestamp(
            script_path.stat().st_mtime,
            tz=timezone.utc,
        )
        return script_dt > ended_dt
    except Exception:
        return False


def select_stages(from_stage=None, to_stage=None):
    keys = [s.key for s in STAGES]

    start = 0
    end = len(STAGES) - 1

    if from_stage:
        if from_stage not in keys:
            raise SystemExit(
                f"Unknown --from-stage {from_stage}. "
                f"Use --list-stages."
            )
        start = keys.index(from_stage)

    if to_stage:
        if to_stage not in keys:
            raise SystemExit(
                f"Unknown --to-stage {to_stage}. "
                f"Use --list-stages."
            )
        end = keys.index(to_stage)

    if start > end:
        raise SystemExit("--from-stage occurs after --to-stage.")

    return STAGES[start:end + 1]


def parse_args():
    p = argparse.ArgumentParser(
        description="Run ATLAS end-to-end with resume/checkpoints."
    )

    p.add_argument("--full", action="store_true")
    p.add_argument("--refresh-data", action="store_true")
    p.add_argument("--force", action="store_true")

    p.add_argument("--from-stage")
    p.add_argument("--to-stage")
    p.add_argument("--rerun-stage")

    p.add_argument("--list-stages", action="store_true")

    p.add_argument(
        "--cmap-poll-minutes",
        type=int,
        default=10,
    )
    p.add_argument(
        "--cmap-max-wait-hours",
        type=int,
        default=12,
    )

    return p.parse_args()


def main():
    args = parse_args()

    if args.list_stages:
        for s in STAGES:
            print(
                f"{s.key:5s}  {s.group:14s}  "
                f"{s.script:45s}  {s.description}"
            )
        return 0

    if args.rerun_stage:
        selected = [s for s in STAGES if s.key == args.rerun_stage]
        if not selected:
            raise SystemExit(
                f"Unknown stage: {args.rerun_stage}. Use --list-stages."
            )
        args.force = True
    else:
        selected = select_stages(args.from_stage, args.to_stage)

    state = load_state()
    state["last_run_started_utc"] = utcnow()
    save_state(state)

    print("=" * 88)
    print("ATLAS FULL PIPELINE AUTOMATION")
    print("=" * 88)
    print(f"Project root: {ROOT}")
    print(f"Python:       {PYTHON}")
    print(f"Stages:       {selected[0].key} -> {selected[-1].key}")
    print(f"Refresh:      {args.refresh_data}")
    print(f"Force:        {args.force}")

    upstream_reran = False

    for stage in selected:
        previous = state["stages"].get(stage.key, {})

        should_skip = False

        # Runtime preflight always runs when it is selected.
        if stage.group == "preflight":
            should_skip = False

        elif not args.force and not args.refresh_data and not upstream_reran:
            previous_success = previous.get("status") == "SUCCESS"
            script_changed = stage_script_changed_since_success(stage, previous)

            if previous_success and not script_changed:
                if stage.checkpoint:
                    should_skip = checkpoint_exists(stage)
                else:
                    should_skip = True

            elif not previous_success and checkpoint_exists(stage):
                # Backward-compatible resume from a legacy checkpoint.
                should_skip = True

        if should_skip:
            print(
                f"\nSKIP [{stage.key}] {stage.description} "
                f"(successful checkpoint exists)"
            )
            continue

        rec = {
            "key": stage.key,
            "script": stage.script,
            "description": stage.description,
            "started_utc": utcnow(),
            "status": "RUNNING",
        }
        state["stages"][stage.key] = rec
        save_state(state)

        try:
            if stage.special == "cmap_poll":
                rc, output, log_path = run_cmap_poll(
                    stage,
                    poll_minutes=max(1, args.cmap_poll_minutes),
                    max_wait_hours=max(1, args.cmap_max_wait_hours),
                )
            else:
                rc, output, log_path = run_command(stage)

        except Exception as exc:
            rec.update({
                "status": "FAILED",
                "ended_utc": utcnow(),
                "error": repr(exc),
            })
            save_state(state)
            print(f"\nFAILED [{stage.key}]: {exc}")
            return 1

        rec.update({
            "ended_utc": utcnow(),
            "return_code": rc,
            "log": str(log_path),
        })

        if rc != 0:
            rec["status"] = "FAILED"
            save_state(state)

            print("\n" + "!" * 88)
            print(
                f"ATLAS STOPPED: stage {stage.key} failed "
                f"with return code {rc}"
            )
            print(f"See: {log_path}")
            print("!" * 88)
            return rc

        rec["status"] = "SUCCESS"
        save_state(state)

        # Executed data/scientific stages invalidate downstream results.
        # Preflight itself does not dirty the scientific pipeline.
        if stage.group != "preflight":
            upstream_reran = True

        print(f"\nSUCCESS [{stage.key}] {stage.description}")

    state["last_run_ended_utc"] = utcnow()
    state["last_run_status"] = "SUCCESS"
    save_state(state)

    print("\n" + "=" * 88)
    print("ATLAS FULL PIPELINE COMPLETE")
    print("=" * 88)
    print(
        ROOT
        / "results"
        / "cmap"
        / "integrated_evidence"
        / "ATLAS_integrated_evidence_matrix.csv"
    )
    print(f"State: {STATE_FILE}")
    print(f"Logs:  {LOG_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
