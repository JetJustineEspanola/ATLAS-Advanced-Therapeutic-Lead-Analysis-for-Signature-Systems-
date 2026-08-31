#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# ATLAS — Unattended Setup + Full Pipeline Launcher
# ============================================================

ATLAS_ROOT="/home/regulus/Documents/ATLAS"
VENV="$ATLAS_ROOT/.venv"
SCRIPTS="$ATLAS_ROOT/scripts"
RESULTS="$ATLAS_ROOT/results"
LOGS="$ATLAS_ROOT/logs"
STATE="$RESULTS/pipeline_state"

echo "============================================================"
echo "ATLAS — UNATTENDED PIPELINE SETUP"
echo "============================================================"
echo "Root: $ATLAS_ROOT"

cd "$ATLAS_ROOT"

# ------------------------------------------------------------
# 1. Validate project
# ------------------------------------------------------------

if [[ ! -d "$VENV" ]]; then
    echo "ERROR: virtual environment not found:"
    echo "$VENV"
    exit 1
fi

if [[ ! -f "$SCRIPTS/run_atlas_full_auto.py" ]]; then
    echo "ERROR: full automation runner not found."
    exit 1
fi

source "$VENV/bin/activate"

PYTHON="$VENV/bin/python"

echo "Python: $($PYTHON --version)"

# ------------------------------------------------------------
# 2. Required directories
# ------------------------------------------------------------

mkdir -p \
    "$RESULTS" \
    "$STATE" \
    "$LOGS" \
    "$LOGS/full_pipeline" \
    "$ATLAS_ROOT/data"

# ------------------------------------------------------------
# 3. Verify metadata tracker
# ------------------------------------------------------------

TRACKER="$SCRIPTS/00aa_dataset_volume_metadata.py"
PATCHER="$SCRIPTS/00ab_patch_runner_data_volume.py"

if [[ ! -f "$TRACKER" ]]; then
    echo "ERROR: metadata tracker missing:"
    echo "$TRACKER"
    exit 1
fi

echo "PASS: dataset-volume metadata tracker installed."

# ------------------------------------------------------------
# 4. Install tracker into full runner if necessary
# ------------------------------------------------------------

if grep -q '00aa_dataset_volume_metadata.py' \
    "$SCRIPTS/run_atlas_full_auto.py"; then

    echo "PASS: 00AA already registered in full runner."

else
    if [[ ! -f "$PATCHER" ]]; then
        echo "ERROR: runner does not contain 00AA and patcher is missing:"
        echo "$PATCHER"
        exit 1
    fi

    echo "Installing dataset-volume stage into runner..."
    "$PYTHON" -u "$PATCHER"
fi

# ------------------------------------------------------------
# 5. Syntax validation
# ------------------------------------------------------------

echo
echo "Checking automation scripts..."

"$PYTHON" -m py_compile \
    "$SCRIPTS/run_atlas_full_auto.py" \
    "$TRACKER"

echo "PASS: syntax checks."

# ------------------------------------------------------------
# 6. Confirm tracker appears in runner
# ------------------------------------------------------------

echo
echo "Checking registered stages..."

STAGES="$("$PYTHON" -u "$SCRIPTS/run_atlas_full_auto.py" --list-stages)"

echo "$STAGES"

if ! grep -q "00aa.*00aa_dataset_volume_metadata.py" <<< "$STAGES"; then
    echo
    echo "ERROR: 00AA metadata stage is not registered."
    exit 1
fi

echo
echo "PASS: metadata tracker registered."

# ------------------------------------------------------------
# 7. Credential check
# ------------------------------------------------------------

if [[ -z "${CLUE_API_KEY:-}" ]]; then
    echo
    echo "ERROR: CLUE_API_KEY is not available."
    echo
    echo "For true unattended execution, export it before running:"
    echo 'export CLUE_API_KEY="YOUR_KEY"'
    echo
    echo "Prefer storing/loading it from a protected environment file"
    echo "rather than hard-coding it into this script."
    exit 1
fi

echo "PASS: CLUE_API_KEY available."

# ------------------------------------------------------------
# 8. Runtime/network/dependency preflight
# ------------------------------------------------------------

echo
echo "============================================================"
echo "ATLAS PREFLIGHT"
echo "============================================================"

"$PYTHON" -u "$SCRIPTS/run_atlas_full_auto.py" \
    --from-stage 00y \
    --to-stage 00y

# ------------------------------------------------------------
# 9. Scientific gate sanity check
# ------------------------------------------------------------

echo
echo "============================================================"
echo "ATLAS SCIENTIFIC GATE"
echo "============================================================"

"$PYTHON" -u "$SCRIPTS/00w_scientific_automation_gate.py"

# ------------------------------------------------------------
# 10. Initial storage snapshot
# ------------------------------------------------------------

echo
echo "============================================================"
echo "INITIAL DATASET/STORAGE SNAPSHOT"
echo "============================================================"

"$PYTHON" -u "$TRACKER"

# ------------------------------------------------------------
# 11. Launch complete unattended pipeline
# ------------------------------------------------------------

echo
echo "============================================================"
echo "STARTING ATLAS FULL UNATTENDED PIPELINE"
echo "============================================================"

"$PYTHON" -u "$SCRIPTS/run_atlas_full_auto.py" \
    --from-stage 00y \
    --to-stage 04u \
    "$@"

# ------------------------------------------------------------
# 12. Final storage snapshot
# ------------------------------------------------------------

echo
echo "============================================================"
echo "FINAL DATASET/STORAGE SNAPSHOT"
echo "============================================================"

"$PYTHON" -u "$TRACKER"

# ------------------------------------------------------------
# 13. Verify final artifact
# ------------------------------------------------------------

FINAL="$RESULTS/cmap/integrated_evidence/ATLAS_integrated_evidence_matrix.csv"

if [[ ! -s "$FINAL" ]]; then
    echo
    echo "ERROR: final ATLAS evidence matrix was not produced."
    exit 1
fi

# ------------------------------------------------------------
# 14. Complete
# ------------------------------------------------------------

echo
echo "============================================================"
echo "ATLAS UNATTENDED RUN COMPLETE"
echo "============================================================"
echo
echo "Final evidence:"
echo "  $FINAL"
echo
echo "Experimental shortlist:"
echo "  $RESULTS/cmap/integrated_evidence/ATLAS_experimental_validation_shortlist.csv"
echo
echo "Dataset/storage metadata:"
echo "  $STATE/dataset_volume_summary.json"
echo "  $STATE/dataset_volume_summary.csv"
echo "  $STATE/dataset_volume_by_accession.csv"
echo "  $STATE/dataset_volume_files.csv"
echo "  $STATE/dataset_volume_summary.txt"
echo
echo "Pipeline state:"
echo "  $STATE/full_pipeline_state.json"
echo
echo "Logs:"
echo "  $LOGS/full_pipeline"
echo "============================================================"
