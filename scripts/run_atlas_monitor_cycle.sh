#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/home/regulus/Documents/ATLAS"
PY="$ROOT/.venv/bin/python"
RUNNER="$ROOT/scripts/run_atlas_full_auto.py"
TRACKER="$ROOT/scripts/00aa_dataset_volume_metadata.py"

STATE="$ROOT/results/pipeline_state/monitor"
LOGS="$ROOT/logs/monitor"

mkdir -p "$STATE" "$LOGS"

STAMP="$(date '+%Y%m%d_%H%M%S')"
LOG="$LOGS/cycle_${STAMP}.log"

exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "ATLAS MONITOR CYCLE"
echo "Started: $(date --iso-8601=seconds)"
echo "============================================================"

cd "$ROOT"

# ------------------------------------------------------------
# Helper: deterministic hash of files that currently exist
# ------------------------------------------------------------

hash_files() {
    local found=0

    for f in "$@"; do
        if [[ -f "$f" ]]; then
            sha256sum "$f"
            found=1
        fi
    done | sort | sha256sum | awk '{print $1}'

    return 0
}

# ------------------------------------------------------------
# 0. Runtime preflight
# ------------------------------------------------------------

echo
echo "[0] Runtime preflight"

"$PY" -u "$RUNNER" \
    --from-stage 00y \
    --to-stage 00y

# ------------------------------------------------------------
# 1. Fingerprint current acquisition/catalog state
# ------------------------------------------------------------

echo
echo "[1] Baseline acquisition fingerprint"

ACQ_FILES=(
    "$ROOT/data/enriched/dataset_candidates_independence_scored.csv"
    "$ROOT/data/enriched/transcriptomic_validation_candidates.csv"
    "$ROOT/data/enriched/dataset_candidates.csv"
    "$ROOT/data/catalog/dataset_candidates.csv"
    "$ROOT/data/catalog/candidates.csv"
)

BEFORE_ACQ="$(hash_files "${ACQ_FILES[@]}")"

echo "Before acquisition hash: $BEFORE_ACQ"

# ------------------------------------------------------------
# 2. Refresh only discovery / metadata / eligibility
# ------------------------------------------------------------

echo
echo "[2] Refreshing online acquisition and metadata"

"$PY" -u "$RUNNER" \
    --from-stage 00a \
    --to-stage 00c4 \
    --refresh-data

# ------------------------------------------------------------
# 3. Always update storage metadata
# ------------------------------------------------------------

echo
echo "[3] Updating 00AA dataset/storage metadata"

"$PY" -u "$TRACKER"

AFTER_ACQ="$(hash_files "${ACQ_FILES[@]}")"

echo "After acquisition hash:  $AFTER_ACQ"

printf '%s\n' "$AFTER_ACQ" > "$STATE/latest_acquisition.sha256"

# ------------------------------------------------------------
# 4. Nothing changed -> end cycle
# ------------------------------------------------------------

if [[ "$BEFORE_ACQ" == "$AFTER_ACQ" ]]; then
    echo
    echo "============================================================"
    echo "NO_NEW_DATA"
    echo "Acquisition metadata is unchanged."
    echo "Skipping validation, CMap, docking, and 04U."
    echo "Finished: $(date --iso-8601=seconds)"
    echo "============================================================"
    exit 0
fi

echo
echo "DATA CHANGE DETECTED."
echo "Recomputing scientific validation."

# ------------------------------------------------------------
# 5. Fingerprint current validated resistance evidence
# ------------------------------------------------------------

RESISTANCE_FILES=(
    "$ROOT/results/external_validation/downstream/validated_resistance_gene_evidence.csv"
    "$ROOT/results/external_validation/downstream/validated_tgfb_module_evidence.csv"
)

BEFORE_RESISTANCE="$(hash_files "${RESISTANCE_FILES[@]}")"

echo "Before resistance hash: $BEFORE_RESISTANCE"

# ------------------------------------------------------------
# 6. Re-run validation through validated evidence export
# ------------------------------------------------------------

echo
echo "[4] Running validation 00W -> 00S"

"$PY" -u "$RUNNER" \
    --from-stage 00w \
    --to-stage 00s \
    --force

# Update volume metadata again because validation may add files
"$PY" -u "$TRACKER"

AFTER_RESISTANCE="$(hash_files "${RESISTANCE_FILES[@]}")"

echo "After resistance hash:  $AFTER_RESISTANCE"

printf '%s\n' "$AFTER_RESISTANCE" \
    > "$STATE/latest_resistance.sha256"

# ------------------------------------------------------------
# 7. Validated resistance evidence unchanged -> reuse downstream
# ------------------------------------------------------------

if [[ "$BEFORE_RESISTANCE" == "$AFTER_RESISTANCE" ]]; then
    echo
    echo "============================================================"
    echo "VALIDATED_SIGNATURE_UNCHANGED"
    echo "New acquisition metadata did not change validated resistance evidence."
    echo "Existing CMap and downstream results are reused."
    echo "Finished: $(date --iso-8601=seconds)"
    echo "============================================================"
    exit 0
fi

# ------------------------------------------------------------
# 8. Resistance signature changed -> downstream recomputation
# ------------------------------------------------------------

echo
echo "VALIDATED RESISTANCE SIGNATURE CHANGED."
echo "Running CMap and downstream stages."

"$PY" -u "$RUNNER" \
    --from-stage 04g \
    --to-stage 04u \
    --force

# ------------------------------------------------------------
# 9. Final storage snapshot
# ------------------------------------------------------------

"$PY" -u "$TRACKER"

echo
echo "============================================================"
echo "ATLAS MONITOR CYCLE COMPLETE"
echo "Final evidence:"
echo "$ROOT/results/cmap/integrated_evidence/ATLAS_integrated_evidence_matrix.csv"
echo
echo "Finished: $(date --iso-8601=seconds)"
echo "============================================================"
