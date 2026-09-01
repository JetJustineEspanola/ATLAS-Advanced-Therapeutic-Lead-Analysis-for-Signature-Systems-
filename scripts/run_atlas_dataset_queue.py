#!/usr/bin/env python3
"""
ATLAS continuous dataset queue worker.

Purpose
-------
Continuously sift candidate GEO datasets until a deadline. Each candidate is
screened one-by-one. If it becomes PRIMARY_VALIDATION, ATLAS re-runs the
validation chain and only re-runs CMap/downstream analysis when the validated
resistance evidence actually changes.

This worker is intentionally conservative:
- it never treats an exploratory dataset as primary evidence;
- it records datasets the current DE driver cannot actually process;
- it does not submit another CMap job unless validated resistance evidence changed;
- it persists queue/history state so the worker can resume after a restart.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import pandas as pd


ROOT = Path("/home/regulus/Documents/ATLAS")
PY = ROOT / ".venv/bin/python"
RUNNER = ROOT / "scripts/run_atlas_full_auto.py"
TRACKER = ROOT / "scripts/00aa_dataset_volume_metadata.py"
ENRICH = ROOT / "scripts/00d_metadata_enrichment.py"

STATE_DIR = ROOT / "results/pipeline_state/dataset_queue"
QUEUE_CSV = STATE_DIR / "dataset_queue.csv"
HISTORY_CSV = STATE_DIR / "dataset_queue_history.csv"
LOCK_FILE = STATE_DIR / "dataset_queue.lock"

INDEPENDENCE = ROOT / "data/enriched/dataset_candidates_independence_scored.csv"
TRANSCRIPTOMIC = ROOT / "data/enriched/transcriptomic_validation_candidates.csv"
DE_SUMMARY = ROOT / "results/external_validation/primary_validation_DE_summary.csv"

RESISTANCE_FILES = [
    ROOT / "results/external_validation/downstream/validated_resistance_gene_evidence.csv",
    ROOT / "results/external_validation/downstream/validated_tgfb_module_evidence.csv",
]

MANILA = ZoneInfo("Asia/Manila")

# Storage safety guard:
# - ATLAS may use up to 40 GB under the project root.
# - Also preserve at least 20 GB free on the filesystem.
MAX_ATLAS_GB = 40.0
MIN_FREE_GB = 20.0

TERMINAL_STATUSES = {
    "COMPLETE_PRIMARY",
    "SCREENED_NOT_PRIMARY",
    "CURRENT_DE_DRIVER_UNSUPPORTED",
    "REJECTED_NON_TRANSCRIPTOMIC",
    "FAILED_PERMANENT",
}

ACTIVE_OTHER_SERVICES = ("atlas-monitor.service", "atlas-full.service")


def now() -> datetime:
    return datetime.now(MANILA)


def stamp() -> str:
    return now().isoformat(timespec="seconds")


def log(msg: str = "") -> None:
    print(msg, flush=True)


def run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    log("$ " + " ".join(map(str, cmd)))
    proc = subprocess.run(
        [str(x) for x in cmd],
        cwd=ROOT,
        env=os.environ.copy(),
        text=True,
    )
    if check and proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)
    return proc


def service_active(name: str) -> bool:
    proc = subprocess.run(
        ["systemctl", "--user", "is-active", "--quiet", name],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return proc.returncode == 0



def bytes_to_gb(n: int) -> float:
    return n / (1024 ** 3)


def directory_size_bytes(path: Path) -> int:
    total = 0
    if not path.exists():
        return 0
    for p in path.rglob("*"):
        try:
            if p.is_file() and not p.is_symlink():
                total += p.stat().st_size
        except (FileNotFoundError, PermissionError):
            continue
    return total


def storage_guard() -> tuple[bool, str]:
    """Return (safe_to_continue, human-readable status)."""
    usage = os.statvfs(ROOT)
    free_bytes = usage.f_bavail * usage.f_frsize
    free_gb = bytes_to_gb(free_bytes)

    atlas_bytes = directory_size_bytes(ROOT)
    atlas_gb = bytes_to_gb(atlas_bytes)

    msg = (
        f"ATLAS storage={atlas_gb:.2f} GB / {MAX_ATLAS_GB:.0f} GB max; "
        f"filesystem free={free_gb:.2f} GB / {MIN_FREE_GB:.0f} GB minimum"
    )

    if atlas_gb >= MAX_ATLAS_GB:
        return False, "STORAGE_LIMIT_REACHED: " + msg

    if free_gb <= MIN_FREE_GB:
        return False, "FREE_SPACE_GUARD_TRIGGERED: " + msg

    return True, msg


def enforce_storage_guard() -> None:
    ok, msg = storage_guard()
    log(msg)
    if not ok:
        run_tracker()
        raise SystemExit(0)


def wait_for_existing_atlas_jobs(deadline: datetime) -> None:
    waiting_logged = False
    while now() < deadline:
        active = [s for s in ACTIVE_OTHER_SERVICES if service_active(s)]
        if not active:
            if waiting_logged:
                log("Existing ATLAS job finished. Dataset queue is taking over.")
            return
        if not waiting_logged:
            log("Existing ATLAS run detected: " + ", ".join(active))
            log("Queue worker will wait and start immediately after it finishes.")
            waiting_logged = True
        time.sleep(15)


def file_hash(paths: list[Path]) -> str:
    h = hashlib.sha256()
    found = False
    for p in sorted(paths, key=lambda x: str(x)):
        if not p.is_file():
            continue
        found = True
        h.update(str(p.relative_to(ROOT)).encode())
        with p.open("rb") as f:
            while True:
                chunk = f.read(1024 * 1024)
                if not chunk:
                    break
                h.update(chunk)
    return h.hexdigest() if found else "MISSING"


def append_history(accession: str, event: str, detail: str = "") -> None:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    new = not HISTORY_CSV.exists()
    with HISTORY_CSV.open("a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["timestamp", "accession", "event", "detail"])
        w.writerow([stamp(), accession, event, detail])


def detect_accession_column(df: pd.DataFrame) -> str:
    for col in ("source_accession", "accession", "dataset_accession", "geo_accession"):
        if col in df.columns:
            return col
    raise RuntimeError("No accession column found in candidate metadata.")


def load_candidate_table() -> pd.DataFrame:
    # Prefer the independence-aware table because 00C4 is authoritative.
    src = INDEPENDENCE if INDEPENDENCE.exists() else TRANSCRIPTOMIC
    if not src.exists():
        raise FileNotFoundError(
            "No candidate table found. Expected "
            f"{INDEPENDENCE} or {TRANSCRIPTOMIC}"
        )

    df = pd.read_csv(src, low_memory=False)
    acc_col = detect_accession_column(df)
    df = df.rename(columns={acc_col: "accession"})
    df["accession"] = df["accession"].astype(str).str.strip()
    df = df[df["accession"].str.match(r"^GSE\d+$", na=False)].copy()

    # One queue item per GEO Series.
    df = df.drop_duplicates("accession", keep="first")

    return df


def existing_processed_accessions() -> set[str]:
    if not DE_SUMMARY.exists():
        return set()
    try:
        df = pd.read_csv(DE_SUMMARY)
    except Exception:
        return set()
    for col in ("accession", "source_accession", "dataset_accession"):
        if col in df.columns:
            return set(df[col].dropna().astype(str).str.strip())
    return set()


def text_from_row(row: pd.Series, names: tuple[str, ...]) -> str:
    vals = []
    for n in names:
        if n in row.index and pd.notna(row[n]):
            vals.append(str(row[n]))
    return " ".join(vals)


def is_obvious_non_transcriptomic(row: pd.Series) -> bool:
    txt = text_from_row(
        row,
        (
            "modality",
            "modality_class",
            "assay",
            "technology",
            "title",
            "summary",
            "independence_aware_category",
        ),
    ).lower()

    non_tx = ("atac", "cut&tag", "cut-tag", "chip-seq", "chip seq", "methyl", "merip")
    tx = ("rna", "expression", "transcript", "microarray", "rna-seq", "rnaseq")

    return any(x in txt for x in non_tx) and not any(x in txt for x in tx)


def numeric_priority(row: pd.Series) -> float:
    for col in (
        "independence_aware_score",
        "independence_score",
        "primary_validation_score",
        "validation_score",
        "eligibility_score",
        "score",
    ):
        if col in row.index:
            try:
                return float(row[col])
            except Exception:
                pass
    return 0.0


def category_of(row: pd.Series) -> str:
    for col in (
        "independence_aware_category",
        "validation_category",
        "category",
        "downstream_category",
    ):
        if col in row.index and pd.notna(row[col]):
            return str(row[col]).strip()
    return ""


def load_queue() -> pd.DataFrame:
    cols = [
        "accession",
        "status",
        "attempts",
        "priority",
        "category",
        "first_seen",
        "started_at",
        "finished_at",
        "last_error",
    ]
    if not QUEUE_CSV.exists():
        return pd.DataFrame(columns=cols)

    # Read queue state as text first. Empty CSV fields otherwise become NaN/float64,
    # which prevents later assignment of ISO timestamp strings under newer pandas.
    q = pd.read_csv(QUEUE_CSV, dtype=str, keep_default_na=False)

    for col in cols:
        if col not in q.columns:
            q[col] = ""

    string_cols = [
        "accession",
        "status",
        "category",
        "first_seen",
        "started_at",
        "finished_at",
        "last_error",
    ]
    for col in string_cols:
        q[col] = q[col].fillna("").astype(str)

    q["attempts"] = (
        pd.to_numeric(q["attempts"], errors="coerce")
        .fillna(0)
        .astype(int)
    )
    q["priority"] = (
        pd.to_numeric(q["priority"], errors="coerce")
        .fillna(0.0)
        .astype(float)
    )

    return q[cols]


def save_queue(q: pd.DataFrame) -> None:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    q.to_csv(QUEUE_CSV, index=False)


def refresh_queue_from_catalog() -> pd.DataFrame:
    candidates = load_candidate_table()
    processed = existing_processed_accessions()
    q = load_queue()
    by_acc = {str(r["accession"]): i for i, r in q.iterrows()}

    for _, row in candidates.iterrows():
        acc = str(row["accession"])
        category = category_of(row)
        priority = numeric_priority(row)

        if acc in by_acc:
            i = by_acc[acc]
            q.at[i, "priority"] = priority
            q.at[i, "category"] = category
            continue

        if acc in processed:
            status = "COMPLETE_PRIMARY"
        elif is_obvious_non_transcriptomic(row):
            status = "REJECTED_NON_TRANSCRIPTOMIC"
        else:
            status = "QUEUED"

        q.loc[len(q)] = {
            "accession": acc,
            "status": status,
            "attempts": 0,
            "priority": priority,
            "category": category,
            "first_seen": stamp(),
            "started_at": "",
            "finished_at": stamp() if status != "QUEUED" else "",
            "last_error": "",
        }
        append_history(acc, "DISCOVERED", f"status={status}; category={category}")

    # Existing PRIMARY_VALIDATION datasets should not be reprocessed endlessly.
    for i, row in q.iterrows():
        if row["accession"] in processed and row["status"] not in ("RUNNING",):
            q.at[i, "status"] = "COMPLETE_PRIMARY"

    q = q.sort_values(
        by=["status", "priority", "accession"],
        key=lambda s: s.map({"QUEUED": 0, "FAILED_RETRYABLE": 1}).fillna(9)
        if s.name == "status"
        else s,
        ascending=[True, False, True],
        kind="stable",
    ).reset_index(drop=True)

    save_queue(q)
    return q


def candidate_row(accession: str) -> pd.Series | None:
    df = load_candidate_table()
    rows = df[df["accession"] == accession]
    if rows.empty:
        return None
    return rows.iloc[0]


def update_queue_row(accession: str, **updates) -> None:
    q = load_queue()
    mask = q["accession"] == accession
    if not mask.any():
        raise RuntimeError(f"Queue row disappeared for {accession}")
    i = q.index[mask][0]
    for k, v in updates.items():
        q.at[i, k] = v
    save_queue(q)


def select_next(max_attempts: int) -> str | None:
    q = load_queue()
    if q.empty:
        return None

    eligible = q[
        ((q["status"] == "QUEUED") | (q["status"] == "FAILED_RETRYABLE"))
        & (q["attempts"] < max_attempts)
    ].copy()

    if eligible.empty:
        return None

    eligible = eligible.sort_values(
        ["priority", "attempts", "first_seen"],
        ascending=[False, True, True],
        kind="stable",
    )
    return str(eligible.iloc[0]["accession"])


def run_tracker() -> None:
    if TRACKER.exists():
        run([PY, "-u", TRACKER])


def refresh_online_catalog() -> None:
    log("")
    log("=== REFRESHING ONLINE DISCOVERY / CATALOG ===")
    run(
        [
            PY,
            "-u",
            RUNNER,
            "--from-stage",
            "00a",
            "--to-stage",
            "00c4",
            "--refresh-data",
        ]
    )
    run_tracker()


def screen_one(accession: str) -> tuple[bool, str]:
    log("")
    log("=" * 80)
    log(f"SCREENING NEXT DATASET: {accession}")
    log("=" * 80)

    # Targeted metadata enrichment.
    run([PY, "-u", ENRICH, "--accessions", accession])

    # Re-run phenotype, relationship and independence scoring from enriched metadata.
    run(
        [
            PY,
            "-u",
            RUNNER,
            "--from-stage",
            "00d1",
            "--to-stage",
            "00c4",
            "--force",
        ]
    )

    run_tracker()

    row = candidate_row(accession)
    if row is None:
        return False, "Candidate disappeared after enrichment"

    category = category_of(row)
    update_queue_row(accession, category=category)

    if category == "PRIMARY_VALIDATION":
        return True, category

    return False, category or "NO_PRIMARY_CATEGORY"


def run_validation_for_primary(accession: str) -> bool:
    before = file_hash(RESISTANCE_FILES)

    log("")
    log(f"{accession} passed PRIMARY_VALIDATION.")
    log("Running 00W -> 00S validation chain.")

    run(
        [
            PY,
            "-u",
            RUNNER,
            "--from-stage",
            "00w",
            "--to-stage",
            "00s",
            "--force",
        ]
    )
    run_tracker()

    # Critical guard: confirm the current DE driver actually included this accession.
    processed = existing_processed_accessions()
    if accession not in processed:
        update_queue_row(
            accession,
            status="CURRENT_DE_DRIVER_UNSUPPORTED",
            finished_at=stamp(),
            last_error=(
                "00C4 classified as PRIMARY_VALIDATION but accession was not present "
                "in primary_validation_DE_summary.csv after 00W->00S."
            ),
        )
        append_history(
            accession,
            "CURRENT_DE_DRIVER_UNSUPPORTED",
            "Not present in primary_validation_DE_summary.csv",
        )
        log(
            f"WARNING: {accession} became PRIMARY_VALIDATION, but the current 00E/00G "
            "driver did not actually include it. It will NOT be counted as validated."
        )
        return False

    after = file_hash(RESISTANCE_FILES)
    changed = before != after

    update_queue_row(
        accession,
        status="COMPLETE_PRIMARY",
        finished_at=stamp(),
        last_error="",
    )
    append_history(accession, "COMPLETE_PRIMARY", f"resistance_changed={changed}")

    if changed:
        log("")
        log("Validated resistance evidence changed.")
        log("Running one new CMap/downstream analysis: 04G -> 04U.")
        run(
            [
                PY,
                "-u",
                RUNNER,
                "--from-stage",
                "04g",
                "--to-stage",
                "04u",
                "--force",
            ]
        )
        run_tracker()
    else:
        log("")
        log("Validated resistance evidence is unchanged.")
        log("Reusing existing CMap / docking / 04U outputs.")

    return True


def process_accession(accession: str, max_attempts: int) -> None:
    q = load_queue()
    row = q[q["accession"] == accession].iloc[0]
    attempts = int(row["attempts"]) + 1

    update_queue_row(
        accession,
        status="RUNNING",
        attempts=attempts,
        started_at=stamp(),
        last_error="",
    )
    append_history(accession, "START", f"attempt={attempts}")

    try:
        primary, category = screen_one(accession)

        if not primary:
            update_queue_row(
                accession,
                status="SCREENED_NOT_PRIMARY",
                finished_at=stamp(),
                last_error="",
            )
            append_history(accession, "SCREENED_NOT_PRIMARY", f"category={category}")
            log(f"{accession}: not PRIMARY_VALIDATION ({category}).")
            log("Moving immediately to the next dataset.")
            return

        run_validation_for_primary(accession)

    except subprocess.CalledProcessError as e:
        status = "FAILED_RETRYABLE" if attempts < max_attempts else "FAILED_PERMANENT"
        msg = f"command failed with exit code {e.returncode}"
        update_queue_row(
            accession,
            status=status,
            finished_at=stamp(),
            last_error=msg,
        )
        append_history(accession, status, msg)
        log(f"{accession}: {msg}; status={status}")

    except Exception as e:
        status = "FAILED_RETRYABLE" if attempts < max_attempts else "FAILED_PERMANENT"
        msg = f"{type(e).__name__}: {e}"
        update_queue_row(
            accession,
            status=status,
            finished_at=stamp(),
            last_error=msg,
        )
        append_history(accession, status, msg)
        log(f"{accession}: {msg}; status={status}")


def parse_deadline(value: str) -> datetime:
    # Accept "2026-09-01 05:00:00" or ISO timestamps.
    value = value.strip()
    try:
        dt = datetime.fromisoformat(value)
    except ValueError as e:
        raise argparse.ArgumentTypeError(str(e))
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=MANILA)
    return dt.astimezone(MANILA)


def acquire_lock():
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    fh = LOCK_FILE.open("w")
    try:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError("Another ATLAS dataset queue worker is already running.")
    fh.write(f"pid={os.getpid()}\nstarted={stamp()}\n")
    fh.flush()
    return fh


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--deadline",
        type=parse_deadline,
        required=True,
        help="Stop time in Asia/Manila if no timezone is supplied, e.g. '2026-09-01 05:00:00'",
    )
    ap.add_argument(
        "--idle-seconds",
        type=int,
        default=300,
        help="Wait before refreshing discovery again when the queue is empty.",
    )
    ap.add_argument(
        "--max-attempts",
        type=int,
        default=2,
        help="Maximum attempts for a dataset that fails with a retryable error.",
    )
    args = ap.parse_args()

    if args.deadline <= now():
        log(f"Deadline already passed: {args.deadline.isoformat()}")
        return 0

    lock_fh = acquire_lock()
    _ = lock_fh  # keep handle alive

    log("=" * 80)
    log("ATLAS CONTINUOUS DATASET QUEUE")
    log("=" * 80)
    log(f"Started:  {stamp()}")
    log(f"Deadline: {args.deadline.isoformat()}")
    log(f"Queue:    {QUEUE_CSV}")
    log(f"History:  {HISTORY_CSV}")
    log("")
    log("Rule: dataset finishes -> immediately screen/process next eligible dataset.")
    log("CMap is rerun only when validated resistance evidence changes.")
    log("=" * 80)

    # Do not collide with the currently running full/monitor job.
    wait_for_existing_atlas_jobs(args.deadline)

    # Preflight once before taking over.
    run(
        [
            PY,
            "-u",
            RUNNER,
            "--from-stage",
            "00y",
            "--to-stage",
            "00y",
        ]
    )

    # Use current catalog first so we do not waste the first cycle.
    refresh_queue_from_catalog()

    while now() < args.deadline:
        accession = select_next(args.max_attempts)

        if accession is not None:
            log("")
            log("=== STORAGE SAFETY CHECK ===")
            enforce_storage_guard()

            process_accession(accession, args.max_attempts)
            refresh_queue_from_catalog()
            # Deliberately no sleep here: immediately take the next dataset.
            continue

        log("")
        log("No unprocessed candidate currently in the queue.")
        log("Refreshing online discovery for additional datasets.")

        try:
            log("")
            log("=== STORAGE SAFETY CHECK BEFORE DISCOVERY REFRESH ===")
            enforce_storage_guard()

            refresh_online_catalog()
            refresh_queue_from_catalog()
        except Exception as e:
            log(f"Discovery refresh failed: {type(e).__name__}: {e}")

        if select_next(args.max_attempts) is not None:
            continue

        remaining = (args.deadline - now()).total_seconds()
        if remaining <= 0:
            break
        sleep_for = min(args.idle_seconds, max(1, int(remaining)))
        log(f"No new eligible dataset yet. Rechecking in {sleep_for} seconds.")
        time.sleep(sleep_for)

    run_tracker()

    log("")
    log("=" * 80)
    log("ATLAS DATASET QUEUE DEADLINE REACHED")
    log("=" * 80)
    log(f"Finished: {stamp()}")
    log(f"Queue state: {QUEUE_CSV}")
    log(f"History:     {HISTORY_CSV}")
    log("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
