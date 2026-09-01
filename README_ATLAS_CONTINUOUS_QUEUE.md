# ATLAS Continuous Dataset Queue — Overnight Run

This adds a continuous dataset-by-dataset worker to ATLAS.

## Behavior

```text
current ATLAS run
      ↓
wait for it to finish
      ↓
build persistent candidate queue
      ↓
next GEO dataset
      ↓
targeted metadata enrichment
      ↓
phenotype / relationship / independence scoring
      ↓
PRIMARY_VALIDATION?
   ├─ no  → record reason → NEXT DATASET immediately
   └─ yes → 00W → 00S validation
                    ↓
          current DE driver actually used it?
             ├─ no  → record unsupported → NEXT
             └─ yes
                    ↓
          validated resistance changed?
             ├─ no  → reuse CMap → NEXT
             └─ yes → 04G → 04U → NEXT
```

The worker stops at **2026-09-01 05:00 Asia/Manila**.

It does not weaken ATLAS scientific gates to increase the dataset count.

## Files

Install:

```text
run_atlas_dataset_queue.py
atlas-dataset-queue.service
setup_atlas_dataset_queue.sh
```

Persistent state:

```text
results/pipeline_state/dataset_queue/dataset_queue.csv
results/pipeline_state/dataset_queue/dataset_queue_history.csv
results/pipeline_state/dataset_queue/dataset_queue.lock
```

Storage metadata continues to be written by `00AA`:

```text
results/pipeline_state/dataset_volume_summary.json
results/pipeline_state/dataset_volume_summary.csv
results/pipeline_state/dataset_volume_by_accession.csv
results/pipeline_state/dataset_volume_files.csv
```

## Install

Assuming the downloaded files are in `~/Downloads`:

```bash
cd /home/regulus/Documents/ATLAS
source .venv/bin/activate

cp ~/Downloads/run_atlas_dataset_queue.py scripts/run_atlas_dataset_queue.py
chmod +x scripts/run_atlas_dataset_queue.py

cp ~/Downloads/atlas-dataset-queue.service \
   ~/.config/systemd/user/atlas-dataset-queue.service

python -m py_compile scripts/run_atlas_dataset_queue.py
systemctl --user daemon-reload
```

Disable the old timers so they do not start overlapping jobs:

```bash
systemctl --user disable --now atlas-monitor.timer 2>/dev/null || true
systemctl --user disable --now atlas-test.timer 2>/dev/null || true
systemctl --user disable --now atlas-full.timer 2>/dev/null || true
```

Do **not** kill the currently active ATLAS analysis if you want it to finish. The queue worker checks for `atlas-monitor.service` and `atlas-full.service`, waits for them to finish, then starts immediately.

## Start

```bash
systemctl --user enable --now atlas-dataset-queue.service
```

Check:

```bash
systemctl --user status atlas-dataset-queue.service
```

Live logs:

```bash
journalctl --user -u atlas-dataset-queue.service -f
```

`Ctrl+C` only exits the log viewer.

## Check queue

```bash
column -s, -t < \
/home/regulus/Documents/ATLAS/results/pipeline_state/dataset_queue/dataset_queue.csv | less -S
```

Quick status counts:

```bash
python - <<'PY'
import pandas as pd

p = "/home/regulus/Documents/ATLAS/results/pipeline_state/dataset_queue/dataset_queue.csv"
df = pd.read_csv(p)
print(df["status"].value_counts(dropna=False).to_string())
PY
```

## Check history

```bash
tail -n 30 \
/home/regulus/Documents/ATLAS/results/pipeline_state/dataset_queue/dataset_queue_history.csv
```

## Stop manually

```bash
systemctl --user stop atlas-dataset-queue.service
```

Disable:

```bash
systemctl --user disable atlas-dataset-queue.service
```

## Important guard

If a newly discovered dataset is classified `PRIMARY_VALIDATION` by `00C4` but the current `00E/00G` implementation does not actually include it in:

```text
results/external_validation/primary_validation_DE_summary.csv
```

the worker records:

```text
CURRENT_DE_DRIVER_UNSUPPORTED
```

and does **not** pretend that dataset contributed to the validated resistance signature. This prevents a false increase in the validation evidence count.

## After 5 AM

The service exits normally at the deadline. Check:

```bash
systemctl --user status atlas-dataset-queue.service
```

and inspect:

```bash
cat \
/home/regulus/Documents/ATLAS/results/pipeline_state/dataset_queue/dataset_queue.csv
```

Run streamlit

source /home/regulus/Documents/ATLAS/.venv/bin/activate
streamlit run ui/app.py
