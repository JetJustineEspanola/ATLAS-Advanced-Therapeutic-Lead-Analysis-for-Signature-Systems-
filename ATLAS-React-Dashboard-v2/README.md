# ATLAS React Dashboard v2

A scalable React + FastAPI research interface for the ATLAS trastuzumab-resistance pipeline.

## What changed in v2

- Built-in **Documentation** page with pipeline, architecture, evidence roles, API reference, interpretation guardrails, operations, troubleshooting, and development notes.
- Functional **Settings** page backed by a persistent server-side JSON configuration.
- **Research Statistics** page with dynamic plots generated from current ATLAS outputs.
- Optional **Developer Mode** with runtime, filesystem, queue-service, queue-state, artifact registry, and API diagnostics.
- Settings now control `ATLAS_ROOT`, auto-refresh, refresh interval, table/chart display limits, developer mode, guardrail visibility, and dense tables.
- Existing ATLAS `.env` one directory above the web project is automatically loaded as a fallback. Secrets are never exposed to React.

## Architecture

```text
React + TypeScript + Vite
          |
          | /api/*
          v
FastAPI dashboard adapter
          |
          v
ATLAS_ROOT
├── data/
├── results/
├── scripts/
└── results/pipeline_state/
```

The UI is intentionally read-focused. It can inspect outputs and diagnostics but does not expose arbitrary shell execution or silently run expensive pipeline stages.

## Quick start

Place/extract the folder inside the existing ATLAS project, for example:

```text
/home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2
```

The backend automatically tries both:

```text
ATLAS-React-Dashboard-v2/.env
../.env
```

So if your main project already has `/home/regulus/Documents/ATLAS/.env`, you do not need to duplicate secrets.

### Terminal 1 — backend

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2
./run_backend.sh
```

Backend:

```text
http://127.0.0.1:8000
```

Health check:

```text
http://127.0.0.1:8000/api/health
```

### Terminal 2 — frontend

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2
./run_frontend.sh
```

Frontend:

```text
http://localhost:5173
```

## Main pages

### Dashboard

Live overview of current evidence outputs, dataset funnel, candidate table, recent artifact updates, and primary-validation cohorts.

### Research Statistics

The backend derives statistics from current files rather than using demo values. The page can display:

- Differential-expression volcano plot
- Exact full-file DEG summary counts
- Dataset evidence-role composition
- Pathway/GSEA NES chart when a compatible file exists
- Integrated candidate score chart when a numeric integrated score is present
- TGF-beta ranked validation chart when a compatible table exists

The chart layer uses graceful unavailable states. If the required output does not exist or its schema is incompatible, the site says so instead of inventing data.

### Documentation

The documentation is inside the app so research users can inspect methodology while looking at results. It covers:

- scientific scope
- software/data architecture
- source databases and evidence roles
- pipeline stage map
- chart definitions
- configuration
- scientific interpretation guardrails
- operations/runbook
- API reference
- developer extension guide
- troubleshooting
- wet-lab handoff principles

### Settings

Saved in:

```text
ATLAS-React-Dashboard-v2/.atlas-dashboard-settings.json
```

This file contains only dashboard preferences, not credentials.

Available settings:

| Setting | Purpose |
|---|---|
| `ATLAS_ROOT` | Select the research project directory. Must exist. |
| Auto refresh | Enable/disable timed dashboard refresh. |
| Refresh seconds | 5–3600 second refresh interval. |
| Table/chart limit | 100–10000 display limit. |
| Developer mode | Unlock read-only operational diagnostics. |
| Scientific guardrails | Keep interpretation warnings visible. |
| Dense tables | Compact large result tables. |

### Developer Mode

Enable it in Settings and save. A **Developer** navigation item appears.

It displays:

- Python/runtime information
- dashboard process PID
- filesystem total/used/free space
- `atlas-dataset-queue.service` state
- dataset queue row/status summary
- important output-file existence, size, and modification times
- API endpoint registry

Developer mode remains read-only.

## API

```text
GET /api/health
GET /api/dashboard
GET /api/datasets
GET /api/signature
GET /api/candidates
GET /api/cmap
GET /api/docking
GET /api/research/statistics
GET /api/settings
PUT /api/settings
GET /api/developer/statistics
```

`GET /api/developer/statistics` returns HTTP 403 unless Developer Mode is enabled.

## Scientific interpretation rules encoded into the UI

- Negative CMap tau/connectivity = transcriptional opposition, not proof of resistance reversal or efficacy.
- Docking = structural-computational evidence, not proof of binding or efficacy.
- PAINS flags = possible assay interference, not toxicity.
- Network/target support = plausibility, not causal proof.
- Regulatory absence is missing evidence and should not automatically become a penalty.
- Mirrors and umbrella datasets should not be double-counted as independent validation.
- Patient datasets provide translational context unless their design directly supports the resistance contrast being claimed.

## Production build

```bash
cd frontend
npm install
npm run build
```

Static assets are written to:

```text
frontend/dist/
```

For deployment, serve the built frontend behind Nginx/Caddy and run FastAPI with a production ASGI process. Keep write/execution privileges separated from the public web process.

## Adding a new research visualization

1. Compute/normalize the required data in the FastAPI layer.
2. Add a typed response shape in `frontend/src/types.ts`.
3. Add an API function in `frontend/src/lib/api.ts`.
4. Build a reusable component under `frontend/src/components/statistics/`.
5. Document the statistic, thresholds, and interpretation limits in the Documentation page.
6. Prefer full-data calculations on the backend and performance-limited rendering in the browser.

## Important design rule

Do not make a visual score look more certain than the underlying science. Keep acquisition quality, validation, perturbational evidence, target evidence, structural evidence, safety/regulatory context, and experimental evidence as distinguishable layers.
