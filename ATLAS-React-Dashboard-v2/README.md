# ATLAS React Dashboard v2

A React + TypeScript frontend and FastAPI backend for exploring local ATLAS research outputs. It includes candidate evidence tables, dataset eligibility, research charts, settings, documentation and optional developer diagnostics.

For research setup, expected input data, pipeline stages, automation and interpretation, see the [main project README](../README.md).

## First-time setup

Commands assume this dashboard is inside `/home/regulus/Documents/ATLAS`; replace the path for your checkout. You need Python with `venv`, Node.js and npm. The local backend was checked with Python 3.12.13; use Node/npm compatible with `frontend/package-lock.json`.

Install the backend into its own environment:

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2
python3 -m venv backend/.venv
backend/.venv/bin/python -m pip install -r backend/requirements.txt
```

Create this directory's `.env` if absent, or edit it without replacing existing settings:

```dotenv
ATLAS_ROOT=/home/regulus/Documents/ATLAS
```

The root contains `data/` and `results/`. Research API keys belong in the parent project's `.env`, not frontend variables.

Install the frontend from its lockfile:

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2/frontend
npm ci
```

## Start development servers

Both servers must run. Start the backend first and keep each terminal open.

### Terminal 1: backend

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2/backend
.venv/bin/python -m uvicorn main:app --host 127.0.0.1 --port 8000 --reload
```

Wait for `Application startup complete.`

### Terminal 2: frontend

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2/frontend
npm run dev
```

Open **http://localhost:5173/** or the URL printed by Vite.

Alternative helpers, run from this dashboard directory in separate terminals, are `bash run_backend.sh` and `bash run_frontend.sh`. The backend helper installs requirements on every run. The frontend helper installs packages only if `node_modules` is absent. Each starts only its own server.

Use `Ctrl+C` in each server terminal to stop it.

## Verify connectivity

```bash
curl --fail http://127.0.0.1:8000/api/health
curl --fail http://127.0.0.1:8000/api/settings
curl --fail http://localhost:5173/api/health
```

The final request checks Vite's proxy. Interactive API documentation is at **http://127.0.0.1:8000/docs**.

If Vite logs `ECONNREFUSED 127.0.0.1:8000`, it cannot reach the backend. Start Terminal 1's command and check its traceback if it exits. Vite startup does not indicate that the API is running. If the API uses another port, update `frontend/vite.config.ts` to match and restart Vite.

## Configuration

| Configuration | Meaning |
|---|---|
| Backend `ATLAS_ROOT` | Initial research data directory. |
| Backend `FRONTEND_ORIGIN` | Extra allowed browser origin for direct cross-origin API requests. Defaults already allow localhost and 127.0.0.1 on port 5173 over HTTP. |
| `frontend/.env.local`: `VITE_API_URL` | Optional API origin without `/api`. Omit for relative requests through Vite. Restart/rebuild after changes. |
| `.atlas-dashboard-settings.json` | Saved UI settings. Its root takes precedence over the initial environment root once the file exists. |

Environment precedence is process variables, this directory's `.env`, then the parent project's `.env`. Existing values are not overwritten. The example `ATLAS_API_HOST` and `ATLAS_API_PORT` variables are currently unused by the launchers; use Uvicorn arguments and matching proxy configuration.

Settings include root, auto refresh, interval (5–3600 seconds), row limit (100–10000), developer mode, scientific guidance and dense tables. Auto refresh controls Dashboard and Research Statistics; generic tables fetch on navigation/query refresh. Developer diagnostics poll every 15 seconds while enabled.

## Pages and data

| Route | Page |
|---|---|
| `/` | Dashboard: metrics, funnel, top candidates and activity |
| `/statistics` | Research Statistics |
| `/datasets` | Dataset eligibility and independence |
| `/signature` | Signature Discovery |
| `/cmap` | CMap Results |
| `/docking` | Docking Results |
| `/candidates` | Final Candidates |
| `/documentation` | In-app documentation |
| `/settings` | Configuration |
| `/developer` | Developer diagnostics; requires developer mode |

Final Candidates reads `results/cmap/integrated_evidence/ATLAS_integrated_evidence_matrix.csv` under the configured root. Top Final Candidates displays its first eight rows in the same order, using final evidence fields. The earlier `results/cmap/final_prioritization/ATLAS_final_candidate_prioritization.csv` is not this page's source.

Datasets prefer the independence-scored CSV in `data/enriched/`; signatures prefer annotated discovery DEGs in `results/differential_expression/`. CMap and docking views can use filename-based source discovery. Inspect displayed source paths when tracing evidence. Generic tables display the first 14 columns and an adjustable row limit; open the source CSV for all columns.

Missing or unreadable files can produce empty views. Starting the app does not generate scientific data. Health `ok: true` means the configured root exists, not that every stage completed.

### Docking score display

The dashboard table labels docking scores in **kcal/mol**. The API's `docking_score` comes from the integrated matrix's `best_affinity_kcal_mol`, not the pose count (`vina_mode_n`) or reference affinity. A dash means the candidate has no score; zero remains a displayed value.

In the artifacts checked on September 6, 2026, sitagliptin has **-9.327 kcal/mol** and clofibrate has **-6.528 kcal/mol**. These values come from the files and will change when evidence is regenerated.

Check the dashboard payload through Vite:

```bash
curl --fail --silent http://localhost:5173/api/dashboard | python3 -c 'import json, sys; rows = json.load(sys.stdin)["top_candidates"]; print(*[(r["name"], r["docking_score"]) for r in rows], sep="\n")'
```

New docking output must be incorporated into the integrated matrix before the dashboard can show it. See the [docking data workflow](../README.md#docking-scores-on-the-dashboard). Restart the backend after code changes if it was started without `--reload`.

## API

| Method | Endpoint |
|---|---|
| GET | `/api/health` |
| GET | `/api/dashboard` |
| GET | `/api/datasets` |
| GET | `/api/signature` |
| GET | `/api/candidates` |
| GET | `/api/cmap` |
| GET | `/api/docking` |
| GET / PUT | `/api/settings` |
| GET | `/api/research/statistics` |
| GET | `/api/developer/statistics` |

Table responses contain `rows`, `columns` and `source_file`. Developer statistics return HTTP 403 until developer mode is enabled. Settings updates validate the root and field values. See `/docs` for schemas.

## Build and checks

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2/frontend
npm run build
```

This performs TypeScript checks and writes `frontend/dist/`. `npm run preview` serves a local build preview; keep the API running and verify `/api/health` through that server. There is no frontend `npm test` script.

Run backend regression checks from this dashboard directory:

```bash
backend/.venv/bin/python -m unittest discover -s backend -p 'test_*.py' -v
```

For deployment, serve `dist/` with an `index.html` fallback for browser routes and proxy `/api/` to the backend, preserving the path. Run Uvicorn without `--reload` under a process supervisor. The API has no authentication; shared deployments need access controls outside this application. See the [deployment guide](../README.md#build-and-deployment).

## Troubleshooting

| Problem | Action |
|---|---|
| Connection refused on port 8000 | Start the backend, confirm startup completes, then retry `/api/health`. |
| Cannot import `main` | Start Uvicorn from the `backend` directory. |
| Missing Python package | Install requirements using `backend/.venv/bin/python -m pip`. |
| Port occupied | Inspect the existing server and reuse it, or stop it in its terminal before restarting. |
| Empty results | Check Settings and `/api/health`, then verify source CSVs exist and are readable. |
| Missing docking score | Inspect `best_affinity_kcal_mol` in the integrated matrix and regenerate downstream evidence if docking results were updated. |
| Pose count appears as score | Restart with the corrected backend mapping, then verify `/api/dashboard` returns the candidate affinity. |
| Root environment change ignored | Update the saved root through Settings. |
| CORS failure | Check `VITE_API_URL` and the allowed frontend origin; default development uses Vite's proxy. |
| HTTP 500 | Inspect the backend traceback. |
| Developer HTTP 403 | Enable developer mode in Settings. |

## Scientific scope

CMap, target/network evidence, docking and structural developability support experimental prioritization. They do not establish resistance reversal or clinical efficacy. The UI preserves unavailable evidence and exposes source paths for checking results against pipeline artifacts.
