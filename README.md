# ATLAS

ATLAS is a computational research project for studying trastuzumab resistance in HER2-positive breast cancer and prioritizing compounds for experimental follow-up. It combines dataset discovery, phenotype and study-independence checks, differential expression, pathway validation, Connectivity Map (CMap) screening, target/network evidence, docking, and structural developability assessment.

The React dashboard displays outputs produced by the Python research pipeline. Its FastAPI backend reads local files under a configured project directory. Starting the dashboard does not start the research pipeline or download research datasets.

## Contents

- [Architecture and repository layout](#architecture-and-repository-layout)
- [Requirements](#requirements)
- [Dashboard setup](#dashboard-setup)
- [Start the backend and frontend](#start-the-backend-and-frontend)
- [Configuration](#configuration)
- [Dashboard pages and data sources](#dashboard-pages-and-data-sources)
- [API reference](#api-reference)
- [Research pipeline setup](#research-pipeline-setup)
- [Run and resume the pipeline](#run-and-resume-the-pipeline)
- [Publication data preparation](#publication-data-preparation)
- [Continuous processing](#continuous-processing)
- [Build and deployment](#build-and-deployment)
- [Development and verification](#development-and-verification)
- [Troubleshooting](#troubleshooting)
- [Data management and interpretation](#data-management-and-interpretation)

## Architecture and repository layout

```text
Browser: React + TypeScript + Tailwind CSS + Recharts
    |
    | /api/* through Vite during development
    v
FastAPI + pandas: ATLAS-React-Dashboard-v2/backend
    |
    | reads files under ATLAS_ROOT
    v
data/ and results/ <--- Python research scripts
```

```text
ATLAS/
├── README.md                         Project setup and operating guide
├── SOP.md                            Research questions
├── requirements.txt                  Base research dependencies
├── .env.example                      CMap credential template
├── atlas_data/                       Shared acquisition/catalog helpers
├── config/dataset_queries.json       Dataset discovery queries
├── scripts/                          Research stages and automation runners
├── data/
│   ├── raw/                          Original inputs, including tx2gene.tsv
│   ├── processed/                    Gene counts and processed inputs
│   ├── discovery/                    Discovered/scored dataset candidates
│   ├── catalog/                      DuckDB metadata catalog
│   ├── enriched/                     Sample metadata and eligibility evidence
│   ├── manifests/                    Acquisition provenance
│   └── validation_expression/        External validation matrices
├── results/
│   ├── qc/                           Discovery quality-control outputs
│   ├── differential_expression/      Discovery DEG tables
│   ├── pathway_analysis/             Discovery pathway outputs
│   ├── external_validation/          Cross-study and mechanism evidence
│   ├── cmap/                         CMap through final integrated evidence
│   └── pipeline_state/               Checkpoints, preflight and queue state
├── logs/                             Pipeline and monitoring logs
└── ATLAS-React-Dashboard-v2/
    ├── README.md                     Dashboard-specific operating guide
    ├── .env.example                  Backend configuration template
    ├── run_backend.sh                Backend development launcher
    ├── run_frontend.sh               Frontend development launcher
    ├── backend/
    │   ├── main.py                   FastAPI routes and configuration
    │   ├── atlas_reader.py           Pipeline file adapter
    │   ├── research_stats.py         Research chart calculations
    │   ├── developer_stats.py        Runtime and queue diagnostics
    │   ├── settings_store.py         Persistent dashboard settings
    │   └── requirements.txt          Dashboard backend dependencies
    └── frontend/
        ├── package.json
        ├── package-lock.json
        ├── vite.config.ts           Development API proxy
        └── src/                     Routes, pages, components and API client
```

Some output directories are created by pipeline stages and may be absent in a fresh checkout. Raw and processed data are excluded by `.gitignore`; a checkout alone is not a complete research dataset.

## Requirements

| Component | Requirement |
|---|---|
| Operating environment | Bash commands below target Linux. Queue/monitor tooling uses Linux facilities such as `fcntl` and user systemd services. |
| Python | Python 3.10 or newer is required by the runtime preflight. The local dashboard environment was checked with Python 3.12.13. Scientific package compatibility also depends on platform and package versions. |
| Frontend | Node.js and npm compatible with the checked-in Vite dependency tree. Local tools were Node.js 26.1.0 and npm 11.14.1 when this guide was written; these are observed versions, not pinned requirements. |
| Storage | Space for the intended research inputs and outputs. Preflight enforces 3 GiB free, a startup floor rather than a total dataset storage estimate. |
| Network | Needed for dependency installation and online research stages. Existing local outputs can be viewed without rerunning those stages. |
| Credentials | A valid `CLUE_API_KEY` for CMap submission. Viewing existing results does not require a key. |

The frontend uses React 19, TypeScript, Vite, React Router, TanStack Query, Tailwind CSS and Recharts. Backend requirements include FastAPI, Uvicorn, pandas, NumPy, DuckDB, python-dotenv, psutil and Pydantic.

## Dashboard setup

Examples use `/home/regulus/Documents/ATLAS`. Replace that path if your checkout is elsewhere.

### 1. Install backend dependencies

```bash
cd /home/regulus/Documents/ATLAS
python3 -m venv ATLAS-React-Dashboard-v2/backend/.venv
ATLAS-React-Dashboard-v2/backend/.venv/bin/python -m pip install -r ATLAS-React-Dashboard-v2/backend/requirements.txt
```

If this environment already exists, reuse it and run only the installation command. It is separate from the root `.venv` used for research scripts.

### 2. Configure the data location

Create `ATLAS-React-Dashboard-v2/.env` if absent, or edit the existing file while preserving its other settings:

```dotenv
ATLAS_ROOT=/home/regulus/Documents/ATLAS
```

This must point to the ATLAS directory containing `data/` and `results/`, not the dashboard subdirectory. See [Configuration](#configuration) for precedence and saved settings.

### 3. Install frontend dependencies

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2/frontend
npm ci
```

`npm ci` installs from `package-lock.json` and replaces an existing `node_modules` directory. Use `npm install` when intentionally changing dependencies and retain the corresponding lockfile changes.

## Start the backend and frontend

Keep both processes running in separate terminals. Run the backend first.

### Terminal 1: backend

After installing dependencies:

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2/backend
.venv/bin/python -m uvicorn main:app --host 127.0.0.1 --port 8000 --reload
```

Wait for `Application startup complete.` The backend listens at **http://127.0.0.1:8000**.

Alternatively, the bundled helper creates the environment if missing, installs requirements on every invocation, and starts Uvicorn with reload on port 8000:

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2
bash run_backend.sh
```

The helper uses `python`, `pip` and `uvicorn` from its shell/activated environment. The explicit `.venv/bin/python -m uvicorn` command selects the interpreter directly and skips repeated installation.

### Terminal 2: frontend

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2/frontend
npm run dev
```

Open the URL printed by Vite, normally **http://localhost:5173/**.

The frontend helper is an alternative:

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2
bash run_frontend.sh
```

It runs `npm install` only when `node_modules` is absent, then starts Vite. It starts only the frontend.

### Check the connection and stop servers

From another terminal:

```bash
curl --fail http://127.0.0.1:8000/api/health
curl --fail http://127.0.0.1:8000/api/settings
curl --fail http://localhost:5173/api/health
```

The first request checks the backend directly; the last checks the Vite proxy. Health includes `ok`, `atlas_root`, `api_version` and `settings_file`. `ok: true` confirms that the configured root exists; it does not certify that all research outputs exist or are valid.

Stop each development server with `Ctrl+C` in its terminal. Starting one server does not start the other.

## Configuration

### Backend environment

The backend loads the dashboard `.env` first, then the project-root `.env`, without overwriting existing process variables. Precedence is therefore process environment, dashboard `.env`, then root `.env`.

| Variable | Purpose |
|---|---|
| `ATLAS_ROOT` | Initial data directory, used before a settings file has been saved. |
| `FRONTEND_ORIGIN` | Additional browser origin allowed by CORS, for example `http://localhost:4173`. Default allowed origins are `http://localhost:5173` and `http://127.0.0.1:5173`. |
| `CLUE_API_KEY` | Research-stage CMap credential; keep in the root `.env` or process environment. Never place it in frontend configuration. |

The dashboard `.env.example` also lists `ATLAS_API_HOST` and `ATLAS_API_PORT`, but the current launcher and Vite configuration do **not** read them. To change ports, set Uvicorn's `--port` and update the proxy target in `frontend/vite.config.ts`, then restart the servers.

### Frontend environment

By default the API client sends relative `/api/*` requests. Vite forwards them to `http://127.0.0.1:8000` during development.

For a separate API origin, create `ATLAS-React-Dashboard-v2/frontend/.env.local`:

```dotenv
VITE_API_URL=http://127.0.0.1:8000
```

Use the origin without `/api`; the client appends endpoint paths. Set backend `FRONTEND_ORIGIN` if the browser origin differs from the defaults. Restart Vite after environment changes and rebuild when changing a production API URL. `VITE_*` values are exposed to the browser and must not contain secrets.

### Settings page

Saving settings writes `ATLAS-React-Dashboard-v2/.atlas-dashboard-settings.json`. Once this file exists, its saved root takes precedence over the initial `ATLAS_ROOT` environment setting. Use **Settings** when switching data directories.

| Setting | Default | Behavior |
|---|---|---|
| Project root | `/home/regulus/Documents/ATLAS` | Must be an existing directory when saved. |
| Auto refresh | Enabled | Controls timed refresh on Dashboard and Research Statistics. |
| Refresh interval | 30 seconds | Allowed range: 5–3600 seconds. |
| Table row limit | 1000 | Allowed range: 100–10000; applied by table views. |
| Developer mode | Disabled | Enables developer diagnostics access. |
| Scientific guardrails | Enabled | Controls interpretation guidance. |
| Dense tables | Disabled | Uses more compact table spacing. |

Generic result tables fetch on navigation/query refresh; they do not use the Dashboard polling interval. Developer statistics poll every 15 seconds while enabled.

## Dashboard pages and data sources

| Route | Page | Purpose |
|---|---|---|
| `/` | Dashboard | Metrics, dataset funnel, top eight candidates, output updates and primary-validation cohorts. |
| `/statistics` | Research Statistics | DEG, dataset, pathway, candidate and TGF-beta visualizations. |
| `/datasets` | Datasets | Eligibility, modality, phenotype confidence and independence. |
| `/signature` | Signature Discovery | Discovery differential-expression table. |
| `/cmap` | CMap Results | Available integrated/CMap evidence. |
| `/docking` | Docking Results | Available structural evidence table. |
| `/candidates` | Final Candidates | Integrated evidence matrix in its stored row order. |
| `/documentation` | Documentation | In-app methodology and operating notes. |
| `/settings` | Settings | Data root and display preferences. |
| `/developer` | Developer | Runtime, services, queue and artifact diagnostics; requires developer mode. |

Core lookup paths are relative to `ATLAS_ROOT`:

| View | Source selection |
|---|---|
| Datasets | First available: `data/enriched/dataset_candidates_independence_scored.csv`, `data/discovery/dataset_candidates_scored.csv`, then `data/discovery/dataset_candidates.csv`. |
| Signature | `results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv`, falling back to `DEGs_resistant_vs_sensitive.csv` in the same directory. |
| Final Candidates | `results/cmap/integrated_evidence/ATLAS_integrated_evidence_matrix.csv`. |
| Top Final Candidates | First eight rows of that same matrix, preserving the Final Candidates page order. |
| CMap | Integrated matrix when present; otherwise a CSV selected by filename keywords. |
| Docking | CSV selected by docking-related filename keywords under `results/` or `data/enriched/`. |

Generic table API responses include `source_file`, `columns` and `rows`, with up to 10000 rows. The UI currently shows at most the first 14 columns; inspect the source CSV for the complete evidence table. Where filename discovery is used, inspect the displayed source path to confirm which artifact was selected.

The `04r` final-prioritization CSV is an earlier artifact. Final Candidates uses the later `04u` integrated matrix, which incorporates subsequent evidence. The dashboard preserves that final ordering and maps compound names, mean tau, validated target symbols, docking affinities, experimental scores and final evidence categories.

Missing outputs normally produce empty views or unavailable charts. The file adapter also returns an empty frame for CSV read failures, so verify readability if an expected table is empty.

### Docking scores on the dashboard

The **Top Final Candidates** table includes **Docking score (kcal/mol)**. The backend maps the candidate's `best_affinity_kcal_mol` from the integrated matrix to the API's `docking_score`. It does not use the number of poses (`vina_mode_n`) or the reference ligand's affinity. Missing values remain a dash (—), and a score of zero is retained.

For the files checked on September 6, 2026, the first two candidates are:

| Compound | Validated target | Docking score (kcal/mol) |
|---|---|---|
| sitagliptin | DPP4 | -9.327 |
| clofibrate | CYP3A4 | -6.528 |

These are a snapshot of the current artifacts, not values hard-coded into the application. Candidates without a docking result still appear in final order with a dash in the score column.

Verify the values delivered to the browser while both servers are running:

```bash
curl --fail --silent http://localhost:5173/api/dashboard | python3 -c 'import json, sys; rows = json.load(sys.stdin)["top_candidates"]; print(*[(r["name"], r["docking_score"]) for r in rows], sep="\n")'
```

Docking results originate in `results/cmap/docking/ATLAS_docking_results.csv` and enter the dashboard through the later integrated matrix. If docking outputs have changed, update the downstream structural assessment and integrated evidence with their prerequisites in place:

```bash
cd /home/regulus/Documents/ATLAS
source .venv/bin/activate
python -u scripts/run_atlas_full_auto.py --from-stage 04t --to-stage 04u --force
```

This recomputes downstream outputs; it does not perform new docking. To intentionally rerun docking as well, start the selected range at `04s`. Refresh the browser after outputs update. With backend `--reload`, Python code edits reload automatically; otherwise restart the backend after changing its field mappings.

## API reference

Interactive documentation: **http://127.0.0.1:8000/docs**. OpenAPI schema: `/openapi.json`.

| Method | Endpoint | Returns / action |
|---|---|---|
| GET | `/api/health` | Root availability, API version and settings file location. |
| GET | `/api/dashboard` | Project metadata, metrics, funnel, candidates, activity, primary validation and warnings. |
| GET | `/api/datasets` | Normalized dataset rows, columns and source path. |
| GET | `/api/signature` | Discovery DEG table. |
| GET | `/api/candidates` | Final integrated candidate table. |
| GET | `/api/cmap` | CMap/integrated evidence table. |
| GET | `/api/docking` | Selected docking evidence table. |
| GET | `/api/settings` | Current settings. |
| PUT | `/api/settings` | Validates and saves settings. Invalid directories return HTTP 400; invalid field values can return HTTP 422. |
| GET | `/api/research/statistics` | Chart data and availability information. |
| GET | `/api/developer/statistics` | Developer diagnostics; HTTP 403 while developer mode is disabled. |

Use the Settings page or Swagger UI to review the schema before submitting an update. The API does not provide endpoints to launch research jobs.

## Research pipeline setup

This setup generates or refreshes research results. Existing output files are enough to use the dashboard.

### Python environment and additional packages

Create or reuse the separate research environment:

```bash
cd /home/regulus/Documents/ATLAS
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The root requirements file is a **base dependency list**, not a complete or pinned environment for every stage. Implemented acquisition, annotation, pathway and docking stages also use these packages:

```bash
python -m pip install requests scipy duckdb python-dotenv pyarrow GEOparse h5py xlrd mygene gseapy
python -m pip install rdkit vina meeko gemmi openbabel-wheel
```

Some chemistry packages require platform-specific wheels or native build dependencies. Docking uses the Vina Python module, Open Babel Python bindings, and Meeko's `mk_prepare_ligand.py` / `mk_prepare_receptor.py` commands. Verify these in the activated research environment:

```bash
python -c 'import pandas, numpy, scipy, requests, duckdb, pydeseq2, gseapy, rdkit, vina; from openbabel import openbabel; print("Core research imports OK")'
command -v mk_prepare_ligand.py
command -v mk_prepare_receptor.py
```

The optional publication script has additional requirements, including PLIP and an `obabel` executable for relevant operations. It is separate from dashboard startup and the main full-runner stage list.

### Credentials and preflight

Create or edit the root `.env` with your own CMap credential:

```dotenv
CLUE_API_KEY=replace_with_your_own_key
```

CMap scripts load this file. Runtime preflight checks the process environment directly, so export the key for that check. In Bash, for your own trusted shell-compatible `.env`:

```bash
set -a
source .env
set +a
python scripts/00y_runtime_preflight.py
```

Preflight checks Python, package imports, stage files, directory writes, free disk space, selected network hosts and discovered credential variables. It writes `results/pipeline_state/runtime_preflight.json` and `.csv`. Passing preflight does not replace study-design, input-file or remote API access validation.

### Discovery input data

`scripts/01_validation.py` currently defines six samples explicitly: TS1–TS3 (Sensitive) and TR1–TR3 (Resistant). Expected files under `data/raw/` are:

```text
GSM9067960_TS1_quant.sf.gz
GSM9067961_TS2_quant.sf.gz
GSM9067962_TS3_quant.sf.gz
GSM9067963_TR1_quant.sf.gz
GSM9067964_TR2_quant.sf.gz
GSM9067965_TR3_quant.sf.gz
tx2gene.tsv
```

`tx2gene.tsv` must be tab-separated and include `Gene stable ID`, `Transcript stable ID` and `Gene name`. Quantification files and transcript mappings must correspond to the intended experiment/reference. Changing studies requires reviewing sample and condition definitions in the discovery scripts; it is not just a file rename.

With those inputs and dependencies prepared:

```bash
python -u scripts/01_validation.py
python -u scripts/02_differential_expression.py
python -u scripts/03_pathway_analysis.py
python -u scripts/04_cmap_analysis.py
```

Stage 01 builds counts and QC; 02 produces DEGs; 03 performs pathway analysis; 04 constructs CMap signatures. Confirm required outputs before proceeding. The full automation runner does not include these four discovery preparation scripts, so they or their existing outputs must be prepared separately.

## Run and resume the pipeline

Run research commands from the ATLAS root with its `.venv` activated. Inspect the authoritative stage list:

```bash
python scripts/run_atlas_full_auto.py --help
python scripts/run_atlas_full_auto.py --list-stages
```

| Group | Work |
|---|---|
| Runtime preflight | Dependency, network, credential and storage checks. |
| Discovery/enrichment | GEO/SRA/GDC/cBioPortal and EBI discovery, cataloging and metadata enrichment. |
| Phenotype/independence | Conservative classification, HER2 confirmation, relationships and scientific gate. |
| External validation | Expression acquisition, matrix/design checks, differential expression and consensus. |
| Mechanism validation | Pathway and TGF-beta audits, evidence synthesis and exports. |
| CMap | Submit, poll, download, parse and prioritize perturbational results. |
| Compound evidence | Identity, regulatory/trial evidence, preliminary safety, targets and networks. |
| Final integration | `04r` prioritization, `04s` docking, `04t` structural assessment, `04u` integrated matrix. |

Start/resume the configured sequence:

```bash
python -u scripts/run_atlas_full_auto.py --full
```

Run a range or rerun one stage with prerequisite artifacts already available:

```bash
python -u scripts/run_atlas_full_auto.py --from-stage 04r --to-stage 04u
python -u scripts/run_atlas_full_auto.py --rerun-stage 04u
```

The full runner uses lowercase stage keys. `scripts/run_atlas_pipeline.py` is a separate downstream runner with `--from` / `--to` and uppercase keys; do not mix their CLI formats.

| Option | Effect |
|---|---|
| `--from-stage` / `--to-stage` | Select an inclusive range. |
| `--rerun-stage` | Force one stage. |
| `--force` | Force selected stages to execute. |
| `--refresh-data` | Bypass checkpoint skipping for the selected sequence; not restricted to download stages. |
| `--cmap-poll-minutes` | Poll interval, default 10 minutes. |
| `--cmap-max-wait-hours` | CMap wait limit, default 12 hours. |

The runner records status and stops on failed commands. It can skip successful unchanged checkpoints; an upstream rerun causes downstream work in the selected sequence to rerun. Review scope before force/refresh operations because stages can submit jobs, download data and replace outputs.

State and logs:

```text
results/pipeline_state/full_pipeline_state.json
results/pipeline_state/runtime_preflight.json
results/pipeline_state/scientific_automation_gate.json
logs/full_pipeline/
```

Metadata enrichment can also run directly after discovery/catalog/initial scoring:

```bash
python -u scripts/00d_metadata_enrichment.py --min-score 40 --sources GEO,SRA --max-datasets 30 --workers 6
python -u scripts/00c2_rescore_after_enrichment.py
```

This writes sample metadata, summaries, rescored candidates and download manifests under `data/enriched/`. Metadata enrichment prepares a download plan; it does not itself download raw FASTQ files.

## Publication data preparation

`scripts/05_publication_remaining_data.py` prepares GO enrichment, STRING interactions, drug–target–resistance network tables, and docking interaction artifacts. It is separate from the full runner and requires these existing inputs:

```text
results/external_validation/three_dataset_strict_core_genes.csv
results/cmap/drug_targets/ATLAS_CMap_drug_target_pairs.csv
results/cmap/network_integration/ATLAS_target_resistance_gene_links.csv
results/cmap/docking/ATLAS_docking_results.csv
```

Inspect the script's configured docking cases and required structural files before running it. Its enrichment/network requests need network access, and relevant structure conversion steps need the `obabel` executable. PLIP supplies detailed interaction analysis; without it, the script records that limitation rather than providing complete interaction results.

With prerequisites prepared, run from the root research environment:

```bash
cd /home/regulus/Documents/ATLAS
source .venv/bin/activate
python -u scripts/05_publication_remaining_data.py
```

The script executes at module level and has no argument parser: importing it or passing `--help` is not a read-only way to inspect it. Read the source instead.

Outputs are written under `results/publication_remaining_data/`:

| Directory/file | Contents |
|---|---|
| `go_enrichment/` | GO enrichment tables. |
| `string_ppi/` | STRING protein interaction nodes and edges. |
| `drug_target_resistance_network/` | Network exports linking compounds, targets and resistance genes. |
| `docking_interactions/` | Prepared complexes and available PLIP interaction outputs. |
| `publication_remaining_data_manifest.json` | Gene counts, analysis settings, prepared complexes and interaction counts. |

Review logs, the manifest and any `PLIP_NOT_INSTALLED.txt` notice before treating the export as complete. These publication files do not replace the integrated matrix used by the dashboard.

## Continuous processing

`scripts/run_atlas_dataset_queue.py` processes datasets with persistent state and a required deadline. Its `ROOT` is currently hard-coded to `/home/regulus/Documents/ATLAS`; review it before using another checkout. Monitor/setup shell scripts also contain machine-specific paths.

Inspect the CLI without starting a worker:

```bash
python scripts/run_atlas_dataset_queue.py --help
```

To start intentionally, supply a future deadline (replace the example text):

```bash
python -u scripts/run_atlas_dataset_queue.py --deadline 'YYYY-MM-DD HH:MM:SS'
```

The deadline defaults to Asia/Manila without an explicit timezone. A past deadline exits immediately. Queue/history/lock files live under `results/pipeline_state/dataset_queue/`. A dataset labeled primary validation must actually appear in the DE driver's outputs before being counted as contributing evidence.

For an already-installed user service:

```bash
systemctl --user status atlas-dataset-queue.service
journalctl --user -u atlas-dataset-queue.service -f
```

The older [continuous queue notes](README_ATLAS_CONTINUOUS_QUEUE.md) describe a historical overnight run and externally supplied service files. Their fixed deadline and Streamlit `ui/app.py` command are not the current dashboard workflow. That UI entry point is not present in this checkout. Review `setup_atlas_unattended.sh` before running it: it configures automation and is not needed to launch React/FastAPI.

## Build and deployment

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2/frontend
npm run build
```

This runs TypeScript build checks and Vite, producing `frontend/dist/`. `npm run preview` serves a local build preview. Keep the backend running for API requests and check `/api/health` through the preview server.

For deployment, serve `dist/` with a fallback to `index.html` for browser routes such as `/candidates`. Proxy `/api/` to the backend while preserving the `/api` path, or build with `VITE_API_URL` and configure `FRONTEND_ORIGIN`.

Run the backend without development reload:

```bash
cd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2/backend
.venv/bin/python -m uvicorn main:app --host 127.0.0.1 --port 8000
```

Use a process supervisor for persistent deployment. The current API has no authentication layer; keep shared deployments behind appropriate access controls. CORS and developer mode are not authentication. The API exposes research paths/diagnostics and saves settings, so local defaults are not a public deployment configuration.

## Development and verification

Paths below are relative to `ATLAS-React-Dashboard-v2/`:

| Change | Main file |
|---|---|
| Frontend routes | `frontend/src/App.tsx` |
| Dashboard layout | `frontend/src/pages/DashboardPage.tsx` |
| Candidate preview | `frontend/src/components/dashboard/TopCandidatesTable.tsx` |
| Client API URLs | `frontend/src/lib/api.ts` |
| Response types | `frontend/src/types.ts` |
| API routes | `backend/main.py` |
| Data selection/mappings | `backend/atlas_reader.py` |
| Chart calculations | `backend/research_stats.py` |

Run candidate regression checks from the project root:

```bash
ATLAS-React-Dashboard-v2/backend/.venv/bin/python -m unittest discover -s ATLAS-React-Dashboard-v2/backend -p 'test_*.py' -v
```

These cover candidate ordering, evidence mappings and missing-data behavior. Use `npm run build` for frontend TypeScript/build validation. The frontend package has no `npm test` script. For an end-to-end local check, verify `/api/health`, `/api/settings` and `/api/dashboard` directly on port 8000 and through Vite on port 5173.

## Troubleshooting

| Symptom | Check and action |
|---|---|
| Vite: `ECONNREFUSED 127.0.0.1:8000` | Start the backend in Terminal 1, wait for startup completion and check `curl --fail http://127.0.0.1:8000/api/health`. Starting Vite alone is insufficient. |
| Backend exits immediately | Read its traceback. Install `backend/requirements.txt` in `backend/.venv`; start with `.venv/bin/python -m uvicorn`. |
| Cannot import module `main` | Start Uvicorn from `ATLAS-React-Dashboard-v2/backend`. |
| `ModuleNotFoundError` | Install packages into the interpreter running that component. Root `.venv` and dashboard `backend/.venv` are separate. |
| Port already in use | Inspect `ss -ltnp '( sport = :8000 or sport = :5173 )'`. Reuse the intended server or stop it in its terminal. Match backend port changes in the Vite proxy. |
| Vite chooses another port | Use the printed URL. With a direct API origin, allow that browser origin using `FRONTEND_ORIGIN`. |
| CORS error | Confirm the exact browser origin including scheme/port. Default development uses relative `/api` requests through Vite. |
| API HTTP 500 | Inspect the backend traceback; a running process can still fail on particular data. |
| Empty dashboard / missing candidates | Check Settings and `/api/health` for the root, then verify source files exist and are readable. Installation does not generate research results. |
| Changing `ATLAS_ROOT` has no effect | Saved `.atlas-dashboard-settings.json` takes precedence. Update the root in Settings. |
| Candidate pages appear different | Confirm both sessions use the same backend/root and refresh them. The dashboard should match the first eight integrated-matrix rows; generic tables do not poll continuously. |
| Docking shows a dash | Check the candidate's `best_affinity_kcal_mol` in the integrated matrix. If a new docking result exists only in the docking CSV, regenerate `04t`–`04u`; otherwise the score is unavailable. |
| Docking shows `9` instead of an affinity | Check that the running backend uses `best_affinity_kcal_mol`, not `vina_mode_n`. Restart it if running without reload, then inspect `/api/dashboard`. |
| Developer API HTTP 403 | Enable developer mode in Settings. |
| User service unavailable | The service may not be installed or user systemd may be unavailable; normal dashboard usage still works. |
| Pipeline preflight failure | Read `runtime_preflight.csv`/`.json` and resolve the reported check before resuming. |
| CMap jobs pending/failed | Inspect stage logs and remote job state. Submission success alone is not completed analysis. |
| Queue exits immediately | Check that `--deadline` is in the future and timezone/root are correct. |

## Data management and interpretation

Keep original inputs, sample provenance and reference versions with analysis records. Back up required `data/`, `results/`, `logs/` and configuration before replacing outputs. Python requirements are not locked, so record the research environment, for example with `python -m pip freeze > results/pipeline_state/research_environment.txt` in the activated research environment after the output directory exists. The frontend lockfile supports repeatable frontend dependency installation.

Evidence layers support experimental prioritization. Negative CMap connectivity indicates transcriptional opposition; docking is structural computational evidence; target/network links support biological plausibility. They do not establish resistance reversal, clinical efficacy or a validated treatment. Structural developability estimates are not a complete experimental ADMET assessment. PAINS alerts are potential assay-interference flags, and absent regulatory records must not be presented as verified regulatory conclusions.

Study mirrors and overlapping cohorts must not be counted as independent validation. Preserve scientific gate results and distinguish mechanistic questions in [SOP.md](SOP.md) from evidence established by a particular run.
