# ATLAS React Dashboard

A scalable React + FastAPI dashboard for the ATLAS trastuzumab-resistance pipeline.

## Included
- React + TypeScript + Vite
- Tailwind CSS theme based on the supplied ATLAS mockup
- React Router for scalable page routing
- TanStack Query for live server-state refresh
- Reusable dashboard components
- FastAPI adapter that reads your existing ATLAS outputs dynamically
- Responsive sidebar/mobile layout and dark mode
- Graceful empty states when a pipeline file is missing

## Configure

```bash
cp .env.example .env
```

Default ATLAS project path:

```text
/home/regulus/Documents/ATLAS
```

## Start backend

```bash
cd backend
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn main:app --reload --port 8000
```

## Start frontend

In another terminal:

```bash
cd frontend
npm install
npm run dev
```

Normally open `http://localhost:5173`.

## Production build

```bash
cd frontend
npm run build
```

The production frontend is written to `frontend/dist/`.

## API

- `GET /api/health`
- `GET /api/dashboard`
- `GET /api/datasets`
- `GET /api/signature`
- `GET /api/candidates`
- `GET /api/cmap`
- `GET /api/docking`

The backend reads current ATLAS CSV/JSON outputs on request, so the UI follows your pipeline rather than depending on hardcoded dashboard numbers.
