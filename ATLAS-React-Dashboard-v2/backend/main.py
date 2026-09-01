from __future__ import annotations

from pathlib import Path
import os

from dotenv import load_dotenv
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

APP_ROOT = Path(__file__).resolve().parents[1]
# Prefer a dashboard-local .env, then the existing ATLAS project .env one level above.
load_dotenv(APP_ROOT / ".env", override=False)
load_dotenv(APP_ROOT.parent / ".env", override=False)

from atlas_reader import AtlasReader  # noqa: E402
from developer_stats import build_developer_statistics  # noqa: E402
from research_stats import build_research_statistics  # noqa: E402
from settings_store import SettingsStore  # noqa: E402

app = FastAPI(
    title="ATLAS Dashboard API",
    version="2.0.0",
    description="Read-focused API adapter and research-statistics layer over current ATLAS pipeline outputs.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",
        "http://127.0.0.1:5173",
    ],
    allow_credentials=True,
    allow_methods=["GET", "PUT", "OPTIONS"],
    allow_headers=["*"],
)

store = SettingsStore(APP_ROOT)


def current_settings():
    data = store.load()
    # Existing ATLAS_ROOT from .env is respected until the dashboard saves an explicit value.
    if not store.path.exists() and os.getenv("ATLAS_ROOT"):
        data["atlas_root"] = os.environ["ATLAS_ROOT"]
    return data


def reader():
    return AtlasReader(root=Path(current_settings()["atlas_root"]).expanduser())


class DashboardSettings(BaseModel):
    atlas_root: str = Field(min_length=1)
    auto_refresh: bool = True
    refresh_seconds: int = Field(default=30, ge=5, le=3600)
    table_row_limit: int = Field(default=1000, ge=100, le=10000)
    developer_mode: bool = False
    show_scientific_guardrails: bool = True
    dense_tables: bool = False


@app.get("/api/health")
def health():
    r = reader()
    return {
        "ok": r.root.exists() and r.root.is_dir(),
        "atlas_root": str(r.root),
        "api_version": "2.0.0",
        "settings_file": str(store.path),
    }


@app.get("/api/settings")
def get_settings():
    return current_settings()


@app.put("/api/settings")
def put_settings(settings: DashboardSettings):
    root = Path(settings.atlas_root).expanduser()
    if not root.exists() or not root.is_dir():
        raise HTTPException(status_code=400, detail=f"ATLAS_ROOT does not exist or is not a directory: {root}")
    saved = store.save(settings.model_dump())
    return {"ok": True, "settings": saved}


@app.get("/api/dashboard")
def dashboard():
    return reader().dashboard()


@app.get("/api/datasets")
def datasets():
    return reader().dataset_rows()


@app.get("/api/signature")
def signature():
    return reader().signature()


@app.get("/api/candidates")
def candidates():
    return reader().candidates()


@app.get("/api/cmap")
def cmap():
    return reader().cmap()


@app.get("/api/docking")
def docking():
    return reader().docking()


@app.get("/api/research/statistics")
def research_statistics():
    settings = current_settings()
    point_limit = min(3000, max(500, int(settings.get("table_row_limit", 1000)) * 2))
    return build_research_statistics(reader(), point_limit=point_limit)


@app.get("/api/developer/statistics")
def developer_statistics():
    settings = current_settings()
    if not settings.get("developer_mode", False):
        raise HTTPException(status_code=403, detail="Developer mode is disabled. Enable it in Settings first.")
    return build_developer_statistics(reader().root)
