from pathlib import Path
from dotenv import load_dotenv
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

load_dotenv(Path(__file__).resolve().parents[1] / ".env")
from atlas_reader import AtlasReader

app = FastAPI(title="ATLAS Dashboard API", version="1.0.0", description="Read-only API adapter over current ATLAS pipeline outputs.")
app.add_middleware(CORSMiddleware, allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"], allow_credentials=True, allow_methods=["GET"], allow_headers=["*"])
reader = AtlasReader.from_env()

@app.get("/api/health")
def health(): return {"ok": reader.root.exists(), "atlas_root": str(reader.root)}

@app.get("/api/dashboard")
def dashboard(): return reader.dashboard()

@app.get("/api/datasets")
def datasets(): return reader.dataset_rows()

@app.get("/api/signature")
def signature(): return reader.signature()

@app.get("/api/candidates")
def candidates(): return reader.candidates()

@app.get("/api/cmap")
def cmap(): return reader.cmap()

@app.get("/api/docking")
def docking(): return reader.docking()
