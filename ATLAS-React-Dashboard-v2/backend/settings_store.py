from __future__ import annotations

import json
from pathlib import Path
from typing import Any

DEFAULTS: dict[str, Any] = {
    "atlas_root": "/home/regulus/Documents/ATLAS",
    "auto_refresh": True,
    "refresh_seconds": 30,
    "table_row_limit": 1000,
    "developer_mode": False,
    "show_scientific_guardrails": True,
    "dense_tables": False,
}


class SettingsStore:
    def __init__(self, app_root: Path):
        self.path = app_root / ".atlas-dashboard-settings.json"

    def load(self) -> dict[str, Any]:
        data = dict(DEFAULTS)
        if self.path.exists():
            try:
                saved = json.loads(self.path.read_text(encoding="utf-8"))
                if isinstance(saved, dict):
                    data.update(saved)
            except Exception:
                pass
        return data

    def save(self, updates: dict[str, Any]) -> dict[str, Any]:
        data = self.load()
        data.update(updates)
        self.path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        return data
