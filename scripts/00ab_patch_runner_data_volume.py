#!/usr/bin/env python3
from __future__ import annotations

from datetime import datetime
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "scripts" / "run_atlas_full_auto.py"
MARKER = "ATLAS_DATA_VOLUME_STAGE_V1"


def main():
    if not TARGET.exists():
        raise SystemExit(f"ERROR: missing {TARGET}")

    text = TARGET.read_text(encoding="utf-8")

    if MARKER in text:
        print("Dataset-volume stage already installed; no changes made.")
        return 0

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = TARGET.with_name(
        f"{TARGET.stem}_pre_data_volume_{stamp}{TARGET.suffix}"
    )
    shutil.copy2(TARGET, backup)

    anchor = (
        '    Stage(\n'
        '        "00f", "00f_primary_validation_matrix_inspection.py", "validation",\n'
    )

    if anchor not in text:
        raise RuntimeError("Could not locate 00F stage anchor.")

    stage = (
        '    # ATLAS_DATA_VOLUME_STAGE_V1\n'
        '    Stage(\n'
        '        "00aa", "00aa_dataset_volume_metadata.py", "metadata",\n'
        '        "Record dataset counts and local/remote storage volume metadata",\n'
        '        "results/pipeline_state/dataset_volume_summary.json",\n'
        '    ),\n\n'
    )

    text = text.replace(anchor, stage + anchor, 1)
    compile(text, str(TARGET), "exec")
    TARGET.write_text(text, encoding="utf-8")

    print("=" * 78)
    print("ATLAS — DATASET VOLUME STAGE INSTALLED")
    print("=" * 78)
    print(f"Patched: {TARGET}")
    print(f"Backup:  {backup}")
    print("Inserted: 00E -> 00AA volume metadata -> 00F")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
