#!/usr/bin/env python3
from __future__ import annotations

from datetime import datetime
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "scripts" / "run_atlas_full_auto.py"
MARKER = "ATLAS_SCIENTIFIC_GATE_STAGE_V1"


def main():
    if not TARGET.exists():
        raise SystemExit(f"ERROR: missing {TARGET}")

    text = TARGET.read_text(encoding="utf-8")

    if MARKER in text:
        print("Full runner already contains the scientific gate; no changes made.")
        return 0

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = TARGET.with_name(
        f"{TARGET.stem}_pre_scientific_gate_{stamp}{TARGET.suffix}"
    )
    shutil.copy2(TARGET, backup)

    anchor = '''    Stage(
        "00e", "00e_primary_validation_expression_fetch.py", "validation",
'''

    if anchor not in text:
        raise RuntimeError("Could not find 00E stage anchor in full runner.")

    gate_stage = '''    # ATLAS_SCIENTIFIC_GATE_STAGE_V1
    Stage(
        "00w", "00w_scientific_automation_gate.py", "phenotype",
        "Scientific automation gate before expression validation",
        "results/pipeline_state/scientific_automation_gate.json",
    ),

'''

    text = text.replace(anchor, gate_stage + anchor, 1)

    TARGET.write_text(text, encoding="utf-8")
    compile(text, str(TARGET), "exec")

    print("=" * 78)
    print("ATLAS — FULL RUNNER SCIENTIFIC GATE PATCH COMPLETE")
    print("=" * 78)
    print(f"Patched: {TARGET}")
    print(f"Backup:  {backup}")
    print("Inserted: 00C4 -> 00W scientific gate -> 00E")
    print("\nNext:")
    print("python -u scripts/run_atlas_full_auto.py --list-stages")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
