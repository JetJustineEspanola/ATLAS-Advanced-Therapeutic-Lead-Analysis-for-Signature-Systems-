#!/usr/bin/env python3
from __future__ import annotations

import re
import shutil
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "scripts" / "run_atlas_full_auto.py"
MARKER = "ATLAS_DEPENDENCY_AWARE_RESUME_V2"


def main():
    if not TARGET.exists():
        raise SystemExit(f"ERROR: missing {TARGET}")

    original = TARGET.read_text(encoding="utf-8")

    if MARKER in original:
        print("Runner already contains dependency-aware resume v2; no changes made.")
        return 0

    text = original

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = TARGET.with_name(
        f"{TARGET.stem}_pre_dependency_resume_v2_{stamp}{TARGET.suffix}"
    )
    shutil.copy2(TARGET, backup)

    # 1. Add 00Y before 00A if absent.
    if '"00y", "00y_runtime_preflight.py"' not in text:
        stage_anchor = re.search(
            r'(?m)^    Stage\(\n'
            r'        "00a",\s*"00a_dataset_discovery\.py",\s*"acquisition",\n',
            text,
        )
        if not stage_anchor:
            raise RuntimeError("Could not locate the 00A stage definition.")

        preflight_stage = (
            '    Stage(\n'
            '        "00y", "00y_runtime_preflight.py", "preflight",\n'
            '        "Runtime, network, dependency, and credential preflight",\n'
            '        "results/pipeline_state/runtime_preflight.json",\n'
            '    ),\n\n'
        )

        pos = stage_anchor.start()
        text = text[:pos] + preflight_stage + text[pos:]

    # 2. Add script-change helper before select_stages().
    if "def stage_script_changed_since_success(" not in text:
        helper = '''
# ATLAS_DEPENDENCY_AWARE_RESUME_V2
def stage_script_changed_since_success(stage: Stage, previous: dict) -> bool:
    """True if this stage script changed after its last successful run."""
    ended = previous.get("ended_utc")
    if not ended:
        return False

    script_path = SCRIPTS / stage.script
    if not script_path.exists():
        return True

    try:
        ended_dt = datetime.fromisoformat(ended.replace("Z", "+00:00"))
        script_dt = datetime.fromtimestamp(
            script_path.stat().st_mtime,
            tz=timezone.utc,
        )
        return script_dt > ended_dt
    except Exception:
        return False


'''
        m = re.search(
            r'(?m)^def select_stages\(from_stage=None, to_stage=None\):',
            text,
        )
        if not m:
            raise RuntimeError("Could not locate select_stages().")
        text = text[:m.start()] + helper + text[m.start():]
    elif MARKER not in text:
        text = text.replace(
            "def stage_script_changed_since_success(",
            "# ATLAS_DEPENDENCY_AWARE_RESUME_V2\n"
            "def stage_script_changed_since_success(",
            1,
        )

    # 3. Add upstream_reran before stage loop.
    if "upstream_reran = False" not in text:
        m = re.search(r'(?m)^    for stage in selected:\s*$', text)
        if not m:
            raise RuntimeError("Could not locate main stage loop.")
        text = (
            text[:m.start()]
            + "    upstream_reran = False\n\n"
            + text[m.start():]
        )

    # 4. Replace skip/resume block.
    start_token = "        should_skip = False\n"
    end_token = "        if should_skip:\n"

    start = text.find(start_token)
    if start < 0:
        raise RuntimeError("Could not locate should_skip initialization.")

    end = text.find(end_token, start)
    if end < 0:
        raise RuntimeError("Could not locate should_skip decision endpoint.")

    end += len(end_token)

    new_skip = '''        should_skip = False

        # Runtime preflight always runs when it is selected.
        if stage.group == "preflight":
            should_skip = False

        elif not args.force and not args.refresh_data and not upstream_reran:
            previous_success = previous.get("status") == "SUCCESS"
            script_changed = stage_script_changed_since_success(stage, previous)

            if previous_success and not script_changed:
                if stage.checkpoint:
                    should_skip = checkpoint_exists(stage)
                else:
                    should_skip = True

            elif not previous_success and checkpoint_exists(stage):
                # Backward-compatible resume from a legacy checkpoint.
                should_skip = True

        if should_skip:
'''

    text = text[:start] + new_skip + text[end:]

    # 5. Dirty downstream only after non-preflight stage succeeds.
    if 'if stage.group != "preflight":' not in text:
        success_token = '        rec["status"] = "SUCCESS"\n'
        pos = text.find(success_token)
        if pos < 0:
            raise RuntimeError("Could not locate successful-stage status update.")

        save_token = "        save_state(state)\n"
        save_pos = text.find(save_token, pos)
        if save_pos < 0:
            raise RuntimeError("Could not locate save_state after stage success.")

        insert_pos = save_pos + len(save_token)

        insertion = (
            "\n"
            "        # Executed data/scientific stages invalidate downstream results.\n"
            "        # Preflight itself does not dirty the scientific pipeline.\n"
            '        if stage.group != "preflight":\n'
            "            upstream_reran = True\n"
        )

        text = text[:insert_pos] + insertion + text[insert_pos:]

    if MARKER not in text:
        text = "# ATLAS_DEPENDENCY_AWARE_RESUME_V2\n" + text

    # Syntax check before writing.
    compile(text, str(TARGET), "exec")
    TARGET.write_text(text, encoding="utf-8")

    print("=" * 88)
    print("ATLAS — FULL RUNNER HARDENING v2 COMPLETE")
    print("=" * 88)
    print(f"Patched: {TARGET}")
    print(f"Backup:  {backup}")
    print()
    print("Behavior:")
    print("- 00Y always runs when selected")
    print("- 00Y does not invalidate downstream outputs")
    print("- --refresh-data reruns the complete selected pipeline")
    print("- any non-preflight rerun forces downstream recomputation")
    print("- changed stage scripts invalidate previous success")
    print("- missing checkpoints invalidate checkpointed success")
    print()
    print("Verify with:")
    print("python -u scripts/run_atlas_full_auto.py --list-stages")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
