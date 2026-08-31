#!/usr/bin/env python3
"""
ATLAS — Full Automation Preflight Audit

Scans the current scripts/ directory and reports:
- missing scripts referenced by run_atlas_full_auto.py
- argparse/CLI parameters
- likely online/network dependencies
- likely manual/curated stages
- output/checkpoint hints
- external binaries/imports referenced
- whether each stage can plausibly run non-interactively

This script is READ-ONLY.

Output:
  results/pipeline_state/automation_preflight_audit.csv
  results/pipeline_state/automation_preflight_audit.txt
"""

from __future__ import annotations

import ast
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
OUTDIR = ROOT / "results" / "pipeline_state"
OUTDIR.mkdir(parents=True, exist_ok=True)

RUNNER = SCRIPTS / "run_atlas_full_auto.py"

OUT_CSV = OUTDIR / "automation_preflight_audit.csv"
OUT_TXT = OUTDIR / "automation_preflight_audit.txt"

NETWORK_TOKENS = [
    "requests.", "urllib", "http://", "https://", "rest.uniprot.org",
    "string-db", "clinicaltrials", "pubchem", "cmap", "clue.io",
    "geo", "ncbi", "gseapy.get_library",
]

MANUAL_TOKENS = [
    "curated", "manual", "confirm_her2", "inspect_dataset_labels",
    "input(", "TODO", "REVIEW_REQUIRED",
]

EXTERNAL_BIN_TOKENS = [
    "vina", "obabel", "openbabel", "docker", "java", "Rscript",
    "wget", "curl", "git ",
]


def get_runner_scripts():
    if not RUNNER.exists():
        return []

    text = RUNNER.read_text(encoding="utf-8", errors="replace")
    # Match Stage("key", "script.py", ...)
    return re.findall(
        r'Stage\(\s*"[^"]+"\s*,\s*"([^"]+\.py)"',
        text,
        flags=re.MULTILINE,
    )


def parse_cli(text):
    args = []
    for m in re.finditer(
        r'add_argument\(\s*["\'](--?[A-Za-z0-9_-]+)["\']',
        text,
    ):
        args.append(m.group(1))
    return sorted(set(args))


def imported_modules(path):
    try:
        tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"))
    except Exception:
        return []

    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for n in node.names:
                mods.add(n.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom) and node.module:
            mods.add(node.module.split(".")[0])
    return sorted(mods)


def classify(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    low = text.lower()

    cli = parse_cli(text)
    imports = imported_modules(path)

    network = any(tok.lower() in low for tok in NETWORK_TOKENS)
    manual = any(tok.lower() in low for tok in MANUAL_TOKENS)
    external = sorted({
        tok for tok in EXTERNAL_BIN_TOKENS
        if tok.lower() in low
    })

    # Heuristic: scripts calling input() are not fully non-interactive.
    interactive = "input(" in text

    # Presence of main guard is useful for orchestration.
    has_main = 'if __name__ == "__main__"' in text

    # Detect argparse required=True.
    required_cli = bool(
        re.search(r"add_argument\([^)]*required\s*=\s*True", text, re.S)
    )

    can_run_default = has_main and not interactive and not required_cli

    notes = []
    if interactive:
        notes.append("USES_INPUT")
    if required_cli:
        notes.append("REQUIRED_CLI_ARG")
    if manual:
        notes.append("MANUAL_OR_CURATED_LOGIC_PRESENT")
    if network:
        notes.append("NETWORK_DEPENDENCY")
    if external:
        notes.append("EXTERNAL_TOOL_REFERENCE")
    if not has_main:
        notes.append("NO_MAIN_GUARD")

    return {
        "script": path.name,
        "exists": True,
        "has_main_guard": has_main,
        "cli_args": " ".join(cli),
        "required_cli_detected": required_cli,
        "network_dependency": network,
        "manual_or_curated_logic": manual,
        "interactive_input": interactive,
        "external_tool_refs": " | ".join(external),
        "imports": " ".join(imports),
        "default_noninteractive_candidate": can_run_default,
        "notes": " | ".join(notes),
    }


def main():
    referenced = get_runner_scripts()

    rows = []
    seen = set()

    for name in referenced:
        path = SCRIPTS / name
        seen.add(name)

        if not path.exists():
            rows.append({
                "script": name,
                "exists": False,
                "has_main_guard": False,
                "cli_args": "",
                "required_cli_detected": False,
                "network_dependency": False,
                "manual_or_curated_logic": False,
                "interactive_input": False,
                "external_tool_refs": "",
                "imports": "",
                "default_noninteractive_candidate": False,
                "notes": "MISSING_SCRIPT_REFERENCED_BY_RUNNER",
            })
        else:
            rows.append(classify(path))

    # Also report scripts not yet orchestrated.
    for path in sorted(SCRIPTS.glob("*.py")):
        if path.name in seen or path.name == RUNNER.name:
            continue

        row = classify(path)
        row["notes"] = (
            (row["notes"] + " | " if row["notes"] else "")
            + "NOT_REFERENCED_BY_FULL_RUNNER"
        )
        rows.append(row)

    fields = [
        "script",
        "exists",
        "has_main_guard",
        "cli_args",
        "required_cli_detected",
        "network_dependency",
        "manual_or_curated_logic",
        "interactive_input",
        "external_tool_refs",
        "imports",
        "default_noninteractive_candidate",
        "notes",
    ]

    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    missing = [r for r in rows if not r["exists"]]
    interactive = [r for r in rows if r["interactive_input"]]
    required = [r for r in rows if r["required_cli_detected"]]
    manual = [r for r in rows if r["manual_or_curated_logic"]]
    network = [r for r in rows if r["network_dependency"]]
    not_runner = [r for r in rows if "NOT_REFERENCED_BY_FULL_RUNNER" in r["notes"]]

    lines = [
        "=" * 88,
        "ATLAS — FULL AUTOMATION PREFLIGHT AUDIT",
        "=" * 88,
        f"Runner: {RUNNER}",
        f"Runner-referenced stages: {len(referenced)}",
        f"Scripts scanned: {len(rows)}",
        "",
        f"Missing runner scripts: {len(missing)}",
        f"Scripts using input(): {len(interactive)}",
        f"Scripts with detected required CLI args: {len(required)}",
        f"Scripts with manual/curated logic: {len(manual)}",
        f"Scripts with network dependencies: {len(network)}",
        f"Scripts not referenced by full runner: {len(not_runner)}",
        "",
    ]

    def section(title, items):
        lines.append(title)
        lines.append("-" * len(title))
        if not items:
            lines.append("None")
        else:
            for r in items:
                lines.append(
                    f"- {r['script']}: {r['notes'] or 'flagged'}"
                )
        lines.append("")

    section("MISSING", missing)
    section("INTERACTIVE", interactive)
    section("REQUIRED CLI", required)
    section("MANUAL / CURATED", manual)
    section("NETWORK-DEPENDENT", network)

    lines.extend([
        "AUTOMATION RULE",
        "---------------",
        "A stage should only be considered fully automatable when it:",
        "1. has no input() prompt,",
        "2. has no unresolved required CLI argument,",
        "3. fails with a nonzero exit code on unsafe/missing inputs,",
        "4. writes a deterministic checkpoint/output,",
        "5. exposes manual-review cases as explicit status rather than silently guessing.",
        "",
        f"CSV: {OUT_CSV}",
    ])

    OUT_TXT.write_text("\n".join(lines), encoding="utf-8")

    print("\n".join(lines))
    print(f"TXT: {OUT_TXT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
