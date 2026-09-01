#!/usr/bin/env python3
"""
ATLAS — 00V Exact Automation Blocker Audit

Uses Python AST instead of keyword heuristics to identify real blockers:
- actual input() calls
- argparse arguments with required=True
- whether scripts execute meaningful top-level code without a main guard
- subprocess/external-tool calls
- network-library usage

Focuses on stages used by run_atlas_full_auto.py.

Output:
  results/pipeline_state/automation_blockers_exact.csv
  results/pipeline_state/automation_blockers_exact.txt
"""

from __future__ import annotations

import ast
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
RUNNER = SCRIPTS / "run_atlas_full_auto.py"
OUTDIR = ROOT / "results" / "pipeline_state"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_CSV = OUTDIR / "automation_blockers_exact.csv"
OUT_TXT = OUTDIR / "automation_blockers_exact.txt"


def runner_scripts():
    if not RUNNER.exists():
        return []
    text = RUNNER.read_text(encoding="utf-8", errors="replace")
    return re.findall(
        r'Stage\(\s*"[^"]+"\s*,\s*"([^"]+\.py)"',
        text,
        flags=re.MULTILINE,
    )


def call_name(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        left = call_name(node.value)
        return f"{left}.{node.attr}" if left else node.attr
    return ""


def has_main_guard(tree):
    for node in tree.body:
        if isinstance(node, ast.If):
            try:
                txt = ast.unparse(node.test)
            except Exception:
                txt = ""
            if "__name__" in txt and "__main__" in txt:
                return True
    return False


def meaningful_top_level(tree):
    """
    True when there is executable top-level code beyond imports, assignments,
    definitions, and the main guard.
    """
    allowed = (
        ast.Import,
        ast.ImportFrom,
        ast.Assign,
        ast.AnnAssign,
        ast.FunctionDef,
        ast.AsyncFunctionDef,
        ast.ClassDef,
    )

    for node in tree.body:
        if isinstance(node, allowed):
            continue

        if isinstance(node, ast.If):
            try:
                txt = ast.unparse(node.test)
            except Exception:
                txt = ""
            if "__name__" in txt and "__main__" in txt:
                continue

        # module docstring
        if (
            isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
        ):
            continue

        return True

    return False


def analyze(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace")

    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        return {
            "script": path.name,
            "syntax_ok": False,
            "input_calls": "",
            "required_cli": "",
            "main_guard": False,
            "top_level_execution": False,
            "network_calls": "",
            "external_calls": "",
            "automation_status": "BLOCKED_SYNTAX",
            "notes": str(exc),
        }

    input_lines = []
    required_cli = []
    network_calls = set()
    external_calls = set()

    network_prefixes = {
        "requests",
        "urllib",
        "httpx",
        "aiohttp",
    }

    external_names = {
        "subprocess.run",
        "subprocess.Popen",
        "subprocess.call",
        "os.system",
    }

    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = call_name(node.func)

            if name == "input":
                input_lines.append(str(getattr(node, "lineno", "?")))

            if name.endswith("add_argument"):
                required = False
                option = ""
                if node.args:
                    try:
                        option = ast.literal_eval(node.args[0])
                    except Exception:
                        option = ast.unparse(node.args[0])

                for kw in node.keywords:
                    if kw.arg == "required":
                        try:
                            required = bool(ast.literal_eval(kw.value))
                        except Exception:
                            pass

                if required:
                    required_cli.append(str(option))

            if any(name == p or name.startswith(p + ".") for p in network_prefixes):
                network_calls.add(name)

            if name in external_names:
                external_calls.add(name)

    main = has_main_guard(tree)
    top_exec = meaningful_top_level(tree)

    blockers = []
    if input_lines:
        blockers.append("INTERACTIVE_INPUT")
    if required_cli:
        blockers.append("REQUIRED_CLI")

    # No main guard is not automatically a blocker if the script intentionally
    # executes top-level code when invoked as `python script.py`.
    if not main and not top_exec:
        blockers.append("NO_ENTRYPOINT")

    status = "READY_NONINTERACTIVE" if not blockers else "BLOCKED"

    notes = []
    if not main and top_exec:
        notes.append("NO_MAIN_GUARD_BUT_HAS_TOP_LEVEL_EXECUTION")
    if network_calls:
        notes.append("NETWORK_DEPENDENT")
    if external_calls:
        notes.append("EXTERNAL_SUBPROCESS")

    return {
        "script": path.name,
        "syntax_ok": True,
        "input_calls": ",".join(input_lines),
        "required_cli": " | ".join(required_cli),
        "main_guard": main,
        "top_level_execution": top_exec,
        "network_calls": " | ".join(sorted(network_calls)),
        "external_calls": " | ".join(sorted(external_calls)),
        "automation_status": status,
        "notes": " | ".join(notes),
    }


def main():
    names = runner_scripts()

    rows = []
    for name in names:
        path = SCRIPTS / name
        if not path.exists():
            rows.append({
                "script": name,
                "syntax_ok": False,
                "input_calls": "",
                "required_cli": "",
                "main_guard": False,
                "top_level_execution": False,
                "network_calls": "",
                "external_calls": "",
                "automation_status": "MISSING",
                "notes": "RUNNER_REFERENCES_MISSING_SCRIPT",
            })
        else:
            rows.append(analyze(path))

    fields = [
        "script",
        "syntax_ok",
        "input_calls",
        "required_cli",
        "main_guard",
        "top_level_execution",
        "network_calls",
        "external_calls",
        "automation_status",
        "notes",
    ]

    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    blockers = [r for r in rows if r["automation_status"] not in {"READY_NONINTERACTIVE"}]

    lines = [
        "=" * 88,
        "ATLAS — EXACT AUTOMATION BLOCKER AUDIT",
        "=" * 88,
        f"Runner stages checked: {len(rows)}",
        f"Actual blockers found: {len(blockers)}",
        "",
        "BLOCKERS",
        "--------",
    ]

    if not blockers:
        lines.append("None")
    else:
        for r in blockers:
            lines.append(
                f"- {r['script']}: status={r['automation_status']} "
                f"input_lines={r['input_calls'] or '-'} "
                f"required_cli={r['required_cli'] or '-'} "
                f"notes={r['notes'] or '-'}"
            )

    lines.extend([
        "",
        "IMPORTANT",
        "---------",
        "Network access, curated rules, or a missing main guard are not treated as blockers by themselves.",
        "A no-main-guard script is acceptable if it has deliberate top-level execution when invoked directly.",
        "",
        f"CSV: {OUT_CSV}",
        f"TXT: {OUT_TXT}",
    ])

    OUT_TXT.write_text("\n".join(lines), encoding="utf-8")

    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
