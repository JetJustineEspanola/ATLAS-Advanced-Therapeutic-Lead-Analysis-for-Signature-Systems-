#!/usr/bin/env python3
"""
ATLAS — 00T1 Focused 04Q Integration Extract

Read-only diagnostic. Finds the actual 04Q script and prints the exact code
around:
- DISCOVERY_DEG
- CONSENSUS_CANDIDATES
- resistance gene loader/selection
- max_resistance_genes
- mapped_resistance
- resistance_genes assignment

This gives enough context to patch 04Q safely without guessing.
"""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
OUT = (
    ROOT
    / "results"
    / "external_validation"
    / "downstream"
    / "04q_focused_integration_extract.txt"
)
OUT.parent.mkdir(parents=True, exist_ok=True)

PATTERNS = [
    r"DISCOVERY_DEG",
    r"CONSENSUS_CANDIDATES",
    r"resistance_genes",
    r"mapped_resistance",
    r"max_resistance_genes",
    r"load.*resistance",
    r"consensus",
    r"fallback",
]


def main():
    candidates = sorted(
        p for p in SCRIPTS.glob("*.py")
        if p.name.lower().startswith("04q")
    )

    if not candidates:
        raise SystemExit("ERROR: No scripts/04q*.py file found.")

    lines_out = []
    lines_out.append("=" * 90)
    lines_out.append("ATLAS — 00T1 FOCUSED 04Q INTEGRATION EXTRACT")
    lines_out.append("=" * 90)

    for path in candidates:
        src = path.read_text(encoding="utf-8", errors="replace").splitlines()

        lines_out.append(f"\nSCRIPT: {path.relative_to(ROOT)}")
        lines_out.append("-" * 90)

        hits = set()
        for i, line in enumerate(src):
            if any(re.search(pat, line, flags=re.I) for pat in PATTERNS):
                for j in range(max(0, i - 12), min(len(src), i + 18)):
                    hits.add(j)

        if not hits:
            lines_out.append("No relevant matches found.")
            continue

        # Merge nearby line indexes into readable blocks.
        ordered = sorted(hits)
        blocks = []
        start = prev = ordered[0]
        for idx in ordered[1:]:
            if idx <= prev + 2:
                prev = idx
            else:
                blocks.append((start, prev))
                start = prev = idx
        blocks.append((start, prev))

        for start, end in blocks:
            lines_out.append(
                f"\n[lines {start + 1}-{end + 1}]"
            )
            for j in range(start, end + 1):
                lines_out.append(f"{j + 1:4d}: {src[j]}")

    text = "\n".join(lines_out)
    OUT.write_text(text, encoding="utf-8")
    print(text)
    print(f"\nOutput: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
