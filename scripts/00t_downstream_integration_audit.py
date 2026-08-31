#!/usr/bin/env python3
"""
ATLAS — 00T Downstream Integration Audit

Purpose
-------
Inspect 04Q / 04R / 04U scripts before replacing the old discovery-only
resistance evidence with the newly validated external-validation layer.

This is read-only. It does not modify any pipeline code or results.

Looks for:
- discovery DEG fallback references
- 03B / 03C references
- network / evidence input files
- gene-symbol joins
- scoring fields that may need to consume validated evidence

Output
------
results/external_validation/downstream/downstream_integration_audit.txt
"""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
OUTDIR = ROOT / "results/external_validation/downstream"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTFILE = OUTDIR / "downstream_integration_audit.txt"

CANDIDATES = [
    "04q",
    "04r",
    "04u",
]

KEYWORDS = [
    "differential_expression",
    "DEGs_resistant_vs_sensitive",
    "03b",
    "03c",
    "consensus",
    "fallback",
    "gene_symbol",
    "gene",
    "network",
    "validated",
    "external_validation",
    "integrated_evidence",
]


def main():
    lines = []
    lines.append("=" * 90)
    lines.append("ATLAS — 00T DOWNSTREAM INTEGRATION AUDIT")
    lines.append("=" * 90)

    matched_scripts = []
    for p in sorted(SCRIPTS.glob("*.py")):
        low = p.name.lower()
        if any(low.startswith(prefix) for prefix in CANDIDATES):
            matched_scripts.append(p)

    if not matched_scripts:
        lines.append("\nNo 04Q/04R/04U Python scripts found under scripts/.")
        OUTFILE.write_text("\n".join(lines), encoding="utf-8")
        print("\n".join(lines))
        print(f"\nOutput: {OUTFILE}")
        return 1

    for path in matched_scripts:
        text = path.read_text(encoding="utf-8", errors="replace").splitlines()

        lines.append("\n" + "-" * 90)
        lines.append(f"SCRIPT: {path.relative_to(ROOT)}")
        lines.append("-" * 90)

        hit_n = 0
        for i, line in enumerate(text, start=1):
            low = line.lower()
            if any(k.lower() in low for k in KEYWORDS):
                hit_n += 1
                start = max(1, i - 2)
                end = min(len(text), i + 2)

                lines.append(f"\n[context around line {i}]")
                for j in range(start, end + 1):
                    marker = ">>" if j == i else "  "
                    lines.append(f"{marker} {j:4d}: {text[j-1]}")

        if hit_n == 0:
            lines.append("\nNo matching integration keywords found.")

    lines.append("\n" + "=" * 90)
    lines.append("NEW VALIDATED EVIDENCE FILES")
    lines.append("=" * 90)
    lines.append(
        "results/external_validation/downstream/validated_resistance_gene_evidence.csv"
    )
    lines.append(
        "results/external_validation/downstream/validated_tgfb_module_evidence.csv"
    )
    lines.append(
        "results/external_validation/downstream/validated_evidence_summary.csv"
    )

    OUTFILE.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nOutput: {OUTFILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
