#!/usr/bin/env python3
"""
ATLAS — 00J1 Discovery ID Diagnostic

Inspects the discovery DE table and consensus signature to determine why
00J reported zero gene overlap.

Outputs:
  results/external_validation/discovery_id_diagnostic.txt
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DISCOVERY = ROOT / "results/differential_expression/DEGs_resistant_vs_sensitive_annotated.csv"
CONSENSUS = ROOT / "results/external_validation/consensus_resistance_signature.csv"
OUT = ROOT / "results/external_validation/discovery_id_diagnostic.txt"


def main():
    d = pd.read_csv(DISCOVERY)
    c = pd.read_csv(CONSENSUS)

    lines = []
    lines.append("ATLAS — 00J1 DISCOVERY ID DIAGNOSTIC")
    lines.append("=" * 78)

    lines.append("\nDiscovery columns:")
    for col in d.columns:
        lines.append(f"- {col}")

    lines.append("\nDiscovery first 10 rows:")
    lines.append(d.head(10).to_string(index=False))

    lines.append("\nConsensus first 20 gene symbols:")
    lines.append(
        ", ".join(c["gene_symbol"].astype(str).head(20).tolist())
    )

    # Check likely gene identifier columns.
    candidate_cols = [
        col for col in d.columns
        if any(k in col.lower() for k in [
            "gene", "symbol", "ensembl", "id", "name"
        ])
    ]

    lines.append("\nCandidate discovery identifier columns:")
    for col in candidate_cols:
        vals = d[col].dropna().astype(str).head(20).tolist()
        lines.append(f"\n[{col}]")
        lines.append(", ".join(vals))

    # Direct overlap checks against each likely column.
    consensus = set(c["gene_symbol"].dropna().astype(str).str.upper().str.strip())
    lines.append("\nDirect overlap by candidate column:")
    for col in candidate_cols:
        vals = set(
            d[col].dropna().astype(str)
            .str.upper().str.strip()
        )
        overlap = len(consensus & vals)
        lines.append(f"{col}: {overlap}")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nOutput: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
