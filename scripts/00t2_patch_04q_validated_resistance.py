#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import re
import shutil
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "scripts/04q_network_integration.py"
MARKER = "ATLAS_04Q_VALIDATED_RESISTANCE_PATCH_V1"

def main():
    if not TARGET.exists():
        raise SystemExit(f"ERROR: missing {TARGET}")

    text = TARGET.read_text(encoding="utf-8")

    if MARKER in text:
        print("04Q is already patched; no changes made.")
        return 0

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = TARGET.with_name(
        f"{TARGET.stem}_pre_validated_resistance_{stamp}{TARGET.suffix}"
    )
    shutil.copy2(TARGET, backup)

    anchor = "CONSENSUS_CANDIDATES = ["
    if anchor not in text:
        raise RuntimeError("Could not find CONSENSUS_CANDIDATES anchor.")

    validated_constant = '''# ATLAS_04Q_VALIDATED_RESISTANCE_PATCH_V1
VALIDATED_RESISTANCE = (
    PROJECT_ROOT
    / "results"
    / "external_validation"
    / "downstream"
    / "validated_resistance_gene_evidence.csv"
)

'''
    text = text.replace(anchor, validated_constant + anchor, 1)

    sig_pattern = re.compile(
        r"(def load_resistance_genes\(\n"
        r"\s*max_genes: int,\n"
        r"\s*padj_cutoff: float,\n"
        r"\s*abs_fc_cutoff: float,\n"
        r"\s*\) -> pd\.DataFrame:\n)"
    )
    m = sig_pattern.search(text)
    if not m:
        raise RuntimeError("Could not locate load_resistance_genes() signature.")

    insertion = '''    # Prefer the externally validated three-dataset strict core.
    # Older 03C/discovery inputs remain fallback only.
    if VALIDATED_RESISTANCE.exists():
        try:
            validated = pd.read_csv(VALIDATED_RESISTANCE)

            if "gene_symbol" not in validated.columns:
                raise ValueError(
                    "Validated resistance file is missing gene_symbol."
                )

            keep_cols = [
                c for c in [
                    "gene_symbol",
                    "evidence_class",
                    "resistance_direction",
                    "validated_evidence_strength",
                    "external_consensus_rank",
                    "minimum_abs_log2FC_across_three",
                    "tgfb_module_role",
                ]
                if c in validated.columns
            ]

            out = validated[keep_cols].copy()
            out["gene_symbol"] = out["gene_symbol"].map(clean_text)
            out = out[out["gene_symbol"].ne("")].drop_duplicates("gene_symbol")

            if "external_consensus_rank" in out.columns:
                out["external_consensus_rank"] = pd.to_numeric(
                    out["external_consensus_rank"],
                    errors="coerce",
                )
                out = out.sort_values(
                    "external_consensus_rank",
                    ascending=True,
                    na_position="last",
                )

            out["resistance_gene_source"] = "THREE_DATASET_STRICT_CORE"

            if len(out) > max_genes:
                print(
                    f"WARNING: validated resistance core contains {len(out):,} genes "
                    f"but --max-resistance-genes={max_genes:,}; truncating by rank.",
                    flush=True,
                )

            if not out.empty:
                print(
                    f"Using validated three-dataset resistance core: "
                    f"{min(len(out), max_genes):,}/{len(out):,} genes",
                    flush=True,
                )
                return out.head(max_genes).reset_index(drop=True)

        except Exception as exc:
            print(
                f"WARNING: could not use validated resistance file "
                f"{VALIDATED_RESISTANCE}: {exc}. "
                f"Falling back to legacy 03C/discovery logic.",
                flush=True,
            )

'''

    text = text[:m.end()] + insertion + text[m.end():]

    if '"discovery_deg": str(DISCOVERY_DEG),' in text:
        text = text.replace(
            '"discovery_deg": str(DISCOVERY_DEG),',
            '"validated_resistance_genes": str(VALIDATED_RESISTANCE),\n'
            '        "discovery_deg_fallback": str(DISCOVERY_DEG),',
            1,
        )

    text = text.replace(
        '"ATLAS discovery/consensus resistance genes",',
        '"ATLAS three-dataset validated resistance genes '
        '(legacy consensus/discovery fallback retained)",',
        1,
    )

    TARGET.write_text(text, encoding="utf-8")
    compile(text, str(TARGET), "exec")

    print("=" * 78)
    print("ATLAS — 00T2 04Q PATCH COMPLETE")
    print("=" * 78)
    print(f"Patched: {TARGET}")
    print(f"Backup:  {backup}")
    print(
        "Primary resistance input: "
        + str(
            ROOT
            / "results"
            / "external_validation"
            / "downstream"
            / "validated_resistance_gene_evidence.csv"
        )
    )
    print("\nNext: rerun 04Q with --max-resistance-genes 242.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
