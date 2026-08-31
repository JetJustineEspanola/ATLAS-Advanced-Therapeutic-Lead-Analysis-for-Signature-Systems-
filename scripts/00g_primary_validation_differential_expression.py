#!/usr/bin/env python3
"""
ATLAS — 00G Primary Validation Differential Expression

Runs DESeq2-style differential expression with PyDESeq2 for:
- GSE121105: curated 6 resistant vs 3 parental samples
- GSE237606: all 40 samples, adjusted for timepoint/stimulus stratum

Outputs:
  results/external_validation/GSE121105_DE.csv
  results/external_validation/GSE237606_DE.csv
  results/external_validation/primary_validation_DE_summary.csv
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/validation_expression"
OUT = ROOT / "results/external_validation"
OUT.mkdir(parents=True, exist_ok=True)

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
except ImportError:
    print("ERROR: PyDESeq2 is not installed.")
    print("Install with: python -m pip install -U pydeseq2")
    raise SystemExit(2)


def read_counts(path: Path):
    if path.name.endswith(".csv.gz"):
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(path, sep="\t")

    gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene_id"})
    df["gene_id"] = df["gene_id"].astype(str)

    # Collapse duplicate gene IDs conservatively by summing counts.
    sample_cols = [c for c in df.columns if c != "gene_id"]
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    df = df.groupby("gene_id", as_index=False)[sample_cols].sum()
    return df


def make_dds(counts, metadata, design):
    # Newer PyDESeq2 API.
    try:
        return DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design=design,
            refit_cooks=True,
            n_cpus=1,
        )
    except TypeError:
        # Compatibility fallback for older releases.
        factors = [x.strip() for x in design.replace("~", "").split("+")]
        return DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design_factors=factors,
            refit_cooks=True,
            n_cpus=1,
        )


def run_deseq(accession, count_path, design_path, design_formula):
    print(f"\n[{accession}]")

    mat = read_counts(count_path)
    design = pd.read_csv(design_path)

    # GSE121105 uses only curated primary contrast.
    if accession == "GSE121105":
        if "included" not in design.columns:
            raise RuntimeError("Missing included column in GSE121105 design")
        design = design[design["included"].astype(str).str.lower().isin(["true", "1"])].copy()

    # Matrix column identifiers are stored in GEO description for both datasets.
    design["count_column"] = design["description"].astype(str)

    missing = [c for c in design["count_column"] if c not in mat.columns]
    if missing:
        raise RuntimeError(
            f"{accession}: {len(missing)} design samples not found in count matrix: {missing[:10]}"
        )

    # Build sample x gene count table expected by PyDESeq2.
    selected = mat[["gene_id"] + list(design["count_column"])].copy()
    counts = selected.set_index("gene_id").T
    counts.index.name = "sample"
    counts = counts.astype(int)

    # Filter very low-count genes before model fitting.
    keep = (counts.sum(axis=0) >= 10) & ((counts > 0).sum(axis=0) >= 2)
    counts = counts.loc[:, keep]

    metadata = pd.DataFrame(index=counts.index)
    design = design.set_index("count_column").loc[counts.index]

    metadata["phenotype"] = design["phenotype"].astype(str).values

    if accession == "GSE237606":
        metadata["stratum"] = (
            design["timepoint"].astype(str) + "_" + design["stimulus"].astype(str)
        ).values

    print(f"  samples: {counts.shape[0]}")
    print(f"  genes after count filter: {counts.shape[1]}")
    print("  phenotype counts:")
    print(metadata["phenotype"].value_counts().to_string())

    if "stratum" in metadata.columns:
        print(f"  strata: {metadata['stratum'].nunique()}")

    dds = make_dds(counts, metadata, design_formula)
    dds.deseq2()

    stat = DeseqStats(
        dds,
        contrast=["phenotype", "RESISTANT", "SENSITIVE_OR_PARENTAL"],
        n_cpus=1,
    )
    stat.summary()

    res = stat.results_df.copy()
    res.index.name = "gene_id"
    res = res.reset_index()

    # Useful derived columns.
    res["direction"] = np.where(
        res["log2FoldChange"] > 0, "UP_IN_RESISTANT",
        np.where(res["log2FoldChange"] < 0, "DOWN_IN_RESISTANT", "NO_CHANGE")
    )
    res["significant_fdr05"] = res["padj"].notna() & (res["padj"] < 0.05)

    res = res.sort_values(
        ["padj", "pvalue"],
        na_position="last"
    )

    out = OUT / f"{accession}_DE.csv"
    res.to_csv(out, index=False)

    sig = res[res["significant_fdr05"]]
    print(f"  FDR<0.05 genes: {len(sig)}")
    if not sig.empty:
        print(f"  up in resistant: {(sig['log2FoldChange'] > 0).sum()}")
        print(f"  down in resistant: {(sig['log2FoldChange'] < 0).sum()}")

    return {
        "accession": accession,
        "samples": counts.shape[0],
        "genes_tested": counts.shape[1],
        "fdr05_genes": len(sig),
        "up_fdr05": int((sig["log2FoldChange"] > 0).sum()),
        "down_fdr05": int((sig["log2FoldChange"] < 0).sum()),
        "design": design_formula,
        "output": str(out),
    }


def main():
    print("=" * 78)
    print("ATLAS — 00G PRIMARY VALIDATION DIFFERENTIAL EXPRESSION")
    print("=" * 78)

    summaries = []

    summaries.append(
        run_deseq(
            "GSE121105",
            DATA / "GSE121105/GSE121105_geneCount.csv.gz",
            DATA / "GSE121105_primary_design.csv",
            "~ phenotype",
        )
    )

    summaries.append(
        run_deseq(
            "GSE237606",
            DATA / "GSE237606/GSE237606_RawCounts.txt.gz",
            DATA / "GSE237606_primary_design.csv",
            "~ stratum + phenotype",
        )
    )

    summary = pd.DataFrame(summaries)
    out = OUT / "primary_validation_DE_summary.csv"
    summary.to_csv(out, index=False)

    print("\n" + "=" * 78)
    print("00G COMPLETE")
    print("=" * 78)
    print(summary.to_string(index=False))
    print(f"\nSummary: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
