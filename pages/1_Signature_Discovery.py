"""
Stage 1: Signature Discovery
Tabs: Validation (QC) | Differential Expression | Pathway Analysis | Gene Prioritization

Reads only from results/ CSVs already produced by
01_validation.py, 02_differential_expression.py, 03_pathway_analysis.py,
0X_gene_prioritization.py.
No computation happens in this app - it is a thin display layer. Filenames
are the stable output contract those scripts write to; every statistic
shown is computed from the live file content at load time, never hardcoded.

Note: CMap (Stage 4) intentionally lives on its own sidebar page, not here -
keeps the sidebar at "4 stages + Results" and gives CMap room for its 5
dedicated tabs (signature construction, cross-signature results, candidate
explorer, pipeline status, provenance) rather than a cramped summary here.
"""

from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

st.set_page_config(page_title="ATLAS - Signature Discovery", layout="wide")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
QC_DIR = PROJECT_ROOT / "results" / "qc"
DE_DIR = PROJECT_ROOT / "results" / "differential_expression"
PATHWAY_DIR = PROJECT_ROOT / "results" / "pathway_analysis"
PRIORITY_DIR = PROJECT_ROOT / "results" / "gene_prioritization"

st.title("Stage 1 - Signature Discovery")


def missing(file_path: Path, run_script: str) -> None:
    st.info(f"`{file_path.name}` not found yet - run `{run_script}` first.")


def last_modified(file_path: Path) -> str:
    return datetime.fromtimestamp(file_path.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")


def provenance_caption(file_path: Path) -> None:
    st.caption(f"Source: `{file_path.name}` \u00b7 Updated {last_modified(file_path)}")


# ============================================================
# Pipeline status strip
# ============================================================
# Each sub-stage's status is derived from whether its key output file
# exists right now - not tracked/cached separately, so it can never drift
# from what the tabs below are actually showing.

status_targets = [
    ("Validation", QC_DIR / "stage1_validation_summary.csv"),
    ("Differential Expression", DE_DIR / "DEGs_resistant_vs_sensitive_annotated.csv"),
    ("Pathway Analysis", PATHWAY_DIR / "targeted_pathway_summary.csv"),
    ("Gene Prioritization", PRIORITY_DIR / "gene_priority.csv"),
]

status_cols = st.columns(len(status_targets))
for col, (label, path) in zip(status_cols, status_targets):
    if path.exists():
        col.markdown(f"\u2705 **{label}**")
        col.caption(f"Updated {last_modified(path)}")
    else:
        col.markdown(f"\u25cb **{label}**")
        col.caption("Not run yet")

st.divider()

tab_validation, tab_deg, tab_pathway, tab_priority = st.tabs(
    ["Validation (QC)", "Differential Expression", "Pathway Analysis", "Gene Prioritization"]
)


# ============================================================
# Tab 1: Validation / QC
# ============================================================
with tab_validation:
    summary_file = QC_DIR / "stage1_validation_summary.csv"
    library_file = QC_DIR / "library_sizes.csv"
    corr_file = QC_DIR / "sample_correlation.csv"
    pca_file = QC_DIR / "pca_coordinates.csv"

    if summary_file.exists():
        provenance_caption(summary_file)
        summary = pd.read_csv(summary_file)
        st.dataframe(summary, width='stretch', hide_index=True)
    else:
        missing(summary_file, "01_validation.py")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Library size per sample")
        if library_file.exists():
            provenance_caption(library_file)
            lib = pd.read_csv(library_file, index_col=0)
            fig = px.bar(
                lib,
                x=lib.index,
                y=lib.columns[0],
                labels={"x": "Sample", lib.columns[0]: "Total estimated counts"},
            )
            st.plotly_chart(fig, width='stretch')
        else:
            missing(library_file, "01_validation.py")

    with col2:
        st.subheader("Sample-to-sample correlation")
        if corr_file.exists():
            provenance_caption(corr_file)
            corr = pd.read_csv(corr_file, index_col=0)
            fig = px.imshow(
                corr,
                text_auto=".3f",
                color_continuous_scale="Blues",
                zmin=0,
                zmax=1,
            )
            st.plotly_chart(fig, width='stretch')
        else:
            missing(corr_file, "01_validation.py")

    st.subheader("PCA")
    if pca_file.exists():
        provenance_caption(pca_file)
        pca_df = pd.read_csv(pca_file, index_col=0)
        fig = px.scatter(
            pca_df,
            x="PC1",
            y="PC2",
            color="Condition",
            text=pca_df.index,
            hover_name=pca_df.index,
        )
        fig.update_traces(textposition="top center", marker=dict(size=12))
        st.plotly_chart(fig, width='stretch')
    else:
        missing(pca_file, "01_validation.py")


# ============================================================
# Tab 2: Differential Expression
# ============================================================
with tab_deg:
    all_degs_file = DE_DIR / "DEGs_resistant_vs_sensitive_annotated.csv"
    sig_degs_file = DE_DIR / "significant_DEGs_annotated.csv"

    if all_degs_file.exists():
        provenance_caption(all_degs_file)
        degs = pd.read_csv(all_degs_file, index_col=0)

        n_total = len(degs)
        n_sig = int((degs["padj"] < 0.05).sum())
        n_up = int(((degs["padj"] < 0.05) & (degs["log2FoldChange"] >= 1)).sum())
        n_down = int(((degs["padj"] < 0.05) & (degs["log2FoldChange"] <= -1)).sum())

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Genes tested", f"{n_total:,}")
        c2.metric("Significant (padj<0.05)", f"{n_sig:,}")
        c3.metric("Up in resistant", f"{n_up:,}")
        c4.metric("Down in resistant", f"{n_down:,}")

        st.subheader("Volcano plot")
        sig_threshold = st.slider(
            "padj significance threshold", 0.001, 0.10, 0.05, step=0.001
        )

        plot_df = degs.copy()
        plot_df["padj_plot"] = plot_df["padj"].clip(lower=1e-300)
        plot_df["neg_log10_padj"] = -np.log10(plot_df["padj_plot"])
        plot_df["significant"] = plot_df["padj"] < sig_threshold

        fig = px.scatter(
            plot_df,
            x="log2FoldChange",
            y="neg_log10_padj",
            color="significant",
            hover_name="Gene name" if "Gene name" in plot_df.columns else None,
            hover_data=["padj", "baseMean"],
            opacity=0.6,
            color_discrete_map={True: "crimson", False: "lightgray"},
        )
        fig.add_hline(y=-np.log10(sig_threshold), line_dash="dash", line_color="gray")
        fig.add_vline(x=1, line_dash="dash", line_color="gray")
        fig.add_vline(x=-1, line_dash="dash", line_color="gray")
        fig.update_layout(yaxis_title="-log10(padj)")
        st.plotly_chart(fig, width='stretch')

        st.subheader("DEG table")
        search = st.text_input("Search gene name")
        table = degs.reset_index()
        if search and "Gene name" in table.columns:
            table = table[table["Gene name"].str.contains(search, case=False, na=False)]
        st.dataframe(
            table.sort_values("padj", na_position="last"),
            width='stretch',
            hide_index=True,
        )
        st.download_button(
            "Download filtered table as CSV",
            table.to_csv(index=False),
            file_name="deg_filtered.csv",
        )
    else:
        missing(all_degs_file, "02_differential_expression.py")


# ============================================================
# Tab 3: Pathway Analysis
# ============================================================
with tab_pathway:
    gsea_file = PATHWAY_DIR / "gsea_results.csv"
    enriched_up_file = PATHWAY_DIR / "enriched_upregulated.csv"
    enriched_down_file = PATHWAY_DIR / "enriched_downregulated.csv"
    targeted_summary_file = PATHWAY_DIR / "targeted_pathway_summary.csv"
    targeted_genes_file = PATHWAY_DIR / "targeted_gene_results.csv"

    st.subheader("Targeted pathway check")
    st.caption("TGF-\u03b2 signaling and the PD-1/PD-L1 axis specifically")

    if targeted_summary_file.exists():
        provenance_caption(targeted_summary_file)
        summary = pd.read_csv(targeted_summary_file)
        cols = st.columns(len(summary))
        for col, (_, row) in zip(cols, summary.iterrows()):
            col.metric(
                row["pathway"],
                f"{row['significant_genes']}/{row['genes_detected']} significant",
            )
            col.caption(row["status"])
    else:
        missing(targeted_summary_file, "03_pathway_analysis.py")

    if targeted_genes_file.exists():
        targeted = pd.read_csv(targeted_genes_file)
        st.dataframe(
            targeted.sort_values("padj", na_position="last"),
            width='stretch',
            hide_index=True,
        )
    else:
        missing(targeted_genes_file, "03_pathway_analysis.py")

    st.divider()
    st.subheader("GSEA (pre-ranked)")

    if gsea_file.exists():
        provenance_caption(gsea_file)
        gsea = pd.read_csv(gsea_file)
        nes_col = "NES" if "NES" in gsea.columns else None
        fdr_col = "FDR q-val" if "FDR q-val" in gsea.columns else None
        term_col = "Term" if "Term" in gsea.columns else gsea.columns[0]

        if nes_col:
            top_terms = gsea.reindex(
                gsea[nes_col].abs().sort_values(ascending=False).index
            ).head(20)
            fig = px.bar(
                top_terms,
                x=nes_col,
                y=term_col,
                orientation="h",
                color=nes_col,
                color_continuous_scale="RdBu_r",
                hover_data=[fdr_col] if fdr_col else None,
            )
            fig.update_layout(yaxis={"categoryorder": "total ascending"})
            st.plotly_chart(fig, width='stretch')

        st.dataframe(gsea, width='stretch', hide_index=True)
    else:
        missing(gsea_file, "03_pathway_analysis.py")

    st.divider()
    st.subheader("Over-representation (Enrichr)")

    enrich_up_tab, enrich_down_tab = st.tabs(["Upregulated genes", "Downregulated genes"])

    def render_enrichment_tab(file_path: Path) -> None:
        if not file_path.exists():
            missing(file_path, "03_pathway_analysis.py")
            return
        provenance_caption(file_path)
        enr = pd.read_csv(file_path)
        score_col = "Combined Score" if "Combined Score" in enr.columns else None
        term_col = "Term" if "Term" in enr.columns else enr.columns[0]
        if score_col:
            top = enr.sort_values(score_col, ascending=False).head(15)
            fig = px.bar(top, x=score_col, y=term_col, orientation="h")
            fig.update_layout(yaxis={"categoryorder": "total ascending"})
            st.plotly_chart(fig, width='stretch')
        st.dataframe(enr, width='stretch', hide_index=True)

    with enrich_up_tab:
        render_enrichment_tab(enriched_up_file)
    with enrich_down_tab:
        render_enrichment_tab(enriched_down_file)


# ============================================================
# Tab 4: Gene Prioritization
# ============================================================
with tab_priority:
    priority_file = PRIORITY_DIR / "gene_priority.csv"

    if priority_file.exists():
        provenance_caption(priority_file)
        priority = pd.read_csv(priority_file, index_col=0)

        tier_counts = priority["tier"].value_counts()
        c1, c2, c3 = st.columns(3)
        c1.metric("Tier 1 (strong)", int(tier_counts.get("Tier 1 - Strong candidate", 0)))
        c2.metric("Tier 2 (moderate)", int(tier_counts.get("Tier 2 - Moderate candidate", 0)))
        c3.metric("Tier 3 (lower)", int(tier_counts.get("Tier 3 - Lower priority", 0)))

        top_n = st.slider("Show top N genes", 5, 50, 20)
        top = priority.sort_values("priority_score", ascending=False).head(top_n)

        label_col = "Gene name" if "Gene name" in top.columns else top.index
        fig = go.Figure(
            go.Bar(
                x=top["priority_score"][::-1],
                y=(top[label_col] if isinstance(label_col, pd.Series) else top.index)[::-1],
                orientation="h",
                marker_color=top["tier"][::-1].map(
                    {
                        "Tier 1 - Strong candidate": "crimson",
                        "Tier 2 - Moderate candidate": "orange",
                        "Tier 3 - Lower priority": "gray",
                    }
                ),
            )
        )
        fig.update_layout(xaxis_title="Priority score", height=max(400, top_n * 25))
        st.plotly_chart(fig, width='stretch')

        st.subheader("Ranked table")
        st.dataframe(top, width='stretch')
        st.caption(
            "priority_score = 0.40 x confidence + 0.35 x effect size + 0.25 x reliability"
        )
    else:
        missing(priority_file, "03_gene_prioritization.py")