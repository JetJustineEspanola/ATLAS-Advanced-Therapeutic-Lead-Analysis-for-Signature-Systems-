"""
Stage 1: Signature Discovery
Tabs: Validation (QC) | Differential Expression | Gene Prioritization

Reads only from results/ CSVs already produced by
01_validation.py, 02_differential_expression.py, 03_gene_prioritization.py.
No computation happens in this app - it is a thin display layer.
"""

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
PRIORITY_DIR = PROJECT_ROOT / "results" / "gene_prioritization"

st.title("Stage 1 - Signature Discovery")

tab_validation, tab_deg, tab_priority = st.tabs(
    ["Validation (QC)", "Differential Expression", "Gene Prioritization"]
)


def missing(file_path: Path, run_script: str) -> None:
    st.info(f"`{file_path.name}` not found yet - run `{run_script}` first.")


# ============================================================
# Tab 1: Validation / QC
# ============================================================
with tab_validation:
    summary_file = QC_DIR / "stage1_validation_summary.csv"
    library_file = QC_DIR / "library_sizes.csv"
    corr_file = QC_DIR / "sample_correlation.csv"
    pca_file = QC_DIR / "pca_coordinates.csv"

    if summary_file.exists():
        summary = pd.read_csv(summary_file)
        st.dataframe(summary, use_container_width=True, hide_index=True)
    else:
        missing(summary_file, "01_validation.py")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Library size per sample")
        if library_file.exists():
            lib = pd.read_csv(library_file, index_col=0)
            fig = px.bar(
                lib,
                x=lib.index,
                y=lib.columns[0],
                labels={"x": "Sample", lib.columns[0]: "Total estimated counts"},
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            missing(library_file, "01_validation.py")

    with col2:
        st.subheader("Sample-to-sample correlation")
        if corr_file.exists():
            corr = pd.read_csv(corr_file, index_col=0)
            fig = px.imshow(
                corr,
                text_auto=".3f",
                color_continuous_scale="Blues",
                zmin=0,
                zmax=1,
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            missing(corr_file, "01_validation.py")

    st.subheader("PCA")
    if pca_file.exists():
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
        st.plotly_chart(fig, use_container_width=True)
    else:
        missing(pca_file, "01_validation.py")


# ============================================================
# Tab 2: Differential Expression
# ============================================================
with tab_deg:
    all_degs_file = DE_DIR / "DEGs_resistant_vs_sensitive_annotated.csv"
    sig_degs_file = DE_DIR / "significant_DEGs_annotated.csv"

    if all_degs_file.exists():
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
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("DEG table")
        search = st.text_input("Search gene name")
        table = degs.reset_index()
        if search and "Gene name" in table.columns:
            table = table[table["Gene name"].str.contains(search, case=False, na=False)]
        st.dataframe(
            table.sort_values("padj", na_position="last"),
            use_container_width=True,
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
# Tab 3: Gene Prioritization
# ============================================================
with tab_priority:
    priority_file = PRIORITY_DIR / "gene_priority.csv"

    if priority_file.exists():
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
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("Ranked table")
        st.dataframe(top, use_container_width=True)
        st.caption(
            "priority_score = 0.40 x confidence + 0.35 x effect size + 0.25 x reliability"
        )
    else:
        missing(priority_file, "03_gene_prioritization.py")
