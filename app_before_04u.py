
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import streamlit as st

from ui.data import (
    ATLASPaths,
    build_candidate_table,
    build_pipeline_status,
    load_external_validation,
    load_resistance_summary,
)
from ui.components import (
    badge,
    evidence_bar,
    evidence_cell,
    info_box,
    metric_card,
    render_candidate_decision,
    render_evidence_trace,
    section_header,
)



def clean_value(value):
    if value is None:
        return "—"
    try:
        if pd.isna(value):
            return "—"
    except Exception:
        pass
    text = str(value).strip()
    return text if text else "—"


def numeric_value(value):
    try:
        if value is None or pd.isna(value):
            return "—"
        x = float(value)
        if abs(x) >= 100:
            return f"{x:,.0f}"
        if abs(x) >= 10:
            return f"{x:.1f}"
        return f"{x:.3g}"
    except Exception:
        return clean_value(value)


st.set_page_config(
    page_title="ATLAS",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

PROJECT_ROOT = Path(__file__).resolve().parent
paths = ATLASPaths(PROJECT_ROOT)

st.markdown(
    """
    <style>
    .block-container {padding-top: 1.4rem; padding-bottom: 4rem; max-width: 1480px;}
    [data-testid="stSidebar"] {min-width: 250px; max-width: 250px;}
    .atlas-title {font-size: 2.2rem; font-weight: 780; letter-spacing: -0.03em; margin-bottom: .15rem;}
    .atlas-subtitle {color: #6b7280; margin-bottom: 1rem;}
    .atlas-card {
        border: 1px solid rgba(120,120,120,.20);
        border-radius: 16px;
        padding: 1rem 1.1rem;
        background: rgba(127,127,127,.035);
        margin-bottom: .75rem;
    }
    .atlas-small {font-size: .84rem; color: #6b7280;}
    .atlas-kpi {font-size: 1.65rem; font-weight: 760; line-height: 1.1;}
    .atlas-label {font-size: .82rem; color: #6b7280; margin-top: .35rem;}
    .atlas-section {font-size: 1.35rem; font-weight: 760; margin-top: .2rem; margin-bottom: .65rem;}
    .atlas-muted {color:#6b7280;}
    .decision {
        border-left: 5px solid #8b5cf6;
        border-radius: 12px;
        padding: .9rem 1rem;
        background: rgba(139,92,246,.08);
        margin: .5rem 0 1rem 0;
    }
    .evidence-row {
        display:flex; align-items:center; justify-content:space-between;
        gap:1rem; padding:.35rem 0;
        border-bottom:1px solid rgba(120,120,120,.12);
    }
    .evidence-name {font-weight:600;}
    .pill {
        display:inline-block; border-radius:999px; padding:.2rem .55rem;
        font-size:.75rem; font-weight:700; border:1px solid rgba(120,120,120,.22);
    }
    </style>
    """,
    unsafe_allow_html=True,
)


@st.cache_data(show_spinner=False)
def candidate_data() -> pd.DataFrame:
    return build_candidate_table(paths)


@st.cache_data(show_spinner=False)
def pipeline_data() -> pd.DataFrame:
    return build_pipeline_status(paths)


@st.cache_data(show_spinner=False)
def external_data():
    return load_external_validation(paths)


@st.cache_data(show_spinner=False)
def biology_data():
    return load_resistance_summary(paths)


candidates = candidate_data()
pipeline = pipeline_data()
external = external_data()
biology = biology_data()


# -------------------------------------------------------------------
# Sidebar
# -------------------------------------------------------------------

with st.sidebar:
    st.markdown("## ATLAS")
    st.caption("Trastuzumab Resistance Decision Support")

    page = st.radio(
        "Navigate",
        [
            "Overview",
            "Candidates",
            "Resistance Biology",
            "Evidence Matrix",
            "Pipeline",
            "Advanced Analysis",
        ],
        label_visibility="collapsed",
    )

    st.divider()

    if not pipeline.empty:
        complete_n = int((pipeline["status"] == "Complete").sum())
        total_n = len(pipeline)
        st.caption("Pipeline progress")
        st.progress(complete_n / total_n if total_n else 0)
        st.caption(f"{complete_n}/{total_n} major stages complete")

    if st.button("Refresh data", use_container_width=True):
        st.cache_data.clear()
        st.rerun()

    st.caption("UI reads pipeline outputs dynamically. No analysis results are hardcoded.")


# -------------------------------------------------------------------
# Header
# -------------------------------------------------------------------

st.markdown('<div class="atlas-title">ATLAS</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="atlas-subtitle">From resistance biology to testable drug candidates</div>',
    unsafe_allow_html=True,
)


# -------------------------------------------------------------------
# Overview
# -------------------------------------------------------------------

if page == "Overview":
    section_header("What matters now", "A decision-first summary of the current ATLAS analysis.")

    deg_n = biology.get("significant_deg_n")
    top_pathway = biology.get("top_pathway_name") or "Pending"
    top_pathway_fdr = biology.get("top_pathway_fdr")
    cmap_n = int((candidates["cmap_evidence"] == "STRONG").sum()) if not candidates.empty else 0
    shortlist_n = int(candidates["recommended_action"].isin(["PRIORITIZE", "PROCEED TO DOCKING"]).sum()) if not candidates.empty else 0

    cols = st.columns(4)
    with cols[0]:
        metric_card("Resistance DEGs", f"{deg_n:,}" if isinstance(deg_n, int) else "—", "Discovery")
    with cols[1]:
        metric_card("Lead pathway", top_pathway, "Resistance biology")
    with cols[2]:
        metric_card("Strong CMap hits", str(cmap_n), "Transcriptomic reversal")
    with cols[3]:
        metric_card("Current shortlist", str(shortlist_n), "Priority decisions")

    st.markdown("")

    left, right = st.columns([1.45, 1])

    with left:
        section_header("Top candidates", "Ranked by the evidence available right now.")

        if candidates.empty:
            info_box(
                "No candidate results found yet",
                "ATLAS will populate this view as 04L/04M/04O/04P/04Q outputs become available.",
            )
        else:
            shown = candidates.head(8)
            for _, row in shown.iterrows():
                c1, c2, c3 = st.columns([2.3, 1, 1.2])
                with c1:
                    st.markdown(f"**#{int(row['ui_rank'])} · {row['pert_iname']}**")
                    reason = row.get("decision_reason") or "Evidence is still accumulating."
                    st.caption(reason)
                with c2:
                    st.markdown(badge(row["overall_evidence"]), unsafe_allow_html=True)
                    st.caption("Evidence")
                with c3:
                    st.markdown(badge(row["recommended_action"]), unsafe_allow_html=True)
                    st.caption("Decision")
                st.divider()

    with right:
        section_header("Resistance story", "The simplest interpretation supported by the current data.")
        if biology.get("top_pathway_name"):
            fdr_text = (
                f"{biology['top_pathway_fdr']:.4g}"
                if isinstance(biology.get("top_pathway_fdr"), (int, float))
                and not pd.isna(biology.get("top_pathway_fdr"))
                else "n/a"
            )
            st.markdown(
                f"""
                <div class="atlas-card">
                <div class="atlas-small">Primary resistance-associated pathway</div>
                <div style="font-size:1.35rem;font-weight:760;margin:.2rem 0 .45rem 0;">
                {biology['top_pathway_name']}
                </div>
                <div><b>FDR:</b> {fdr_text}</div>
                <div class="atlas-small" style="margin-top:.65rem;">
                Association with resistance does not establish causality.
                </div>
                </div>
                """,
                unsafe_allow_html=True,
            )
        else:
            info_box("Pathway results pending", "The UI will summarize Stage 03 outputs when present.")

        if external and external.get("dataset_count") is not None:
            st.markdown(
                f"""
                <div class="atlas-card">
                <div class="atlas-small">External validation</div>
                <div style="font-size:1.35rem;font-weight:760;margin:.2rem 0 .45rem 0;">
                {external.get('successful_count', 0)} / {external.get('dataset_count', 0)} datasets analyzed
                </div>
                <div class="atlas-small">
                Statistical and directional evidence are kept separate when replication is limited.
                </div>
                </div>
                """,
                unsafe_allow_html=True,
            )

    section_header("What should happen next", "ATLAS translates evidence into an explicit next action.")

    if candidates.empty:
        st.info("Run the downstream candidate stages and refresh this page.")
    else:
        action_counts = candidates["recommended_action"].value_counts()
        a, b, c, d = st.columns(4)
        with a:
            metric_card("Prioritize", str(int(action_counts.get("PRIORITIZE", 0))), "Keep at the front")
        with b:
            metric_card("Dock next", str(int(action_counts.get("PROCEED TO DOCKING", 0))), "Mechanistically supported")
        with c:
            metric_card("Manual review", str(int(action_counts.get("MANUAL REVIEW", 0))), "Needs expert judgment")
        with d:
            metric_card("Deprioritize", str(int(action_counts.get("DEPRIORITIZE", 0))), "Weak/high-risk evidence")


# -------------------------------------------------------------------
# Candidates
# -------------------------------------------------------------------

elif page == "Candidates":
    section_header(
        "Candidate Explorer",
        "Choose a compound and see the evidence that affects the decision — not every chart at once.",
    )

    if candidates.empty:
        st.warning("No candidate-stage outputs are available yet.")
        st.stop()

    c1, c2, c3 = st.columns([2.4, 1, 1])
    with c1:
        selected_name = st.selectbox(
            "Candidate",
            candidates["pert_iname"].tolist(),
            index=0,
        )
    with c2:
        tier_filter = st.selectbox("CMap tier", ["All"] + sorted(candidates["priority_tier"].dropna().astype(str).unique().tolist()))
    with c3:
        action_filter = st.selectbox(
            "Decision",
            ["All"] + sorted(candidates["recommended_action"].dropna().astype(str).unique().tolist()),
        )

    row = candidates[candidates["pert_iname"] == selected_name].iloc[0]

    top1, top2, top3 = st.columns([1.5, 1, 1])
    with top1:
        st.markdown(f"## {row['pert_iname']}")
        st.caption(f"Current UI rank #{int(row['ui_rank'])}")
    with top2:
        st.markdown(badge(row["overall_evidence"]), unsafe_allow_html=True)
        st.caption("Overall evidence")
    with top3:
        st.markdown(badge(row["recommended_action"]), unsafe_allow_html=True)
        st.caption("Recommended action")

    render_candidate_decision(row)

    tabs = st.tabs(["Overview", "Transcriptomics", "Safety", "Targets", "Network", "Regulatory", "Evidence trace"])

    with tabs[0]:
        section_header("Evidence profile")
        evidence_bar("CMap reversal", row["cmap_evidence"])
        evidence_bar("External validation", row["external_evidence"])
        evidence_bar("Safety", row["safety_evidence"])
        evidence_bar("Target evidence", row["target_evidence"])
        evidence_bar("Network support", row["network_evidence"])
        evidence_bar("Regulatory", row["regulatory_evidence"])

    with tabs[1]:
        section_header("Transcriptomic reversal")
        a, b, c = st.columns(3)
        with a:
            metric_card("CMap tier", clean_value(row.get("priority_tier")), "Consensus tier")
        with b:
            metric_card("Mean tau", numeric_value(row.get("mean_tau")), "Negative = opposing signature")
        with c:
            metric_card("Strong negatives", numeric_value(row.get("n_strong_negative")), "Across queried signatures")
        st.caption(
            "Negative CMap tau indicates transcriptional opposition to the resistance signature; it is not proof of efficacy or resistance reversal."
        )

    with tabs[2]:
        section_header("Preliminary safety screen")
        a, b, c = st.columns(3)
        with a:
            metric_card("Safety decision", clean_value(row.get("safety_screening_recommendation")), "04O")
        with b:
            metric_card("Risk score", numeric_value(row.get("safety_risk_score")), "Rule-based screening")
        with c:
            metric_card("PAINS", clean_value(row.get("pains_flag")), "Assay-interference alert")
        reason = clean_value(row.get("safety_screening_reasons"))
        if reason != "—":
            st.info(reason)
        st.caption("PAINS and PubChem safety signals are screening evidence, not clinical toxicity conclusions.")

    with tabs[3]:
        section_header("Target evidence")
        a, b, c = st.columns(3)
        with a:
            metric_card("Target support", clean_value(row.get("target_support_category")), "04P")
        with b:
            metric_card("ChEMBL targets", numeric_value(row.get("chembl_target_n")), "Activity-supported")
        with c:
            metric_card("Max pChEMBL", numeric_value(row.get("chembl_max_pchembl")), "Activity strength")
        targets = clean_value(row.get("chembl_target_ids"))
        if targets != "—":
            st.code(targets, language=None)

    with tabs[4]:
        section_header("Resistance-network support")
        a, b, c = st.columns(3)
        with a:
            metric_card("Best target", clean_value(row.get("best_network_target")), "04Q")
        with b:
            metric_card("Network support", clean_value(row.get("best_network_support_category")), "STRING")
        with c:
            metric_card("Linked genes", numeric_value(row.get("best_target_linked_resistance_gene_n")), "Resistance genes")
        linked = clean_value(row.get("best_target_linked_resistance_genes"))
        if linked != "—":
            st.code(linked, language=None)
        st.caption("STRING edges are functional associations and may not represent direct physical binding.")

    with tabs[5]:
        section_header("Regulatory / clinical evidence")
        a, b, c = st.columns(3)
        with a:
            metric_card("Regulatory evidence", clean_value(row.get("regulatory_evidence_category")), "04N")
        with b:
            metric_card("FDA application", clean_value(row.get("fda_application_record_found")), "Evidence only")
        with c:
            metric_card("Clinical trial", clean_value(row.get("clinical_trial_record_found")), "Evidence only")
        st.caption("A regulatory database hit does not automatically mean FDA approval for HER2+ breast cancer or trastuzumab resistance.")

    with tabs[6]:
        render_evidence_trace(row)


# -------------------------------------------------------------------
# Resistance Biology
# -------------------------------------------------------------------

elif page == "Resistance Biology":
    section_header(
        "Resistance Biology",
        "Explain the biological story first; expose detailed plots only when they answer a specific question.",
    )

    if not biology:
        st.warning("No discovery/pathway outputs found.")
        st.stop()

    c1, c2, c3 = st.columns(3)
    with c1:
        metric_card(
            "Significant DEGs",
            f"{biology.get('significant_deg_n'):,}" if isinstance(biology.get("significant_deg_n"), int) else "—",
            "|log2FC| ≥ 1, FDR < 0.05",
        )
    with c2:
        metric_card("Top pathway", biology.get("top_pathway_name") or "—", "Strongest detected pathway")
    with c3:
        metric_card(
            "Pathway FDR",
            numeric_value(biology.get("top_pathway_fdr")),
            "Multiple-testing adjusted",
        )

    if biology.get("top_pathway_name"):
        st.markdown(
            f"""
            <div class="decision">
            <b>Interpretation</b><br>
            <span class="atlas-muted">
            {biology.get('interpretation', 'This pathway is associated with the resistant transcriptional state.')}
            </span>
            </div>
            """,
            unsafe_allow_html=True,
        )

    section_header("Key resistance-associated genes")
    key_genes = biology.get("key_genes", [])
    if key_genes:
        gene_df = pd.DataFrame(key_genes)
        st.dataframe(
            gene_df,
            use_container_width=True,
            hide_index=True,
        )
    else:
        st.caption("Key-gene summary will appear when DEG/pathway outputs can be resolved.")

    if external:
        section_header("External validation")
        ext_df = external.get("summary")
        if isinstance(ext_df, pd.DataFrame) and not ext_df.empty:
            cols = [
                c for c in [
                    "dataset",
                    "validation_evidence_type",
                    "overlapping_genes",
                    "direction_agreement_fraction",
                    "informative_effect_direction_agreement_fraction",
                    "atlas_significant_direction_agreement_fraction",
                    "both_significant_direction_agreement_fraction",
                ]
                if c in ext_df.columns
            ]
            st.dataframe(ext_df[cols], use_container_width=True, hide_index=True)
            st.caption(
                "Whole-transcriptome direction agreement can be noisy; informative-effect and ATLAS-significant subsets are more decision-relevant."
            )


# -------------------------------------------------------------------
# Evidence Matrix
# -------------------------------------------------------------------

elif page == "Evidence Matrix":
    section_header(
        "Integrated Evidence Matrix",
        "Scan across evidence layers to see where a candidate is strong, weak, risky, or still pending.",
    )

    if candidates.empty:
        st.warning("No candidate data available.")
        st.stop()

    top_n = st.slider("Candidates shown", 5, min(50, len(candidates)), min(20, len(candidates)))

    matrix = candidates.head(top_n)[
        [
            "ui_rank",
            "pert_iname",
            "cmap_evidence",
            "external_evidence",
            "safety_evidence",
            "target_evidence",
            "network_evidence",
            "regulatory_evidence",
            "recommended_action",
        ]
    ].copy()

    matrix.columns = [
        "Rank",
        "Candidate",
        "CMap",
        "External",
        "Safety",
        "Target",
        "Network",
        "Regulatory",
        "Decision",
    ]

    st.dataframe(
        matrix,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Rank": st.column_config.NumberColumn(width="small"),
            "Candidate": st.column_config.TextColumn(width="medium"),
        },
    )

    st.caption(
        "This matrix intentionally avoids a fake probability score. Each column preserves the meaning of its own evidence layer."
    )


# -------------------------------------------------------------------
# Pipeline
# -------------------------------------------------------------------

elif page == "Pipeline":
    section_header(
        "Research Pipeline",
        "Pipeline stage names belong here as project status — not as the primary user navigation.",
    )

    if pipeline.empty:
        st.warning("Pipeline status could not be determined.")
        st.stop()

    for _, r in pipeline.iterrows():
        c1, c2, c3 = st.columns([1.5, 2.8, 1])
        with c1:
            st.markdown(f"**{r['stage']}**")
        with c2:
            st.write(r["label"])
            st.caption(r["detail"])
        with c3:
            st.markdown(badge(r["status"]), unsafe_allow_html=True)
        st.divider()


# -------------------------------------------------------------------
# Advanced
# -------------------------------------------------------------------

elif page == "Advanced Analysis":
    section_header(
        "Advanced Analysis",
        "Raw tables and detailed outputs live here so they do not overwhelm the decision workflow.",
    )

    dataset_options = paths.available_output_tables()

    if not dataset_options:
        st.info("No pipeline output tables were found.")
        st.stop()

    choice = st.selectbox("Output table", list(dataset_options.keys()))
    selected_path = dataset_options[choice]

    st.caption(str(selected_path.relative_to(PROJECT_ROOT)))

    try:
        df = pd.read_csv(selected_path)
        st.write(f"{len(df):,} rows × {len(df.columns):,} columns")
        st.dataframe(df, use_container_width=True, hide_index=True)
    except Exception as exc:
        st.error(f"Could not read this output: {exc}")

