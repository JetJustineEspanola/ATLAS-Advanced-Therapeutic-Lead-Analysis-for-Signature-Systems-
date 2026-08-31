from pathlib import Path
import streamlit as st
import pandas as pd

from ui.data import (
    PROJECT_ROOT,
    load_all_data,
    safe_col,
    numeric,
    clean,
)
from ui.components import (
    page_header,
    metric_card,
    evidence_badge,
    decision_badge,
    candidate_header,
    evidence_row,
    empty_state,
)

st.set_page_config(
    page_title="ATLAS",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    .block-container { padding-top: 1.4rem; padding-bottom: 3rem; max-width: 1500px; }
    [data-testid="stSidebar"] { min-width: 250px; max-width: 250px; }
    .atlas-kicker { color:#888; font-size:.78rem; letter-spacing:.08em; text-transform:uppercase; margin-bottom:.25rem; }
    .atlas-title { font-size:2rem; font-weight:750; margin-bottom:.15rem; }
    .atlas-subtitle { color:#888; margin-bottom:1.1rem; }
    .atlas-card { border:1px solid rgba(128,128,128,.22); border-radius:14px; padding:1rem 1.05rem; margin-bottom:.8rem; }
    .atlas-metric { font-size:1.65rem; font-weight:750; line-height:1.1; }
    .atlas-label { color:#888; font-size:.8rem; margin-top:.3rem; }
    .atlas-badge { display:inline-block; padding:.22rem .48rem; border-radius:999px; border:1px solid rgba(128,128,128,.28); font-size:.75rem; font-weight:650; margin-right:.3rem; margin-bottom:.25rem; }
    .atlas-muted { color:#888; }
    .atlas-divider { height:1px; background:rgba(128,128,128,.18); margin:1rem 0; }
    </style>
    """,
    unsafe_allow_html=True,
)

data = load_all_data()
matrix = data["matrix"]
shortlist = data["shortlist"]
summary = data["summary"]
metadata = data["metadata"]
docking = data["docking"]
biology = data["biology"]

with st.sidebar:
    st.markdown("## ATLAS")
    st.caption("HER2+ trastuzumab resistance")
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
    if not matrix.empty:
        st.caption(f"{len(matrix)} integrated candidates")
    st.caption("Decision-support interface")
    st.caption("Computational evidence ≠ clinical efficacy")

if page == "Overview":
    page_header(
        "ATLAS",
        "Integrated computational prioritization for trastuzumab resistance in HER2-positive breast cancer.",
        kicker="Decision Support",
    )

    if matrix.empty:
        empty_state(
            "04U outputs not found",
            "Run Stage 04U first. The UI expects results/cmap/integrated_evidence/ATLAS_integrated_evidence_matrix.csv.",
        )
        st.stop()

    top = matrix.sort_values("04u_rank").iloc[0]
    n_priority = int((matrix["final_evidence_category"] == "PRIORITY_EXPERIMENTAL_VALIDATION").sum())
    n_caution = int((matrix["final_evidence_category"] == "EXPERIMENTAL_VALIDATION_WITH_CAUTION").sum())
    n_secondary = int((matrix["final_evidence_category"] == "SECONDARY_EXPERIMENTAL_CANDIDATE").sum())

    c1, c2, c3, c4 = st.columns(4)
    with c1:
        metric_card("Top candidate", clean(top.get("pert_iname")) or "—")
    with c2:
        metric_card("Priority validation", n_priority)
    with c3:
        metric_card("Validation with caution", n_caution)
    with c4:
        metric_card("Secondary candidates", n_secondary)

    st.markdown("### Recommended next action")

    if not shortlist.empty:
        lead = shortlist.sort_values("experimental_shortlist_rank").iloc[0]
        with st.container(border=True):
            left, right = st.columns([2, 1])
            with left:
                candidate_header(
                    clean(lead.get("pert_iname")),
                    clean(lead.get("final_evidence_category")),
                    f"Integrated rank #{int(numeric(lead.get('04u_rank'), 1))}",
                )
                st.write(clean(lead.get("final_evidence_reason")))
                evidence_row(lead)
            with right:
                st.metric(
                    "Priority score",
                    f"{numeric(lead.get('experimental_priority_score'), 0):.3f}",
                )
                st.caption("Transparent ranking score, not probability.")
                target = clean(lead.get("best_network_target"))
                if target:
                    st.metric("Network target", target)

    st.markdown("### Experimental validation shortlist")
    if shortlist.empty:
        st.info("No experimental validation shortlist is available.")
    else:
        cols = [
            "experimental_shortlist_rank",
            "pert_iname",
            "best_network_target",
            "experimental_priority_score",
            "final_evidence_category",
        ]
        show = shortlist[[c for c in cols if c in shortlist.columns]].copy()
        st.dataframe(show, use_container_width=True, hide_index=True)

    st.markdown("### Evidence distribution")
    counts = (
        matrix["final_evidence_category"]
        .value_counts()
        .rename_axis("category")
        .reset_index(name="candidates")
    )
    st.dataframe(counts, use_container_width=True, hide_index=True)

elif page == "Candidates":
    page_header(
        "Candidates",
        "Inspect one candidate at a time and trace the evidence behind its recommendation.",
        kicker="Decision Review",
    )

    if matrix.empty:
        empty_state("No integrated candidates", "Stage 04U output is missing.")
        st.stop()

    matrix = matrix.sort_values(["04u_rank", "priority_rank"], na_position="last")
    names = matrix["pert_iname"].dropna().astype(str).tolist()
    selected = st.selectbox("Candidate", names, index=0)
    row = matrix[matrix["pert_iname"].astype(str) == selected].iloc[0]

    candidate_header(
        selected,
        clean(row.get("final_evidence_category")),
        f"04U rank #{int(numeric(row.get('04u_rank'), 0))}",
    )

    c1, c2, c3, c4 = st.columns(4)
    with c1:
        metric_card("Priority score", f"{numeric(row.get('experimental_priority_score'), 0):.3f}")
    with c2:
        metric_card("CMap", clean(row.get("cmap_evidence_badge")) or "—")
    with c3:
        metric_card("Safety", clean(row.get("safety_evidence_badge")) or "—")
    with c4:
        metric_card("Docking", clean(row.get("docking_evidence_badge")) or "—")

    st.markdown("### Evidence trace")
    evidence_row(row)

    a, b = st.columns(2)
    with a:
        st.markdown("#### Biological support")
        st.write(f"**Target support:** {clean(row.get('target_evidence_badge')) or '—'}")
        st.write(f"**Network support:** {clean(row.get('network_evidence_badge')) or '—'}")
        st.write(f"**Best network target:** {clean(row.get('best_network_target')) or '—'}")
        st.write(f"**UniProt:** {clean(row.get('validated_uniprot_accession')) or '—'}")
        st.write(
            f"**Linked resistance genes:** {int(numeric(row.get('best_target_linked_resistance_gene_n'), 0))}"
        )

    with b:
        st.markdown("#### Safety / developability")
        st.write(f"**Safety recommendation:** {clean(row.get('safety_screening_recommendation')) or '—'}")
        st.write(f"**Developability:** {clean(row.get('structural_developability_category')) or '—'}")
        st.write(f"**QED:** {numeric(row.get('qed'), float('nan')):.3f}" if pd.notna(numeric(row.get('qed'), float('nan'))) else "**QED:** —")
        st.write(f"**ESOL logS:** {numeric(row.get('esol_logS'), float('nan')):.3f}" if pd.notna(numeric(row.get('esol_logS'), float('nan'))) else "**ESOL logS:** —")
        st.write(f"**PAINS alerts:** {int(numeric(row.get('pains_alert_n_04t'), 0))}")
        st.write(f"**Brenk alerts:** {int(numeric(row.get('brenk_alert_n'), 0))}")

    st.markdown("#### Docking")
    if clean(row.get("docking_status")) == "COMPLETED":
        d1, d2, d3 = st.columns(3)
        with d1:
            st.metric("Best affinity", f"{numeric(row.get('best_affinity_kcal_mol'), 0):.3f} kcal/mol")
        with d2:
            rmsd = numeric(row.get("reference_redocking_rmsd_A"), float("nan"))
            st.metric("Redocking RMSD", f"{rmsd:.3f} Å" if pd.notna(rmsd) else "Unavailable")
        with d3:
            st.metric("Protocol", clean(row.get("protocol_validation")) or "—")
    else:
        st.info(clean(row.get("docking_status")) or "No validated docking evidence.")

    st.markdown("#### Interpretation")
    st.write(clean(row.get("final_evidence_reason")) or "No final interpretation available.")
    st.caption(
        "Negative CMap tau, network support, docking, and structural screening are hypothesis-generating evidence and do not prove reversal of trastuzumab resistance."
    )

elif page == "Resistance Biology":
    page_header(
        "Resistance Biology",
        "Core discovery signals that define the resistance context used downstream.",
        kicker="Biology",
    )

    if biology:
        for title, value, note in biology:
            with st.container(border=True):
                st.markdown(f"#### {title}")
                st.markdown(f"**{value}**")
                if note:
                    st.caption(note)
    else:
        st.info("No curated biology summary was found. This page remains intentionally compact.")

elif page == "Evidence Matrix":
    page_header(
        "Evidence Matrix",
        "Compare candidates across the evidence layers used in final prioritization.",
        kicker="Integrated Evidence",
    )

    if matrix.empty:
        empty_state("No matrix available", "Run Stage 04U first.")
        st.stop()

    categories = sorted(matrix["final_evidence_category"].dropna().astype(str).unique())
    selected_categories = st.multiselect(
        "Final evidence category",
        categories,
        default=categories,
    )

    filtered = matrix[
        matrix["final_evidence_category"].isin(selected_categories)
    ].copy()

    display_cols = [
        "04u_rank",
        "pert_iname",
        "cmap_evidence_badge",
        "safety_evidence_badge",
        "target_evidence_badge",
        "network_evidence_badge",
        "developability_evidence_badge",
        "docking_evidence_badge",
        "strong_evidence_layers_n",
        "supportive_evidence_layers_n",
        "caution_layers_n",
        "experimental_priority_score",
        "final_evidence_category",
    ]
    st.dataframe(
        filtered[[c for c in display_cols if c in filtered.columns]],
        use_container_width=True,
        hide_index=True,
    )

elif page == "Pipeline":
    page_header(
        "Pipeline",
        "Stage-level provenance and completion status.",
        kicker="Reproducibility",
    )

    stages = [
        ("01", "Data validation / QC", True),
        ("02", "Differential expression", True),
        ("03", "Pathway analysis", True),
        ("03B", "External dataset validation", False),
        ("03C", "Consensus resistance signature", False),
        ("04K", "CMap parsing", True),
        ("04L", "CMap prioritization", True),
        ("04M", "Compound identity", True),
        ("04N", "Regulatory / clinical evidence", data["regulatory_exists"]),
        ("04O", "Safety / promiscuity / PAINS", data["safety_exists"]),
        ("04P", "Drug-target annotation", data["targets_exists"]),
        ("04Q", "PPI / resistance-network integration", data["network_exists"]),
        ("04R", "Final candidate prioritization", data["final_exists"]),
        ("04S", "Target-supported docking", data["docking_exists"]),
        ("04T", "ADMET / structural developability", data["admet_exists"]),
        ("04U", "Integrated evidence matrix", data["integrated_exists"]),
    ]

    for code, label, complete in stages:
        with st.container(border=True):
            c1, c2, c3 = st.columns([1, 5, 2])
            c1.markdown(f"**{code}**")
            c2.write(label)
            c3.markdown("✅ Complete" if complete else "◻ Pending / optional")

    st.caption(
        "03B/03C are displayed as optional/backfill here because the current downstream pipeline can proceed without them. Update this page once those stages are rerun."
    )

elif page == "Advanced Analysis":
    page_header(
        "Advanced Analysis",
        "Raw tables for audit, debugging, and methods review.",
        kicker="Advanced",
    )

    tabs = st.tabs(["04U Matrix", "Experimental shortlist", "04S Docking", "Metadata"])

    with tabs[0]:
        if matrix.empty:
            st.info("No 04U matrix.")
        else:
            st.dataframe(matrix, use_container_width=True, hide_index=True)

    with tabs[1]:
        if shortlist.empty:
            st.info("No shortlist.")
        else:
            st.dataframe(shortlist, use_container_width=True, hide_index=True)

    with tabs[2]:
        if docking.empty:
            st.info("No docking output.")
        else:
            st.dataframe(docking, use_container_width=True, hide_index=True)

    with tabs[3]:
        if metadata:
            st.json(metadata)
        else:
            st.info("No 04U metadata JSON found.")
