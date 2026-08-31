
from __future__ import annotations

import html

import pandas as pd
import streamlit as st


def _safe(value) -> str:
    if value is None:
        return "—"
    try:
        if pd.isna(value):
            return "—"
    except Exception:
        pass
    text = str(value).strip()
    return html.escape(text) if text else "—"


def badge(label: str) -> str:
    label = _safe(label)
    return f'<span class="pill">{label}</span>'


def section_header(title: str, subtitle: str | None = None):
    st.markdown(f'<div class="atlas-section">{html.escape(title)}</div>', unsafe_allow_html=True)
    if subtitle:
        st.markdown(f'<div class="atlas-small">{html.escape(subtitle)}</div>', unsafe_allow_html=True)


def metric_card(label: str, value: str, caption: str = ""):
    st.markdown(
        f"""
        <div class="atlas-card">
            <div class="atlas-kpi">{_safe(value)}</div>
            <div class="atlas-label">{html.escape(label)}</div>
            <div class="atlas-small" style="margin-top:.35rem;">{html.escape(caption)}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def info_box(title: str, body: str):
    st.markdown(
        f"""
        <div class="atlas-card">
        <b>{html.escape(title)}</b>
        <div class="atlas-small" style="margin-top:.35rem;">{html.escape(body)}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def evidence_bar(name: str, level: str):
    st.markdown(
        f"""
        <div class="evidence-row">
            <span class="evidence-name">{html.escape(name)}</span>
            {badge(level)}
        </div>
        """,
        unsafe_allow_html=True,
    )


def evidence_cell(level: str) -> str:
    return _safe(level)


def render_candidate_decision(row: pd.Series):
    action = _safe(row.get("recommended_action"))
    reason = _safe(row.get("decision_reason"))
    st.markdown(
        f"""
        <div class="decision">
        <b>ATLAS recommendation: {action}</b><br>
        <span class="atlas-muted">{reason}</span>
        </div>
        """,
        unsafe_allow_html=True,
    )


def render_evidence_trace(row: pd.Series):
    st.markdown("### Why is this candidate here?")

    steps = [
        ("1. CMap reversal", row.get("cmap_evidence"), "Does the compound oppose the resistance signature?"),
        ("2. Safety screen", row.get("safety_evidence"), "Are there major preliminary safety or assay-liability flags?"),
        ("3. Target evidence", row.get("target_evidence"), "Are pharmacological targets supported by PubChem/ChEMBL evidence?"),
        ("4. Resistance network", row.get("network_evidence"), "Do those targets connect to ATLAS resistance genes?"),
        ("5. Regulatory evidence", row.get("regulatory_evidence"), "What regulatory/clinical evidence exists?"),
    ]

    for title, level, question in steps:
        c1, c2 = st.columns([3, 1])
        with c1:
            st.markdown(f"**{title}**")
            st.caption(question)
        with c2:
            st.markdown(badge(level), unsafe_allow_html=True)
        st.divider()

    st.markdown(
        """
        **Interpretation rule:** ATLAS does not treat any single layer as proof.
        A candidate becomes more compelling when independent evidence layers agree.
        """
    )
