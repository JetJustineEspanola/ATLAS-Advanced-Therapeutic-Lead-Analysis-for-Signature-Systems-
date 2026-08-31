import streamlit as st

def _esc(x):
    return (
        str(x)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )

def page_header(title, subtitle, kicker=""):
    if kicker:
        st.markdown(f'<div class="atlas-kicker">{_esc(kicker)}</div>', unsafe_allow_html=True)
    st.markdown(f'<div class="atlas-title">{_esc(title)}</div>', unsafe_allow_html=True)
    st.markdown(f'<div class="atlas-subtitle">{_esc(subtitle)}</div>', unsafe_allow_html=True)

def metric_card(label, value):
    st.markdown(
        f"""
        <div class="atlas-card">
          <div class="atlas-metric">{_esc(value)}</div>
          <div class="atlas-label">{_esc(label)}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

def evidence_badge(text):
    text = text or "—"
    st.markdown(
        f'<span class="atlas-badge">{_esc(text)}</span>',
        unsafe_allow_html=True,
    )

def decision_badge(text):
    evidence_badge(text)

def candidate_header(name, decision, subtitle=""):
    st.markdown(f"## {name}")
    if subtitle:
        st.caption(subtitle)
    decision_badge(decision)

def evidence_row(row):
    cols = st.columns(6)
    labels = [
        ("CMap", row.get("cmap_evidence_badge", row.get("cmap_label", ""))),
        ("Safety", row.get("safety_evidence_badge", row.get("safety_label", ""))),
        ("Target", row.get("target_evidence_badge", row.get("target_label", ""))),
        ("Network", row.get("network_evidence_badge", row.get("network_label", ""))),
        ("Developability", row.get("developability_evidence_badge", row.get("structural_developability_category", ""))),
        ("Docking", row.get("docking_evidence_badge", row.get("protocol_validation", ""))),
    ]
    for col, (label, value) in zip(cols, labels):
        with col:
            st.caption(label)
            evidence_badge(value)

def empty_state(title, body):
    st.info(f"**{title}**\n\n{body}")
