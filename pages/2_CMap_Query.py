"""
Stage 2: CMap / L1000 Analysis
Tabs: Signature Construction | Cross-Signature Results | Candidate Explorer |
      Pipeline Status | Data Provenance

This page is a pure viewer. It reads only from results/cmap/ and never
hard-codes compound names, tau values, job IDs, signature names, thresholds,
or counts - all of that is discovered from the files at runtime, since the
pipeline (04_cmap_analysis.py -> 04g..04l) gets rerun and those values change
between runs.

Primary source of truth for candidates: results/cmap/prioritized/
ATLAS_CMap_prioritized_all.csv (produced by 04l_cmap_prioritize.py). It is a
strict superset of the earlier parsed/cross_signature_tau.csv - same columns
plus tier/rank/consistency score - and 04l re-derives + audits the summary
stats from the raw tau values rather than trusting them blindly, so this page
reads that file instead of the earlier parsed/ one to avoid ever showing two
disagreeing "top candidate" lists.
"""

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

st.set_page_config(page_title="ATLAS - CMap Analysis", layout="wide")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CMAP_DIR = PROJECT_ROOT / "results" / "cmap"
JOBS_DIR = CMAP_DIR / "jobs"
PARSED_DIR = CMAP_DIR / "parsed"
PRIORITIZED_DIR = CMAP_DIR / "prioritized"

SIGNATURE_MANIFEST_FILE = CMAP_DIR / "signature_manifest.csv"
SIGNATURE_DEFINITIONS_FILE = CMAP_DIR / "signature_definitions.csv"
ENTREZ_SUMMARY_FILE = CMAP_DIR / "cmap_entrez_conversion_summary.csv"
JOB_MANIFEST_FILE = JOBS_DIR / "cmap_job_manifest.json"
JOB_STATUS_FILE = JOBS_DIR / "cmap_job_status_all.json"
DOWNLOAD_MANIFEST_FILE = JOBS_DIR / "cmap_download_manifest.json"
PRIORITIZED_ALL_FILE = PRIORITIZED_DIR / "ATLAS_CMap_prioritized_all.csv"
COMPOUND_TAU_LONG_FILE = PARSED_DIR / "ATLAS_CMap_compound_tau_long.csv"

# Non-signature columns known to exist in prioritized_all.csv - used only as
# a last-resort fallback for detecting signature columns if signature_manifest
# isn't available. See the "Assumptions" note in the chat for why this is the
# one deliberately-fragile fallback path.
KNOWN_NON_SIGNATURE_COLUMNS = {
    "priority_rank", "priority_tier_number", "priority_tier",
    "cross_signature_consistency_score", "pert_id", "pert_iname",
    "n_signatures", "n_negative", "n_strong_negative",
    "mean_tau", "median_tau", "minimum_tau",
    "strong_negative_threshold", "consensus_rank",
}


# ============================================================
# Data loading layer
# ============================================================

@st.cache_data
def load_csv(path_str: str):
    path = Path(path_str)
    if not path.exists():
        return None
    return pd.read_csv(path)


@st.cache_data
def load_json(path_str: str):
    path = Path(path_str)
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def file_info(path: Path) -> dict:
    if not path.exists():
        return {"file": path.name, "exists": False}
    stat = path.stat()
    return {
        "file": path.name,
        "exists": True,
        "path": str(path),
        "modified": datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
        "size_kb": round(stat.st_size / 1024, 1),
    }


def missing(file_path: Path, produced_by: str) -> None:
    st.info(f"`{file_path.name}` not found yet - produced by `{produced_by}`.")


def detect_signature_columns(prioritized: pd.DataFrame, manifest) -> list:
    """
    Which prioritized_all.csv columns are per-signature tau columns.
    Preferred: intersect signature_manifest's signature names with the
    actual columns present (handles signatures that were built but never
    submitted to CMap, so have no tau column).
    Fallback: anything not in the known fixed-column set.
    """
    if manifest is not None and "signature" in manifest.columns:
        detected = [s for s in manifest["signature"].tolist() if s in prioritized.columns]
        if detected:
            return detected
    return [c for c in prioritized.columns if c not in KNOWN_NON_SIGNATURE_COLUMNS]


# ============================================================
# Header
# ============================================================

st.title("Stage 2 - CMap / L1000 Analysis")

st.warning(
    "CMap reversal identifies compounds whose transcriptional signatures "
    "oppose the queried resistance signature. This is a computational "
    "prioritization result and does not establish drug efficacy, clinical "
    "suitability, safety, FDA approval, or mechanism of action."
)

col_spacer, col_refresh = st.columns([6, 1])
with col_refresh:
    if st.button("Refresh data"):
        st.cache_data.clear()
        st.rerun()

manifest = load_csv(str(SIGNATURE_MANIFEST_FILE))
definitions = load_csv(str(SIGNATURE_DEFINITIONS_FILE))
entrez_summary = load_csv(str(ENTREZ_SUMMARY_FILE))
prioritized = load_csv(str(PRIORITIZED_ALL_FILE))
compound_long = load_csv(str(COMPOUND_TAU_LONG_FILE))

# ------------------------------------------------------------
# Top-line overview metrics - all derived, nothing hardcoded
# ------------------------------------------------------------
overview_cols = st.columns(4)

n_signatures_built = len(manifest) if manifest is not None else None
n_compounds = len(prioritized) if prioritized is not None else None
n_tier1 = (
    int((prioritized["priority_tier_number"] == 1).sum())
    if prioritized is not None and "priority_tier_number" in prioritized.columns
    else None
)
n_submitted_signatures = (
    len(detect_signature_columns(prioritized, manifest)) if prioritized is not None else None
)

overview_cols[0].metric("Signatures constructed", n_signatures_built if n_signatures_built is not None else "-")
overview_cols[1].metric("Signatures queried in CMap", n_submitted_signatures if n_submitted_signatures is not None else "-")
overview_cols[2].metric("Compounds scored", n_compounds if n_compounds is not None else "-")
overview_cols[3].metric("Tier 1 (strong consensus)", n_tier1 if n_tier1 is not None else "-")

tab_sig, tab_cross, tab_explorer, tab_jobs, tab_provenance = st.tabs(
    [
        "Signature Construction",
        "Cross-Signature Results",
        "Candidate Explorer",
        "Pipeline Status",
        "Data Provenance",
    ]
)


# ============================================================
# Tab 1: Signature Construction
# ============================================================
with tab_sig:
    if manifest is not None and definitions is not None:
        merged = manifest.merge(definitions, on="signature", how="left")
        display_cols = [c for c in ["signature", "description", "cmap_role", "up_count", "down_count"] if c in merged.columns]
        st.dataframe(merged[display_cols], width="stretch", hide_index=True)
    elif manifest is not None:
        st.dataframe(manifest, width="stretch", hide_index=True)
    else:
        missing(SIGNATURE_MANIFEST_FILE, "04_cmap_analysis.py")

    st.divider()
    st.subheader("Entrez ID mapping (required for CMap submission)")

    if entrez_summary is not None:
        st.dataframe(entrez_summary, width="stretch", hide_index=True)

        plot_df = entrez_summary.copy()
        if {"mapped_up", "original_up", "mapped_down", "original_down"}.issubset(plot_df.columns):
            plot_df["UP mapped %"] = (plot_df["mapped_up"] / plot_df["original_up"].replace(0, np.nan) * 100).round(1)
            plot_df["DOWN mapped %"] = (plot_df["mapped_down"] / plot_df["original_down"].replace(0, np.nan) * 100).round(1)
            fig = px.bar(
                plot_df,
                x="signature",
                y=["UP mapped %", "DOWN mapped %"],
                barmode="group",
                labels={"value": "% mapped to Entrez", "variable": ""},
            )
            st.plotly_chart(fig, width="stretch")
    else:
        missing(ENTREZ_SUMMARY_FILE, "04_cmap_analysis.py")


# ============================================================
# Tab 2: Cross-Signature Results
# ============================================================
with tab_cross:
    if prioritized is None:
        missing(PRIORITIZED_ALL_FILE, "04l_cmap_prioritize.py")
    else:
        sig_cols = detect_signature_columns(prioritized, manifest)

        threshold = (
            float(prioritized["strong_negative_threshold"].iloc[0])
            if "strong_negative_threshold" in prioritized.columns and len(prioritized) > 0
            else None
        )

        st.caption(
            f"Strong-negative tau threshold in this run: {threshold:.1f}"
            if threshold is not None
            else "Strong-negative tau threshold not found in the data."
        )

        stat_cols = st.columns(4)
        if "n_negative" in prioritized.columns:
            stat_cols[0].metric("Compounds with any negative tau", int((prioritized["n_negative"] > 0).sum()))
        if "n_strong_negative" in prioritized.columns:
            stat_cols[1].metric("Compounds with strong negative tau", int((prioritized["n_strong_negative"] > 0).sum()))
        if "mean_tau" in prioritized.columns:
            stat_cols[2].metric("Median of mean_tau", round(float(prioritized["mean_tau"].median()), 1))
        if "minimum_tau" in prioritized.columns:
            stat_cols[3].metric("Most negative single tau", round(float(prioritized["minimum_tau"].min()), 1))

        st.divider()
        st.subheader("Tier distribution")
        if "priority_tier" in prioritized.columns:
            tier_counts = prioritized["priority_tier"].value_counts().reset_index()
            tier_counts.columns = ["tier", "count"]
            # order by leading digit in the label, not the label string itself,
            # so a future wording change doesn't break the sort
            tier_counts["order"] = tier_counts["tier"].str.extract(r"(\d+)").astype(float)
            tier_counts = tier_counts.sort_values("order")
            fig = px.bar(tier_counts, x="tier", y="count", color="tier")
            st.plotly_chart(fig, width="stretch")

        st.divider()
        st.subheader("Tau distribution across signatures")
        if sig_cols:
            melted = prioritized.melt(
                id_vars=["pert_iname"] if "pert_iname" in prioritized.columns else None,
                value_vars=sig_cols,
                var_name="signature",
                value_name="tau",
            ).dropna(subset=["tau"])
            fig = px.histogram(melted, x="tau", color="signature", barmode="overlay", opacity=0.6, nbins=60)
            if threshold is not None:
                fig.add_vline(x=threshold, line_dash="dash", line_color="red")
            st.plotly_chart(fig, width="stretch")
        else:
            st.info("No per-signature tau columns detected.")

        st.divider()
        st.subheader("Cross-signature heatmap (top candidates)")
        if sig_cols and "priority_rank" in prioritized.columns:
            top_n_heatmap = st.slider("Number of top candidates to show", 5, 100, 30, key="heatmap_n")
            top_for_heatmap = prioritized.sort_values("priority_rank").head(top_n_heatmap)
            label_col = "pert_iname" if "pert_iname" in top_for_heatmap.columns else top_for_heatmap.index
            fig = px.imshow(
                top_for_heatmap[sig_cols].to_numpy(),
                x=sig_cols,
                y=top_for_heatmap[label_col] if isinstance(label_col, pd.Series) else top_for_heatmap.index,
                color_continuous_scale="RdBu",
                aspect="auto",
                labels={"color": "tau"},
            )
            st.plotly_chart(fig, width="stretch")

        st.divider()
        st.subheader("Top candidates")
        if "priority_rank" in prioritized.columns and "mean_tau" in prioritized.columns:
            top_n_bar = st.slider("Number of candidates to show", 5, 50, 20, key="bar_n")
            top_for_bar = prioritized.sort_values("priority_rank").head(top_n_bar)
            label_col = "pert_iname" if "pert_iname" in top_for_bar.columns else top_for_bar.index
            fig = go.Figure(
                go.Bar(
                    x=top_for_bar["mean_tau"][::-1],
                    y=(top_for_bar[label_col] if isinstance(label_col, pd.Series) else top_for_bar.index)[::-1],
                    orientation="h",
                    marker_color=top_for_bar["priority_tier_number"][::-1] if "priority_tier_number" in top_for_bar.columns else None,
                )
            )
            fig.update_layout(xaxis_title="mean_tau", height=max(400, top_n_bar * 25))
            st.plotly_chart(fig, width="stretch")

        st.divider()
        st.subheader("Signature agreement")
        if len(sig_cols) >= 2:
            c1, c2 = st.columns(2)
            sig_x = c1.selectbox("Signature (x-axis)", sig_cols, index=0)
            sig_y = c2.selectbox("Signature (y-axis)", sig_cols, index=1)
            fig = px.scatter(
                prioritized,
                x=sig_x,
                y=sig_y,
                color="priority_tier" if "priority_tier" in prioritized.columns else None,
                hover_name="pert_iname" if "pert_iname" in prioritized.columns else None,
            )
            st.plotly_chart(fig, width="stretch")
        else:
            st.info("Need at least 2 queried signatures to compare agreement.")


# ============================================================
# Tab 3: Candidate Explorer
# ============================================================
with tab_explorer:
    if prioritized is None:
        missing(PRIORITIZED_ALL_FILE, "04l_cmap_prioritize.py")
    else:
        sig_cols = detect_signature_columns(prioritized, manifest)

        f1, f2, f3 = st.columns(3)
        search = f1.text_input("Search compound name or pert_id")

        tier_options = (
            sorted(prioritized["priority_tier"].dropna().unique().tolist())
            if "priority_tier" in prioritized.columns
            else []
        )
        selected_tiers = f2.multiselect("Tier", tier_options, default=tier_options)

        sort_options = [c for c in ["priority_rank", "mean_tau", "median_tau", "minimum_tau"] if c in prioritized.columns]
        sort_by = f3.selectbox("Sort by", sort_options) if sort_options else None

        filtered = prioritized.copy()
        if search:
            name_cols = [c for c in ["pert_iname", "pert_id"] if c in filtered.columns]
            if name_cols:
                mask = np.zeros(len(filtered), dtype=bool)
                for c in name_cols:
                    mask |= filtered[c].astype(str).str.contains(search, case=False, na=False)
                filtered = filtered[mask]
        if selected_tiers and "priority_tier" in filtered.columns:
            filtered = filtered[filtered["priority_tier"].isin(selected_tiers)]

        f4, f5 = st.columns(2)
        if "n_negative" in filtered.columns:
            min_neg = f4.slider("Minimum negative signatures", 0, int(prioritized["n_negative"].max()), 0)
            filtered = filtered[filtered["n_negative"] >= min_neg]
        if "n_strong_negative" in filtered.columns:
            min_strong = f5.slider("Minimum strong-negative signatures", 0, int(prioritized["n_strong_negative"].max()), 0)
            filtered = filtered[filtered["n_strong_negative"] >= min_strong]

        if sort_by:
            ascending = sort_by == "priority_rank"
            filtered = filtered.sort_values(sort_by, ascending=ascending)

        st.caption(f"{len(filtered):,} of {len(prioritized):,} compounds match the current filters")
        st.dataframe(filtered, width="stretch", hide_index=True)
        st.download_button("Download filtered table as CSV", filtered.to_csv(index=False), file_name="cmap_filtered.csv")

        st.divider()
        st.subheader("Compound detail")
        if "pert_iname" in filtered.columns and len(filtered) > 0:
            chosen = st.selectbox("Select a compound", filtered["pert_iname"].tolist())
            row = filtered[filtered["pert_iname"] == chosen].iloc[0]
            detail_cols = st.columns(3)
            shown = 0
            for label, key in [
                ("Perturbagen ID", "pert_id"),
                ("Tier", "priority_tier"),
                ("Priority rank", "priority_rank"),
                ("Mean tau", "mean_tau"),
                ("Median tau", "median_tau"),
                ("Minimum tau", "minimum_tau"),
                ("Signatures negative", "n_negative"),
                ("Signatures strong-negative", "n_strong_negative"),
            ]:
                if key in row.index:
                    detail_cols[shown % 3].metric(label, row[key])
                    shown += 1
            if sig_cols:
                st.write("Per-signature tau:")
                st.dataframe(row[sig_cols].to_frame(name="tau"), width="stretch")


# ============================================================
# Tab 4: Pipeline Status
# ============================================================
with tab_jobs:
    job_manifest = load_json(str(JOB_MANIFEST_FILE))
    job_status = load_json(str(JOB_STATUS_FILE))
    download_manifest = load_json(str(DOWNLOAD_MANIFEST_FILE))

    st.subheader("Submitted jobs")
    if job_manifest is not None:
        jobs = job_manifest.get("existing_completed_queries", []) + job_manifest.get("new_submissions", [])
        if jobs:
            st.dataframe(pd.json_normalize(jobs), width="stretch", hide_index=True)
        else:
            st.info("Job manifest has no job entries yet.")
    else:
        missing(JOB_MANIFEST_FILE, "04g_cmap_submit_all.py")

    st.divider()
    st.subheader("Job status")
    if job_status is not None:
        jobs_list = job_status.get("jobs", [])
        if jobs_list:
            status_df = pd.json_normalize(jobs_list)
            st.dataframe(status_df, width="stretch", hide_index=True)

            active_statuses = {"pending", "submitted", "running", "queued"}
            failed_statuses = {"failed", "error", "cancelled"}
            status_col = status_df["status"] if "status" in status_df.columns else pd.Series(dtype=str)
            c1, c2, c3 = st.columns(3)
            c1.metric("Completed", int((status_col == "completed").sum()))
            c2.metric("Active", int(status_col.isin(active_statuses).sum()))
            c3.metric("Failed", int(status_col.isin(failed_statuses).sum()))
        else:
            st.info("No jobs found in the status file.")
    else:
        missing(JOB_STATUS_FILE, "04h_cmap_check_all_jobs.py")

    st.divider()
    st.subheader("Downloads")
    if download_manifest is not None:
        records = download_manifest if isinstance(download_manifest, list) else download_manifest.get("jobs", [])
        if records:
            st.dataframe(pd.json_normalize(records), width="stretch", hide_index=True)
        else:
            st.info("No download records found.")
    else:
        missing(DOWNLOAD_MANIFEST_FILE, "04j_cmap_download_all.py")


# ============================================================
# Tab 5: Data Provenance
# ============================================================
with tab_provenance:
    st.caption("Which files are currently being visualized on this page")

    provenance_targets = [
        (SIGNATURE_MANIFEST_FILE, manifest),
        (SIGNATURE_DEFINITIONS_FILE, definitions),
        (ENTREZ_SUMMARY_FILE, entrez_summary),
        (PRIORITIZED_ALL_FILE, prioritized),
        (COMPOUND_TAU_LONG_FILE, compound_long),
        (JOB_MANIFEST_FILE, None),
        (JOB_STATUS_FILE, None),
        (DOWNLOAD_MANIFEST_FILE, None),
    ]

    rows = []
    for path, df in provenance_targets:
        info = file_info(path)
        if df is not None:
            info["rows"] = len(df)
            info["columns"] = len(df.columns)
        rows.append(info)

    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)