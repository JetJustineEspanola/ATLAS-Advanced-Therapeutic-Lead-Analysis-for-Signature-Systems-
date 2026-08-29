"""
Final Results - combined view across all 4 stages.
Placeholder until Stages 2-4 exist; will eventually show the final
composite drug_priority_score ranking (same weighted-percentile pattern
as the Stage 1 gene prioritization, extended with CMap/docking/ADMET scores).
"""

import streamlit as st

st.set_page_config(page_title="ATLAS - Results", layout="wide")

st.title("Final Results")
st.info(
    "Once all 4 stages have run, this page will combine CMap connectivity score, "
    "docking binding energy, and ADMET desirability into one final ranked "
    "candidate list - the same weighted-scoring approach used for gene "
    "prioritization in Stage 1, extended across the full pipeline."
)
