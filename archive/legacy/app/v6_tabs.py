"""
v6_tabs.py  -  BrainSafe AI v6
Rendering functions for the 3 new v6 tabs:
  render_neurotransmitter_panel, render_brain_region_panel,
  render_dose_response_panel, render_v6_tabs
Extracted from patch_app_v6.NEW_FUNCTIONS.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

# ── BrainSafe AI v6 module imports (guarded) ─────────────────────────────────
# Defines _V6_MODULES_OK / _V6_ERR_MSG used by the panel functions below, and
# brings in the mapper/dose-response/scorer symbols those panels call.
# Ported from patch_app_v6.NEW_IMPORTS — this block was dropped when the panel
# functions were extracted into v6_tabs.py, leaving _V6_MODULES_OK undefined.
try:
    from neurotransmitter_mapper import compute_nt_effects, format_nt_summary
    from brain_region_mapper import compute_region_scores, format_region_summary
    from dose_response import (compute_dose_response, compute_dose_risk,
                                all_dimensions_dose_response)
    from scorer import neuro_score, neuro_score_breakdown, disease_nps
    _V6_MODULES_OK = True
    _V6_ERR_MSG = ""
except ImportError as _v6_err:
    _V6_MODULES_OK = False
    _V6_ERR_MSG = str(_v6_err)



# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# BrainSafe AI v6 â€” New feature panels
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

def render_neurotransmitter_panel(entry: dict, patient_mode: bool = False) -> None:
    """
    Render the Neurotransmitter Effect Panel.
    Shows whether a compound raises / lowers / modulates each of the
    5 key neurotransmitters relevant to NDDs.
    """
    import streamlit as st

    if not _V6_MODULES_OK:
        st.warning(f"NT module unavailable: {_V6_ERR_MSG}")
        return

    st.markdown("### Neurotransmitter Effects")
    st.caption(
        "Based on known enzyme and transporter interactions. "
        "Arrows show predicted net effect direction at research doses."
    )

    # Build enzyme list from entry data
    enzymes = []
    enzyme_data = entry.get("enzyme_data", [])
    if not enzyme_data:
        # Build from class templates if enzyme_data not present
        comp_type = entry.get("compound_type", "general").lower()
        templates = entry.get("enzyme_actions", [])
        for t in templates:
            enzymes.append({
                "name":     t.get("enzyme", t.get("name", "")),
                "action":   t.get("action", "Modulation"),
                "strength": t.get("strength", "Weak"),
            })

    if not enzymes and entry.get("mechanisms"):
        for m in entry.get("mechanisms", []):
            enzymes.append({
                "name":     m.get("target_name", m.get("target", "")),
                "action":   m.get("action_type", m.get("action", "Modulation")),
                "strength": "Moderate",
            })

    effects = compute_nt_effects(enzymes)
    rows    = format_nt_summary(effects, patient_mode=patient_mode)

    if not rows:
        st.info(
            "No specific neurotransmitter interactions found for this compound in the database. "
            "This may indicate a compound without well-characterised enzyme targets, "
            "not necessarily low activity."
        )
        return

    NT_COLORS = {
        "Dopamine":             "#7C3AED",
        "Serotonin":            "#059669",
        "Acetylcholine":        "#2563EB",
        "GABA (calming signal)":"#D97706",
        "Glutamate (excitatory signal)": "#DC2626",
    }

    for row in rows:
        nt   = row["neurotransmitter"]
        arr  = row["arrow"]
        dir_ = row["direction"] if not patient_mode else row["direction"]
        str_ = row["strength"]
        mech = row["mechanism"]
        color = NT_COLORS.get(nt, "#64748B")

        arrow_colors = {"â†‘": "green", "â†“": "red", "âŸ·": "orange", "â†’": "gray"}
        arrow_color = arrow_colors.get(arr, "gray")

        with st.container():
            col1, col2, col3 = st.columns([1, 3, 4])
            with col1:
                st.markdown(
                    f"<span style='font-size:28px;color:{arrow_color}'>{arr}</span>",
                    unsafe_allow_html=True
                )
            with col2:
                st.markdown(
                    f"**<span style='color:{color}'>{nt}</span>**  "
                    f"<span style='font-size:11px;background:#F1F5F9;"
                    f"padding:2px 6px;border-radius:10px'>{str_}</span>",
                    unsafe_allow_html=True
                )
            with col3:
                if patient_mode:
                    st.caption(dir_)
                else:
                    st.caption(mech[:120] + ("..." if len(mech) > 120 else ""))

    if patient_mode:
        st.info(
            "âš ï¸  This information is for educational purposes only. "
            "Neurotransmitter effects depend on dose, formulation, individual variation, "
            "and concurrent medications. Always consult your neurologist or pharmacist."
        )


def render_brain_region_panel(dimension_scores: dict, patient_mode: bool = False) -> None:
    """
    Render the Brain Region Specificity panel.
    Shows which brain regions are most affected based on the compound's
    7-axis neuroprotective profile.
    """
    import streamlit as st
    import plotly.graph_objects as go

    if not _V6_MODULES_OK:
        st.warning(f"Brain region module unavailable: {_V6_ERR_MSG}")
        return

    st.markdown("### Brain Region Activity Map")
    st.caption(
        "Predicted regional activity based on the compound's 7-axis profile and "
        "known neuroanatomical associations. Higher score = stronger predicted activity in that region."
    )

    regions = compute_region_scores(dimension_scores)
    rows    = format_region_summary(regions)

    if not rows:
        st.info("No region data available.")
        return

    # Bar chart of region scores
    region_names  = [r["region"] for r in rows[:8]]
    region_scores = [r["score"] for r in rows[:8]]
    region_colors = [
        "#16A34A" if r["score"] >= 65 else "#D97706" if r["score"] >= 35 else "#DC2626"
        for r in rows[:8]
    ]

    fig = go.Figure(go.Bar(
        x=region_scores,
        y=region_names,
        orientation="h",
        marker_color=region_colors,
        text=[f"{s:.0f}" for s in region_scores],
        textposition="outside",
    ))
    fig.update_layout(
        xaxis_title="Regional activity score (0-100)",
        yaxis=dict(autorange="reversed"),
        height=350,
        margin=dict(l=10, r=60, t=20, b=40),
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(size=12),
        xaxis=dict(range=[0, 110]),
    )
    st.plotly_chart(fig, use_container_width=True)

    # NDD disease relevance callout
    top_regions = rows[:3]
    diseases_affected = set()
    for r in top_regions:
        diseases_affected.update(r.get("diseases", []))

    dis_display = {
        "alzheimers": "Alzheimer's disease",
        "parkinsons":  "Parkinson's disease",
        "als":         "ALS",
        "huntingtons": "Huntington's disease",
    }
    if diseases_affected:
        dis_names = [dis_display.get(d, d) for d in sorted(diseases_affected)]
        st.success(
            f"**Primary disease relevance:** {', '.join(dis_names)}  "
            f"(based on top activated brain regions)"
        )

    # Mechanism explanations
    with st.expander("Show region-mechanism details"):
        for r in rows[:6]:
            st.markdown(f"**{r['region']}** (score: {r['score']:.0f})")
            st.caption(r["mechanism"])


def render_dose_response_panel(
    compound_name: str,
    compound_type: str,
    mw: float,
    dimension_scores: dict,
    patient_mode: bool = False
) -> None:
    """
    Render the Dose-Response panel.
    Interactive dose slider showing predicted activity across all 7 dimensions.
    """
    import streamlit as st
    import plotly.graph_objects as go
    import numpy as np

    if not _V6_MODULES_OK:
        st.warning(f"Dose-response module unavailable: {_V6_ERR_MSG}")
        return

    st.markdown("### Dose-Response Analysis")
    st.caption(
        "4-parameter logistic (4PL) model. Class-based ECâ‚…â‚€ estimates â€” "
        "uncertainty is Ã—3 to Ã—10 fold. Use compound-specific ICâ‚…â‚€ data for precision."
    )

    DIMENSION_COLS = [
        "antioxidant", "anti_inflammatory", "mitochondrial_support",
        "aggregation_modulation", "cognitive_enhancement",
        "neurogenesis", "synaptic_plasticity",
    ]
    DIM_COLORS = {
        "antioxidant":            "#F59E0B",
        "anti_inflammatory":      "#10B981",
        "mitochondrial_support":  "#6366F1",
        "aggregation_modulation": "#EC4899",
        "cognitive_enhancement":  "#3B82F6",
        "neurogenesis":           "#8B5CF6",
        "synaptic_plasticity":    "#14B8A6",
    }

    selected_dim = st.selectbox(
        "Select dimension to model:",
        options=DIMENSION_COLS,
        format_func=lambda x: x.replace("_", " ").title(),
        index=0,
    )

    dose_val = st.slider(
        "Oral dose (mg/day):" if not patient_mode else "Daily dose (mg):",
        min_value=1,
        max_value=2000,
        value=200,
        step=10,
    )

    # Compute curve
    curve = compute_dose_response(
        compound_name, compound_type, mw,
        base_nps=dimension_scores.get(selected_dim, 50.0),
        dimension=selected_dim,
    )

    risk = compute_dose_risk(dose_val, compound_type)

    # Plot
    fig = go.Figure()

    # Main curve
    fig.add_trace(go.Scatter(
        x=curve.x_doses,
        y=curve.y_responses,
        mode="lines",
        line=dict(color=DIM_COLORS.get(selected_dim, "#3B82F6"), width=3),
        name=selected_dim.replace("_", " ").title(),
    ))

    # Current dose marker
    dose_response_at_slider = float(np.interp(dose_val, curve.x_doses, curve.y_responses))
    fig.add_vline(x=dose_val, line_dash="dash", line_color="gray", line_width=1.5)
    fig.add_trace(go.Scatter(
        x=[dose_val], y=[dose_response_at_slider],
        mode="markers+text",
        marker=dict(size=12, color=DIM_COLORS.get(selected_dim, "#3B82F6"),
                    line=dict(width=2, color="white")),
        text=[f"  {dose_response_at_slider:.1f}"],
        textposition="top right",
        showlegend=False,
    ))

    # EC50 marker
    fig.add_vline(x=curve.ec50_mg, line_dash="dot", line_color="#DC2626", line_width=1,
                  annotation_text=f"ECâ‚…â‚€ â‰ˆ{curve.ec50_mg:.0f}mg",
                  annotation_font_color="#DC2626",
                  annotation_font_size=11)

    fig.update_layout(
        xaxis_title="Oral dose (mg/day)" if not patient_mode else "Daily dose (mg)",
        xaxis_type="log",
        yaxis_title="Predicted activity (0-100)",
        yaxis_range=[0, 105],
        height=350,
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(size=12),
        margin=dict(l=10, r=10, t=30, b=50),
        legend=dict(orientation="h", y=-0.25),
    )
    st.plotly_chart(fig, use_container_width=True)

    # Metrics
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("At selected dose", f"{dose_response_at_slider:.1f}/100")
    with col2:
        st.metric("ECâ‚…â‚€ (oral equivalent)", f"~{curve.ec50_mg:.0f} mg")
    with col3:
        st.metric("Maximum activity (plateau)", f"{curve.plateau:.1f}/100")

    # Risk badge
    risk_colors = {"Low": "success", "Moderate": "warning", "High": "error"}
    risk_fns = {"Low": st.success, "Moderate": st.warning, "High": st.error}
    risk_fn = risk_fns.get(risk.risk_level, st.info)
    risk_fn(
        f"**Safety at {dose_val} mg/day:** {risk.risk_level} risk â€” {risk.recommendation}"
    )
    if risk.concerns:
        for concern in risk.concerns:
            st.caption(f"âš ï¸ {concern}")

    # Disclaimer
    st.caption(
        "âš ï¸ DISCLAIMER: Dose estimates are computational predictions based on "
        "class-level pharmacokinetic assumptions. Not clinical dosing guidance. "
        "Consult a healthcare professional before supplementation."
    )

    with st.expander("Technical note (for researchers)"):
        st.markdown(f"""
**Model:** 4-Parameter Logistic (4PL)  
`y = Bottom + (Top - Bottom) / (1 + (EC50/x)^Hill)`

**Compound class:** {compound_type}  
**Molecular weight:** {mw:.1f} Da  
**ECâ‚…â‚€ estimate:** {curve.ec50_mg:.1f} mg oral (class-based)  
**Uncertainty:** Â±{curve.compound_type}-class fold-range (see dose_response.py)  
**Note:** {curve.note}
        """)


def render_v6_tabs(
    entry: dict,
    dimension_scores: dict,
    patient_mode: bool = False,
) -> None:
    """
    Master function: renders all 3 new v6 tabs in a Streamlit tab group.
    Call this after the existing radar chart section in app.py.

    Usage in app.py (add after existing score display):
        from patch_app_v6 import render_v6_tabs
        render_v6_tabs(entry, dimension_scores, patient_mode=patient_mode)
    """
    import streamlit as st

    compound_name = entry.get("name", "Unknown")
    compound_type = entry.get("compound_type", "general")
    mw            = float(entry.get("mw", 350.0))

    tab_nt, tab_regions, tab_dose = st.tabs([
        "ðŸ§  Neurotransmitters",
        "ðŸ“ Brain Regions",
        "ðŸ“ˆ Dose-Response",
    ])

    with tab_nt:
        render_neurotransmitter_panel(entry, patient_mode=patient_mode)

    with tab_regions:
        render_brain_region_panel(dimension_scores, patient_mode=patient_mode)

    with tab_dose:
        render_dose_response_panel(
            compound_name, compound_type, mw,
            dimension_scores, patient_mode=patient_mode
        )


# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# Auto-patcher
