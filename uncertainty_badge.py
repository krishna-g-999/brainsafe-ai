"""
uncertainty_badge.py  —  BrainSafe AI SAI-Net
Renders a per-dimension uncertainty mini-bar in the Streamlit UI.
Import and call render_uncertainty_table(data) inside any Streamlit tab.
"""

def uncertainty_html(uncertainty: dict, scores: dict) -> str:
    """
    Returns an HTML table showing each dimension's predicted score
    alongside its uncertainty (std) as a coloured confidence band.
    """
    SCORE_LABELS = {
        "antioxidant":           "Antioxidant",
        "anti_inflammatory":     "Anti-Inflammatory",
        "mitochondrial_support": "Mito. Support",
        "aggregation_modulation":"Aggreg. Mod.",
        "cognitive_enhancement": "Cogn. Enhancement",
        "neurogenesis":          "Neurogenesis",
        "synaptic_plasticity":   "Synaptic Plasticity",
    }
    rows = ""
    for dim, label in SCORE_LABELS.items():
        score = float(scores.get(dim, 5.0))
        std   = float(uncertainty.get(dim, 0.0))
        pct   = score / 10 * 100
        # Confidence colour: green < 1.0 std, amber 1-2, red >2
        conf_col = ("#15803D" if std < 1.0 else
                    "#B45309" if std < 2.0 else "#B91C1C")
        conf_lbl = ("High" if std < 1.0 else "Medium" if std < 2.0 else "Low")
        rows += f"""
        <tr>
          <td style="width:36%;font-weight:500;font-size:0.84rem">{label}</td>
          <td style="width:28%">
            <div style="background:#EEF1F6;border-radius:4px;height:8px">
              <div style="width:{pct:.0f}%;background:#1B6B8A;height:8px;border-radius:4px"></div>
            </div>
          </td>
          <td style="width:12%;font-weight:700;font-size:0.84rem;text-align:right">{score}</td>
          <td style="width:12%;font-size:0.78rem;color:{conf_col};font-weight:700">{conf_lbl}</td>
          <td style="width:12%;font-size:0.78rem;color:#94A3B8">±{std}</td>
        </tr>"""

    return f"""
    <table style="width:100%;border-collapse:collapse">
      <thead>
        <tr style="border-bottom:2px solid #E2E8F0">
          <th style="text-align:left;font-size:0.72rem;color:#94A3B8;padding:5px 6px">Dimension</th>
          <th style="font-size:0.72rem;color:#94A3B8;padding:5px 6px">Profile</th>
          <th style="text-align:right;font-size:0.72rem;color:#94A3B8;padding:5px 6px">Score</th>
          <th style="font-size:0.72rem;color:#94A3B8;padding:5px 6px">Confidence</th>
          <th style="font-size:0.72rem;color:#94A3B8;padding:5px 6px">Std</th>
        </tr>
      </thead>
      <tbody>{rows}</tbody>
    </table>"""
