"""
scorer.py  –  BrainSafe AI v6
NPS and utility scoring functions.

SCIENTIFIC FIX (v5 → v6):
  The v5 formula used only 4 of 7 dimensions (antioxidant×3, anti-inflam×3,
  mitochondrial×2, aggregation×2) and excluded cognitive, neurogenesis, and
  synaptic plasticity entirely. This contradicted the manuscript Methods 2.3
  and produced the same NPS for compounds with wildly different cognitive
  or neurogenesis profiles.

  v6 uses all 7 dimensions with disease-biology-grounded weights:
    Antioxidant          ×3  (ROS is a shared driver across all 4 NDDs)
    Anti-inflammatory    ×3  (neuroinflammation universal in AD/PD/ALS/HD)
    Mitochondrial        ×2  (complex I dysfunction central to PD/ALS)
    Aggregation          ×2  (amyloid/tau/TDP-43/HTT aggregation in all 4)
    Cognitive            ×2  (cholinergic/BDNF; especially AD/HD)
    Neurogenesis         ×1  (BDNF/Wnt; secondary mechanism, important for AD)
    Synaptic plasticity  ×1  (LTP/AMPA modulation; tertiary, supportive role)

  Total denominator = (3+3+2+2+2+1+1) × max_dim_score
  Normalised to 0–100. All 7 dimensions must be on the SAME 0–100 scale.

  Matches: Figure 1 weighted composite (updated), manuscript Methods 2.3 (updated).
"""


# ---------------------------------------------------------------------------
# Dimension weights (scientifically justified for NDD neuroprotection)
# ---------------------------------------------------------------------------
DIM_WEIGHTS: dict[str, int] = {
    "antioxidant":            3,
    "anti_inflammatory":      3,
    "mitochondrial_support":  2,
    "aggregation_modulation": 2,
    "cognitive_enhancement":  2,
    "neurogenesis":           1,
    "synaptic_plasticity":    1,
}

_WEIGHT_SUM = sum(DIM_WEIGHTS.values())   # = 14
_MAX_DIM    = 100                          # all dimensions scored 0–100
_DENOMINATOR = _WEIGHT_SUM * _MAX_DIM     # = 1400

SCORE_MAX: float = 100.0  # canonical max score, re-exported here for convenience


def neuro_score(data: dict) -> float:
    """
    Compute the Neuroprotective Score (NPS, 0–100) from 7 mechanistic
    dimension scores (each 0–100).

    Parameters
    ----------
    data : dict
        Must contain keys for all 7 dimensions (see DIM_WEIGHTS).
        Any missing key is treated as 0.

    Returns
    -------
    float  –  NPS in [0, 100], rounded to 1 decimal place.
    """
    raw = sum(
        DIM_WEIGHTS[dim] * float(data.get(dim, 0.0))
        for dim in DIM_WEIGHTS
    )
    return round(min(100.0, (raw / _DENOMINATOR) * 100.0), 1)


def neuro_score_breakdown(data: dict) -> dict:
    """
    Return per-dimension weighted contributions alongside the final NPS.
    Useful for the radar chart and the 'what drives the score' explanation.

    Returns
    -------
    dict with keys:
        nps          – overall NPS (0–100)
        contributions – {dim: weighted_pct_of_max} for each dimension
        dominant_dim  – dimension with the highest weighted contribution
    """
    contributions = {
        dim: round(
            (DIM_WEIGHTS[dim] * float(data.get(dim, 0.0))) / _DENOMINATOR * 100,
            2
        )
        for dim in DIM_WEIGHTS
    }
    nps = round(sum(contributions.values()), 1)
    dominant = max(contributions, key=contributions.__getitem__)
    return {
        "nps":          min(100.0, nps),
        "contributions": contributions,
        "dominant_dim":  dominant,
    }


# ---------------------------------------------------------------------------
# Disease-specific NPS — applies disease-relevant dimension re-weighting
# ---------------------------------------------------------------------------
DISEASE_WEIGHTS: dict[str, dict[str, int]] = {
    "alzheimers": {
        "antioxidant": 3, "anti_inflammatory": 3,
        "mitochondrial_support": 2, "aggregation_modulation": 3,
        "cognitive_enhancement": 3, "neurogenesis": 2, "synaptic_plasticity": 2,
    },
    "parkinsons": {
        "antioxidant": 3, "anti_inflammatory": 3,
        "mitochondrial_support": 3, "aggregation_modulation": 2,
        "cognitive_enhancement": 1, "neurogenesis": 1, "synaptic_plasticity": 1,
    },
    "als": {
        "antioxidant": 3, "anti_inflammatory": 3,
        "mitochondrial_support": 3, "aggregation_modulation": 2,
        "cognitive_enhancement": 1, "neurogenesis": 1, "synaptic_plasticity": 1,
    },
    "huntingtons": {
        "antioxidant": 2, "anti_inflammatory": 3,
        "mitochondrial_support": 2, "aggregation_modulation": 3,
        "cognitive_enhancement": 2, "neurogenesis": 2, "synaptic_plasticity": 2,
    },
}


def disease_nps(data: dict, disease: str) -> float:
    """
    Compute a disease-specific NPS using disease-appropriate dimension weights.
    Falls back to the composite NPS if disease is not recognised.
    """
    weights = DISEASE_WEIGHTS.get(disease.lower(), DIM_WEIGHTS)
    denom   = sum(weights.values()) * _MAX_DIM
    raw     = sum(weights[d] * float(data.get(d, 0.0)) for d in weights)
    return round(min(100.0, (raw / denom) * 100.0), 1)


# ---------------------------------------------------------------------------
# Colour helpers (UI rendering)
# ---------------------------------------------------------------------------
def score_color(score: float) -> str:
    if score >= 70:
        return "green"
    elif score >= 40:
        return "orange"
    return "red"


def score_label(score: float) -> str:
    if score >= 70:
        return "Strong"
    elif score >= 40:
        return "Moderate"
    return "Limited"


def bbb_color(bbb: str) -> str:
    return {"High": "green", "Medium": "orange",
            "Low-Med": "orange", "Low": "red"}.get(bbb, "gray")


def disease_color(level: str) -> str:
    return {"High": "green", "Med": "orange", "Low": "red"}.get(level, "gray")


# ---------------------------------------------------------------------------
# Confidence badge helpers
# ---------------------------------------------------------------------------
def confidence_label(max_tanimoto: float, tier: str | None = None) -> str:
    """
    Return human-readable confidence tier for a predicted compound.

    Parameters
    ----------
    max_tanimoto : float  –  maximum Tanimoto similarity to any training compound
    tier         : str   –  'gold', 'silver', or None (live prediction)
    """
    if tier == "gold":
        return "High (literature-curated)"
    if tier == "silver":
        return "Medium (pseudo-labelled)"
    if max_tanimoto >= 0.70:
        return "High (structurally close to training data)"
    if max_tanimoto >= 0.40:
        return "Medium"
    if max_tanimoto >= 0.30:
        return "Low-Medium (near applicability boundary)"
    return "Low (outside applicability domain)"


def confidence_color(max_tanimoto: float, tier: str | None = None) -> str:
    if tier == "gold":
        return "green"
    if tier == "silver":
        return "orange"
    if max_tanimoto >= 0.50:
        return "green"
    if max_tanimoto >= 0.30:
        return "orange"
    return "red"
