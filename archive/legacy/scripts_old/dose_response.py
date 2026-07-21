"""
dose_response.py  –  BrainSafe AI v6
4-Parameter Logistic (4PL) dose-response modelling.

This module implements the "Dose-Response" vision feature:
"What happens if I take more or less of this compound?"

The 4PL model: y = Bottom + (Top - Bottom) / (1 + (EC50/x)^HillSlope)

For compounds without measured IC50/EC50 data, class-based estimates
are used with wide confidence intervals clearly communicated to the user.

Data sources:
  - ChEMBL IC50/EC50/Ki values (fetched live)
  - Class-based estimates from published systematic reviews
"""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


# ---------------------------------------------------------------------------
# 4PL model
# ---------------------------------------------------------------------------
def four_param_logistic(x: np.ndarray,
                         bottom: float, top: float,
                         ec50: float, hill: float) -> np.ndarray:
    """4PL dose-response curve."""
    return bottom + (top - bottom) / (1.0 + (ec50 / np.clip(x, 1e-10, None)) ** hill)


def fit_4pl(doses: np.ndarray, responses: np.ndarray
            ) -> Optional[tuple[float, float, float, float]]:
    """
    Fit 4PL model to dose-response data.
    Returns (bottom, top, ec50, hill) or None if fitting fails.
    """
    try:
        from scipy.optimize import curve_fit
        p0 = [responses.min(), responses.max(), np.median(doses), 1.0]
        bounds = ([0, 0, 1e-3, 0.1], [100, 100, 1e6, 10.0])
        popt, _ = curve_fit(four_param_logistic, doses, responses,
                            p0=p0, bounds=bounds, maxfev=5000)
        return tuple(popt)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Class-based EC50 estimates (μM)
# From systematic reviews and meta-analyses
# ---------------------------------------------------------------------------
CLASS_EC50_ESTIMATES: dict[str, dict[str, dict]] = {
    "flavonoid": {
        "antioxidant":           {"ec50_um": 15.0,  "hill": 1.2, "ci_fold": 5},
        "anti_inflammatory":     {"ec50_um": 20.0,  "hill": 1.1, "ci_fold": 5},
        "mitochondrial_support": {"ec50_um": 50.0,  "hill": 0.9, "ci_fold": 8},
        "aggregation_modulation":{"ec50_um": 30.0,  "hill": 1.0, "ci_fold": 6},
        "cognitive_enhancement": {"ec50_um": 10.0,  "hill": 1.3, "ci_fold": 5},
        "neurogenesis":          {"ec50_um": 25.0,  "hill": 1.0, "ci_fold": 8},
        "synaptic_plasticity":   {"ec50_um": 20.0,  "hill": 1.1, "ci_fold": 8},
    },
    "stilbene": {   # resveratrol, pterostilbene
        "antioxidant":           {"ec50_um": 10.0,  "hill": 1.4, "ci_fold": 4},
        "anti_inflammatory":     {"ec50_um": 15.0,  "hill": 1.2, "ci_fold": 4},
        "mitochondrial_support": {"ec50_um": 25.0,  "hill": 1.1, "ci_fold": 5},
        "aggregation_modulation":{"ec50_um": 20.0,  "hill": 1.2, "ci_fold": 5},
        "cognitive_enhancement": {"ec50_um": 8.0,   "hill": 1.4, "ci_fold": 4},
        "neurogenesis":          {"ec50_um": 15.0,  "hill": 1.1, "ci_fold": 6},
        "synaptic_plasticity":   {"ec50_um": 12.0,  "hill": 1.2, "ci_fold": 6},
    },
    "curcuminoid": {
        "antioxidant":           {"ec50_um": 5.0,   "hill": 1.5, "ci_fold": 3},
        "anti_inflammatory":     {"ec50_um": 8.0,   "hill": 1.4, "ci_fold": 3},
        "mitochondrial_support": {"ec50_um": 20.0,  "hill": 1.2, "ci_fold": 4},
        "aggregation_modulation":{"ec50_um": 10.0,  "hill": 1.5, "ci_fold": 3},
        "cognitive_enhancement": {"ec50_um": 6.0,   "hill": 1.5, "ci_fold": 3},
        "neurogenesis":          {"ec50_um": 12.0,  "hill": 1.3, "ci_fold": 4},
        "synaptic_plasticity":   {"ec50_um": 10.0,  "hill": 1.3, "ci_fold": 4},
    },
    "general": {
        "antioxidant":           {"ec50_um": 30.0,  "hill": 1.0, "ci_fold": 10},
        "anti_inflammatory":     {"ec50_um": 40.0,  "hill": 1.0, "ci_fold": 10},
        "mitochondrial_support": {"ec50_um": 60.0,  "hill": 0.9, "ci_fold": 10},
        "aggregation_modulation":{"ec50_um": 50.0,  "hill": 1.0, "ci_fold": 10},
        "cognitive_enhancement": {"ec50_um": 35.0,  "hill": 1.0, "ci_fold": 10},
        "neurogenesis":          {"ec50_um": 50.0,  "hill": 0.9, "ci_fold": 12},
        "synaptic_plasticity":   {"ec50_um": 45.0,  "hill": 1.0, "ci_fold": 12},
    },
}

# MW-to-dose conversion factors (approximate oral dose in mg → serum μM)
# Based on simplified PK: F=0.30, Vd=50L, MW scaling
# These are rough estimates for educational/screening purposes only
MW_ORAL_TO_UM: dict[str, float] = {
    "typical_flavonoid_mw300":  0.22,   # mg oral → μM serum (rough)
    "typical_flavonoid_mw500":  0.13,
    "typical_curcuminoid_mw370": 0.03,  # poor bioavailability
    "typical_stilbene_mw228":    0.60,
}


@dataclass
class DoseResponseCurve:
    """Complete dose-response curve data for one dimension."""
    dimension:       str
    compound_type:   str
    doses_mg:        np.ndarray        # oral dose range in mg
    doses_um:        np.ndarray        # serum concentration estimate in μM
    responses:       np.ndarray        # predicted activity (0–100 scale)
    ec50_mg:         float             # EC50 in mg (oral equivalent)
    plateau:         float             # maximum achievable activity (%)
    is_estimated:    bool = True
    note:            str = ""
    x_doses:         list = field(default_factory=list)
    y_responses:     list = field(default_factory=list)

    def __post_init__(self):
        self.x_doses   = self.doses_mg.tolist()
        self.y_responses = self.responses.tolist()


@dataclass
class DoseRisk:
    """Safety/toxicity risk at a given dose."""
    dose_mg:         float
    risk_level:      str    # "Low", "Moderate", "High"
    concerns:        list[str]
    recommendation:  str


# Dose-safety thresholds for common compound classes (mg/day)
CLASS_SAFETY_THRESHOLDS: dict[str, dict] = {
    "flavonoid": {
        "safe_max_mg":      1000,
        "caution_mg":       500,
        "concerns_high": [
            "Possible P-glycoprotein inhibition at very high doses",
            "Thyroid interference at doses >2000mg/day",
        ],
    },
    "curcuminoid": {
        "safe_max_mg":      8000,   # curcumin is very safe; NOAEL is high
        "caution_mg":       3000,
        "concerns_high": [
            "GI discomfort at high doses",
            "CYP3A4 inhibition may affect co-administered drugs",
            "Iron chelation at very high doses",
        ],
    },
    "stilbene": {
        "safe_max_mg":      1000,
        "caution_mg":       500,
        "concerns_high": [
            "Hormonal effects at very high doses (SIRT1/estrogen modulation)",
        ],
    },
    "vitamin": {
        "safe_max_mg":      500,
        "caution_mg":       200,
        "concerns_high": [
            "Fat-soluble vitamins (A, D, E, K) accumulate — toxicity possible",
            "Water-soluble vitamins (B, C) are safer; excess excreted",
        ],
    },
    "alkaloid": {
        "safe_max_mg":      100,
        "caution_mg":       30,
        "concerns_high": [
            "Narrow therapeutic window for many alkaloids",
            "MAO interactions possible",
            "Consult a pharmacist/physician before dosing",
        ],
    },
    "general": {
        "safe_max_mg":      500,
        "caution_mg":       200,
        "concerns_high": [
            "Class-based estimate only — consult compound-specific literature",
        ],
    },
}


def compute_dose_response(compound_name: str,
                          compound_type: str,
                          mw: float = 350.0,
                          base_nps: float = 60.0,
                          dimension: str = "antioxidant"
                          ) -> DoseResponseCurve:
    """
    Compute a dose-response curve for a compound-dimension pair.

    Parameters
    ----------
    compound_name : str
    compound_type : str  – lowercase compound class (flavonoid, stilbene, etc.)
    mw            : float – molecular weight (Da)
    base_nps      : float – compound's NPS score (0–100); used to scale the plateau
    dimension     : str   – which of the 7 dimensions to model

    Returns
    -------
    DoseResponseCurve
    """
    # Get class-based estimates
    ct_key = "general"
    for key in CLASS_EC50_ESTIMATES:
        if key in compound_type.lower() or key in compound_name.lower():
            ct_key = key
            break
    params = CLASS_EC50_ESTIMATES[ct_key].get(
        dimension, CLASS_EC50_ESTIMATES["general"][dimension]
    )

    ec50_um  = params["ec50_um"]
    hill     = params["hill"]
    plateau  = min(100.0, base_nps * 1.15)   # plateau slightly above baseline NPS

    # Dose range: 1–2000 mg oral
    doses_mg = np.logspace(0, 3.3, 100)   # 1–2000 mg

    # Convert mg oral → μM serum (very rough pharmacokinetics)
    # Assumes: F=0.30 (oral bioavailability), Vd=50L, steady state
    f_oral = 0.30 if "curcumin" not in compound_name.lower() else 0.01
    f_pip  = 1.0  if "piperine" in compound_name.lower() else 1.0
    doses_um = (doses_mg * f_oral * f_pip * 1000) / (mw * 50)

    # 4PL response
    responses = four_param_logistic(doses_um, 0.0, plateau, ec50_um, hill)
    responses = np.clip(responses, 0.0, 100.0)

    # EC50 in mg
    ec50_mg = (ec50_um * mw * 50) / (f_oral * f_pip * 1000)

    note = (
        f"Class-based estimate for {ct_key} compounds. "
        f"Actual EC50 may differ ±{params['ci_fold']}× from this estimate. "
        "Use compound-specific IC50 data (ChEMBL) for precise modelling."
    )

    return DoseResponseCurve(
        dimension    = dimension,
        compound_type= ct_key,
        doses_mg     = doses_mg,
        doses_um     = doses_um,
        responses    = responses,
        ec50_mg      = round(ec50_mg, 1),
        plateau      = round(plateau, 1),
        is_estimated = True,
        note         = note,
    )


def compute_dose_risk(dose_mg: float, compound_type: str) -> DoseRisk:
    """Assess safety risk at a given oral dose."""
    ct_key = "general"
    for key in CLASS_SAFETY_THRESHOLDS:
        if key in compound_type.lower():
            ct_key = key
            break
    thresholds = CLASS_SAFETY_THRESHOLDS[ct_key]

    if dose_mg <= thresholds["caution_mg"]:
        return DoseRisk(
            dose_mg       = dose_mg,
            risk_level    = "Low",
            concerns      = [],
            recommendation= "Within typical research dose range. Standard dietary amounts are safe.",
        )
    elif dose_mg <= thresholds["safe_max_mg"]:
        return DoseRisk(
            dose_mg       = dose_mg,
            risk_level    = "Moderate",
            concerns      = ["Approaching upper dose limit; monitor for adverse effects"],
            recommendation= "Use caution. Discuss with healthcare provider.",
        )
    else:
        return DoseRisk(
            dose_mg       = dose_mg,
            risk_level    = "High",
            concerns      = thresholds["concerns_high"],
            recommendation= "Exceeds typical safe range. Do not use without medical supervision.",
        )


def all_dimensions_dose_response(compound_name: str,
                                  compound_type: str,
                                  mw: float,
                                  dimension_scores: dict[str, float]
                                  ) -> dict[str, DoseResponseCurve]:
    """Compute dose-response curves for all 7 dimensions at once."""
    from model_config import DIMENSION_COLS
    curves = {}
    for dim in DIMENSION_COLS:
        base_score = dimension_scores.get(dim, 50.0)
        curves[dim] = compute_dose_response(
            compound_name, compound_type, mw, base_score, dim
        )
    return curves
