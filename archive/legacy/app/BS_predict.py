"""
BS_predict.py — BrainSafe AI (BS)
Genuine, scaffold-CV-validated inference for the ONE structure-predictable
endpoint (antioxidant). Uses the saved interpretable Ridge-on-descriptors model.

Scientific status (see BS_MODEL_CARD.md section 10, BS_antioxidant_report.json):
  * antioxidant: GENUINE — scaffold-CV R2 = 0.267, 95% CI [0.182, 0.346], on
    human/literature labels, structure-only features (no circular disease counts).
  * the other 6 neuroprotection dimensions are NOT structure-predictable
    (validated R2 <= 0) and are deliberately NOT predicted here.
"""
from __future__ import annotations
import os
import numpy as np

_MODEL = None
_SPEC = None
_DIR = os.path.dirname(os.path.abspath(__file__))
# Prefer the MEASURED DPPH model (ChEMBL radical-scavenging assays); fall back to the
# curated descriptor model only if the measured artifact is unavailable.
_MEASURED_PATH = os.path.join(_DIR, "models_genuine", "antioxidant_measured_dpph.joblib")
_CURATED_PATH = os.path.join(_DIR, "models_genuine", "antioxidant_genuine_ridge.joblib")

# Validated performance of the MEASURED model (models_genuine/antioxidant_measured_meta.json).
ANTIOXIDANT_VALIDATION = {
    "scaffold_cv_r2": 0.43,
    "rmse": 0.60,
    "spearman": 0.636,
    "method": "RF+ExtraTrees+HistGB ensemble on Morgan-1024 + 24 RDKit descriptors",
    "target": "DPPH radical-scavenging pIC50 (measured, ChEMBL)",
    "validation": "scaffold GroupKFold(5) on 2,862 measured compounds; temporal R2~0 "
                  "(pooled cross-lab DPPH protocols do not generalise across time)",
}


def _load():
    """Load the measured DPPH ensemble if present, else the curated fallback.
    Returns (model_obj, kind) where kind is 'measured' or 'curated'."""
    global _MODEL
    if _MODEL is None:
        import joblib
        if os.path.exists(_MEASURED_PATH):
            _MODEL = (joblib.load(_MEASURED_PATH), "measured")
        elif os.path.exists(_CURATED_PATH):
            _MODEL = (joblib.load(_CURATED_PATH), "curated")
        else:
            _MODEL = (None, None)
    return _MODEL


def antioxidant_available() -> bool:
    return os.path.exists(_MEASURED_PATH) or os.path.exists(_CURATED_PATH)


def predict_antioxidant(smiles: str | None) -> dict | None:
    """
    Predict antioxidant activity from SMILES.
    Measured model: predicts DPPH radical-scavenging pIC50 (and a 0-100 display scale,
    pIC50 3->0, 7->100). Returns None if no model / no SMILES / unparseable structure.
    """
    mdl, kind = _load()
    if mdl is None or not smiles or str(smiles).strip().lower() in ("nan", "n/a", "none", ""):
        return None
    try:
        from BS_predictive_model import morgan, descriptors
        if kind == "measured":
            X = np.hstack([morgan([str(smiles)]), descriptors([str(smiles)])])
            pic50 = float(np.mean([m.predict(X)[0] for m in mdl["models"]]))
            score = float(np.clip((pic50 - 3.0) / 4.0 * 100.0, 0.0, 100.0))  # 3..7 -> 0..100
            return {"antioxidant": round(score, 1), "pic50": round(pic50, 2), **ANTIOXIDANT_VALIDATION}
        else:  # curated fallback (0-100 directly)
            X = descriptors([str(smiles)])
            val = float(np.clip(mdl.predict(X)[0], 0.0, 100.0))
            return {"antioxidant": round(val, 1), "scaffold_cv_r2": 0.25,
                    "method": "curated Ridge (fallback)", "spearman": 0.47}
    except Exception:
        return None
