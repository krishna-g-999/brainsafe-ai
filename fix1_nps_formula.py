#!/usr/bin/env python3
"""
FIX 1: Unified NPS formula — replace _nps() in ml_expander.py
Run from: ~/brainsafe_ai/
"""

import re

with open("ml_expander.py") as f:
    src = f.read()

# Old broken 4-dim formula (only antioxidant, anti_inflammatory,
# mitochondrial_support, aggregation_modulation)
OLD = """def _nps(pred: dict) -> float:
    raw = (pred.get("antioxidant", 5) * 3
        + pred.get("anti_inflammatory", 5) * 3
        + pred.get("mitochondrial_support", 5) * 2
        + pred.get("aggregation_modulation", 5) * 2)
    return min(100.0, raw)"""

# New correct 7-dim normalized formula (max raw = 130)
NEW = """# ── NPS_WEIGHTS ─────────────────────────────────────────────────────────────
# Antioxidant      w=3  (23.1%)  upstream ROS driver — Lin & Beal, Nature 2006
# Anti-Inflam.     w=3  (23.1%)  upstream NF-kB/TNF-a — Heneka et al., Lancet Neurol 2015
# Mitoch. Sup.     w=2  (15.4%)  ETC dysfunction amplifier — Lin & Beal 2006
# Aggregation Mod. w=2  (15.4%)  protein misfolding — Hardy & Selkoe, Science 2002
# Cognitive Enh.   w=1  ( 7.7%)  downstream functional endpoint
# Neurogenesis     w=1  ( 7.7%)  BDNF/TrkB mediated
# Synaptic Plas.   w=1  ( 7.7%)  LTP/BDNF mediated
# Max raw = 10*(3+3+2+2+1+1+1) = 130  → normalized to 0-100
NPS_WEIGHTS = {
    "antioxidant":           3,
    "anti_inflammatory":     3,
    "mitochondrial_support": 2,
    "aggregation_modulation":2,
    "cognitive_enhancement": 1,
    "neurogenesis":          1,
    "synaptic_plasticity":   1,
}
NPS_MAX_RAW = sum(10 * w for w in NPS_WEIGHTS.values())  # = 130

def _nps(pred: dict) -> float:
    """
    Neuroprotective Score (NPS) — 7-dimension weighted formula.
    NPS = (sum of score_i * weight_i) / 130 * 100  ∈ [0, 100]
    Weights are literature-justified; see NPS_WEIGHTS above.
    Reference: BrainSafe AI manuscript (2026), Methods §2.3
    """
    raw = sum(pred.get(col, 5.0) * w for col, w in NPS_WEIGHTS.items())
    return round(min(100.0, raw / NPS_MAX_RAW * 100.0), 1)"""

if OLD in src:
    src = src.replace(OLD, NEW)
    with open("ml_expander.py", "w") as f:
        f.write(src)
    print("✅ FIX 1 APPLIED: _nps() updated to full 7-dim normalized formula (max=130)")
    print("   NPS formula: sum(score_i * weight_i) / 130 * 100")
else:
    print("⚠️  Pattern not found — check indentation in ml_expander.py manually")
    print("   Add this function body replacing the existing _nps():")
    print(NEW)
