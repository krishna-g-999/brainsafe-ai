#!/usr/bin/env python3
"""
FIX 2: Ensure scorer.py NPS formula matches ml_expander.py (both must use same 7-dim formula).
Run from: ~/brainsafe_ai/
"""

import re

with open("scorer.py") as f:
    src = f.read()

# Common patterns found in scorer.py that might use partial NPS
# Pattern A: inline arithmetic
PATTERN_A = re.compile(
    r"""(antioxidant.*?\*\s*3.*?anti_inflammatory.*?\*\s*3.*?mitochondrial_support.*?\*\s*2.*?aggregation_modulation.*?\*\s*2)""",
    re.DOTALL
)

# Canonical function to inject if needed
CANONICAL_NPS = """
# ── CANONICAL NPS FUNCTION (must match ml_expander.py) ───────────────────────
NPS_WEIGHTS = {
    "antioxidant":           3,
    "anti_inflammatory":     3,
    "mitochondrial_support": 2,
    "aggregation_modulation":2,
    "cognitive_enhancement": 1,
    "neurogenesis":          1,
    "synaptic_plasticity":   1,
}
NPS_MAX_RAW = 130  # 10 * (3+3+2+2+1+1+1)

def compute_nps(compound: dict) -> float:
    """
    NPS = sum(score_i * w_i) / 130 * 100  ∈ [0, 100]
    Literature basis:
      Antioxidant/Anti-Inflam (w=3 each): Lin & Beal 2006; Heneka et al. 2015
      Mitochon./Aggregation (w=2 each):   Lin & Beal 2006; Hardy & Selkoe 2002
      Cogn./Neuro./Synap. (w=1 each):     Wager et al. ACS Chem Neurosci 2010
    """
    raw = sum(float(compound.get(col, 5.0)) * w for col, w in NPS_WEIGHTS.items())
    return round(min(100.0, raw / NPS_MAX_RAW * 100.0), 1)
"""

if "compute_nps" in src or "def nps" in src.lower() or "calc_nps" in src:
    print("⚠️  scorer.py has its own NPS function — check manually that it matches:")
    print("   NPS = sum(score_i * w_i) / 130 * 100")
    print("   Weights: antioxidant=3, anti_inflammatory=3, mitochondrial=2,")
    print("            aggregation=2, cognitive=1, neurogenesis=1, synaptic=1")
else:
    # Inject canonical function after imports block
    insert_after = "from pathlib import Path"
    if insert_after in src:
        src = src.replace(insert_after, insert_after + CANONICAL_NPS, 1)
        with open("scorer.py", "w") as f:
            f.write(src)
        print("✅ FIX 2 APPLIED: compute_nps() injected into scorer.py")
    else:
        print("⚠️  Could not locate injection point in scorer.py — add manually:")
        print(CANONICAL_NPS)
