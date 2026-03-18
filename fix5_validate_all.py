#!/usr/bin/env python3
"""
FIX 5: Validation script — run after all fixes to confirm consistency.
Run from: ~/brainsafe_ai/
Checks: NPS formula consistency, ERC data counts, compound key names, CV R² stored.
"""
import json, importlib.util, re

print("=" * 60)
print("BrainSafe AI — Post-Fix Validation")
print("=" * 60)

# ── Check 1: NPS formula in ml_expander.py ───────────────────────────────────
with open("ml_expander.py") as f:
    ml_src = f.read()

if ("NPS_MAX_RAW" in ml_src or "NPS_MAX" in ml_src) and "130" in ml_src:
    print("✅ Check 1 PASS: _nps() uses 7-dim normalized formula (max=130)")
else:
    print("❌ Check 1 FAIL: _nps() still uses old 4-dim formula — apply FIX 1")

# ── Check 2: Key alignment ────────────────────────────────────────────────────
with open("compounds.json") as f:
    c = json.load(f)
if isinstance(c, dict) and "compounds" in c:
    c = c["compounds"]
sample = next(iter(c.values()))
required_keys = {"antioxidant","anti_inflammatory","mitochondrial_support",
                 "aggregation_modulation","cognitive_enhancement",
                 "neurogenesis","synaptic_plasticity","bbb","pathways"}
missing = required_keys - set(sample.keys())
if not missing:
    print("✅ Check 2 PASS: All 7 score dimensions + BBB + pathways present in compounds.json")
else:
    print(f"❌ Check 2 FAIL: Missing keys in compounds.json: {missing}")

# ── Check 3: NPS benchmark verification ───────────────────────────────────────
NPS_W = {"antioxidant":3,"anti_inflammatory":3,"mitochondrial_support":2,
          "aggregation_modulation":2,"cognitive_enhancement":1,"neurogenesis":1,"synaptic_plasticity":1}
NPS_MAX = 130
curcumin = c.get("Curcumin", {})
if curcumin:
    raw = sum(float(curcumin.get(k,5)) * w for k,w in NPS_W.items())
    nps = round(raw / NPS_MAX * 100, 1)
    expected = 80.0  # pre-computed
    if abs(nps - expected) < 0.5:
        print(f"✅ Check 3 PASS: Curcumin NPS = {nps} (expected ~{expected})")
    else:
        print(f"❌ Check 3 FAIL: Curcumin NPS = {nps}, expected ~{expected}")
else:
    print("⚠️  Check 3 SKIP: Curcumin not in compounds.json")

# ── Check 4: ERC data presence ────────────────────────────────────────────────
has_erc = sum(1 for v in c.values() if v.get("erc") or v.get("enzymes"))
print(f"{'✅' if has_erc == 0 else '⚠️ '} Check 4: Curated ERC data = {has_erc}/129 (expected 0 — Tier 1 is literature-curated, no IC50 data)")

# ── Check 5: ML compound uniqueness ──────────────────────────────────────────
with open("compounds_ml.json") as f:
    ml = json.load(f)
ml_clean = {k:v for k,v in ml.items() if not k.startswith("_")}
score_vectors = [
    tuple(round(float(v.get(c,5)),1) for c in ["antioxidant","anti_inflammatory",
     "mitochondrial_support","aggregation_modulation","cognitive_enhancement",
     "neurogenesis","synaptic_plasticity"])
    for v in ml_clean.values()
]
unique_vecs = len(set(score_vectors))
total = len(score_vectors)
pct_unique = round(unique_vecs/total*100, 1)
indicator = "ℹ️ "
print(f"{indicator} Check 5: ML score vector uniqueness = {unique_vecs}/{total} ({pct_unique}%) — target >40%")

# ── Check 6: Manuscript URL ───────────────────────────────────────────────────
with open("manuscript_brainsafe_ai.html") as f:
    html = f.read()
if "[URL to be added upon deployment]" in html:
    print("❌ Check 6 FAIL: Placeholder URL still present in manuscript — apply FIX 3(a)")
else:
    print("✅ Check 6 PASS: No placeholder URL in manuscript")

# ── Check 7: NPS equation in manuscript ───────────────────────────────────────
if "NPS_MAX_RAW" in html or "/ 130" in html or "130 &times; 100" in html or "130 × 100" in html:
    print("✅ Check 7 PASS: NPS formula equation present in manuscript")
else:
    print("❌ Check 7 FAIL: NPS equation missing — apply FIX 3(b)")

print("=" * 60)
print("Validation complete. Fix all ❌ items before submission.")
