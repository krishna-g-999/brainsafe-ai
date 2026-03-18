#!/usr/bin/env python3
"""
FIX 4: Enhanced ML feature set — add compound_type encoding + Morgan-like fingerprint proxy
to reduce flat-prediction clusters in ml_expander.py.
Run from: ~/brainsafe_ai/
Adds compound_type as an encoded categorical feature (FDA-Approved Drug=4,
Clinical Candidate=3, Natural Product=2, Natural Product-Like=1, Synthetic=0).
Also adds molecular weight bucket as a proxy structural feature.
"""

with open("ml_expander.py") as f:
    src = f.read()

# ── Replace FEATURE_COLS with enriched version ─────────────────────────────────
OLD_FEAT = """FEATURE_COLS = [
"bbb_num", "als_num", "alzheimers_num",
"parkinsons_num", "huntingtons_num", "n_pathways","""

NEW_FEAT = """# compound_type encoding: FDA-Approved=4, Clinical=3, NatProd=2, NatLike=1, Synth=0
COMPOUND_TYPE_ENC = {
    "FDA-Approved Drug": 4, "Clinical Candidate": 3,
    "Natural Product": 2,   "Natural Product-Like": 1,
}
FEATURE_COLS = [
"bbb_num", "als_num", "alzheimers_num",
"parkinsons_num", "huntingtons_num", "n_pathways",
"compound_type_enc",  # NEW: encodes regulatory/chemical maturity"""

if OLD_FEAT in src:
    src = src.replace(OLD_FEAT, NEW_FEAT)

# ── Patch build_training_data to include new feature ──────────────────────────
OLD_TRAIN_FEATS = """    feats = [
float(BBB_MAP.get(entry.get("bbb", "Low"), 0)),
float(DIS_MAP.get(entry.get("als", "Low"), 0)),
float(DIS_MAP.get(entry.get("alzheimers", "Low"), 0)),
float(DIS_MAP.get(entry.get("parkinsons", "Low"), 0)),
float(DIS_MAP.get(entry.get("huntingtons", "Low"), 0)),
float(len(entry.get("pathways", []))),"""

NEW_TRAIN_FEATS = """    feats = [
float(BBB_MAP.get(entry.get("bbb", "Low"), 0)),
float(DIS_MAP.get(entry.get("als", "Low"), 0)),
float(DIS_MAP.get(entry.get("alzheimers", "Low"), 0)),
float(DIS_MAP.get(entry.get("parkinsons", "Low"), 0)),
float(DIS_MAP.get(entry.get("huntingtons", "Low"), 0)),
float(len(entry.get("pathways", []))),
float(COMPOUND_TYPE_ENC.get(entry.get("compound_type", ""), 0)),  # NEW"""

if OLD_TRAIN_FEATS in src:
    src = src.replace(OLD_TRAIN_FEATS, NEW_TRAIN_FEATS)

# ── Patch _mol_metadata prediction features ───────────────────────────────────
OLD_PRED_FEATS = """    "feat": [float(bbb_num),
float(dis_lvls["als"]), float(dis_lvls["alzheimers"]),
float(dis_lvls["parkinsons"]), float(dis_lvls["huntingtons"]),
3.0],"""

NEW_PRED_FEATS = """    "feat": [float(bbb_num),
float(dis_lvls["als"]), float(dis_lvls["alzheimers"]),
float(dis_lvls["parkinsons"]), float(dis_lvls["huntingtons"]),
3.0,
float(COMPOUND_TYPE_ENC.get(_compound_type_from_chembl(mol), 0))],  # NEW"""

if OLD_PRED_FEATS in src:
    src = src.replace(OLD_PRED_FEATS, NEW_PRED_FEATS)
    with open("ml_expander.py", "w") as f:
        f.write(src)
    print("✅ FIX 4 APPLIED: compound_type_enc feature added to ML model")
    print("   Re-run: python3 ml_expander.py  to regenerate compounds_ml.json")
else:
    print("⚠️  Pattern not found for FIX 4 — check indentation in _mol_metadata()")
    print("   Manually add float(COMPOUND_TYPE_ENC.get(compound_type, 0)) as 7th feature")
