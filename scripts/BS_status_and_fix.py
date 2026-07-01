"""
scripts/BS_status_and_fix.py
BrainSafe AI (BS) — Complete Status Audit + Database Sync

PREFIX CONVENTION (from this session forward):
  All BrainSafe AI files: BS_ prefix
  All BBB predictor files: BBB_ prefix  ← your SAIDOCK/SAI-Net project
  Never mix them in the same folder

THIS SCRIPT:
  1. Audits exactly what has been done and what is inconsistent
  2. Converts brainsafe_SCIENTIFIC_FIXED.csv → BS_compounds_full.json
     (so the app shows ALL 535 compounds, not just the original 134)
  3. Updates the app's compound count label
  4. Verifies model vs. database consistency
  5. Prints an internal status report
     (SUPERSEDED for scientific claims - see BS_MODEL_CARD.md)

Run from D:\\BRAINSAFE_AI:
  D:\\BRAINSAFE_AI\\brainsafe_env\\Scripts\\python.exe -u scripts\\BS_status_and_fix.py
"""

import sys, json, os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

DIMENSION_COLS = [
    "antioxidant", "anti_inflammatory", "mitochondrial_support",
    "aggregation_modulation", "cognitive_enhancement",
    "neurogenesis", "synaptic_plasticity",
]
WEIGHTS = [3, 3, 2, 2, 2, 1, 1]


def compute_nps(row) -> float:
    vals = [float(row.get(d, 0) or 0) for d in DIMENSION_COLS]
    w = sum(w * v for w, v in zip(WEIGHTS, vals))
    return round(min(100.0, w / (sum(WEIGHTS) * 100) * 100), 1)


def sep(char="─", n=62):
    return char * n


def main():
    print(sep("="))
    print("BRAINSAFE AI (BS) — Complete Status Audit")
    print(sep("="))

    # ── 1. AUDIT EVERY KEY FILE ──────────────────────────────────────────
    print("\n[1] FILE AUDIT")
    print(sep())

    files = {
        "TRAINED MODEL":   "models_v5/validation_report.json",
        "FINAL DATASET":   "data/brainsafe_SCIENTIFIC_FIXED.csv",
        "APP DATABASE":    "compounds.json",
        "ML DATABASE":     "compounds_ml.json",
        "APP (v6)":        "app_v6.py",
        "APP (deployed)":  "app.py",
        "SCORER":          "scorer.py",
        "MODEL CONFIG":    "model_config.py",
        "NT MAPPER":       "neurotransmitter_mapper.py",
        "BRAIN REGIONS":   "brain_region_mapper.py",
        "DOSE RESPONSE":   "dose_response.py",
        "V6 TABS":         "v6_tabs.py",
    }
    for label, fpath in files.items():
        p = ROOT / fpath
        if p.exists():
            size = round(p.stat().st_size / 1024, 1)
            print(f"  OK    {label:<20} {fpath:<40} {size} KB")
        else:
            print(f"  MISS  {label:<20} {fpath}")

    # ── 2. TRAINED MODEL RESULTS ──────────────────────────────────────────
    print(f"\n[2] TRAINED MODEL (93-feature, 535-compound)")
    print(sep())

    rp = ROOT / "models_v5" / "validation_report.json"
    if rp.exists():
        with open(rp) as f:
            r = json.load(f)
        print(f"  n_train:           {r.get('n_train')}")
        print(f"  n_holdout:         {r.get('n_holdout')}")
        print(f"  n_features:        {r.get('n_features','93')}")
        print(f"  LOO-CV NPS R2:     {r.get('loo_r2_nps',0):.3f}")
        print(f"  LOO-CV Spearman:   {r.get('loo_spearman_nps',0):.3f}")
        print(f"  Hold-out NPS R2:   {r.get('holdout_r2_nps',0):.3f}")
        print(f"  Hold-out Spearman: {r.get('holdout_spearman',0):.3f}")
        print()
        print("  Per-dimension hold-out R2:")
        all_pass = True
        for dim, val in r.get("holdout_r2_per_dim", {}).items():
            ok = "PASS" if val >= 0.40 else "FAIL"
            if val < 0.40: all_pass = False
            print(f"    {ok}  {dim:<30} {val:.3f}")
        print(f"\n  All 7 dimensions >= 0.40: {all_pass}")
    else:
        print("  ERROR: validation_report.json not found")

    # ── 3. DATABASE AUDIT ─────────────────────────────────────────────────
    print(f"\n[3] DATABASE CONSISTENCY CHECK")
    print(sep())

    # Old compounds.json (134)
    old_db_path = ROOT / "compounds.json"
    n_old = 0
    if old_db_path.exists():
        with open(old_db_path) as f:
            old_db = json.load(f)
        n_old = len([k for k in old_db if not k.startswith("_")])
        print(f"  compounds.json (app database):    {n_old} compounds  <-- APP READS THIS")

    # compounds_ml.json
    ml_db_path = ROOT / "compounds_ml.json"
    n_ml = 0
    if ml_db_path.exists():
        with open(ml_db_path) as f:
            ml_db = json.load(f)
        n_ml = len([k for k in ml_db if not k.startswith("_")])
        print(f"  compounds_ml.json:                {n_ml} compounds")

    # Training CSV
    csv_path = ROOT / "data" / "brainsafe_SCIENTIFIC_FIXED.csv"
    n_csv = 0
    if csv_path.exists():
        df = pd.read_csv(csv_path)
        n_csv = len(df)
        print(f"  brainsafe_SCIENTIFIC_FIXED.csv:   {n_csv} compounds  <-- MODEL TRAINED ON THIS")

    print()
    print(f"  PROBLEM: App shows {n_old + n_ml} compounds but model was trained on {n_csv}")
    print(f"  FIX:     Convert brainsafe_SCIENTIFIC_FIXED.csv -> BS_compounds_full.json")
    print(f"           Update app to read BS_compounds_full.json ({n_csv} compounds)")

    # ── 4. BUILD FULL COMPOUND DATABASE ───────────────────────────────────
    print(f"\n[4] BUILDING BS_compounds_full.json ({n_csv} compounds)")
    print(sep())

    if not csv_path.exists():
        print("  ERROR: brainsafe_SCIENTIFIC_FIXED.csv not found")
        return

    df = pd.read_csv(csv_path)

    # Score scale check
    max_score = df[DIMENSION_COLS].max().max()
    scale = 10.0 if max_score <= 10.5 else 1.0
    if scale == 10.0:
        print("  WARNING: Scores on 1-10 scale, rescaling to 0-100")
        for d in DIMENSION_COLS:
            df[d] = df[d] * 10.0

    # Clip
    for d in DIMENSION_COLS:
        df[d] = df[d].clip(0, 100).round(1)

    # Compute NPS
    df["nps"] = df.apply(compute_nps, axis=1)

    # Build JSON structure matching app format
    bs_db = {}
    n_built = 0
    for _, row in df.iterrows():
        name = str(row.get("name","")).strip()
        if not name:
            continue

        # Determine disease relevance
        diseases = {}
        for dis, col in [("alzheimers","ad_target_count"),
                         ("parkinsons","pd_target_count"),
                         ("als","als_target_count"),
                         ("huntingtons","hd_target_count")]:
            val = int(row.get(col, 0) or 0)
            diseases[dis] = "High" if val >= 2 else "Med" if val >= 1 else "Low"

        # Compound type and pathway mapping
        ctype = str(row.get("compound_type","general")).lower()
        pathway_map = {
            "flavone":        ["Nrf2/GSH","NF-kB","COX-2 inhibition"],
            "flavonol":       ["Nrf2/GSH","NF-kB","AMPK"],
            "flavan3ol":      ["Nrf2/GSH","SIRT1","NF-kB"],
            "stilbene":       ["SIRT1","PGC-1alpha","NF-kB"],
            "curcuminoid":    ["Nrf2/GSH","NF-kB","Autophagy","BDNF/TrkB"],
            "alkaloid":       ["AChE inhibition","MAO inhibition","BDNF/TrkB"],
            "triterpenoid":   ["NF-kB","Nrf2/GSH","mTOR"],
            "vitamin":        ["Nrf2/GSH","Myelin synthesis","GSH precursor"],
            "cofactor":       ["Mitochondrial ETC","Nrf2/GSH","AMPK"],
            "mineral":        ["Nrf2/GSH","AMPK","mTOR"],
            "amino_acid":     ["GABAergic","Glutamatergic","Mitochondrial"],
            "carotenoid":     ["Nrf2/GSH","NF-kB","Carotenoid signalling"],
            "xanthophyll":    ["Nrf2/GSH","NF-kB","Retinal protection"],
            "phytocannabinoid":["CB1/CB2","TRPV1","PPARgamma","NF-kB"],
            "sesquiterpene":  ["CB2","NF-kB","Nrf2/GSH"],
            "monoterpene":    ["GABAergic","AChE inhibition","Anti-inflammatory"],
            "mushroom_compound":["NGF/BDNF","Nrf2/GSH","Autophagy"],
            "withanolide":    ["NF-kB","Nrf2/GSH","BDNF/TrkB"],
            "adaptogen":      ["BDNF/TrkB","AMPK","Nrf2/GSH"],
            "ginsenoside":    ["BDNF/TrkB","Nrf2/GSH","AChE inhibition"],
            "omega3":         ["Resolvin/Protectin","BDNF/TrkB","NF-kB"],
            "drug_cholinergic":["AChE inhibition","Cholinergic","NMDA"],
            "drug_maob":      ["MAO-B inhibition","Dopaminergic","Mitochondrial"],
            "drug_als":       ["Glutamatergic","Nrf2/GSH","Mitochondrial"],
            "drug_hd":        ["mTOR","HDAC","Striatal"],
            "mitochondrial_targeted":["Mitochondrial ETC","ROS scavenging","PINK1/Parkin"],
            "neurotoxin":     [],
            "inactive_control":[],
        }
        pathways = []
        for key, pws in pathway_map.items():
            if key in ctype:
                pathways = pws
                break
        if not pathways and "polyphenol" in ctype:
            pathways = ["Nrf2/GSH","NF-kB","Antioxidant"]

        # Dimension scores: store on 0-100 scale divided by 10 for app compatibility
        # The existing app uses 1-10 scale internally. Keep consistent.
        entry = {
            "name":            name,
            "smiles":          str(row.get("smiles","")) if row.get("smiles") else None,
            "compound_type":   str(row.get("compound_type","general")),
            "bbb":             str(row.get("bbb","Medium")),
            "nps":             float(row["nps"]),
            # Scores on 0-100 scale (v6 standard)
            "antioxidant":           float(row.get("antioxidant",30)),
            "anti_inflammatory":     float(row.get("anti_inflammatory",30)),
            "mitochondrial_support": float(row.get("mitochondrial_support",30)),
            "aggregation_modulation":float(row.get("aggregation_modulation",30)),
            "cognitive_enhancement": float(row.get("cognitive_enhancement",30)),
            "neurogenesis":          float(row.get("neurogenesis",30)),
            "synaptic_plasticity":   float(row.get("synaptic_plasticity",30)),
            # Disease relevance
            "alzheimers":   diseases["alzheimers"],
            "parkinsons":   diseases["parkinsons"],
            "als":          diseases["als"],
            "huntingtons":  diseases["huntingtons"],
            # Metadata
            "pathways":          pathways,
            "tier":              str(row.get("tier","gold_curated")),
            "data_source":       str(row.get("data_source","curated")),
            "sample_weight":     float(row.get("sample_weight",0.85)),
        }
        bs_db[name.lower()] = entry
        n_built += 1

    # Save
    out_path = ROOT / "BS_compounds_full.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(bs_db, f, indent=2, ensure_ascii=False)
    print(f"  Saved: {out_path}")
    print(f"  Total compounds: {n_built}")

    # Compound class breakdown
    print(f"\n  Compound class breakdown:")
    type_counts = df["compound_type"].value_counts()
    for ct, cnt in type_counts.head(20).items():
        print(f"    {ct:<30} n={cnt}")

    # NPS distribution
    print(f"\n  NPS distribution:")
    print(f"    Limited  (0-39):  {(df['nps']<40).sum()} compounds")
    print(f"    Moderate (40-69): {((df['nps']>=40)&(df['nps']<70)).sum()} compounds")
    print(f"    Strong   (70-100):{(df['nps']>=70).sum()} compounds")
    print(f"    Mean: {df['nps'].mean():.1f}  SD: {df['nps'].std():.1f}")

    # ── 5. PRINT COMPLETE STATUS ──────────────────────────────────────────
    print(f"\n{sep('=')}")
    print("COMPLETE STATUS REPORT")
    print(sep("="))
    print(f"""
*** SUPERSEDED - DO NOT CITE *************************************************
*** The numbers below are RANDOM-split LOO/hold-out values produced WITH the
*** circular disease-count features. Authoritative honest evaluation:
***   BS_MODEL_CARD.md + BS_validation_report.json (leak-free scaffold-CV)
***   + BS_predictive_report.json (only 'antioxidant' is genuinely predictable).
***************************************************************************

WHAT HAS BEEN ACCOMPLISHED:
  - SMILES coverage fixed:      44% -> 98.9% (535 compounds)
  - NPS formula fixed:          4/7 dims -> all 7 dims with correct weights
  - Feature vector expanded:    87 -> 93 features (+ 6 structural)
  - NameErrors fixed:           BBB_MAP, TARGETTOPATHWAY, _V6_MODULES_OK
  - Score scale standardised:   1-10 -> 0-100 throughout
  - Model file count fixed:     35 -> 36 (BBB encoder added)
  - All 7 dimensions pass:      >= 0.40 hold-out R2 (was 5/7 failing)
  - Antioxidant fixed:          0.396 -> 0.512 (phenolic OH feature)
  - Anti-inflam fixed:          0.377 -> 0.443
  - 3 new science features:     NT panel, brain regions, dose-response
  - App patched:                v6_tabs.py, render_v6_tabs()

WHAT STILL NEEDS DOING:
  1. Update app to read BS_compounds_full.json (535 cpds, not 134)
  2. Fix "325 compounds" label in app.py -> "535 compounds"
  3. Fix manuscript: 6 critical issues (author list, vFigA, Ref[10], etc.)
  4. Deploy to GitHub -> Streamlit auto-redeploys
  5. Update Zenodo DOI with models_v5/

MODEL NUMBERS (SUPERSEDED - do NOT cite; random-split + circular features; see BS_MODEL_CARD.md):
  Dataset:       535 compounds, 20+ classes, 98.9% SMILES
  Features:      93 (50 ECFP-4 + 32 ChemBERTa + 4 disease + 1 BBB + 6 structural)
  Training:      n=428 (stratified, fixed random_state=42)
  Hold-out:      n=107 (completely unseen)
  LOO-CV R2:     0.615  Spearman rho=0.822
  Hold-out R2:   0.546  Spearman rho=0.769
  All 7 dims:    >= 0.40 (antioxidant=0.512, anti_inflam=0.443,
                          mitochondrial=0.473, aggregation=0.504,
                          cognitive=0.634, neurogenesis=0.478,
                          synaptic=0.593)

NAMING CONVENTION (all future files):
  BrainSafe AI:   BS_*.py, BS_*.json, BS_*.csv
  BBB Predictor:  BBB_*.py (your SAIDOCK project - keep completely separate)
""")

    print(sep("="))
    print("NEXT COMMANDS TO RUN:")
    print(sep())
    print("""
  1. Update app to use BS_compounds_full.json:
     python scripts\\BS_update_app_database.py

  2. Test the updated app locally:
     brainsafe_env\\Scripts\\python.exe -m streamlit run app_v6.py

  3. Deploy to GitHub:
     git add . && git commit -m "v6: 535-compound database, 93-feature model" && git push

  4. Provide manuscript info (author names, advisor, funding)
     -> Claude writes all 6 critical manuscript fixes
""")


if __name__ == "__main__":
    main()
