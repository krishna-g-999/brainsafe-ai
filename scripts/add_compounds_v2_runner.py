"""
scripts/add_compounds_v2_runner.py
Adds the 360 additional compounds from additional_compounds_v2.py
to the existing master dataset and retrains.

ONLY RUN THIS AFTER master_expansion_pipeline.py completed successfully
and hold-out NPS R² >= 0.65.

Run from D:\BRAINSAFE_AI:
  D:\BRAINSAFE_AI\brainsafe_env\Scripts\python.exe -u scripts\add_compounds_v2_runner.py 2>&1 | Tee-Object logs\add_v2.log

Time: ~25 min (SMILES retrieval + retraining)
"""

import sys, json, time, logging, subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
LOG_DIR = ROOT / "logs"
LOG_DIR.mkdir(exist_ok=True)

# Force UTF-8 on stdout so Unicode log lines don't crash on Windows cp1252
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOG_DIR / "add_v2.log", encoding="utf-8", mode="w"),
    ]
)
log = logging.getLogger(__name__)

from additional_compounds_v2 import ADDITIONAL_COMPOUNDS_V2
from master_expansion_data import CLASS_PRIORS, NEGATIVE_CONTROL_SCORES, DIMENSION_COLS

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "BrainSafe-AI-v6-PhD/1.0"})


def validate_smiles(smi) -> bool:
    if not smi or str(smi).strip().lower() in ("nan","none","n/a",""):
        return False
    try:
        from rdkit import Chem
        return Chem.MolFromSmiles(str(smi).strip()) is not None
    except Exception:
        return False


def canonical(smi: str) -> str:
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smi)
        return Chem.MolToSmiles(mol, canonical=True) if mol else smi
    except Exception:
        return smi


def fetch_pubchem(name: str) -> str:
    try:
        url = (f"{PUBCHEM_URL}/compound/name/{requests.utils.quote(name)}"
               "/property/IsomericSMILES,CanonicalSMILES/JSON")
        r = SESSION.get(url, timeout=15)
        if r.status_code == 200:
            props = r.json().get("PropertyTable",{}).get("Properties",[])
            if props:
                return props[0].get("IsomericSMILES") or props[0].get("CanonicalSMILES","")
    except Exception:
        pass
    return ""


def estimate_bbb(smiles: str) -> str:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return "Medium"
        mw   = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd  = Descriptors.NumHDonors(mol)
        if mw <= 360 and 1.0 <= logp <= 3.0 and tpsa <= 60 and hbd <= 1:
            return "High"
        if mw <= 450 and 0.0 <= logp <= 4.0 and tpsa <= 90 and hbd <= 3:
            return "Medium"
        if mw <= 500 and tpsa <= 120:
            return "Low-Med"
        return "Low"
    except Exception:
        return "Medium"


def main():
    log.info("=" * 62)
    log.info(f"BrainSafe AI v6 — Adding {len(ADDITIONAL_COMPOUNDS_V2)} new compounds")
    log.info("=" * 62)

    # Load existing master dataset
    for fname in ["brainsafe_master.csv", "brainsafe_fixed.csv",
                  "brainsafe_training_set_325.csv"]:
        p = ROOT / "data" / fname
        if p.exists():
            df_existing = pd.read_csv(p)
            log.info(f"Loaded {len(df_existing)} from {fname}")
            break
    else:
        log.error("No existing dataset found. Run master_expansion_pipeline.py first.")
        sys.exit(1)

    # Build dedup sets
    exist_names = set(df_existing["name"].str.lower().str.strip())
    exist_smi = set()
    for s in df_existing["smiles"].fillna(""):
        if validate_smiles(s):
            exist_smi.add(canonical(s))

    WEIGHTS = [3,3,2,2,2,1,1]

    new_rows = []
    skipped = 0

    log.info(f"\nProcessing {len(ADDITIONAL_COMPOUNDS_V2)} compounds...")
    for i, c in enumerate(ADDITIONAL_COMPOUNDS_V2):
        name   = c["name"]
        smiles = c.get("smiles","")
        ctype  = c.get("type","general")
        diseases = c.get("diseases",[])

        # Dedup
        if name.lower() in exist_names:
            skipped += 1
            continue
        if validate_smiles(smiles) and canonical(smiles) in exist_smi:
            skipped += 1
            continue

        # Validate/fix SMILES
        if not validate_smiles(smiles):
            smiles = fetch_pubchem(name)
            time.sleep(0.35)
            if not validate_smiles(smiles):
                log.info(f"  [{i+1}] ✗ {name:<35} SMILES not found")
            else:
                log.info(f"  [{i+1}] ✓ {name:<35} retrieved from PubChem")
        else:
            log.info(f"  [{i+1}] ✓ {name:<35} SMILES ok")

        exist_names.add(name.lower())
        if validate_smiles(smiles):
            exist_smi.add(canonical(smiles))

        bbb = estimate_bbb(smiles) if validate_smiles(smiles) else "Medium"

        dis_map = {
            "ad_target_count":  int("alzheimers" in diseases),
            "pd_target_count":  int("parkinsons"  in diseases),
            "als_target_count": int("als"         in diseases),
            "hd_target_count":  int("huntingtons" in diseases),
        }

        row = {"name":name,"smiles":smiles,"compound_type":ctype,
               "bbb":bbb,"tier":"gold_curated","sample_weight":0.85,
               "data_source":"curated_v2", **dis_map}

        # Scores
        if name in NEGATIVE_CONTROL_SCORES:
            row.update(NEGATIVE_CONTROL_SCORES[name])
        else:
            prior_key = None
            for key in CLASS_PRIORS:
                if key in ctype.lower():
                    prior_key = key
                    break
            prior = CLASS_PRIORS.get(prior_key, CLASS_PRIORS["general"])
            for dim in DIMENSION_COLS:
                row[dim] = prior[dim]

        new_rows.append(row)

    log.info(f"\nNew unique compounds: {len(new_rows)}  (skipped {skipped} duplicates)")

    df_new = pd.DataFrame(new_rows)
    df_all = pd.concat([df_existing, df_new], ignore_index=True)

    # Fix scores
    for dim in DIMENSION_COLS:
        if dim in df_all.columns:
            df_all[dim] = pd.to_numeric(df_all[dim], errors="coerce").clip(0,100)
            if df_all[dim].isna().any():
                gm = df_all[dim].median()
                df_all[dim].fillna(gm if not np.isnan(gm) else 30.0, inplace=True)

    # NPS
    df_all["nps"] = (
        df_all[DIMENSION_COLS].fillna(0).mul(WEIGHTS).sum(axis=1)
        / (sum(WEIGHTS)*100) * 100
    ).clip(0,100).round(1)

    # SMILES coverage
    valid = df_all["smiles"].fillna("").apply(validate_smiles)
    pct = round(100*valid.sum()/len(df_all),1)

    log.info(f"\n── DATASET SUMMARY ────────────────────────────────────────")
    log.info(f"  Total compounds: {len(df_all)}")
    log.info(f"  SMILES coverage: {valid.sum()}/{len(df_all)} ({pct}%)")
    log.info(f"  NPS: mean={df_all['nps'].mean():.1f} ± {df_all['nps'].std():.1f}  "
             f"range={df_all['nps'].min():.1f}–{df_all['nps'].max():.1f}")
    if "compound_type" in df_all.columns:
        for ct, cnt in df_all["compound_type"].value_counts().head(15).items():
            log.info(f"  {ct:<25} n={cnt}")

    # Save
    out = ROOT / "data" / "brainsafe_master_v2.csv"
    df_all.to_csv(out, index=False)

    # Update SMILES list
    smiles_list = df_all[valid]["smiles"].tolist()
    with open(ROOT/"models_v5"/"training_smiles.json","w") as f:
        json.dump(smiles_list, f)

    log.info(f"\n  ✓ Saved: {out}")
    log.info(f"  ✓ Training SMILES: {len(smiles_list)}")

    # Retrain
    log.info(f"\n{'='*62}")
    log.info(f"RETRAINING on {len(df_all)}-compound master v2 dataset")
    log.info(f"{'='*62}\n")

    python = sys.executable
    cmd = [python, "-u", str(ROOT/"ml_v5_training.py"),
           "--data", str(out), "--out", str(ROOT/"models_v5"/"")]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace", cwd=str(ROOT)
    )
    for line in proc.stdout:
        print(line, end="", flush=True)
    proc.wait()

    if proc.returncode == 0:
        rp = ROOT/"models_v5"/"validation_report.json"
        if rp.exists():
            with open(rp) as f:
                r = json.load(f)
            log.info(f"\n{'='*62}")
            log.info("FINAL RESULTS — PASTE THIS BLOCK TO CLAUDE")
            log.info(f"{'='*62}")
            log.info(f"  Total compounds:      {len(df_all)}")
            log.info(f"  n_train:              {r.get('n_train')}")
            log.info(f"  n_holdout:            {r.get('n_holdout')}")
            log.info(f"  LOO-CV  NPS R²:       {r.get('loo_r2_nps',0):.3f}")
            log.info(f"  Hold-out NPS R²:      {r.get('holdout_r2_nps',0):.3f}")
            log.info(f"  Hold-out Spearman ρ:  {r.get('holdout_spearman',0):.3f}")
            log.info(f"\n  Per-dimension hold-out R²:")
            for dim, val in r.get("holdout_r2_per_dim",{}).items():
                ok = "✓" if val >= 0.40 else "✗"
                log.info(f"    {ok} {dim:<30} {val:.3f}")
            log.info(f"{'='*62}")
    else:
        log.error(f"Training failed (exit {proc.returncode}) — paste log to Claude")


if __name__ == "__main__":
    main()
