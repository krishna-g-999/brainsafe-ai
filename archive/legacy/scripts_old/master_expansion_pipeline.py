"""
scripts/master_expansion_pipeline.py
BrainSafe AI v6 — Master Dataset Expansion + SMILES Fix + Retrain

ONE SCRIPT. ONE COMMAND. EVERYTHING.

Does in sequence:
  Step 1. Load existing 325-compound gold standard
  Step 2. Retrieve missing SMILES for all 201 compounds (PubChem + ChEMBL)
  Step 3. Add 420 new compounds (master_expansion_data.py) with validated SMILES
  Step 4. Deduplicate by canonical SMILES
  Step 5. Fix all NaN dimension scores (negatives fixed, class priors for gaps)
  Step 6. Compute NPS for all compounds
  Step 7. Validate dataset (NAR quality checks)
  Step 8. Save data/brainsafe_master.csv  (~750-850 compounds)
  Step 9. Retrain full 4-model ensemble on expanded set

Expected final dataset:  ~800 compounds, ≥92% SMILES, NPS R² ≥ 0.80

Run from D:\BRAINSAFE_AI:
  D:\BRAINSAFE_AI\brainsafe_env\Scripts\python.exe -u scripts\master_expansion_pipeline.py 2>&1 | Tee-Object logs\master_pipeline.log

Time estimate:
  SMILES retrieval: 20-30 min (rate-limited API)
  Training:         25-45 min (ChemBERTa + LOO-CV)
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

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOG_DIR / "master_pipeline.log",
                            encoding="utf-8", mode="w"),
    ]
)
log = logging.getLogger(__name__)

# Import compound library
from master_expansion_data import (
    MASTER_NEW_COMPOUNDS, NEGATIVE_CONTROL_SCORES,
    CLASS_PRIORS, DIMENSION_COLS
)

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL_URL  = "https://www.ebi.ac.uk/chembl/api/data"
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "BrainSafe-AI-v6-PhD/1.0"})


# ─── UTILITIES ────────────────────────────────────────────────────────────────

def validate_smiles(smiles) -> bool:
    if not smiles or str(smiles).strip().lower() in ("nan","none","n/a",""):
        return False
    try:
        from rdkit import Chem
        return Chem.MolFromSmiles(str(smiles).strip()) is not None
    except Exception:
        return False


def canonical_smiles(smiles: str) -> str:
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol, canonical=True) if mol else smiles
    except Exception:
        return smiles


def get_mw(smiles: str) -> float:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        mol = Chem.MolFromSmiles(smiles)
        return round(Descriptors.MolWt(mol), 1) if mol else 350.0
    except Exception:
        return 350.0


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


def api_get(url, params=None, retries=3) -> dict:
    for i in range(retries):
        try:
            r = SESSION.get(url, params=params, timeout=15)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return {}
        except Exception:
            pass
        time.sleep(1.5 ** i)
    return {}


# ─── STEP 2: SMILES RETRIEVAL ─────────────────────────────────────────────────

def fetch_smiles_pubchem(name: str) -> str:
    data = api_get(
        f"{PUBCHEM_URL}/compound/name/{requests.utils.quote(name)}"
        "/property/IsomericSMILES,CanonicalSMILES/JSON"
    )
    props = data.get("PropertyTable", {}).get("Properties", [])
    if props:
        return props[0].get("IsomericSMILES") or props[0].get("CanonicalSMILES", "")
    return ""


def fetch_smiles_chembl(name: str) -> str:
    # Exact match
    data = api_get(f"{CHEMBL_URL}/molecule.json",
                   {"pref_name__iexact": name, "limit": 1, "format": "json"})
    mols = data.get("molecules", [])
    if mols and mols[0] is not None:
        s = mols[0].get("molecule_structures") or {}
        s = s.get("canonical_smiles", "") if isinstance(s, dict) else ""
        if s:
            return s
    # Fuzzy match
    first_word = name.split()[0]
    data = api_get(f"{CHEMBL_URL}/molecule.json",
                   {"pref_name__icontains": first_word, "limit": 5, "format": "json"})
    for mol in (data.get("molecules", []) or []):
        if mol is None:
            continue
        pref = (mol.get("pref_name") or "").lower()
        if name.lower() in pref or pref in name.lower():
            structs = mol.get("molecule_structures") or {}
            s = structs.get("canonical_smiles", "") if isinstance(structs, dict) else ""
            if s:
                return s
    return ""


def retrieve_smiles_for_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    missing = df["smiles"].fillna("").apply(lambda s: not validate_smiles(s))
    n_miss  = missing.sum()
    log.info(f"  Compounds needing SMILES: {n_miss}/{len(df)}")
    retrieved, failed = 0, []

    for idx, row in df[missing].iterrows():
        name = str(row.get("name", "")).strip()
        if not name:
            continue

        smi = fetch_smiles_pubchem(name)
        time.sleep(0.35)
        if not validate_smiles(smi):
            smi = fetch_smiles_chembl(name)
            time.sleep(0.25)

        if validate_smiles(smi):
            df.at[idx, "smiles"] = smi
            retrieved += 1
            log.info(f"    ✓ {name:<35} {smi[:45]}...")
        else:
            failed.append(name)
            log.info(f"    ✗ {name:<35} not found")

    log.info(f"  Retrieved: {retrieved}  |  Still missing: {len(failed)}")
    return df


# ─── STEP 4: DEDUPLICATION ────────────────────────────────────────────────────

def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    before = len(df)
    # Name dedup (case-insensitive)
    df["_name_lower"] = df["name"].str.lower().str.strip()
    df = df.drop_duplicates(subset="_name_lower", keep="first")
    # SMILES dedup (canonical)
    valid_smi = df["smiles"].fillna("").apply(validate_smiles)
    df.loc[valid_smi, "_canonical"] = df.loc[valid_smi, "smiles"].apply(canonical_smiles)
    dup_smi = df["_canonical"].notna() & df["_canonical"].duplicated()
    df = df[~dup_smi]
    df.drop(columns=["_name_lower", "_canonical"], errors="ignore", inplace=True)
    df.reset_index(drop=True, inplace=True)
    log.info(f"  Removed {before - len(df)} duplicates → {len(df)} unique compounds")
    return df


# ─── STEP 5: NaN SCORE FIX ────────────────────────────────────────────────────

def fix_scores(df: pd.DataFrame) -> pd.DataFrame:
    # Apply named negative control scores
    for name, scores in NEGATIVE_CONTROL_SCORES.items():
        mask = df["name"].str.lower() == name.lower()
        if mask.any():
            for dim, val in scores.items():
                if dim in df.columns:
                    df.loc[mask, dim] = val
            log.info(f"    Fixed negative: {name}")

    # Apply class priors for NaN dimensions
    if "compound_type" in df.columns:
        for comp_type, prior in CLASS_PRIORS.items():
            mask = df["compound_type"].str.lower().str.contains(
                comp_type, na=False, regex=False
            )
            if not mask.any():
                continue
            for dim in DIMENSION_COLS:
                if dim not in df.columns:
                    continue
                nan_mask = mask & df[dim].isna()
                if nan_mask.any():
                    df.loc[nan_mask, dim] = prior[dim]

    # Last-resort: global median for any remaining NaN
    for dim in DIMENSION_COLS:
        if dim in df.columns and df[dim].isna().any():
            gm = df[dim].median()
            df[dim].fillna(gm if not np.isnan(gm) else 30.0, inplace=True)

    # Clip to 0-100
    for dim in DIMENSION_COLS:
        if dim in df.columns:
            df[dim] = pd.to_numeric(df[dim], errors="coerce").clip(0.0, 100.0)

    remaining = df[DIMENSION_COLS].isna().sum().sum()
    log.info(f"    NaN after fix: {remaining}")
    return df


# ─── STEP 6: COMPUTE NPS ──────────────────────────────────────────────────────

def compute_nps(df: pd.DataFrame) -> pd.DataFrame:
    WEIGHTS = [3, 3, 2, 2, 2, 1, 1]
    df["nps"] = (
        df[DIMENSION_COLS].fillna(0).mul(WEIGHTS).sum(axis=1)
        / (sum(WEIGHTS) * 100) * 100
    ).clip(0, 100).round(1)
    return df


# ─── STEP 7: VALIDATION ───────────────────────────────────────────────────────

def validate_dataset(df: pd.DataFrame):
    n = len(df)
    valid_smi = df["smiles"].fillna("").apply(validate_smiles)
    pct_smi = round(100 * valid_smi.sum() / n, 1)
    log.info(f"\n── DATASET VALIDATION ─────────────────────────────────")
    log.info(f"  Total compounds:      {n}")
    log.info(f"  SMILES coverage:      {valid_smi.sum()}/{n} ({pct_smi}%)")
    log.info(f"  NPS mean ± sd:        {df['nps'].mean():.1f} ± {df['nps'].std():.1f}")
    log.info(f"  NPS range:            {df['nps'].min():.1f} – {df['nps'].max():.1f}")

    log.info(f"\n  Score distributions (must all be 0-100 scale):")
    for dim in DIMENSION_COLS:
        col = df[dim].dropna()
        ok = col.max() > 10
        log.info(f"    {'✓' if ok else '✗'} {dim:<30} mean={col.mean():.1f}  "
                 f"max={col.max():.1f}  n={len(col)}")

    log.info(f"\n  Negative control check (must be < 20):")
    negs = df[df["compound_type"].isin(["neurotoxin","inactive_control"])]
    for _, r in negs.iterrows():
        mean_s = np.nanmean([r.get(d, 0) for d in DIMENSION_COLS])
        status = "✓" if mean_s < 20 else "✗ TOO HIGH"
        log.info(f"    {status}  {r['name']:<25} mean={mean_s:.1f}")

    log.info(f"\n  Compound class distribution:")
    if "compound_type" in df.columns:
        for ct, cnt in df["compound_type"].value_counts().items():
            log.info(f"    {ct:<25} n={cnt}")

    if pct_smi < 80:
        log.warning(f"\n  ⚠️  SMILES coverage {pct_smi}% is below 80%.")
        log.warning(f"  Model performance will be limited.")
        log.warning(f"  Consider manual SMILES curation for the {n - valid_smi.sum()} missing.")
    else:
        log.info(f"\n  ✓ Dataset quality: GOOD ({pct_smi}% SMILES coverage)")

    return pct_smi


# ─── STEP 9: RETRAIN ──────────────────────────────────────────────────────────

def retrain(data_path: str):
    python = sys.executable
    cmd = [python, "-u", str(ROOT / "ml_v5_training.py"),
           "--data", data_path, "--out", str(ROOT / "models_v5" / "")]
    log.info(f"\n{'='*62}")
    log.info("STARTING TRAINING ON EXPANDED DATASET")
    log.info(f"{'='*62}\n")

    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace", cwd=str(ROOT)
    )
    for line in proc.stdout:
        print(line, end="", flush=True)
    proc.wait()

    if proc.returncode != 0:
        log.error(f"Training failed (exit code {proc.returncode})")
        log.error("Paste the full output above to Claude chat.")
        sys.exit(1)

    # Print final report
    rp = ROOT / "models_v5" / "validation_report.json"
    if rp.exists():
        with open(rp) as f:
            r = json.load(f)
        log.info(f"\n{'='*62}")
        log.info("FINAL TRAINING RESULTS — PASTE THIS TO CLAUDE")
        log.info(f"{'='*62}")
        log.info(f"  n_train:              {r.get('n_train')}")
        log.info(f"  n_holdout:            {r.get('n_holdout')}")
        log.info(f"  n_silver:             {r.get('n_silver',0)}")
        log.info(f"  LOO-CV  NPS R²:       {r.get('loo_r2_nps',0):.3f}")
        log.info(f"  Hold-out NPS R²:      {r.get('holdout_r2_nps',0):.3f}")
        log.info(f"  Hold-out Spearman ρ:  {r.get('holdout_spearman',0):.3f}")
        log.info(f"  n_features:           {r.get('n_features')}")
        log.info(f"\n  Per-dimension hold-out R²:")
        for dim, val in r.get("holdout_r2_per_dim", {}).items():
            ok = "✓" if val >= 0.40 else "✗"
            log.info(f"    {ok} {dim:<30} {val:.3f}")
        log.info(f"{'='*62}")


# ─── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    log.info("=" * 62)
    log.info("BrainSafe AI v6 — Master Expansion Pipeline")
    log.info(f"Target: ~850 compounds, ≥92% SMILES, NPS R²≥0.80")
    log.info("=" * 62)

    # ── Step 1: Load existing 325 ───────────────────────────────────────
    log.info("\n[1/9] Loading existing gold standard...")
    for fname in ["brainsafe_training_set_325.csv",
                  "brainsafe_expanded.csv",
                  "brainsafe_training_set.csv"]:
        p = ROOT / "data" / fname
        if p.exists():
            df_existing = pd.read_csv(p)
            log.info(f"  Loaded {len(df_existing)} from {fname}")
            break
    else:
        log.error("No existing CSV found. Run generate_training_data.py first.")
        sys.exit(1)

    # Ensure smiles column exists
    if "smiles" not in df_existing.columns:
        df_existing["smiles"] = ""

    # ── Step 2: Retrieve missing SMILES for existing 325 ────────────────
    log.info(f"\n[2/9] Retrieving SMILES for existing compounds...")
    df_existing = retrieve_smiles_for_dataframe(df_existing)

    # ── Step 3: Add new compounds ────────────────────────────────────────
    log.info(f"\n[3/9] Adding {len(MASTER_NEW_COMPOUNDS)} new compounds...")

    # Build set of existing names + canonical SMILES for deduplication
    exist_names = set(df_existing["name"].str.lower().str.strip())
    exist_smi   = set()
    for s in df_existing["smiles"].dropna():
        if validate_smiles(s):
            exist_smi.add(canonical_smiles(s))

    new_rows = []
    skipped  = 0
    for c in MASTER_NEW_COMPOUNDS:
        name  = c["name"]
        smiles = c.get("smiles", "")
        ctype  = c.get("type", "general")
        diseases = c.get("diseases", [])

        # Dedup check
        if name.lower() in exist_names:
            skipped += 1
            continue
        if validate_smiles(smiles):
            csmi = canonical_smiles(smiles)
            if csmi in exist_smi:
                skipped += 1
                continue
            exist_smi.add(csmi)
        exist_names.add(name.lower())

        # Disease target counts
        dis_map = {"ad_target_count":  int("alzheimers" in diseases),
                   "pd_target_count":  int("parkinsons"  in diseases),
                   "als_target_count": int("als"         in diseases),
                   "hd_target_count":  int("huntingtons" in diseases)}

        # BBB from SMILES
        bbb = estimate_bbb(smiles) if validate_smiles(smiles) else "Medium"

        row = {"name": name, "smiles": smiles, "compound_type": ctype,
               "bbb": bbb, "tier": "gold_curated", "sample_weight": 0.85,
               "data_source": "curated_smiles", **dis_map}

        # Check for named negatives first
        if name in NEGATIVE_CONTROL_SCORES:
            row.update(NEGATIVE_CONTROL_SCORES[name])
        else:
            # Get class priors as defaults (ChEMBL scoring happens in training)
            prior_key = None
            for key in CLASS_PRIORS:
                if key in ctype.lower():
                    prior_key = key
                    break
            if prior_key:
                for dim in DIMENSION_COLS:
                    row[dim] = CLASS_PRIORS[prior_key][dim]
            else:
                for dim in DIMENSION_COLS:
                    row[dim] = CLASS_PRIORS["general"][dim]

        new_rows.append(row)

    log.info(f"  New compounds added: {len(new_rows)}  (skipped {skipped} duplicates)")

    df_new = pd.DataFrame(new_rows)
    df_all = pd.concat([df_existing, df_new], ignore_index=True)

    # ── Step 4: Deduplicate ──────────────────────────────────────────────
    log.info(f"\n[4/9] Deduplicating...")
    df_all = deduplicate(df_all)

    # ── Step 5: Fix NaN scores ───────────────────────────────────────────
    log.info(f"\n[5/9] Fixing NaN dimension scores...")
    df_all = fix_scores(df_all)

    # ── Step 6: Compute NPS ──────────────────────────────────────────────
    log.info(f"\n[6/9] Computing NPS for all compounds...")
    df_all = compute_nps(df_all)

    # ── Step 7: Retrieve SMILES for new additions still missing ──────────
    log.info(f"\n[7/9] Retrieving SMILES for new additions...")
    df_all = retrieve_smiles_for_dataframe(df_all)

    # ── Step 8: Save + validate ──────────────────────────────────────────
    log.info(f"\n[8/9] Saving and validating...")
    pct_smi = validate_dataset(df_all)

    out_path = ROOT / "data" / "brainsafe_master.csv"
    df_all.to_csv(out_path, index=False)
    log.info(f"\n  ✓ Master dataset saved: {out_path}")
    log.info(f"  ✓ Total compounds: {len(df_all)}")

    # Save SMILES for applicability domain
    valid_smi_list = df_all[
        df_all["smiles"].fillna("").apply(validate_smiles)
    ]["smiles"].tolist()
    with open(ROOT / "models_v5" / "training_smiles.json", "w") as f:
        json.dump(valid_smi_list, f)
    log.info(f"  ✓ Training SMILES: {len(valid_smi_list)}")

    # ── Step 9: Retrain ──────────────────────────────────────────────────
    log.info(f"\n[9/9] Training on {len(df_all)}-compound master dataset...")
    retrain(str(out_path))

    log.info("\n" + "=" * 62)
    log.info("MASTER PIPELINE COMPLETE")
    log.info("Paste the FINAL TRAINING RESULTS block above to Claude.")
    log.info("=" * 62)


if __name__ == "__main__":
    main()
