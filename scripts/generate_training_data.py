"""
scripts/generate_training_data.py
Converts existing compounds.json + smiles_cache.json into a properly
formatted CSV for ml_v5_training.py.

Run from D:\\BRAINSAFE_AI:
    python scripts/generate_training_data.py

Output: data/brainsafe_training_set.csv
"""

import json
import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path

# Add parent directory to path so we can import model_config
sys.path.insert(0, str(Path(__file__).parent.parent))

ROOT = Path(__file__).parent.parent
DATA_DIR = ROOT / "data"
DATA_DIR.mkdir(exist_ok=True)

DIMENSION_COLS = [
    "antioxidant", "anti_inflammatory", "mitochondrial_support",
    "aggregation_modulation", "cognitive_enhancement",
    "neurogenesis", "synaptic_plasticity",
]

DIS_MAP = {"Low": 0, "Med": 1, "High": 2}
BBB_MAP = {"Low": 0, "Low-Med": 1, "Medium": 2, "High": 3}


def load_smiles_cache(root: Path) -> dict:
    cache_path = root / "smiles_cache.json"
    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)
    return {}


def detect_score_scale(values: list) -> float:
    """Return 10.0 if scores appear to be on 1–10 scale, else 1.0."""
    vals = [v for v in values if v is not None and not np.isnan(v)]
    if not vals:
        return 1.0
    return 10.0 if max(vals) <= 10.5 else 1.0


def compounds_json_to_df(compounds: dict, smiles_cache: dict) -> pd.DataFrame:
    rows = []
    for name, entry in compounds.items():
        if name.startswith("_"):
            continue

        # SMILES: try entry first, then smiles_cache
        smiles = entry.get("smiles", "") or smiles_cache.get(name, "")
        smiles = str(smiles).strip() if smiles else ""
        if smiles.lower() in ("n/a", "na", "none", "null", ""):
            smiles = ""

        # Extract dimension scores
        raw_scores = {}
        for col in DIMENSION_COLS:
            val = entry.get(col)
            if val is None:
                val = entry.get(col.replace("_", ""))  # try without underscore
            try:
                raw_scores[col] = float(val) if val is not None else np.nan
            except (TypeError, ValueError):
                raw_scores[col] = np.nan

        # Auto-detect scale and normalise to 0–100
        scale = detect_score_scale(list(raw_scores.values()))
        scores_100 = {
            col: round(float(np.clip(v * scale, 0.0, 100.0)), 1)
            if not np.isnan(v) else np.nan
            for col, v in raw_scores.items()
        }

        # Disease relevance (High/Med/Low → numeric count proxy)
        def dis_count(key):
            val = entry.get(key, "Low")
            if isinstance(val, (int, float)):
                return int(val)
            return DIS_MAP.get(str(val), 0)

        row = {
            "name":            name,
            "smiles":          smiles,
            "bbb":             entry.get("bbb", "Low"),
            "compound_type":   entry.get("compound_type", "Unknown"),
            "data_source":     entry.get("data_source", "curated"),
            "ad_target_count": dis_count("alzheimers"),
            "pd_target_count": dis_count("parkinsons"),
            "als_target_count":dis_count("als"),
            "hd_target_count": dis_count("huntingtons"),
            **scores_100,
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    return df


def main():
    print("=" * 60)
    print("BrainSafe AI v6 — Training Data Generator")
    print("=" * 60)

    # Load SMILES cache
    smiles_cache = load_smiles_cache(ROOT)
    print(f"  SMILES cache: {len(smiles_cache)} entries")

    # Load curated compounds
    curated_path = ROOT / "compounds.json"
    if not curated_path.exists():
        print(f"ERROR: {curated_path} not found. Run from D:\\BRAINSAFE_AI")
        sys.exit(1)
    with open(curated_path) as f:
        curated = json.load(f)
    print(f"  Curated compounds: {len(curated)}")

    df_curated = compounds_json_to_df(curated, smiles_cache)
    df_curated["tier"] = "gold"

    # Load ML compounds (if present) — these become silver candidates
    ml_path = ROOT / "compounds_ml.json"
    df_silver = None
    if ml_path.exists():
        with open(ml_path) as f:
            ml_raw = json.load(f)
        ml_raw.pop("_ml_metadata", None)
        df_silver = compounds_json_to_df(ml_raw, smiles_cache)
        df_silver["tier"] = "silver"
        print(f"  ML (silver) compounds: {len(df_silver)}")

    # ── Validation: check score distributions ─────────────────────────────
    print("\n  Score validation (curated set):")
    for col in DIMENSION_COLS:
        col_data = df_curated[col].dropna()
        if len(col_data) > 0:
            print(f"    {col:<30} n={len(col_data):3d}  "
                  f"min={col_data.min():.1f}  mean={col_data.mean():.1f}  "
                  f"max={col_data.max():.1f}")
        else:
            print(f"    {col:<30} *** NO DATA ***")

    # ── Check for missing dimensions ──────────────────────────────────────
    missing_any = df_curated[DIMENSION_COLS].isnull().any(axis=1).sum()
    if missing_any > 0:
        print(f"\n  WARNING: {missing_any} compounds have at least one missing dimension score.")
        print("  These will be imputed with column medians during training.")
        # Impute with column medians for training
        for col in DIMENSION_COLS:
            median_val = df_curated[col].median()
            df_curated[col].fillna(median_val, inplace=True)

    # ── SMILES coverage ───────────────────────────────────────────────────
    n_smiles = (df_curated["smiles"] != "").sum()
    print(f"\n  SMILES coverage: {n_smiles}/{len(df_curated)} "
          f"({100*n_smiles/len(df_curated):.1f}%)")

    # ── Save gold-standard training set ──────────────────────────────────
    out_path = DATA_DIR / "brainsafe_training_set.csv"
    df_curated.to_csv(out_path, index=False)
    print(f"\n  ✓ Gold training set saved: {out_path}")
    print(f"    Shape: {df_curated.shape}")

    # ── Save silver candidates ────────────────────────────────────────────
    if df_silver is not None:
        silver_path = DATA_DIR / "silver_candidates.csv"
        df_silver.to_csv(silver_path, index=False)
        print(f"  ✓ Silver candidates saved: {silver_path}")
        print(f"    Shape: {df_silver.shape}")

    # ── Save training SMILES for applicability domain ─────────────────────
    smiles_list = df_curated[df_curated["smiles"] != ""]["smiles"].tolist()
    smiles_out = ROOT / "models_v5" / "training_smiles.json"
    smiles_out.parent.mkdir(exist_ok=True)
    with open(smiles_out, "w") as f:
        json.dump(smiles_list, f)
    print(f"  ✓ Training SMILES saved: {smiles_out} ({len(smiles_list)} entries)")

    print("\n  DONE. Ready to run:")
    print("  python ml_v5_training.py --data data/brainsafe_training_set.csv --out models_v5/")


if __name__ == "__main__":
    main()
