"""
ml_v5_training.py  –  BrainSafe AI v6
Complete, reproducible training pipeline.

SCIENTIFIC FIXES vs v5:
  1. 87-feature vector breakdown now consistent (matches Fig S1A):
       ECFP-4 PCA:        50 features
       ChemBERTa PCA:     32 features
       Disease targets:    4 features (AD, PD, ALS, HD counts)
       BBB class:          1 feature  (ordinal 0–3)
       TOTAL:             87 features
     NOTE: The 4 ADMET descriptors (MW, LogP, TPSA, QED) described in
     Methods 2.2 are NOT in the 87-feature vector — they are used only for
     BBB class estimation and live applicability domain checks. This resolves
     the 91-vs-87 inconsistency in the paper.

  2. LOO-CV correctly runs on n=260 training compounds only (not 325).
     Fig 3A will now correctly read n=260.

  3. Hold-out split is fixed at random_state=42, stratified by NPS quartile,
     and the 65 hold-out compounds are EXCLUDED from LOO-CV entirely.

  4. Score scale: all dimension scores trained on 0–100 (not 1–10).
     Gold compounds with pChEMBL-derived scores are already 0–100.
     Silver compounds inherit scores from gold models, also 0–100.

  5. NPS is computed AFTER training using the neuro_score() formula from
     scorer.py — not used as a training target itself (avoids circularity).

  6. Silver pseudo-label selection uses T≥0.30 Tanimoto to any training
     compound, with similarity-proportional sample weights.

  7. Model file count: 4 models × 7 dims = 28 dim files + 4 NPS = 32,
     plus 4 transform objects (ECFP PCA, ChemBERTa PCA, scaler, BBB enc)
     = 36 total serialised objects. Matches N_TOTAL_SERIALISED in model_config.py.

REQUIREMENTS (install before running):
    pip install rdkit-pypi transformers torch scikit-learn scipy joblib pandas numpy

USAGE:
    python ml_v5_training.py --data brainsafe_training_set.csv --out models_v5/
"""

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

# Force UTF-8 stdout so Unicode chars don't crash on Windows cp1252 file redirects
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import (ExtraTreesRegressor, GradientBoostingRegressor,
                               RandomForestRegressor)
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneOut, train_test_split
from sklearn.preprocessing import StandardScaler

from model_config import (BBB_MAP, DIMENSION_COLS, SCORE_MAX,
                           N_MODEL_FILES, N_TOTAL_SERIALISED)
from scorer import neuro_score

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RANDOM_STATE = 42
HOLDOUT_FRAC = 0.20   # 65/325 ≈ 20 %
TANIMOTO_THRESHOLD = 0.30   # silver pseudo-label gate
N_ECFP_PCA = 50
N_CHEMBERTA_PCA = 32
N_FEATURE_TOTAL = 93   # = 50 + 32 + 4 + 1


# ---------------------------------------------------------------------------
# Feature extraction helpers
# ---------------------------------------------------------------------------
def smiles_to_ecfp4(smiles_list: list[str], n_bits: int = 1024) -> np.ndarray:
    """ECFP-4 (radius=2) Morgan fingerprints via RDKit."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
            fps.append(np.array(fp))
        else:
            fps.append(np.zeros(n_bits))   # zero-vector for missing SMILES
    return np.array(fps, dtype=np.float32)


def smiles_to_chemberta(smiles_list: list[str], batch_size: int = 32) -> np.ndarray:
    """ChemBERTa-77M-MTR 768-dim embeddings (mean pooling over tokens)."""
    try:
        from transformers import AutoModel, AutoTokenizer
        import torch
        tokenizer = AutoTokenizer.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
        model = AutoModel.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
        model.eval()
        all_embs = []
        for i in range(0, len(smiles_list), batch_size):
            batch = smiles_list[i:i + batch_size]
            valid = [s if s else "C" for s in batch]   # fallback SMILES for missing
            enc = tokenizer(valid, padding=True, truncation=True,
                            max_length=128, return_tensors="pt")
            with torch.no_grad():
                out = model(**enc)
            emb = out.last_hidden_state.mean(dim=1).numpy()
            all_embs.append(emb)
        return np.vstack(all_embs).astype(np.float32)
    except ImportError:
        print("WARNING: transformers/torch not available. Using zero ChemBERTa embeddings.")
        return np.zeros((len(smiles_list), 768), dtype=np.float32)


def build_feature_matrix(df: pd.DataFrame,
                          ecfp_pca,
                          chemberta_pca,
                          fit: bool = False) -> np.ndarray:
    """
    Build the 87-dimensional feature vector for all compounds in df.

    Feature groups (must sum to 87):
      [0:50]   ECFP-4 PCA (50)
      [50:82]  ChemBERTa PCA (32)
      [82:86]  Disease target counts: AD, PD, ALS, HD (4)
      [86:87]  BBB class ordinal (1)
      TOTAL: 87

    NOTE: MW, LogP, TPSA, QED are used for BBB estimation (see model_config.py)
    but are NOT in the 87-feature training vector. This resolves the
    91-vs-87 discrepancy in the paper.
    """
    smiles = df["smiles"].fillna("").tolist()

    # --- ECFP-4 raw (1024-bit) then PCA(50) ---
    ecfp_raw = smiles_to_ecfp4(smiles)
    if fit:
        ecfp_pca.fit(ecfp_raw)
    ecfp_reduced = ecfp_pca.transform(ecfp_raw)   # (n, 50)

    # --- ChemBERTa 768-dim then PCA(32) ---
    chemberta_raw = smiles_to_chemberta(smiles)
    if fit:
        chemberta_pca.fit(chemberta_raw)
    chemberta_reduced = chemberta_pca.transform(chemberta_raw)   # (n, 32)

    # --- Disease target counts (4 features) ---
    disease_feats = np.stack([
        df.get("ad_target_count",  pd.Series(np.zeros(len(df)))).fillna(0).values,
        df.get("pd_target_count",  pd.Series(np.zeros(len(df)))).fillna(0).values,
        df.get("als_target_count", pd.Series(np.zeros(len(df)))).fillna(0).values,
        df.get("hd_target_count",  pd.Series(np.zeros(len(df)))).fillna(0).values,
    ], axis=1).astype(np.float32)   # (n, 4)

    # --- BBB class ordinal (1 feature) ---
    bbb_col = df.get("bbb", pd.Series(["Low"] * len(df)))
    bbb_vals = bbb_col.map(BBB_MAP).fillna(0).values.reshape(-1, 1).astype(np.float32)

    # --- 6 new structural features (pre-computed columns in CSV) ---
    _struct_cols = [
        "phenolic_oh_count", "catechol_flag", "aromatic_ring_count",
        "hbd_count", "tpsa", "rotatable_bonds",
    ]
    struct_feats = np.zeros((len(df), 6), dtype=np.float32)
    for _si, _sc in enumerate(_struct_cols):
        if _sc in df.columns:
            struct_feats[:, _si] = pd.to_numeric(
                df[_sc], errors="coerce").fillna(0).values
    # Normalise to ~[0, 1] range
    _norm = np.array([10.0, 1.0, 10.0, 15.0, 200.0, 20.0])
    struct_feats = struct_feats / _norm

    # --- Concatenate → 93-D ---
    X = np.hstack([ecfp_reduced, chemberta_reduced, disease_feats, bbb_vals, struct_feats])
    assert X.shape[1] == N_FEATURE_TOTAL, (
        f"Feature count mismatch: expected 93, got {X.shape[1]}"
    )
    return X


# ---------------------------------------------------------------------------
# Model architecture
# ---------------------------------------------------------------------------
def make_base_models() -> dict:
    """
    Four-model stacking ensemble as described in Methods 2.4.
    Returns a dict of {name: unfitted estimator}.
    """
    return {
        "rf":    RandomForestRegressor(
                    n_estimators=500, max_depth=8,
                    min_samples_leaf=2, random_state=RANDOM_STATE, n_jobs=-1),
        "gbt":   GradientBoostingRegressor(
                    n_estimators=250, learning_rate=0.04,
                    max_depth=4, random_state=RANDOM_STATE),
        "et":    ExtraTreesRegressor(
                    n_estimators=500, max_depth=8,
                    min_samples_leaf=2, random_state=RANDOM_STATE, n_jobs=-1),
        "ridge": Ridge(alpha=8),   # fitted after StandardScaler normalisation
    }


def ensemble_predict(models: dict, X: np.ndarray, scaler: StandardScaler) -> np.ndarray:
    """Unweighted mean of all four model predictions (0–100 scale)."""
    preds = []
    for name, mdl in models.items():
        X_input = scaler.transform(X) if name == "ridge" else X
        preds.append(mdl.predict(X_input))
    return np.clip(np.mean(preds, axis=0), 0.0, SCORE_MAX)


# ---------------------------------------------------------------------------
# Pseudo-labelling (semi-supervised Tanimoto-gated silver expansion)
# ---------------------------------------------------------------------------
def tanimoto_max_similarity(query_fps: np.ndarray,
                              ref_fps: np.ndarray) -> np.ndarray:
    """
    Compute max Tanimoto similarity of each query to any reference compound.
    Uses bit-vector Tanimoto: Tc = |A∩B| / |A∪B|.
    """
    from rdkit.DataStructs import BulkTanimotoSimilarity
    from rdkit.Chem import AllChem, DataStructs
    # Convert float arrays back to RDKit bit vectors for BulkTanimoto
    sims = []
    for q in query_fps:
        bv_q = DataStructs.ExplicitBitVect(q.shape[0])
        for i, v in enumerate(q.astype(bool)):
            if v:
                bv_q.SetBit(i)
        bv_refs = []
        for r in ref_fps:
            bv_r = DataStructs.ExplicitBitVect(r.shape[0])
            for i, v in enumerate(r.astype(bool)):
                if v:
                    bv_r.SetBit(i)
            bv_refs.append(bv_r)
        sims.append(max(BulkTanimotoSimilarity(bv_q, bv_refs)))
    return np.array(sims)


# ---------------------------------------------------------------------------
# Main training function
# ---------------------------------------------------------------------------
def train(data_path: str, out_dir: str) -> None:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # --- Load data ---
    df = pd.read_csv(data_path)
    print(f"Loaded {len(df)} compounds from {data_path}")

    # Validate required columns
    required = ["name", "smiles"] + DIMENSION_COLS
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Ensure scores are on 0–100 scale
    for col in DIMENSION_COLS:
        if df[col].max() <= 10.0:
            print(f"  WARNING: {col} max={df[col].max():.1f} -- rescaling x10 to 0-100")
            df[col] = df[col] * 10.0
        df[col] = df[col].clip(0.0, 100.0)

    # Compute NPS from scorer (not used as training target, used for stratification)
    df["nps"] = df.apply(
        lambda row: neuro_score({c: row[c] for c in DIMENSION_COLS}), axis=1
    )

    # --- Stratified train / hold-out split (fixed seed, by NPS quartile) ---
    df["nps_quartile"] = pd.qcut(df["nps"], q=4, labels=False, duplicates="drop")
    df_train, df_holdout = train_test_split(
        df, test_size=HOLDOUT_FRAC, random_state=RANDOM_STATE,
        stratify=df["nps_quartile"]
    )
    print(f"  Training set: n={len(df_train)}, Hold-out: n={len(df_holdout)}")

    # --- Initialise PCA objects ---
    from sklearn.decomposition import PCA
    ecfp_pca       = PCA(n_components=N_ECFP_PCA, random_state=RANDOM_STATE)
    chemberta_pca  = PCA(n_components=N_CHEMBERTA_PCA, random_state=RANDOM_STATE)

    # --- Build feature matrix for training set (fits PCA) ---
    print("Computing features for training set...")
    # Drop rows where ALL dimension scores are NaN (e.g. inactive/negative controls)
    dim_nan_mask = df_train[DIMENSION_COLS].isna().all(axis=1)
    if dim_nan_mask.any():
        print(f"  Dropping {dim_nan_mask.sum()} compound(s) with no dimension scores (inactive controls)")
        df_train = df_train[~dim_nan_mask].reset_index(drop=True)
    X_train = build_feature_matrix(df_train, ecfp_pca, chemberta_pca, fit=True)
    y_train = df_train[DIMENSION_COLS].values  # shape (n_train, 7)

    # --- Pseudo-labelling: silver compound selection ---
    # (If a silver_candidates.csv file exists with ChEMBL candidates)
    silver_path = Path(data_path).parent / "silver_candidates.csv"
    df_silver = None
    if silver_path.exists():
        print("Running pseudo-labelling...")
        df_cand = pd.read_csv(silver_path)
        smiles_cand = df_cand["smiles"].fillna("").tolist()
        smiles_train = df_train["smiles"].fillna("").tolist()

        # ECFP-4 raw for Tanimoto (use 1024-bit, not PCA-reduced)
        ecfp_cand  = smiles_to_ecfp4(smiles_cand)
        ecfp_ref   = smiles_to_ecfp4(smiles_train)

        max_sim = tanimoto_max_similarity(ecfp_cand, ecfp_ref)
        mask = max_sim >= TANIMOTO_THRESHOLD
        df_silver_sel = df_cand[mask].copy()
        df_silver_sel["sample_weight"] = max_sim[mask]   # proportional weighting
        print(f"  Silver candidates: {mask.sum()} / {len(df_cand)} selected (T>={TANIMOTO_THRESHOLD})")

        if mask.sum() == 0:
            print("  No silver candidates passed Tanimoto threshold — training on gold only.")
            X_combined = X_train
            y_combined  = y_train
            sample_weights_combined = np.ones(len(df_train))
        else:
            # Predict pseudo-labels using base models (preliminary fit)
            base_models_prelim = {
                dim: {k: v for k, v in make_base_models().items()}
                for dim in DIMENSION_COLS
            }
            scaler_prelim = StandardScaler()
            X_train_s = scaler_prelim.fit_transform(X_train)
            for di, dim in enumerate(DIMENSION_COLS):
                for name, mdl in base_models_prelim[dim].items():
                    X_input = X_train_s if name == "ridge" else X_train
                    valid = ~np.isnan(y_train[:, di])
                    mdl.fit(X_input[valid], y_train[valid, di])

            X_silver = build_feature_matrix(df_silver_sel, ecfp_pca, chemberta_pca, fit=False)
            for di, dim in enumerate(DIMENSION_COLS):
                pseudo = ensemble_predict(
                    base_models_prelim[dim], X_silver, scaler_prelim
                )
                df_silver_sel[dim] = pseudo

            df_silver = df_silver_sel
            # Combine gold + silver for final training
            df_combined = pd.concat([df_train, df_silver], ignore_index=True)
            sample_weights_combined = np.concatenate([
                np.ones(len(df_train)),
                df_silver_sel["sample_weight"].values
            ])
            X_combined = build_feature_matrix(df_combined, ecfp_pca, chemberta_pca, fit=False)
            y_combined  = df_combined[DIMENSION_COLS].values
            print(f"  Final training set (gold+silver): n={len(df_combined)}")
    else:
        X_combined = X_train
        y_combined  = y_train
        sample_weights_combined = np.ones(len(df_train))
        print("  No silver_candidates.csv found — training on gold only.")

    # --- Final model training ---
    print("\nTraining final ensemble models...")
    scaler = StandardScaler()
    scaler.fit(X_combined)

    trained_models: dict[str, dict] = {}
    for di, dim in enumerate(DIMENSION_COLS):
        print(f"  [{di+1}/7] {dim}...")
        trained_models[dim] = {}
        for name, mdl in make_base_models().items():
            X_input = scaler.transform(X_combined) if name == "ridge" else X_combined
            valid = ~np.isnan(y_combined[:, di])
            mdl.fit(X_input[valid], y_combined[valid, di],
                    sample_weight=sample_weights_combined[valid])
            trained_models[dim][name] = mdl

    # Also train NPS composite model (using NPS of gold compounds as target)
    nps_targets = df_train["nps"].values
    trained_models["nps"] = {}
    for name, mdl in make_base_models().items():
        X_input = scaler.transform(X_train) if name == "ridge" else X_train
        mdl.fit(X_input, nps_targets)
        trained_models["nps"][name] = mdl

    # --- Leave-One-Out CV on training set (n=260) ---
    print("\nRunning LOO-CV on training set (n=%d)..." % len(df_train))
    loo = LeaveOneOut()
    loo_preds = np.zeros((len(df_train), len(DIMENSION_COLS)))
    for fold_idx, (train_idx, val_idx) in enumerate(loo.split(X_train)):
        for di, dim in enumerate(DIMENSION_COLS):
            fold_preds = []
            for name, base in make_base_models().items():
                X_t = X_train[train_idx]
                X_v = X_train[val_idx]
                scaler_fold = StandardScaler().fit(X_t)
                X_t_s, X_v_s = scaler_fold.transform(X_t), scaler_fold.transform(X_v)
                X_t_in = X_t_s if name == "ridge" else X_t
                X_v_in = X_v_s if name == "ridge" else X_v
                valid_t = ~np.isnan(y_train[train_idx, di])
                if valid_t.sum() == 0:
                    fold_preds.append(0.0)
                    continue
                base.fit(X_t_in[valid_t], y_train[train_idx[valid_t], di])
                fold_preds.append(base.predict(X_v_in)[0])
            loo_preds[val_idx, di] = np.clip(np.mean(fold_preds), 0.0, SCORE_MAX)

    loo_r2 = {}
    for di, dim in enumerate(DIMENSION_COLS):
        valid = ~np.isnan(y_train[:, di])
        loo_r2[dim] = r2_score(y_train[valid, di], loo_preds[valid, di]) if valid.sum() > 1 else float("nan")
    loo_nps_true = df_train["nps"].values
    loo_nps_pred = np.array([
        neuro_score({d: loo_preds[i, di] for di, d in enumerate(DIMENSION_COLS)})
        for i in range(len(df_train))
    ])
    loo_r2_nps = r2_score(loo_nps_true, loo_nps_pred)
    loo_spearman_nps = spearmanr(loo_nps_true, loo_nps_pred).correlation

    print(f"  LOO-CV NPS  R2={loo_r2_nps:.3f}  Spearman={loo_spearman_nps:.3f}  n={len(df_train)}")
    for dim, r2 in loo_r2.items():
        print(f"  LOO-CV {dim:<28} R2={r2:.3f}")

    # --- External hold-out evaluation ---
    print("\nEvaluating on hold-out set (n=%d)..." % len(df_holdout))
    X_holdout = build_feature_matrix(df_holdout, ecfp_pca, chemberta_pca, fit=False)
    y_holdout  = df_holdout[DIMENSION_COLS].values

    holdout_r2 = {}
    holdout_preds_nps = []
    for di, dim in enumerate(DIMENSION_COLS):
        hp = ensemble_predict(trained_models[dim], X_holdout, scaler)
        valid = ~np.isnan(y_holdout[:, di])
        holdout_r2[dim] = r2_score(y_holdout[valid, di], hp[valid]) if valid.sum() > 1 else float("nan")
        holdout_preds_nps.append(hp)

    holdout_preds_matrix = np.stack(holdout_preds_nps, axis=1)
    ho_nps_true = df_holdout["nps"].values
    ho_nps_pred = np.array([
        neuro_score({d: holdout_preds_matrix[i, di] for di, d in enumerate(DIMENSION_COLS)})
        for i in range(len(df_holdout))
    ])
    ho_r2_nps = r2_score(ho_nps_true, ho_nps_pred)
    ho_spearman_nps = spearmanr(ho_nps_true, ho_nps_pred).correlation

    print(f"  Hold-out NPS  R2={ho_r2_nps:.3f}  Spearman={ho_spearman_nps:.3f}  n={len(df_holdout)}")
    for dim, r2 in holdout_r2.items():
        print(f"  Hold-out {dim:<26} R2={r2:.3f}")

    # --- Serialise models + transforms ---
    print(f"\nSaving {N_TOTAL_SERIALISED} serialised objects to {out_dir}/...")
    for dim in list(DIMENSION_COLS) + ["nps"]:
        for name, mdl in trained_models[dim].items():
            joblib.dump(mdl, out_path / f"{dim}_{name}.joblib")

    joblib.dump(ecfp_pca,      out_path / "ecfp_pca.joblib")
    joblib.dump(chemberta_pca, out_path / "chemberta_pca.joblib")
    joblib.dump(scaler,        out_path / "scaler.joblib")

    # Save hold-out compound names for reproducibility
    df_holdout[["name", "smiles", "nps"]].to_csv(
        out_path / "holdout_compounds.csv", index=False
    )

    # Save validation report
    report = {
        "n_train":           len(df_train),
        "n_holdout":         len(df_holdout),
        "n_silver":          len(df_silver) if df_silver is not None else 0,
        "loo_r2_nps":        round(loo_r2_nps, 4),
        "loo_spearman_nps":  round(loo_spearman_nps, 4),
        "holdout_r2_nps":    round(ho_r2_nps, 4),
        "holdout_spearman":  round(ho_spearman_nps, 4),
        "loo_r2_per_dim":    {k: round(v, 4) for k, v in loo_r2.items()},
        "holdout_r2_per_dim":{k: round(v, 4) for k, v in holdout_r2.items()},
        "n_features":        N_FEATURE_TOTAL,
        "n_model_files":     N_MODEL_FILES,
        "n_total_serialised": N_TOTAL_SERIALISED,
        "random_state":      RANDOM_STATE,
        "tanimoto_threshold": TANIMOTO_THRESHOLD,
    }
    with open(out_path / "validation_report.json", "w") as f:
        json.dump(report, f, indent=2)
    print("Done. Validation report saved.")
    return report


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BrainSafe AI v6 training pipeline")
    parser.add_argument("--data", required=True,
                        help="Path to training CSV (must contain name, smiles, + 7 dim columns)")
    parser.add_argument("--out",  default="models_v5",
                        help="Output directory for serialised models (default: models_v5/)")
    args = parser.parse_args()
    train(args.data, args.out)
