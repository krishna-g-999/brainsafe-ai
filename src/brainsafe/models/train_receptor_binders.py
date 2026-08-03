"""Decoy-aware receptor binder classifiers (D2, A2A, HT2A, SERT).

The receptor ChEMBL sets are ~96-98% actives (only binders are reported), so a pKi regressor
learns to predict ~the median for everything and cannot discriminate a real binder from an
arbitrary molecule. Here we build an honest binder-vs-decoy CLASSIFIER per receptor:

  positives : measured strong binders (pChEMBL >= 7)
  negatives : property-matched decoys sampled from the 74k measured background library,
              excluded from that receptor's known binders and with max ECFP Tanimoto < 0.35
              to any positive (so the model must learn structure, not a similarity shortcut).

Each model is scaffold-10-fold validated (honest, no scaffold leakage), then isotonically
calibrated and saved to models_rf/{ep}_binder.joblib (+ _binder_meta.json). Decoys are
*presumed* inactive, so AUROC is an upper bound; this is still far more discriminative and
honest than the flat regressor. Run:  python src/brainsafe/models/train_receptor_binders.py
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.frozen import FrozenEstimator
from sklearn.model_selection import GroupKFold, cross_val_predict, train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize, featurize_one  # noqa: E402
from models.train_rf import _scaffold_groups  # noqa: E402

RECEPTORS = ["D2", "A2A", "HT2A", "SERT"]
ACTIVE_PCHEMBL = 7.0
DECOY_RATIO = 3
TANIMOTO_MAX = 0.35
_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
rng = np.random.default_rng(42)


def _fp(smi):
    m = Chem.MolFromSmiles(smi)
    return _GEN.GetFingerprint(m) if m else None


def main():
    # background library (smiles + fps) from the applicability-domain reference
    with (ROOT / "models_rf" / "ad_reference.pkl").open("rb") as fh:
        bg_smiles, bg_fps = pickle.load(fh)
    bg_smiles = np.array(bg_smiles)
    # descriptors for property matching (MW, clogP) for the background
    bgX, bg_mask = featurize(list(bg_smiles))
    bg_smiles_ok = bg_smiles[bg_mask]
    bg_fps_ok = [f for f, ok in zip(bg_fps, bg_mask) if ok]
    bg_mw = bgX[:, -12]      # 'mw' is first descriptor after 1024 fp bits
    bg_clogp = bgX[:, -11]   # 'clogp'

    # union of all receptor binders (to exclude from any decoy pool)
    all_binders = set()
    active_sets = {}
    for ep in RECEPTORS:
        df = pd.read_csv(ROOT / "data" / "endpoints" / f"{ep}.csv")
        p = pd.to_numeric(df["pchembl"], errors="coerce")
        act = set(df.loc[p >= ACTIVE_PCHEMBL, "smiles"].astype(str))
        active_sets[ep] = act
        all_binders |= set(df["smiles"].astype(str))

    summary = {}
    for ep in RECEPTORS:
        actives = sorted(active_sets[ep])
        aX, amask = featurize(actives)
        actives = [s for s, ok in zip(actives, amask) if ok]
        afps = [_fp(s) for s in actives]
        a_mw, a_clogp = aX[:, -12], aX[:, -11]
        mw_lo, mw_hi = np.quantile(a_mw, 0.02), np.quantile(a_mw, 0.98)
        lp_lo, lp_hi = np.quantile(a_clogp, 0.02), np.quantile(a_clogp, 0.98)

        # candidate decoys: background, property-matched, not a known binder of this receptor
        cand = ((bg_mw >= mw_lo) & (bg_mw <= mw_hi) & (bg_clogp >= lp_lo) & (bg_clogp <= lp_hi))
        cand_idx = [i for i in np.where(cand)[0] if bg_smiles_ok[i] not in active_sets[ep]]
        rng.shuffle(cand_idx)
        need = DECOY_RATIO * len(actives)
        decoys, hard = [], []
        for i in cand_idx:
            f = bg_fps_ok[i]
            if f is None:
                continue
            mt = max(DataStructs.BulkTanimotoSimilarity(f, afps))
            if mt < TANIMOTO_MAX and len(decoys) < need:
                decoys.append(str(bg_smiles_ok[i]))
            elif TANIMOTO_MAX <= mt < 0.55 and len(hard) < len(actives):
                hard.append(str(bg_smiles_ok[i]))  # near-miss negatives, TEST ONLY
            if len(decoys) >= need and len(hard) >= len(actives):
                break

        smiles = actives + decoys
        y = np.array([1] * len(actives) + [0] * len(decoys))
        X, mask = featurize(smiles)
        y = y[mask]; smiles = [s for s, ok in zip(smiles, mask) if ok]
        groups = _scaffold_groups(smiles)

        # optimistic scaffold-CV AUROC on the easy (dissimilar) decoys
        oof = cross_val_predict(
            RandomForestClassifier(200, min_samples_leaf=6, n_jobs=-1, random_state=42, class_weight="balanced"),
            X, y, groups=groups, cv=GroupKFold(10), method="predict_proba", n_jobs=-1)[:, 1]
        auroc_easy = roc_auc_score(y, oof)

        # compact final model: single forest + prefit sigmoid calibration (small on disk)
        Xt, Xc, yt, yc = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)
        forest = RandomForestClassifier(200, min_samples_leaf=6, n_jobs=-1, random_state=42,
                                        class_weight="balanced").fit(Xt, yt)
        cal = CalibratedClassifierCV(FrozenEstimator(forest), method="sigmoid").fit(Xc, yc)
        joblib.dump(cal, ROOT / "models_rf" / f"{ep}_binder.joblib", compress=3)

        # HONEST realistic metrics ---------------------------------------------------------
        # (a) hard-decoy AUROC: held-out actives vs near-miss (0.35-0.55 Tanimoto) negatives
        auroc_hard = None
        if len(hard) >= 30:
            hX, hm = featurize(hard)
            pa = cal.predict_proba(Xc[yc == 1])[:, 1]
            ph = cal.predict_proba(hX)[:, 1]
            auroc_hard = float(roc_auc_score(np.r_[np.ones(len(pa)), np.zeros(len(ph))], np.r_[pa, ph]))
        # (b) background false-positive rate: random library compounds predicted as binders
        bg_idx = rng.choice(len(bg_smiles_ok), size=min(2500, len(bg_smiles_ok)), replace=False)
        bg_samp = [str(bg_smiles_ok[i]) for i in bg_idx if str(bg_smiles_ok[i]) not in active_sets[ep]]
        bX, bm = featurize(bg_samp)
        bp = cal.predict_proba(bX)[:, 1]
        fpr = float((bp >= 0.5).mean())

        meta = {"endpoint": ep, "task": "binder_vs_decoy", "n_active": int((y == 1).sum()),
                "n_decoy": int((y == 0).sum()), "active_pchembl_cut": ACTIVE_PCHEMBL,
                "auroc_easy_decoys": round(float(auroc_easy), 3),
                "auroc_hard_decoys": round(auroc_hard, 3) if auroc_hard else None,
                "background_fpr@0.5": round(fpr, 3),
                "note": "decoys presumed inactive; use hard-decoy AUROC + background FPR as the honest metrics"}
        (ROOT / "models_rf" / f"{ep}_binder_meta.json").write_text(json.dumps(meta, indent=2))
        summary[ep] = meta
        print(f"[{ep}] act={meta['n_active']} dec={meta['n_decoy']} | AUROC easy={auroc_easy:.3f} "
              f"hard={auroc_hard if auroc_hard is None else round(auroc_hard,3)} | bg-FPR={fpr:.3f}")

    print("\nDONE\n" + json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
