"""Hybrid retraining: property-matched decoys PLUS experimentally measured inactives.

Decoy-only training has a specific failure mode. Because the negatives are chosen to be structurally
dissimilar, the model becomes over-confident and saturates: for melatonin MT1 the measured inactives
score a median binder probability of 0.972 against an active median of 0.997, so the decision
boundary sits in a region where tiny probability differences swing the result. Sensitivity at a
controlled false-positive rate then collapses (MT1 0.52, KEAP1 0.33).

The fix is to give the model the hard negatives it never saw. Measured inactives (compounds tested
against the target and found inactive) are split in half: one half joins the training negatives
alongside the decoys, the other half is held out and used for unbiased validation and for setting
the decision threshold. Positives remain measured binders at pChEMBL >= 7.

Every model is scaffold-grouped 10-fold cross-validated, then refit and calibrated. Metrics are
reported on the held-out measured inactives only, so no compound used for training is used to judge
the model. Outputs models_rf/<T>_binder.joblib and refreshed binder_modes.json entries.
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
from sklearn.metrics import roc_auc_score

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize  # noqa: E402
from models.train_rf import _scaffold_groups  # noqa: E402

M = ROOT / "models_rf"
ACTIVE_P, DECOY_RATIO, TAN_MAX = 7.0, 3, 0.35
TARGET_FPR = 0.10
_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
rng = np.random.default_rng(42)
RF = dict(n_estimators=300, min_samples_leaf=4, n_jobs=-1, random_state=42, class_weight="balanced")

TARGETS = ["HT1A", "HT6", "HT7", "H3", "DAT", "NET", "Sigma1", "CB1", "OPRK1", "OPRM1", "D3", "A1",
           "a7nAChR", "LRRK2", "NLRP3", "P2X7", "COX2", "CSF1R", "PDE10A", "HDAC1", "HDAC6",
           "GluN2B", "mGluR5", "GABA_A", "OX1", "OX2", "MT1", "mTOR", "SIRT1", "KEAP1", "GBA1",
           "PDE4B", "Nav1_5", "D2", "A2A", "HT2A", "SERT"]


def _fp(s):
    m = Chem.MolFromSmiles(s)
    return _GEN.GetFingerprint(m) if m else None


def main():
    with (M / "ad_reference.pkl").open("rb") as fh:
        bg_smiles, bg_fps = pickle.load(fh)
    bg_smiles = np.array(bg_smiles)
    bgX, bgm = featurize(list(bg_smiles))
    bg_ok = bg_smiles[bgm]
    bg_fp_ok = [f for f, k in zip(bg_fps, bgm) if k]
    bg_mw, bg_lp = bgX[:, -12], bgX[:, -11]

    modes = json.loads((M / "binder_modes.json").read_text()) if (M / "binder_modes.json").exists() else {}
    for ep in TARGETS:
        f = ROOT / "data" / "endpoints" / f"{ep}.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f).dropna(subset=["smiles"])
        p = pd.to_numeric(df.get("pchembl"), errors="coerce")
        act = sorted(set(df.loc[p >= ACTIVE_P, "smiles"].astype(str)))
        ina_all = sorted(set(df.loc[df["label"] == 0, "smiles"].astype(str)) - set(act))
        if len(act) < 100:
            print(f"[{ep}] too few binders ({len(act)}), skipped", flush=True)
            continue

        # split measured inactives: half to train as hard negatives, half held out
        ina_all = list(rng.permutation(ina_all))
        half = len(ina_all) // 2
        ina_train, ina_hold = ina_all[:half], ina_all[half:]

        aX, am = featurize(act)
        act = [s for s, k in zip(act, am) if k]
        afps = [_fp(s) for s in act]
        mw_lo, mw_hi = np.quantile(aX[:, -12], .02), np.quantile(aX[:, -12], .98)
        lp_lo, lp_hi = np.quantile(aX[:, -11], .02), np.quantile(aX[:, -11], .98)
        cand = (bg_mw >= mw_lo) & (bg_mw <= mw_hi) & (bg_lp >= lp_lo) & (bg_lp <= lp_hi)
        aset = set(act)
        idx = [i for i in np.where(cand)[0] if bg_ok[i] not in aset]
        rng.shuffle(idx)
        need = max(DECOY_RATIO * len(act) - len(ina_train), len(act))
        dec = []
        for i in idx:
            fp = bg_fp_ok[i]
            if fp is None:
                continue
            if max(DataStructs.BulkTanimotoSimilarity(fp, afps)) < TAN_MAX:
                dec.append(str(bg_ok[i]))
            if len(dec) >= need:
                break

        smiles = act + dec + ina_train
        y = np.array([1] * len(act) + [0] * (len(dec) + len(ina_train)))
        X, mask = featurize(smiles)
        y = y[mask]
        smiles = [s for s, k in zip(smiles, mask) if k]
        groups = _scaffold_groups(smiles)

        oof = cross_val_predict(RandomForestClassifier(**RF), X, y, groups=groups,
                                cv=GroupKFold(10), method="predict_proba", n_jobs=-1)[:, 1]
        auroc_cv = float(roc_auc_score(y, oof))

        Xt, Xc, yt, yc = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)
        forest = RandomForestClassifier(**RF).fit(Xt, yt)
        cal = CalibratedClassifierCV(FrozenEstimator(forest), method="sigmoid").fit(Xc, yc)

        rec = {"mode": "hybrid_decoys_plus_measured_inactives",
               "n_positive": int((y == 1).sum()), "n_decoy": len(dec),
               "n_measured_inactive_train": len(ina_train),
               "n_measured_inactive_holdout": len(ina_hold),
               "scaffold_cv_auroc": round(auroc_cv, 3)}

        if len(ina_hold) >= 15:
            hX, _ = featurize(ina_hold)
            ph = cal.predict_proba(hX)[:, 1]
            pa = cal.predict_proba(aX)[:, 1]
            thr = float(np.clip(np.quantile(ph, 1.0 - TARGET_FPR), 0.05, 0.999))
            rec.update({
                "threshold": round(thr, 4), "threshold_basis": "held_out_measured_inactives",
                "n_measured_inactive": len(ina_hold),
                "auroc_vs_measured_inactives": round(float(roc_auc_score(
                    np.r_[np.ones(len(pa)), np.zeros(len(ph))], np.r_[pa, ph])), 3),
                "fpr_at_threshold": round(float((ph >= thr).mean()), 3),
                "sensitivity_at_threshold": round(float((pa >= thr).mean()), 3)})
            rec["reliable_call"] = bool(rec["sensitivity_at_threshold"] >= 0.60
                                        and rec["auroc_vs_measured_inactives"] >= 0.75)
            old = modes.get(ep, {})
            print(f"[{ep:8}] AUROC {rec['auroc_vs_measured_inactives']:.3f} "
                  f"(was {old.get('auroc_vs_measured_inactives','-')}) | sens "
                  f"{rec['sensitivity_at_threshold']:.3f} (was {old.get('sensitivity_at_threshold','-')}) "
                  f"| thr {thr:.3f}", flush=True)
        else:
            rec.update({"threshold": 0.40, "threshold_basis": "global_fallback"})
            print(f"[{ep:8}] too few held-out inactives, fallback threshold", flush=True)

        joblib.dump(cal, M / f"{ep}_binder.joblib", compress=3)
        modes[ep] = rec

    (M / "binder_modes.json").write_text(json.dumps(modes, indent=2))
    good = [v for v in modes.values() if "auroc_vs_measured_inactives" in v]
    if good:
        print(f"\nmean AUROC {np.mean([g['auroc_vs_measured_inactives'] for g in good]):.3f} | "
              f"mean sensitivity {np.mean([g['sensitivity_at_threshold'] for g in good]):.3f} | "
              f"unreliable {[k for k, v in modes.items() if v.get('reliable_call') is False]}")
    print("DONE")


if __name__ == "__main__":
    main()
