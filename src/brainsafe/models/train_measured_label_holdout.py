"""Retrain the remaining measured-label models with a held-out inactive split.

SIRT1 and Nav1.5 have too few binders at pChEMBL >= 7 for the hybrid binder scheme, so they are
trained on the measured activity label instead (active at pChEMBL >= 6, inactive below 5). In the
earlier pass these models were evaluated on the full set of measured inactives, including the ones
used to fit them, which inflates their apparent performance. Here half the measured inactives are
withheld from training and used for validation and threshold calibration, matching the protocol used
for the hybrid binder panel so all numbers in the paper are comparable and none are leaked.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.frozen import FrozenEstimator
from sklearn.model_selection import GroupKFold, cross_val_predict, train_test_split
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize  # noqa: E402
from models.train_rf import _scaffold_groups  # noqa: E402

M = ROOT / "models_rf"
TARGETS = ["SIRT1", "Nav1_5"]
TARGET_FPR = 0.10
RF = dict(n_estimators=300, min_samples_leaf=4, n_jobs=-1, random_state=42, class_weight="balanced")
rng = np.random.default_rng(42)


def main():
    modes = json.loads((M / "binder_modes.json").read_text())
    for ep in TARGETS:
        f = ROOT / "data" / "endpoints" / f"{ep}.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f).dropna(subset=["smiles", "label"])
        pos = df.loc[df.label == 1, "smiles"].astype(str).tolist()
        neg_all = list(rng.permutation(df.loc[df.label == 0, "smiles"].astype(str).tolist()))
        half = len(neg_all) // 2
        neg_train, neg_hold = neg_all[:half], neg_all[half:]

        smiles = pos + neg_train
        y = np.array([1] * len(pos) + [0] * len(neg_train))
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

        hX, _ = featurize(neg_hold)
        pX, _ = featurize(pos)
        ph = cal.predict_proba(hX)[:, 1]
        pa = cal.predict_proba(pX)[:, 1]
        thr = float(np.clip(np.quantile(ph, 1.0 - TARGET_FPR), 0.05, 0.999))
        rec = {"mode": "measured_labels_holdout",
               "n_positive": int(len(pos)), "n_measured_inactive_train": len(neg_train),
               "n_measured_inactive": len(neg_hold),
               "scaffold_cv_auroc": round(auroc_cv, 3),
               "threshold": round(thr, 4), "threshold_basis": "held_out_measured_inactives",
               "auroc_vs_measured_inactives": round(float(roc_auc_score(
                   np.r_[np.ones(len(pa)), np.zeros(len(ph))], np.r_[pa, ph])), 3),
               "fpr_at_threshold": round(float((ph >= thr).mean()), 3),
               "sensitivity_at_threshold": round(float((pa >= thr).mean()), 3)}
        rec["reliable_call"] = bool(rec["sensitivity_at_threshold"] >= 0.60
                                    and rec["auroc_vs_measured_inactives"] >= 0.75)
        old = modes.get(ep, {})
        print(f"[{ep:8}] AUROC {rec['auroc_vs_measured_inactives']:.3f} "
              f"(previously reported {old.get('auroc_vs_measured_inactives','-')}, leaked) | "
              f"sens {rec['sensitivity_at_threshold']:.3f} | thr {thr:.3f} | "
              f"held-out inactives {len(neg_hold)}", flush=True)
        joblib.dump(cal, M / f"{ep}_binder.joblib", compress=3)
        modes[ep] = rec

    (M / "binder_modes.json").write_text(json.dumps(modes, indent=2))
    good = [v for v in modes.values() if "auroc_vs_measured_inactives" in v]
    print(f"\npanel mean AUROC {np.mean([g['auroc_vs_measured_inactives'] for g in good]):.3f} | "
          f"mean sensitivity {np.mean([g['sensitivity_at_threshold'] for g in good]):.3f} | n={len(good)}")
    print("DONE")


if __name__ == "__main__":
    main()
