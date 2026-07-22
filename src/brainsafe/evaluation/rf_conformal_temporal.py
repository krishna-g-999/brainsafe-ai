"""Conformal prediction and temporal validation for the random-forest models.

Two analyses the manuscript body needs, computed directly on the deployed random forests:

  Conformal prediction (classification). Inductive Mondrian (class-conditional) conformal at a 10%
  significance level. Data are split train / calibration / test; the nonconformity of a compound is
  1 - predicted probability of its true class; a per-class threshold is taken from the calibration
  quantile; the test prediction set contains every class whose nonconformity is below its threshold.
  Empirical coverage should be near the target 0.90, with average set size reported.

  Temporal validation. For endpoints that carry a document year, the model is trained on compounds up
  to the 75th-percentile year and tested on strictly later compounds, measuring generalisation to
  chemistry that appeared after training (AUROC for classifiers, R-squared for regressors).

Outputs: results/tables/rf_conformal.csv, results/tables/rf_temporal.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import RDLogger
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score, r2_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from features.featurize import featurize  # noqa: E402
from models.train_rf import (CLASSIFICATION, REGRESSION, RF_COMMON, SEED, _load)  # noqa: E402

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[3]
EPS = 0.10  # 90% target coverage


def conformal(endpoint):
    df = _load(endpoint).dropna(subset=["smiles", "label"]).reset_index(drop=True)
    X, mask = featurize(df["smiles"].tolist())
    y = df.loc[mask, "label"].to_numpy().astype(int)
    rng = np.random.RandomState(SEED)
    idx = rng.permutation(len(y))
    n = len(y); a, b = int(0.6 * n), int(0.8 * n)
    tr, cal, te = idx[:a], idx[a:b], idx[b:]
    model = RandomForestClassifier(class_weight="balanced", **RF_COMMON).fit(X[tr], y[tr])
    proba = model.predict_proba(X)
    classes = list(model.classes_)
    # Mondrian: per-class nonconformity threshold from calibration compounds of that class
    thr = {}
    for c in classes:
        ci = classes.index(c)
        cal_c = cal[y[cal] == c]
        if len(cal_c) == 0:
            thr[c] = 1.0; continue
        alpha = 1.0 - proba[cal_c, ci]
        k = int(np.ceil((len(cal_c) + 1) * (1 - EPS)))
        thr[c] = np.sort(alpha)[min(k, len(cal_c)) - 1]
    covered, sizes = [], []
    for i in te:
        pset = [c for c in classes if (1.0 - proba[i, classes.index(c)]) <= thr[c]]
        covered.append(y[i] in pset)
        sizes.append(len(pset))
    return {"endpoint": endpoint, "n_test": len(te), "target_coverage": 1 - EPS,
            "empirical_coverage": round(float(np.mean(covered)), 3),
            "avg_set_size": round(float(np.mean(sizes)), 3)}


def temporal(endpoint, task):
    target = "label" if task == "classification" else ("y" if endpoint == "antioxidant_DPPH" else "pchembl")
    df = _load(endpoint)
    if "year" not in df.columns:
        return None
    df = df.dropna(subset=["smiles", target, "year"]).reset_index(drop=True)
    if len(df) < 200:
        return None
    cutoff = np.percentile(df["year"], 75)
    if (df["year"] > cutoff).sum() < 30:
        return None
    X, mask = featurize(df["smiles"].tolist())
    df = df.loc[mask].reset_index(drop=True)
    y = df[target].to_numpy()
    trm = (df["year"] <= cutoff).to_numpy(); tem = (df["year"] > cutoff).to_numpy()
    if task == "classification":
        y = y.astype(int)
        if len(set(y[tem])) < 2:
            return None
        m = RandomForestClassifier(class_weight="balanced", **RF_COMMON).fit(X[trm], y[trm])
        score = roc_auc_score(y[tem], m.predict_proba(X[tem])[:, 1]); metric = "auroc"
    else:
        y = y.astype(float)
        m = RandomForestRegressor(**RF_COMMON).fit(X[trm], y[trm])
        score = r2_score(y[tem], m.predict(X[tem])); metric = "r2"
    return {"endpoint": endpoint, "task": task, "cutoff_year": int(cutoff),
            "n_train": int(trm.sum()), "n_test": int(tem.sum()),
            "metric": metric, "score": round(float(score), 3)}


def main():
    conf = [conformal(ep) for ep in CLASSIFICATION]
    pd.DataFrame(conf).to_csv(ROOT / "results" / "tables" / "rf_conformal.csv", index=False)
    print("=== conformal (target coverage 0.90) ===")
    print(pd.DataFrame(conf).to_string(index=False))

    temp = []
    for ep in CLASSIFICATION:
        r = temporal(ep, "classification")
        if r:
            temp.append(r)
    for ep in list(REGRESSION):
        r = temporal(ep, "regression")
        if r:
            temp.append(r)
    pd.DataFrame(temp).to_csv(ROOT / "results" / "tables" / "rf_temporal.csv", index=False)
    print("\n=== temporal (train past, test future) ===")
    print(pd.DataFrame(temp).to_string(index=False))
    print("\nwrote results/tables/rf_conformal.csv, rf_temporal.csv")


if __name__ == "__main__":
    main()
