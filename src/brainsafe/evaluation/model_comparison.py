"""Does a different estimator beat the random forest? RF vs XGBoost vs histogram gradient boosting.

Same features, same folds, three estimators. Random and scaffold-grouped 5-fold (the scaffold column
is the one that matters). This answers, with numbers rather than assertion, whether switching to
gradient boosting would materially improve the models.

Output: results/tables/model_comparison.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import (HistGradientBoostingClassifier, HistGradientBoostingRegressor,
                              RandomForestClassifier, RandomForestRegressor)
from sklearn.model_selection import GroupKFold, KFold, StratifiedKFold
from sklearn.metrics import roc_auc_score, r2_score
from xgboost import XGBClassifier, XGBRegressor

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from features.featurize import featurize  # noqa: E402
from models.train_rf import CLASSIFICATION, REGRESSION, SEED, _load, _scaffold_groups  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
N_SPLITS = 5


def classifiers(pos_weight):
    return {
        "RandomForest": RandomForestClassifier(n_estimators=300, min_samples_leaf=2,
                                               class_weight="balanced", n_jobs=-1, random_state=SEED),
        "XGBoost": XGBClassifier(n_estimators=400, max_depth=6, learning_rate=0.05, subsample=0.8,
                                 colsample_bytree=0.8, tree_method="hist", eval_metric="logloss",
                                 scale_pos_weight=pos_weight, n_jobs=-1, random_state=SEED),
        "HistGradientBoosting": HistGradientBoostingClassifier(max_iter=400, learning_rate=0.06,
                                                               random_state=SEED),
    }


def regressors():
    return {
        "RandomForest": RandomForestRegressor(n_estimators=300, min_samples_leaf=2, n_jobs=-1,
                                              random_state=SEED),
        "XGBoost": XGBRegressor(n_estimators=400, max_depth=6, learning_rate=0.05, subsample=0.8,
                                colsample_bytree=0.8, tree_method="hist", n_jobs=-1, random_state=SEED),
        "HistGradientBoosting": HistGradientBoostingRegressor(max_iter=400, learning_rate=0.06,
                                                              random_state=SEED),
    }


def cv_score(task, model, X, y, groups, split):
    splitter = (GroupKFold(N_SPLITS) if split == "scaffold" else
                (StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED) if task == "classification"
                 else KFold(N_SPLITS, shuffle=True, random_state=SEED)))
    it = splitter.split(X, y, groups) if split == "scaffold" else splitter.split(X, y)
    scores = []
    for tr, te in it:
        m = model.__class__(**model.get_params())
        m.fit(X[tr], y[tr])
        if task == "classification":
            p = m.predict_proba(X[te])[:, 1]
            scores.append(roc_auc_score(y[te], p) if len(set(y[te])) > 1 else np.nan)
        else:
            scores.append(r2_score(y[te], m.predict(X[te])))
    return float(np.nanmean(scores)), float(np.nanstd(scores, ddof=1))


def main():
    rows = []
    for ep in CLASSIFICATION + list(REGRESSION):
        task = "classification" if ep in CLASSIFICATION else "regression"
        target = "label" if task == "classification" else REGRESSION[ep]
        df = _load(ep).dropna(subset=["smiles", target]).reset_index(drop=True)
        X, mask = featurize(df["smiles"].tolist())
        df = df.loc[mask].reset_index(drop=True)
        y = df[target].to_numpy()
        y = y.astype(int) if task == "classification" else y.astype(float)
        groups = _scaffold_groups(df["smiles"].tolist())
        pos_weight = (len(y) - y.sum()) / max(y.sum(), 1) if task == "classification" else None
        models = classifiers(pos_weight) if task == "classification" else regressors()
        print(f"\n[{ep}] {task}, {len(y)} compounds")
        for name, model in models.items():
            for split in ("random", "scaffold"):
                mean, sd = cv_score(task, model, X, y, groups, split)
                rows.append({"endpoint": ep, "task": task, "model": name, "split": split,
                             "metric": "roc_auc" if task == "classification" else "r2",
                             "mean": round(mean, 4), "sd": round(sd, 4)})
            sc = [r for r in rows if r["endpoint"] == ep and r["model"] == name and r["split"] == "scaffold"][0]
            print(f"    {name:22s} scaffold {sc['metric']} = {sc['mean']:.3f}")
    out = pd.DataFrame(rows)
    (ROOT / "results" / "tables").mkdir(parents=True, exist_ok=True)
    out.to_csv(ROOT / "results" / "tables" / "model_comparison.csv", index=False)

    # headline: mean scaffold metric per model
    scaf = out[out.split == "scaffold"]
    print("\n=== mean scaffold metric by model ===")
    for task in ("classification", "regression"):
        t = scaf[scaf.task == task]
        print(f"  {task}:")
        for name in ("RandomForest", "XGBoost", "HistGradientBoosting"):
            print(f"    {name:22s} {t[t.model == name]['mean'].mean():.4f}")
    print("wrote results/tables/model_comparison.csv")


if __name__ == "__main__":
    main()
