"""Reproduce the published random 10-fold CV, then repeat it with duplicates collapsed.

Arm A reproduces train_rf.py exactly (same seed, same splitter, same RF hyper-parameters) to
establish that the reported number regenerates. Arm B is identical except that rows sharing an
identical feature vector are collapsed to one row before splitting, so no compound the model
cannot distinguish can appear in both train and test. The difference between A and B is the
part of the reported score that comes from duplicate rows rather than from generalisation.
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import RDLogger
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

RDLogger.DisableLog("rdApp.*")
ROOT = Path(r"D:\BRAINSAFE_AI")
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize  # noqa: E402

SEED = 42
N_SPLITS = 10
RF_COMMON = dict(n_estimators=300, min_samples_leaf=2, n_jobs=-1, random_state=SEED)

ENDPOINTS = ["BBB", "AChE", "MAO_A", "hERG"]


def cv_auroc(X, y):
    skf = StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED)
    aucs = []
    for tr, te in skf.split(X, y):
        m = RandomForestClassifier(class_weight="balanced", **RF_COMMON)
        m.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs))


published = pd.read_csv(ROOT / "results" / "tables" / "rf_cv_summary.csv")
out = []
for ep in ENDPOINTS:
    df = pd.read_csv(ROOT / "data" / "endpoints" / f"{ep}.csv")
    X, mask = featurize(df["smiles"].astype(str).tolist())
    y = df.loc[mask, "label"].to_numpy().astype(int)

    a_mean, a_sd = cv_auroc(X, y)

    groups = defaultdict(list)
    for i, vec in enumerate(X):
        groups[vec.tobytes()].append(i)
    keep, dropped_conflict = [], 0
    for idxs in groups.values():
        ys = {y[i] for i in idxs}
        if len(ys) > 1:
            dropped_conflict += 1
            continue          # contradictory duplicate group: unusable, exclude entirely
        keep.append(idxs[0])
    keep = np.array(sorted(keep))
    Xd, yd = X[keep], y[keep]
    b_mean, b_sd = cv_auroc(Xd, yd)

    row = published[(published.endpoint == ep) & (published.split == "random")]
    pub = float(row["roc_auc_mean"].iloc[0]) if len(row) else float("nan")
    pub_n = int(row["n"].iloc[0]) if len(row) else -1
    out.append({
        "endpoint": ep,
        "published_n": pub_n, "published_auroc": pub,
        "reproduced_n": len(y), "reproduced_auroc": round(a_mean, 4), "reproduced_sd": round(a_sd, 4),
        "delta_vs_published": round(a_mean - pub, 4),
        "dedup_n": len(yd), "rows_removed": len(y) - len(yd),
        "conflicting_groups_dropped": dropped_conflict,
        "dedup_auroc": round(b_mean, 4), "dedup_sd": round(b_sd, 4),
        "inflation_from_duplicates": round(a_mean - b_mean, 4),
    })
    print(out[-1], flush=True)

res = pd.DataFrame(out)
res.to_csv(Path(__file__).with_name("repro_cv.csv"), index=False)
print()
print(res.to_string(index=False))
