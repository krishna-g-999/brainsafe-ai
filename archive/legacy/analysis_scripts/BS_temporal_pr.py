"""
BS_temporal_pr.py — temporal (time-split) validation + PR-curve / threshold-sensitivity.

  * TEMPORAL: per ChEMBL target, assign each compound its earliest document_year, train on
    compounds reported <= 75th-percentile year, test on the most recent ~25%. This is the
    gold-standard "will it work on future compounds" test (computational, no wet-lab).
  * PR / THRESHOLD: precision, recall, F1 at thresholds 0.3/0.5/0.7 and Youden-J, plus
    PR-AUC, on a held-out set — directly addresses imbalanced-endpoint reporting.
  * BBB (B3DB, undated) gets a scaffold-holdout PR/threshold report instead of temporal.
Out: BS_temporal_pr_report.json
"""
import os, glob, json, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from rdkit import Chem
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import (roc_auc_score, average_precision_score, roc_curve,
                             precision_recall_fscore_support)
from BS_predictive_model import morgan, descriptors, scaffold, bvs
from BS_train_endpoints import models


def canon_y(df):
    has_y = "year" in df.columns
    rows = []
    for _, r in df.iterrows():
        m = Chem.MolFromSmiles(str(r["smiles"]))
        if m is None:
            continue
        try: ik = Chem.MolToInchiKey(m)
        except Exception: ik = Chem.MolToSmiles(m)
        yr = float(r["year"]) if (has_y and pd.notna(r["year"])) else None
        rows.append((Chem.MolToSmiles(m), int(r["label"]), ik, yr))
    return pd.DataFrame(rows, columns=["smiles", "label", "ik", "year"]).drop_duplicates("ik").reset_index(drop=True)


def fit_predict(Xtr, ytr, Xte):
    ps = []
    for m in models().values():
        m.fit(Xtr, ytr); ps.append(m.predict_proba(Xte)[:, 1])
    return np.mean(ps, axis=0)


def thresholds(y, p):
    out = {}
    for t in (0.3, 0.5, 0.7):
        pr, rc, f1, _ = precision_recall_fscore_support(y, (p >= t).astype(int), average="binary", zero_division=0)
        out[f"@{t}"] = {"precision": round(pr, 3), "recall": round(rc, 3), "f1": round(f1, 3)}
    fpr, tpr, thr = roc_curve(y, p); j = int(np.argmax(tpr - fpr)); yt = float(thr[j])
    pr, rc, f1, _ = precision_recall_fscore_support(y, (p >= yt).astype(int), average="binary", zero_division=0)
    out["@youden"] = {"threshold": round(yt, 3), "precision": round(pr, 3), "recall": round(rc, 3), "f1": round(f1, 3)}
    return out


def evaluate(name, df):
    d = canon_y(df)
    if d["label"].nunique() < 2 or len(d) < 120:
        return None
    smi = d["smiles"].tolist(); y = d["label"].values.astype(int)
    X = np.hstack([morgan(smi), descriptors(smi)])
    res = {"n": len(d), "pos_rate": round(float(y.mean()), 3)}

    # --- temporal split ---
    dated = int(d["year"].notna().sum())
    if dated >= 0.6 * len(d) and d["year"].nunique() > 3:
        cut = int(d["year"].quantile(0.75))
        te = (d["year"] > cut).values; tr = (~te) & d["year"].notna().values
        if te.sum() >= 40 and len(np.unique(y[te])) == 2 and len(np.unique(y[tr])) == 2:
            p = fit_predict(X[tr], y[tr], X[te])
            res["temporal"] = {"cutoff_year": cut, "n_train": int(tr.sum()), "n_test": int(te.sum()),
                               "test_years": f">{cut}", "pos_rate_test": round(float(y[te].mean()), 3),
                               "auroc": round(float(roc_auc_score(y[te], p)), 3),
                               "pr_auc": round(float(average_precision_score(y[te], p)), 3),
                               "thresholds": thresholds(y[te], p)}

    # --- scaffold holdout (PR/threshold; only split available for undated BBB) ---
    scaf = np.array([scaffold(s) for s in smi])
    tr, te = next(GroupShuffleSplit(1, test_size=0.25, random_state=42).split(X, y, groups=scaf))
    if len(np.unique(y[te])) == 2:
        p = fit_predict(X[tr], y[tr], X[te])
        res["scaffold_holdout"] = {"n_test": int(len(te)), "pos_rate_test": round(float(y[te].mean()), 3),
                                   "auroc": round(float(roc_auc_score(y[te], p)), 3),
                                   "pr_auc": round(float(average_precision_score(y[te], p)), 3),
                                   "thresholds": thresholds(y[te], p)}
    return res


def main():
    rep = {}
    for f in sorted(glob.glob("data/endpoints/*.csv")):
        nm = os.path.basename(f).replace(".csv", "")
        if nm.startswith("_"):
            continue
        r = evaluate(nm, pd.read_csv(f))
        if r:
            rep[nm] = r
            t = r.get("temporal", {})
            print(f"  [{nm:6}] temporal AUROC={t.get('auroc','-')} (train<= {t.get('cutoff_year','-')}, "
                  f"test n={t.get('n_test','-')}) | scaffold-holdout AUROC={r.get('scaffold_holdout',{}).get('auroc','-')}")
    json.dump(rep, open("BS_temporal_pr_report.json", "w"), indent=2)
    print("Saved BS_temporal_pr_report.json")


if __name__ == "__main__":
    main()
