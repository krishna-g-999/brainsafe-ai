"""
BS_cv_comparison.py
-------------------
Answers the presentation question: "did you use Random Forest, and why not 10-fold?"
For each of the eight deployed classification endpoints we compute, on the SAME data:
  - scaffold GroupKFold(5) ensemble AUROC  (the deployed, honest number; from the OOF cache)
  - random StratifiedKFold(10) ensemble AUROC
  - random StratifiedKFold(10) Random-Forest-only AUROC
  - random StratifiedKFold(5)  ensemble AUROC
This shows (a) Random Forest IS used (it is one of the three ensemble members and is
reported here on its own), and (b) random k-fold inflates AUROC relative to the
scaffold split, because analogue-dense ChEMBL series leak between random folds.
Nothing is fabricated. Output: BS_cv_comparison.json ; supplementary/STable15_cv_comparison.csv
"""
import os, json
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from BS_predictive_model import morgan, descriptors
from BS_train_endpoints import canon as clf_canon, models as clf_models

EPS = ["BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_B", "MAO_A", "hERG"]
disp = {"MAO_A": "MAO-A", "MAO_B": "MAO-B", "GSK3B": "GSK-3β"}
cache = {k: np.load("figures/_oof_cache.npz", allow_pickle=True)[k]
         for k in np.load("figures/_oof_cache.npz", allow_pickle=True).files}

def feats(smi):
    return np.hstack([np.asarray(morgan(smi)), np.asarray(descriptors(smi))]).astype(np.float32)

def cv_auroc(X, y, k, members):
    oof = np.zeros(len(y))
    for tr, te in StratifiedKFold(k, shuffle=True, random_state=42).split(X, y):
        oof[te] = np.mean([m.fit(X[tr], y[tr]).predict_proba(X[te])[:, 1] for m in members], axis=0)
    return roc_auc_score(y, oof)

rows = []
for ep in EPS:
    d = clf_canon(pd.read_csv(f"data/endpoints/{ep}.csv"))
    smi = d["smiles"].tolist(); y = d["label"].values.astype(int); X = feats(smi)
    scaf = roc_auc_score(cache[f"clf_{ep}_y"].astype(int), cache[f"clf_{ep}_p"].astype(float))  # deployed
    ens10 = cv_auroc(X, y, 10, list(clf_models().values()))
    rf10 = cv_auroc(X, y, 10, [clf_models()["rf"]])
    ens5 = cv_auroc(X, y, 5, list(clf_models().values()))
    rows.append({"endpoint": disp.get(ep, ep), "n": int(len(y)),
                 "scaffold5_ensemble": round(scaf, 3),
                 "random10_ensemble": round(ens10, 3),
                 "random10_RFonly": round(rf10, 3),
                 "random5_ensemble": round(ens5, 3),
                 "inflation_random10_minus_scaffold": round(ens10 - scaf, 3)})
    print(f"{ep:6} scaffold5-ens {scaf:.3f} | random10-ens {ens10:.3f} | random10-RF {rf10:.3f} | "
          f"random5-ens {ens5:.3f} | inflation +{ens10-scaf:.3f}", flush=True)

df = pd.DataFrame(rows)
df.to_csv("supplementary/STable15_cv_comparison.csv", index=False)
summ = {"mean_scaffold5_ensemble": round(df["scaffold5_ensemble"].mean(), 3),
        "mean_random10_ensemble": round(df["random10_ensemble"].mean(), 3),
        "mean_random10_RFonly": round(df["random10_RFonly"].mean(), 3),
        "mean_inflation_random10_vs_scaffold": round(df["inflation_random10_minus_scaffold"].mean(), 3)}
json.dump({"per_endpoint": rows, "summary": summ}, open("BS_cv_comparison.json", "w"), indent=2)
print("\nSUMMARY:", json.dumps(summ))
print("Wrote BS_cv_comparison.json + STable15_cv_comparison.csv")
