"""
BS_significance.py
------------------
Reviewer point (2): the ensemble-vs-baseline AUROC deltas need a significance test,
not just point estimates. For each of the eight deployed classification endpoints we
pair the deployed ensemble's out-of-fold scores (from the figures OOF cache, built with
scaffold GroupKFold(5)) with kNN-Tanimoto and logistic-regression scores computed on the
IDENTICAL folds, then run:
  - DeLong's test for two correlated ROC AUCs (analytic p-value), and
  - a paired bootstrap (2,000 resamples) for the delta 95% CI and a one-sided p-value.

GroupKFold is deterministic, so re-loading the same data + scaffolds reproduces the exact
folds behind the cached ensemble OOF; we assert label alignment before pairing.
Nothing is fabricated: every number is computed here from the measured data.

Outputs: BS_significance_report.json ; supplementary/STable14_significance.csv
"""
import os, json
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from scipy import stats
from rdkit import DataStructs
from BS_predictive_model import morgan, descriptors, scaffold, bvs
from BS_train_endpoints import canon as clf_canon

EPS = ["BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_B", "MAO_A", "hERG"]
CACHE = "figures/_oof_cache.npz"
cache = {k: np.load(CACHE, allow_pickle=True)[k] for k in np.load(CACHE, allow_pickle=True).files}

# ---------- fast DeLong (Sun & Xu 2014) for two correlated AUCs ----------
def _midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x); T = np.zeros(N); i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]: j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1; i = j
    T2 = np.empty(N); T2[J] = T; return T2

def delong(y, p1, p2):
    order = (-y).argsort(kind="mergesort"); m = int(y.sum())
    preds = np.vstack((p1, p2))[:, order]; n = preds.shape[1] - m
    pos, neg = preds[:, :m], preds[:, m:]
    tx = np.vstack([_midrank(pos[r]) for r in range(2)])
    ty = np.vstack([_midrank(neg[r]) for r in range(2)])
    tz = np.vstack([_midrank(preds[r]) for r in range(2)])
    aucs = tz[:, :m].sum(1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n; v10 = 1.0 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n
    L = np.array([[1.0, -1.0]]); var = float((L @ cov @ L.T).item())
    z = (aucs[0] - aucs[1]) / np.sqrt(var) if var > 0 else 0.0
    return float(aucs[0]), float(aucs[1]), z, float(2 * stats.norm.sf(abs(z)))

def paired_boot(y, p1, p2, B=2000, seed=42):
    rng = np.random.default_rng(seed); n = len(y); d = []
    for _ in range(B):
        s = rng.integers(0, n, n); ys = y[s]
        if ys.sum() == 0 or ys.sum() == n: continue
        d.append(roc_auc_score(ys, p1[s]) - roc_auc_score(ys, p2[s]))
    d = np.array(d)
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), float((d <= 0).mean())

# kNN + logistic OOF replicate BS_figures_v2.py exactly (the source of the reported
# baseline table): similarity-weighted top-5 Tanimoto kNN via RDKit; balanced logistic.
def baselines_oof(X, bv, y, g):
    knn = np.zeros(len(y)); lr = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, groups=g):
        bt = [bv[i] for i in tr]
        for i in te:
            sims = np.array(DataStructs.BulkTanimotoSimilarity(bv[i], bt))
            idx = np.argsort(sims)[::-1][:5]; w = sims[idx]
            knn[i] = (w * y[tr][idx]).sum() / w.sum() if w.sum() > 0 else y[tr].mean()
        sc = StandardScaler().fit(X[tr])
        lr[te] = LogisticRegression(max_iter=2000, class_weight="balanced").fit(
            sc.transform(X[tr]), y[tr]).predict_proba(sc.transform(X[te]))[:, 1]
    return knn, lr

rep = {}; rows = []
for ep in EPS:
    y_cache = cache[f"clf_{ep}_y"].astype(int); p_ens = cache[f"clf_{ep}_p"].astype(float)
    d = clf_canon(pd.read_csv(f"data/endpoints/{ep}.csv"))
    smi = d["smiles"].tolist(); y = d["label"].values.astype(int)
    assert len(y) == len(y_cache) and (y == y_cache).all(), f"{ep}: data/cache misalignment"
    X = np.hstack([np.asarray(morgan(smi)), np.asarray(descriptors(smi))]).astype(np.float32)
    g = np.array([scaffold(s) for s in smi]); bv = bvs(smi)
    p_knn, p_lr = baselines_oof(X, bv, y, g)
    a_ens, a_knn, z_k, p_k = delong(y, p_ens, p_knn)
    lo_k, hi_k, pb_k = paired_boot(y, p_ens, p_knn)
    _, a_lr, z_l, p_l = delong(y, p_ens, p_lr)
    lo_l, hi_l, pb_l = paired_boot(y, p_ens, p_lr)
    rep[ep] = {
        "n": int(len(y)), "AUROC_ensemble": round(a_ens, 3),
        "AUROC_kNN": round(a_knn, 3), "delta_vs_kNN": round(a_ens - a_knn, 3),
        "delong_p_vs_kNN": p_k, "boot_delta95_vs_kNN": [round(lo_k, 3), round(hi_k, 3)],
        "AUROC_logistic": round(a_lr, 3), "delta_vs_LR": round(a_ens - a_lr, 3),
        "delong_p_vs_LR": p_l, "boot_delta95_vs_LR": [round(lo_l, 3), round(hi_l, 3)],
        "sig_vs_kNN_0.05": bool(p_k < 0.05), "sig_vs_LR_0.05": bool(p_l < 0.05),
    }
    rows.append({"endpoint": ep, "AUROC_ensemble": round(a_ens, 3), "AUROC_kNN": round(a_knn, 3),
                 "delta_vs_kNN": round(a_ens - a_knn, 3), "delta95_CI": f"[{lo_k:.3f}, {hi_k:.3f}]",
                 "DeLong_p": f"{p_k:.4f}", "significant_0.05": p_k < 0.05})
    print(f"{ep:6} ens {a_ens:.3f} vs kNN {a_knn:.3f}  d={a_ens-a_knn:+.3f}  DeLong p={p_k:.4f}  "
          f"boot95=[{lo_k:+.3f},{hi_k:+.3f}]  {'SIG' if p_k<0.05 else 'n.s.'}", flush=True)

n_sig = sum(rep[e]["sig_vs_kNN_0.05"] for e in EPS)
rep["_summary"] = {"n_endpoints": len(EPS), "n_significant_vs_kNN_0.05": n_sig,
                   "not_significant_vs_kNN": [e for e in EPS if not rep[e]["sig_vs_kNN_0.05"]]}
json.dump(rep, open("BS_significance_report.json", "w"), indent=2)
pd.DataFrame(rows).to_csv("supplementary/STable14_significance.csv", index=False)
print(f"\n{n_sig}/{len(EPS)} endpoints significant vs kNN at p<0.05; not significant: "
      f"{rep['_summary']['not_significant_vs_kNN']}")
print("Wrote BS_significance_report.json + STable14_significance.csv")
