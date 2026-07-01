"""
BS_external_validation.py — publication-grade, no-compromise validation.

For every measured endpoint:
  1. SCAFFOLD GroupKFold(5) out-of-fold (OOF) probabilities (same ensemble as deployed).
  2. SIMILARITY-BINNED AUROC: performance as a function of each test compound's max
     Tanimoto to its training fold -> the honest novel-chemotype generalisation curve.
  3. STRICT LEAVE-CLUSTER-OUT: LeaderPicker sphere-exclusion clusters (distance 0.4);
     GroupShuffleSplit holds out whole clusters -> AUROC on structurally novel clusters.
  4. MONDRIAN (class-conditional) INDUCTIVE CONFORMAL PREDICTION calibrated on OOF;
     empirically validated coverage at 90%; per-class calibration saved for the engine.

Outputs: BS_external_validation_report.json, models_brain/<ep>_conformal.json
"""
import os, glob, json, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.SimDivFilters.rdSimDivPickers import LeaderPicker
from sklearn.model_selection import GroupKFold, GroupShuffleSplit
from sklearn.metrics import roc_auc_score
from BS_predictive_model import morgan, descriptors, scaffold, bvs
from BS_train_endpoints import models, canon

SEED = 42
OUT_REPORT = "BS_external_validation_report.json"


def ensemble_oof_proba(X, y, groups, splitter):
    oof = np.full(len(y), np.nan)
    for tr, te in splitter:
        ps = []
        for m in models().values():
            m.fit(X[tr], y[tr]); ps.append(m.predict_proba(X[te])[:, 1])
        oof[te] = np.mean(ps, axis=0)
    return oof


def leader_clusters(fps, dist=0.4):
    lp = LeaderPicker()
    leaders = list(lp.LazyBitVectorPick(fps, len(fps), dist))
    # assign each compound to the most similar leader
    lead_fps = [fps[i] for i in leaders]
    cl = np.zeros(len(fps), dtype=int)
    for i, fp in enumerate(fps):
        sims = DataStructs.BulkTanimotoSimilarity(fp, lead_fps)
        cl[i] = int(np.argmax(sims))
    return cl, len(leaders)


def sim_binned_auroc(y, oof, maxsim):
    bins = [(0.0, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01)]
    res = {}
    for lo, hi in bins:
        m = (maxsim >= lo) & (maxsim < hi)
        key = f"T[{lo:.1f}-{hi:.1f})"
        if m.sum() > 20 and len(np.unique(y[m])) == 2:
            res[key] = {"n": int(m.sum()), "auroc": round(float(roc_auc_score(y[m], oof[m])), 3)}
        else:
            res[key] = {"n": int(m.sum()), "auroc": None}
    return res


def mondrian_conformal(y, p, eps=0.10):
    """Class-conditional inductive conformal. Split OOF 50/50 into calib/test;
    nonconformity = 1 - prob(true class). Report coverage + mean set size at 1-eps."""
    rng = np.random.default_rng(SEED)
    idx = rng.permutation(len(y)); half = len(y) // 2
    cal, tst = idx[:half], idx[half:]
    A1 = sorted(1 - p[cal][y[cal] == 1])   # active calib nonconf
    A0 = sorted(p[cal][y[cal] == 0])       # inactive calib nonconf
    def pval(prob, cls):
        a = (1 - prob) if cls == 1 else prob
        ref = A1 if cls == 1 else A0
        if not ref: return 1.0
        ge = sum(1 for x in ref if x >= a)
        return (ge + 1) / (len(ref) + 1)
    covered, sizes = 0, []
    for i in tst:
        ps = {c: pval(p[i], c) for c in (0, 1)}
        pset = [c for c in (0, 1) if ps[c] > eps]
        sizes.append(len(pset))
        if y[i] in pset: covered += 1
    return {"target_coverage": 1 - eps,
            "empirical_coverage": round(covered / len(tst), 3),
            "mean_set_size": round(float(np.mean(sizes)), 2),
            "frac_singletons": round(float(np.mean([s == 1 for s in sizes])), 3),
            # full-data calibration scores for deployment
            "calib_active_nonconf": [round(float(x), 4) for x in sorted(1 - p[y == 1])],
            "calib_inactive_nonconf": [round(float(x), 4) for x in sorted(p[y == 0])]}


def evaluate(name, df):
    df = canon(df)
    if df["label"].nunique() < 2 or len(df) < 100:
        print(f"  [{name}] insufficient"); return None
    smi = df["smiles"].tolist(); y = df["label"].values.astype(int)
    X = np.hstack([morgan(smi), descriptors(smi)]); bv = bvs(smi)
    scaf = np.array([scaffold(s) for s in smi])

    # 1) scaffold OOF + max sim to train fold
    gkf = list(GroupKFold(5).split(X, groups=scaf))
    oof = ensemble_oof_proba(X, y, scaf, gkf)
    maxsim = np.zeros(len(y))
    for tr, te in gkf:
        bt = [bv[i] for i in tr if bv[i] is not None]
        for i in te:
            maxsim[i] = max(DataStructs.BulkTanimotoSimilarity(bv[i], bt)) if (bv[i] is not None and bt) else 0.0
    auroc_scaffold = round(float(roc_auc_score(y, oof)), 3)
    binned = sim_binned_auroc(y, oof, maxsim)

    # 2) strict leave-cluster-out (sphere-exclusion clusters)
    cl, n_lead = leader_clusters([f for f in bv], dist=0.4)
    gss = list(GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=SEED).split(X, y, groups=cl))
    tr, te = gss[0]
    auroc_cluster = None; cluster_leak = None
    if len(np.unique(y[te])) == 2:
        ps = []
        for m in models().values():
            m.fit(X[tr], y[tr]); ps.append(m.predict_proba(X[te])[:, 1])
        pte = np.mean(ps, axis=0)
        auroc_cluster = round(float(roc_auc_score(y[te], pte)), 3)
        bt = [bv[i] for i in tr if bv[i] is not None]
        ms = [max(DataStructs.BulkTanimotoSimilarity(bv[i], bt)) for i in te if bv[i] is not None]
        cluster_leak = round(float(np.median(ms)), 3)

    # 3) conformal on scaffold-OOF
    conf = mondrian_conformal(y, oof)
    json.dump({"calib_active_nonconf": conf.pop("calib_active_nonconf"),
               "calib_inactive_nonconf": conf.pop("calib_inactive_nonconf"),
               "eps": 0.10},
              open(f"models_brain/{name}_conformal.json", "w"))

    print(f"  [{name:6}] scaffold AUROC={auroc_scaffold} | cluster-split AUROC={auroc_cluster} "
          f"(median T={cluster_leak}, {n_lead} clusters) | conformal cov={conf['empirical_coverage']} "
          f"setsize={conf['mean_set_size']}")
    print(f"           similarity-binned: " +
          " ".join(f"{k}={v['auroc']}(n{v['n']})" for k, v in binned.items()))
    return {"n": len(y), "auroc_scaffold_cv": auroc_scaffold,
            "auroc_strict_cluster_split": auroc_cluster, "cluster_split_median_tanimoto": cluster_leak,
            "n_clusters": n_lead, "similarity_binned_auroc": binned, "conformal": conf}


def main():
    rep = {}
    for f in sorted(glob.glob("data/endpoints/*.csv")):
        nm = os.path.basename(f).replace(".csv", "")
        if nm.startswith("_"): continue
        r = evaluate(nm, pd.read_csv(f))
        if r: rep[nm] = r
    json.dump(rep, open(OUT_REPORT, "w"), indent=2)
    print("\nSaved", OUT_REPORT, "and per-endpoint conformal calibrations.")


if __name__ == "__main__":
    main()
