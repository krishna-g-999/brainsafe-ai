"""
BS_assay_sensitivity.py -- single-assay sensitivity for threat-to-validity #1.
Using the standard_type-tagged cache (data/_chembl_cache/<name>_std.json), retrain scaffold-CV
on (a) the deployed POOLED set and (b) the dominant SINGLE assay type only (IC50), for the most
heterogeneous endpoints. If AUROC is stable, pooling on the pChEMBL scale is empirically justified.

Output: appends to BS_flaw_fixes.json (key 'assay_type_sensitivity');
        supplementary/STable12_assay_sensitivity.csv
"""
import os, json, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from rdkit import Chem
from sklearn.ensemble import (RandomForestClassifier, ExtraTreesClassifier,
                              HistGradientBoostingClassifier)
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from BS_predictive_model import morgan, descriptors, scaffold
SEED = 42
OUT = "supplementary"

# GSK3B is the most heterogeneous (IC50 0.49); MAO_B strongly IC50 (0.89); hERG large IC50 (0.85)
ENDPOINTS = ["GSK3B", "MAO_B", "hERG"]

def models():
    return {"rf": RandomForestClassifier(n_estimators=300, min_samples_leaf=2,
                                         class_weight="balanced_subsample", n_jobs=-1, random_state=SEED),
            "et": ExtraTreesClassifier(n_estimators=300, min_samples_leaf=2,
                                       class_weight="balanced_subsample", n_jobs=-1, random_state=SEED),
            "hgb": HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06, random_state=SEED)}

def compound_labels(rows, only_type=None):
    df = pd.DataFrame(rows)
    if only_type:
        df = df[df["std"] == only_type]
    rec = {}
    for s, v in zip(df["smiles"], df["pchembl"]):
        m = Chem.MolFromSmiles(str(s))
        if m is None:
            continue
        try:
            ik = Chem.MolToInchiKey(m)
        except Exception:
            ik = Chem.MolToSmiles(m)
        rec.setdefault(ik, [Chem.MolToSmiles(m), []])[1].append(float(v))
    smis = [x[0] for x in rec.values()]
    med = np.array([np.median(x[1]) for x in rec.values()])
    y = np.full(len(med), -1); y[med >= 6.0] = 1; y[med < 5.0] = 0
    mask = y >= 0
    return [s for s, k in zip(smis, mask) if k], y[mask].astype(int)

def auroc(smis, y):
    if len(set(y)) < 2 or min(int(y.sum()), int((1 - y).sum())) < 20:
        return None, len(y)
    X = np.hstack([morgan(smis), descriptors(smis)])
    g = np.array([scaffold(s) for s in smis])
    gkf = GroupKFold(min(5, len(set(g)))); oof = np.zeros(len(y))
    for tr, te in gkf.split(X, groups=g):
        ps = [m.fit(X[tr], y[tr]).predict_proba(X[te])[:, 1] for m in models().values()]
        oof[te] = np.mean(ps, axis=0)
    return round(float(roc_auc_score(y, oof)), 3), len(y)

rows12 = []
for ep in ENDPOINTS:
    rows = json.load(open(f"data/_chembl_cache/{ep}_std.json"))
    sp, yp = compound_labels(rows)                 # pooled (all assay types)
    si, yi = compound_labels(rows, only_type="IC50")  # single assay type
    ap, np_ = auroc(sp, yp)
    ai, ni = auroc(si, yi)
    delta = round(ap - ai, 3) if (ap and ai) else None
    rows12.append({"endpoint": ep, "pooled_n": np_, "pooled_AUROC": ap,
                   "IC50only_n": ni, "IC50only_AUROC": ai, "AUROC_delta_pooled_minus_IC50": delta})
    print(f"  [{ep}] pooled AUROC={ap} (n={np_}) | IC50-only AUROC={ai} (n={ni}) | Δ={delta}")

pd.DataFrame(rows12).to_csv(f"{OUT}/STable12_assay_sensitivity.csv", index=False)
try:
    b = json.load(open("BS_flaw_fixes.json"))
except Exception:
    b = {}
b["assay_type_sensitivity"] = rows12
b["assay_type_sensitivity_note"] = ("Scaffold-CV AUROC of the deployed ensemble trained on the "
    "pooled multi-assay set vs the dominant single assay type (IC50). Small |Δ| indicates pooling on "
    "the standardised pChEMBL scale does not materially distort discrimination.")
json.dump(b, open("BS_flaw_fixes.json", "w"), indent=2)
print("Wrote STable12 + updated BS_flaw_fixes.json")
