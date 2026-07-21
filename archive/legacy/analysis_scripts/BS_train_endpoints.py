"""
BS_train_endpoints.py — train genuine, scaffold-validated classifiers for each
measured brain-relevant endpoint (BBB + CNS targets).

Rigor:
  * RDKit canonicalisation + InChIKey dedup.
  * Features: Morgan-1024 + 24 RDKit descriptors (structure only).
  * Honest split: scaffold GroupKFold(5) on Bemis-Murcko generic scaffolds;
    out-of-fold probabilities -> AUROC, PR-AUC, balanced accuracy, MCC.
    (AUROC ~0.5 == no signal; class imbalance handled by AUROC/MCC + balanced models.)
  * Applicability domain: max Tanimoto to a 2000-sample of training compounds.
  * Final model trained on all data; saved with metrics, operating threshold
    (Youden's J on OOF) and AD fingerprints.

Out: models_brain/<endpoint>.joblib, models_brain/<endpoint>_meta.json,
     models_brain/endpoints_report.json
"""
import os, json, glob, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, joblib
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.ensemble import (RandomForestClassifier, ExtraTreesClassifier,
                              HistGradientBoostingClassifier)
from sklearn.model_selection import GroupKFold
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             balanced_accuracy_score, matthews_corrcoef, roc_curve,
                             brier_score_loss)
from sklearn.isotonic import IsotonicRegression
from BS_predictive_model import morgan, descriptors, scaffold, bvs

OUT = "models_brain"; os.makedirs(OUT, exist_ok=True)
SEED = 42
MEANING = {"BBB":"blood-brain-barrier penetration","AChE":"Alzheimer's / cognition (AChE inhib.)",
           "BChE":"Alzheimer's / cholinergic (BuChE inhib.)",
           "BACE1":"Alzheimer's / amyloid (BACE1 inhib.)","MAO_B":"Parkinson's / dopamine (MAO-B inhib.)",
           "MAO_A":"mood / depression (MAO-A inhib.)","GSK3B":"tau / neuroprotection (GSK-3b inhib.)",
           "D2":"Parkinson's / psychosis (dopamine D2)","A2A":"Parkinson's (adenosine A2A)",
           "HT2A":"mood / psychosis (5-HT2A)","SERT":"depression (serotonin transporter)",
           "hERG":"SAFETY: cardiotoxicity (hERG block)"}


def canon(df):
    has_val = "pchembl" in df.columns
    vals = df["pchembl"].tolist() if has_val else [None]*len(df)
    out=[]
    for s,l,v in zip(df["smiles"],df["label"],vals):
        m=Chem.MolFromSmiles(str(s))
        if m is None: continue
        try: ik=Chem.MolToInchiKey(m)
        except Exception: ik=Chem.MolToSmiles(m)
        out.append((Chem.MolToSmiles(m), int(l), ik, (float(v) if v is not None else None)))
    d=pd.DataFrame(out,columns=["smiles","label","ik","value"]).drop_duplicates("ik")
    return d.reset_index(drop=True)


def models():
    return {
        "rf": RandomForestClassifier(n_estimators=300, min_samples_leaf=2,
                                     class_weight="balanced_subsample", n_jobs=-1, random_state=SEED),
        "et": ExtraTreesClassifier(n_estimators=300, min_samples_leaf=2,
                                   class_weight="balanced_subsample", n_jobs=-1, random_state=SEED),
        "hgb": HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06,
                                              random_state=SEED),
    }


def proba(mdls, Xtr, ytr, Xte):
    ps=[]
    for name,m in mdls.items():
        m.fit(Xtr,ytr); ps.append(m.predict_proba(Xte)[:,1])
    return np.mean(ps,axis=0)


def evaluate(name, df):
    df=canon(df)
    if df["label"].nunique()<2 or len(df)<80:
        print(f"  [{name}] insufficient"); return None
    smi=df["smiles"].tolist(); y=df["label"].values.astype(int)
    X=np.hstack([morgan(smi), descriptors(smi)])
    groups=np.array([scaffold(s) for s in smi])
    bv=bvs(smi)
    gkf=GroupKFold(min(5,len(set(groups))))
    oof=np.zeros(len(y)); cross=[]
    for tr,te in gkf.split(X,groups=groups):
        oof[te]=proba(models(),X[tr],y[tr],X[te])
        bt=[bv[i] for i in tr if bv[i] is not None]
        for i in te:
            if bv[i] is not None and bt: cross.append(max(DataStructs.BulkTanimotoSimilarity(bv[i],bt)))
    auroc=roc_auc_score(y,oof); prauc=average_precision_score(y,oof)
    # Isotonic calibration on scaffold-CV out-of-fold scores -> trustworthy probabilities.
    iso=IsotonicRegression(out_of_bounds="clip").fit(oof,y)
    oof_cal=iso.predict(oof)
    brier=brier_score_loss(y,oof_cal)
    fpr,tpr,thr=roc_curve(y,oof_cal); j=int(np.argmax(tpr-fpr)); cut=float(thr[j])
    pred=(oof_cal>=cut).astype(int)
    ba=balanced_accuracy_score(y,pred); mcc=matthews_corrcoef(y,pred)
    print(f"  [{name:6}] n={len(y):5d} pos={int(y.sum()):5d} | AUROC={auroc:.3f} PR-AUC={prauc:.3f} "
          f"BA={ba:.3f} MCC={mcc:.3f} Brier={brier:.3f} | scaffold-CV, leak median T={np.median(cross):.2f}")
    # final ensemble on all data
    final=models(); [m.fit(X,y) for m in final.values()]
    # Evidence index: every valid training compound (fp + smiles + label + measured value)
    ev=[i for i in range(len(bv)) if bv[i] is not None]
    has_val = "value" in df.columns
    evidence={"fps":[bv[i] for i in ev], "smiles":[smi[i] for i in ev],
              "label":[int(y[i]) for i in ev],
              "value":[(float(df["value"].iloc[i]) if has_val and pd.notna(df["value"].iloc[i]) else None) for i in ev]}
    joblib.dump({"models":final,"calibrator":iso,"evidence":evidence,
                 "threshold":cut,"n_bits":1024}, f"{OUT}/{name}.joblib")
    meta={"endpoint":name,"meaning":MEANING.get(name,name),"n":len(y),"pos":int(y.sum()),
          "pos_rate":round(float(y.mean()),3),"auroc":round(auroc,3),"pr_auc":round(prauc,3),
          "balanced_acc":round(ba,3),"mcc":round(mcc,3),"brier":round(float(brier),3),
          "calibrated":True,"threshold":round(cut,3),
          "cv":"scaffold GroupKFold(5)","leak_median_tanimoto":round(float(np.median(cross)),3),
          "source":"B3DB (measured)" if name=="BBB" else "ChEMBL measured pChEMBL (>=6 active, <5 inactive)"}
    json.dump(meta,open(f"{OUT}/{name}_meta.json","w"),indent=2)
    return meta


def main():
    rep={}
    for f in sorted(glob.glob("data/endpoints/*.csv")):
        name=os.path.basename(f).replace(".csv","")
        if name.startswith("_"): continue
        meta=evaluate(name, pd.read_csv(f))
        if meta: rep[name]=meta
    json.dump(rep,open(f"{OUT}/endpoints_report.json","w"),indent=2)
    print("\nSaved", len(rep), "genuine endpoint models to", OUT)
    print(f"\n{'endpoint':7} {'AUROC':>6} {'PR-AUC':>7} {'MCC':>6} {'n':>6}  meaning")
    for k,m in rep.items():
        print(f"{k:7} {m['auroc']:6.3f} {m['pr_auc']:7.3f} {m['mcc']:6.3f} {m['n']:6d}  {m['meaning']}")


if __name__=="__main__":
    main()
