"""
BS_predictive_model.py — BrainSafe AI (BS)  Tier-1/2 genuine predictor build.

GOAL: determine, honestly, whether neuroprotection dimensions are predictable
FROM CHEMICAL STRUCTURE ALONE, and build + save genuine models only where a
real signal exists.

SCIENTIFIC LOCKS (no shortcuts):
  * LABELS: human/literature-curated ONLY (data_source != 'chembl_ml_predicted')
    -> removes the 35% ML pseudo-labels.
  * FEATURES: STRUCTURE ONLY (Morgan-1024 bits + RDKit physchem descriptor panel).
    The disease_target_count features are EXCLUDED (they are curation-coupled to
    the labels -> circular).  This is the honest test of structure -> activity.
  * SPLIT: scaffold GroupKFold(5) on Bemis-Murcko generic scaffolds; every
    transform fit inside the fold; leakage audited.
  * UNCERTAINTY: 2000x bootstrap 95% CI on out-of-fold R2.
  * BASELINES: mean predictor + Tanimoto k-NN. A dimension is declared
    GENUINELY PREDICTABLE only if ensemble R2 CI lower-bound > 0 AND it beats
    both baselines.
  * FINAL MODELS saved (models_genuine/) ONLY for predictable dimensions.

Run:  python BS_predictive_model.py
Out:  BS_predictive_report.json (+ models_genuine/ for predictable dims)
"""
from __future__ import annotations
import os, sys, json, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np, pandas as pd, joblib
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.ensemble import (RandomForestRegressor, ExtraTreesRegressor,
                              GradientBoostingRegressor)
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.Lipinski import NumHDonors, NumHAcceptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from scorer import neuro_score

DIMS = ["antioxidant","anti_inflammatory","mitochondrial_support","aggregation_modulation",
        "cognitive_enhancement","neurogenesis","synaptic_plasticity"]
HUMAN = {"curated","curated_v2","curated_smiles","literature_derived"}
DATA = "data/brainsafe_SCIENTIFIC_FIXED.csv"
SEED = 42
NBITS = 1024

_DESCS = [
    ("MolWt", Descriptors.MolWt), ("MolLogP", Crippen.MolLogP), ("TPSA", Descriptors.TPSA),
    ("HBD", NumHDonors), ("HBA", NumHAcceptors),
    ("RotB", rdMolDescriptors.CalcNumRotatableBonds),
    ("AromRings", rdMolDescriptors.CalcNumAromaticRings),
    ("AliphRings", rdMolDescriptors.CalcNumAliphaticRings),
    ("SatRings", rdMolDescriptors.CalcNumSaturatedRings),
    ("RingCount", rdMolDescriptors.CalcNumRings),
    ("FracCSP3", rdMolDescriptors.CalcFractionCSP3),
    ("Heteroatoms", lambda m: rdMolDescriptors.CalcNumHeteroatoms(m)),
    ("HeavyAtoms", lambda m: m.GetNumHeavyAtoms()),
    ("QED", lambda m: __import__("rdkit.Chem.QED", fromlist=["qed"]).qed(m)),
    ("MolMR", Crippen.MolMR), ("BertzCT", Descriptors.BertzCT),
    ("LabuteASA", Descriptors.LabuteASA), ("NHOH", Descriptors.NHOHCount),
    ("NO", Descriptors.NOCount),
    ("AromHeterocyc", rdMolDescriptors.CalcNumAromaticHeterocycles),
    ("FormalCharge", Chem.GetFormalCharge),
    ("ValenceE", Descriptors.NumValenceElectrons),
    # antioxidant-relevant structural motifs
    ("PhenolicOH", lambda m: sum(1 for a in m.GetAtoms() if a.GetAtomicNum()==8
                                 and a.GetTotalNumHs()>=1
                                 and any(n.GetIsAromatic() for n in a.GetNeighbors()))),
    ("Catechol", lambda m: 1 if m.HasSubstructMatch(Chem.MolFromSmarts("c1c(O)c(O)ccc1")) else 0),
]


def morgan(smiles):
    X = np.zeros((len(smiles), NBITS), dtype=np.float32)
    for i,s in enumerate(smiles):
        m = Chem.MolFromSmiles(s) if s else None
        if m: X[i] = np.array(AllChem.GetMorganFingerprintAsBitVect(m,2,NBITS), dtype=np.float32)
    return X

def bvs(smiles):
    out=[]
    for s in smiles:
        m=Chem.MolFromSmiles(s) if s else None
        out.append(AllChem.GetMorganFingerprintAsBitVect(m,2,NBITS) if m else None)
    return out

def descriptors(smiles):
    X=np.zeros((len(smiles),len(_DESCS)),dtype=np.float32)
    for i,s in enumerate(smiles):
        m=Chem.MolFromSmiles(s) if s else None
        if m is None: continue
        for j,(_,fn) in enumerate(_DESCS):
            try: X[i,j]=float(fn(m))
            except Exception: X[i,j]=0.0
    return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

def scaffold(s):
    m=Chem.MolFromSmiles(s) if s else None
    if m is None: return "INVALID"
    try:
        g=MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m))
        c=Chem.MolToSmiles(g); return c if c else "ACYCLIC"
    except Exception: return "ACYCLIC"

def ens_models():
    return {
        "rf":  RandomForestRegressor(n_estimators=400, max_depth=10, min_samples_leaf=2,
                                     random_state=SEED, n_jobs=-1),
        "et":  ExtraTreesRegressor(n_estimators=400, max_depth=10, min_samples_leaf=2,
                                   random_state=SEED, n_jobs=-1),
        "gbt": GradientBoostingRegressor(n_estimators=150, learning_rate=0.05, max_depth=3,
                                         random_state=SEED),
        "ridge": Ridge(alpha=10.0),
    }

def predict_ensemble(Xtr,ytr,Xte):
    """Per-dimension 4-model ensemble (GBT is single-output)."""
    sc=StandardScaler().fit(Xtr); Xtr_s,Xte_s=sc.transform(Xtr),sc.transform(Xte)
    ytr=np.atleast_2d(ytr); ytr=ytr if ytr.shape[0]==Xtr.shape[0] else ytr.T
    out=np.zeros((Xte.shape[0],ytr.shape[1]),dtype=float)
    for di in range(ytr.shape[1]):
        preds=[]
        for name,mdl in ens_models().items():
            if name=="ridge": mdl.fit(Xtr_s,ytr[:,di]); preds.append(mdl.predict(Xte_s))
            else: mdl.fit(Xtr,ytr[:,di]); preds.append(mdl.predict(Xte))
        out[:,di]=np.mean(preds,axis=0)
    return np.clip(out,0,100) if out.shape[1]>1 else np.clip(out[:,0],0,100)

def predict_knn(bvtr,ytr,bvte,k=5):
    out=np.zeros(len(bvte))
    for i,q in enumerate(bvte):
        if q is None: out[i]=ytr.mean(); continue
        sims=np.array(DataStructs.BulkTanimotoSimilarity(q,bvtr))
        idx=np.argsort(sims)[::-1][:k]; w=sims[idx]
        out[i]=ytr[idx].mean() if w.sum()==0 else float((w*ytr[idx]).sum()/w.sum())
    return np.clip(out,0,100)

def boot_ci(yt,yp,B=2000):
    n=len(yt); rng=np.random.default_rng(SEED); r2s=[]
    for _ in range(B):
        idx=rng.integers(0,n,n)
        if np.var(yt[idx])==0: continue
        r2s.append(r2_score(yt[idx],yp[idx]))
    lo,hi=np.percentile(r2s,[2.5,97.5]); return round(float(lo),3),round(float(hi),3)


def main():
    df=pd.read_csv(DATA)
    for c in DIMS:
        if df[c].max()<=10: df[c]=df[c]*10
        df[c]=df[c].clip(0,100)
    df=df[df["smiles"].notna() & (df["smiles"].astype(str).str.lower()!="nan")]
    df=df[df["data_source"].isin(HUMAN)]                 # quarantine ML pseudo-labels
    df=df[df[DIMS].notna().all(axis=1)].reset_index(drop=True)
    smi=df["smiles"].astype(str).tolist()
    print(f"Human-curated, structure-only, scaffold-CV genuine-predictor test")
    print(f"  n={len(df)} compounds | features = Morgan-{NBITS} + {len(_DESCS)} RDKit descriptors (NO disease counts)")

    Xm, Xd, bv = morgan(smi), descriptors(smi), bvs(smi)
    X = np.hstack([Xm, Xd])
    groups = np.array([scaffold(s) for s in smi]);
    print(f"  unique generic scaffolds: {len(set(groups))}")
    y = df[DIMS].values.astype(float)
    gkf = GroupKFold(5)

    # leakage audit
    cross=[]
    for tr,te in gkf.split(df,groups=groups):
        bt=[bv[i] for i in tr if bv[i] is not None]
        for i in te:
            if bv[i] is not None and bt: cross.append(max(DataStructs.BulkTanimotoSimilarity(bv[i],bt)))
    print(f"  leakage audit test->train Tanimoto: median={np.median(cross):.2f} frac>=0.7={np.mean(np.array(cross)>=0.7):.2f}")

    oof={m:np.zeros_like(y) for m in ["mean","knn","ensemble"]}
    for tr,te in gkf.split(df,groups=groups):
        for di in range(7):
            ytr=y[tr,di]
            oof["mean"][te,di]=ytr.mean()
            oof["knn"][te,di]=predict_knn([bv[i] for i in tr],ytr,[bv[i] for i in te])
        oof["ensemble"][te]=predict_ensemble(X[tr], y[tr], X[te])

    print(f"\n  {'dimension':24} {'ens R2':>7} {'95% CI':>16} {'kNN R2':>7} {'mean':>6}  genuine?")
    results={}; predictable=[]
    for di,d in enumerate(DIMS):
        yt=y[:,di]
        r2e=r2_score(yt,oof["ensemble"][:,di]); lo,hi=boot_ci(yt,oof["ensemble"][:,di])
        r2k=r2_score(yt,oof["knn"][:,di]); r2m=r2_score(yt,oof["mean"][:,di])
        rho=spearmanr(yt,oof["ensemble"][:,di]).correlation
        mae=mean_absolute_error(yt,oof["ensemble"][:,di])
        genuine = (lo>0) and (r2e>r2k) and (r2e>r2m)
        if genuine: predictable.append(d)
        results[d]={"ensemble_r2":round(r2e,3),"ci":[lo,hi],"spearman":round(float(rho),3),
                    "mae":round(mae,2),"knn_r2":round(r2k,3),"mean_r2":round(r2m,3),"genuine":bool(genuine)}
        print(f"  {d:24} {r2e:7.3f} {('['+str(lo)+','+str(hi)+']'):>16} {r2k:7.3f} {r2m:6.3f}  {'YES' if genuine else 'no'}")

    # NPS from ensemble dims
    nps_t=np.array([neuro_score({d:y[i,di] for di,d in enumerate(DIMS)}) for i in range(len(y))])
    nps_p=np.array([neuro_score({d:oof["ensemble"][i,di] for di,d in enumerate(DIMS)}) for i in range(len(y))])
    nps_r2=r2_score(nps_t,nps_p); nlo,nhi=boot_ci(nps_t,nps_p)
    print(f"\n  NPS (aggregate): ensemble R2={nps_r2:.3f} CI=[{nlo},{nhi}] Spearman={spearmanr(nps_t,nps_p).correlation:.3f}")

    # Save genuine final models (trained on ALL human data) for predictable dims
    saved=[]
    if predictable:
        os.makedirs("models_genuine",exist_ok=True)
        sc=StandardScaler().fit(X); joblib.dump(sc,"models_genuine/scaler.joblib")
        for d in predictable:
            di=DIMS.index(d);
            for name,mdl in ens_models().items():
                Xin=sc.transform(X) if name=="ridge" else X
                mdl.fit(Xin,y[:,di]); joblib.dump(mdl,f"models_genuine/{d}_{name}.joblib")
            saved.append(d)
        joblib.dump({"n_bits":NBITS,"descriptors":[n for n,_ in _DESCS]}, "models_genuine/feature_spec.joblib")

    report={
        "design":"structure-only (Morgan-%d + %d descriptors), human-curated labels only, "
                 "scaffold GroupKFold(5), 2000x bootstrap CI; disease-count features EXCLUDED."%(NBITS,len(_DESCS)),
        "n_compounds":len(df),"n_scaffolds":int(len(set(groups))),
        "leakage_median_tanimoto":round(float(np.median(cross)),3),
        "per_dimension":results,
        "nps":{"ensemble_r2":round(nps_r2,3),"ci":[nlo,nhi]},
        "genuinely_predictable_dimensions":predictable,
        "models_saved":saved,
        "verdict":("GENUINE structural signal in: "+", ".join(predictable)) if predictable
                  else "No dimension is genuinely predictable from structure with current labels.",
    }
    json.dump(report,open("BS_predictive_report.json","w"),indent=2)
    print("\nVERDICT:",report["verdict"])
    print("Saved BS_predictive_report.json", "+ models_genuine/" if saved else "")


if __name__=="__main__":
    main()
