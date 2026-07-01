"""
BS_antioxidant_select.py — Tier-2 model/representation selection for the ONE
genuinely structure-predictable endpoint (antioxidant). Honest scaffold-CV,
bootstrap CI, multiple representations and model families. Saves the winner.
"""
import os, sys, json, warnings
os.environ.setdefault("HF_HUB_OFFLINE","1"); os.environ.setdefault("TRANSFORMERS_OFFLINE","1")
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
if hasattr(sys.stdout,"reconfigure"): sys.stdout.reconfigure(encoding="utf-8",errors="replace")
import numpy as np, pandas as pd, joblib
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, GradientBoostingRegressor
from sklearn.svm import SVR
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error
import BS_predictive_model as P   # reuse featurizers/scaffold

SEED=42; DIM="antioxidant"
df=pd.read_csv(P.DATA)
if df[DIM].max()<=10: df[DIM]=df[DIM]*10
df[DIM]=df[DIM].clip(0,100)
df=df[df["smiles"].notna() & (df["smiles"].astype(str).str.lower()!="nan")]
df=df[df["data_source"].isin(P.HUMAN)]
df=df[df[DIM].notna()].reset_index(drop=True)
smi=df["smiles"].astype(str).tolist(); y=df[DIM].values.astype(float)
groups=np.array([P.scaffold(s) for s in smi])
print(f"antioxidant | n={len(df)} | scaffolds={len(set(groups))} | label unique vals={len(np.unique(y))}")

Xm=P.morgan(smi); Xd=P.descriptors(smi)
def chemberta(smiles,batch=32):
    from transformers import AutoModel,AutoTokenizer; import torch
    tok=AutoTokenizer.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
    mdl=AutoModel.from_pretrained("seyonec/ChemBERTa-zinc-base-v1"); mdl.eval(); E=[]
    for i in range(0,len(smiles),batch):
        ch=[s if s else "C" for s in smiles[i:i+batch]]
        enc=tok(ch,padding=True,truncation=True,max_length=128,return_tensors="pt")
        with torch.no_grad(): out=mdl(**enc)
        E.append(out.last_hidden_state.mean(1).numpy().astype(np.float32))
    return np.vstack(E)
Xcb=chemberta(smi)

FS={"descriptors":Xd,"morgan":Xm,"morgan+desc":np.hstack([Xm,Xd]),
    "morgan+desc+chemberta":np.hstack([Xm,Xd,Xcb])}

def model_set():
    return {"ridge":Ridge(alpha=10.0),
            "rf":RandomForestRegressor(400,max_depth=10,min_samples_leaf=2,random_state=SEED,n_jobs=-1),
            "et":ExtraTreesRegressor(400,max_depth=10,min_samples_leaf=2,random_state=SEED,n_jobs=-1),
            "gbt":GradientBoostingRegressor(n_estimators=200,learning_rate=0.05,max_depth=3,random_state=SEED),
            "svr":SVR(C=10,gamma="scale")}

def boot(yt,yp,B=2000):
    rng=np.random.default_rng(SEED); n=len(yt); r=[]
    for _ in range(B):
        idx=rng.integers(0,n,n)
        if np.var(yt[idx])>0: r.append(r2_score(yt[idx],yp[idx]))
    return round(float(np.percentile(r,2.5)),3),round(float(np.percentile(r,97.5)),3)

def cv_eval(X, use_ensemble=True, single=None):
    gkf=GroupKFold(5); oof=np.zeros_like(y)
    for tr,te in gkf.split(X,groups=groups):
        sc=StandardScaler().fit(X[tr]); Xtr_s,Xte_s=sc.transform(X[tr]),sc.transform(X[te])
        preds=[]
        for name,mdl in model_set().items():
            if single and name!=single: continue
            Xtr_in,Xte_in=(Xtr_s,Xte_s) if name in("ridge","svr") else (X[tr],X[te])
            mdl.fit(Xtr_in,y[tr]); preds.append(mdl.predict(Xte_in))
        oof[te]=np.clip(np.mean(preds,axis=0),0,100)
    return r2_score(y,oof), boot(y,oof), spearmanr(y,oof).correlation, mean_absolute_error(y,oof)

print(f"\n{'representation':24} {'model':10} {'R2':>7} {'95% CI':>16} {'rho':>6} {'MAE':>6}")
best=None
for fsname,X in FS.items():
    r2,ci,rho,mae=cv_eval(X,use_ensemble=True)
    print(f"{fsname:24} {'ensemble':10} {r2:7.3f} {('['+str(ci[0])+','+str(ci[1])+']'):>16} {rho:6.3f} {mae:6.2f}")
    if best is None or r2>best[1]: best=(fsname,r2,ci,rho,mae,X)
# also per-model on the best representation
print("\nper-model on best representation (%s):"%best[0])
for m in model_set():
    r2,ci,rho,mae=cv_eval(best[5],single=m)
    print(f"  {m:10} R2={r2:.3f} CI=[{ci[0]},{ci[1]}] rho={rho:.3f}")

print(f"\nBEST: {best[0]} ensemble R2={best[1]:.3f} CI={best[2]} Spearman={best[3]:.3f} MAE={best[4]:.2f}")
json.dump({"endpoint":DIM,"n":len(df),"label_unique_values":int(len(np.unique(y))),
           "best_representation":best[0],"r2":round(best[1],3),"ci":best[2],
           "spearman":round(float(best[3]),3),"mae":round(best[4],2),
           "note":"Honest scaffold-CV; coarse curated labels (few unique values) cap achievable R2; "
                  "continuous assay endpoints (ORAC/DPPH/TEAC) would improve this."},
          open("BS_antioxidant_report.json","w"),indent=2)
print("Saved BS_antioxidant_report.json")
