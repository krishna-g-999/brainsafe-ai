"""
Train a MEASURED antioxidant model (DPPH pIC50) to replace the weak curated R2~0.25 model,
and CROSS-CHECK the old curated model against the measured data.
"""
import os, json, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, joblib
from rdkit import Chem
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_squared_error
from BS_predictive_model import morgan, descriptors, scaffold
SEED = 42

df = pd.read_csv("data/endpoints_reg/antioxidant_dpph.csv")
rows = []
for _, r in df.iterrows():
    m = Chem.MolFromSmiles(str(r["smiles"]))
    if m is None: continue
    try: ik = Chem.MolToInchiKey(m)
    except Exception: ik = Chem.MolToSmiles(m)
    yr = float(r["year"]) if pd.notna(r.get("year")) else None
    rows.append((Chem.MolToSmiles(m), float(r["y"]), ik, yr))
d = pd.DataFrame(rows, columns=["smiles", "y", "ik", "year"]).drop_duplicates("ik").reset_index(drop=True)
smi = d["smiles"].tolist(); y = d["y"].values
Xd = descriptors(smi); X = np.hstack([morgan(smi), Xd])
scaf = np.array([scaffold(s) for s in smi])
print(f"Measured DPPH antioxidant: n={len(d)} pIC50 {y.min():.1f}-{y.max():.1f} (std {y.std():.2f})")

def ens(): return [RandomForestRegressor(300, min_samples_leaf=2, n_jobs=-1, random_state=SEED),
                   ExtraTreesRegressor(300, min_samples_leaf=2, n_jobs=-1, random_state=SEED),
                   HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06, random_state=SEED)]
def pred(Xtr, ytr, Xte): return np.mean([m.fit(Xtr, ytr).predict(Xte) for m in ens()], axis=0)

# scaffold-CV (morgan+desc) and descriptor-only
oof = np.zeros(len(y)); oofd = np.zeros(len(y))
for tr, te in GroupKFold(5).split(X, groups=scaf):
    oof[te] = pred(X[tr], y[tr], X[te])
    oofd[te] = pred(Xd[tr], y[tr], Xd[te])
r2 = r2_score(y, oof); rmse = mean_squared_error(y, oof) ** 0.5; rho = spearmanr(y, oof).correlation
r2d = r2_score(y, oofd)
print(f"  scaffold-CV (morgan+desc) R2={r2:.3f} RMSE={rmse:.2f} rho={rho:.3f} | (desc-only R2={r2d:.3f})")

# temporal
temporal = None
if d["year"].notna().sum() >= 0.5 * len(d) and d["year"].nunique() > 3:
    cut = int(d["year"].quantile(0.75)); te = (d["year"] > cut).values; tr = (~te) & d["year"].notna().values
    if te.sum() >= 30:
        p = pred(X[tr], y[tr], X[te])
        temporal = {"cutoff_year": cut, "n_test": int(te.sum()),
                    "r2": round(float(r2_score(y[te], p)), 3), "rmse": round(float(mean_squared_error(y[te], p) ** 0.5), 3)}
        print(f"  temporal: train<= {cut}, n_test={te.sum()} R2={temporal['r2']} RMSE={temporal['rmse']}")

# CROSS-CHECK old curated model vs measured DPPH
cc = None
try:
    old = joblib.load("models_genuine/antioxidant_genuine_ridge.joblib")
    old_pred = old.predict(Xd)             # curated 0-100 score
    cc = round(float(spearmanr(old_pred, y).correlation), 3)
    print(f"  CROSS-CHECK: old curated 0-100 score vs measured DPPH pIC50 -> Spearman rho={cc}")
except Exception as e:
    print("  cross-check failed:", e)

# save measured model (morgan+desc ensemble)
os.makedirs("models_genuine", exist_ok=True)
finals = ens(); [m.fit(X, y) for m in finals]
joblib.dump({"models": finals, "n_bits": 1024, "task": "regression_pIC50"},
            "models_genuine/antioxidant_measured_dpph.joblib")
meta = {"endpoint": "antioxidant_DPPH_measured", "n": len(d), "task": "regression (DPPH pIC50)",
        "scaffold_cv_r2": round(float(r2), 3), "rmse": round(float(rmse), 3),
        "spearman": round(float(rho), 3), "temporal": temporal,
        "crosscheck_curated_vs_measured_spearman": cc,
        "source": "ChEMBL DPPH radical-scavenging assays (measured IC50/EC50 -> pIC50)",
        "vs_old_curated_r2": 0.25}
json.dump(meta, open("models_genuine/antioxidant_measured_meta.json", "w"), indent=2)
print("Saved models_genuine/antioxidant_measured_dpph.joblib + meta")
