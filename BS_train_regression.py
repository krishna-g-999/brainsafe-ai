"""
BS_train_regression.py — fix the receptor endpoints (D2, A2A, 5-HT2A, SERT) properly.

Binary active/inactive QSAR was ill-posed (96-98% actives). The scientifically correct
frame is POTENCY REGRESSION on measured pChEMBL (-log10 IC50/Ki, molar). We train a
RandomForest+ExtraTrees+HistGB regression ensemble and validate honestly:
  * scaffold GroupKFold(5) out-of-fold R2 / RMSE / Spearman
  * temporal split (train<=75th-pct year, test newest ~25%) R2
Models saved to models_brain_reg/<ep>.joblib (+ meta).
"""
import os, glob, json, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, joblib
from rdkit import Chem
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_squared_error
from BS_predictive_model import morgan, descriptors, scaffold

OUT = "models_brain_reg"; os.makedirs(OUT, exist_ok=True)
SEED = 42
TARGETS = {"D2": "Parkinson's / psychosis (dopamine D2) — pKi/pIC50",
           "A2A": "Parkinson's (adenosine A2A) — pKi/pIC50",
           "HT2A": "mood / psychosis (5-HT2A) — pKi/pIC50",
           "SERT": "depression (serotonin transporter) — pKi/pIC50"}


def canon(df):
    rows = []
    for _, r in df.iterrows():
        m = Chem.MolFromSmiles(str(r["smiles"]))
        if m is None or pd.isna(r.get("pchembl")):
            continue
        try: ik = Chem.MolToInchiKey(m)
        except Exception: ik = Chem.MolToSmiles(m)
        yr = float(r["year"]) if ("year" in df.columns and pd.notna(r["year"])) else None
        rows.append((Chem.MolToSmiles(m), float(r["pchembl"]), ik, yr))
    return pd.DataFrame(rows, columns=["smiles", "y", "ik", "year"]).drop_duplicates("ik").reset_index(drop=True)


def reg_ensemble():
    return [RandomForestRegressor(n_estimators=300, min_samples_leaf=2, n_jobs=-1, random_state=SEED),
            ExtraTreesRegressor(n_estimators=300, min_samples_leaf=2, n_jobs=-1, random_state=SEED),
            HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06, random_state=SEED)]


def predict(models, Xtr, ytr, Xte):
    return np.mean([m.fit(Xtr, ytr).predict(Xte) for m in models], axis=0)


def evaluate(name, df):
    d = canon(df)
    if len(d) < 150:
        return None
    smi = d["smiles"].tolist(); y = d["y"].values
    X = np.hstack([morgan(smi), descriptors(smi)])
    scaf = np.array([scaffold(s) for s in smi])

    # scaffold-CV OOF regression
    oof = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, groups=scaf):
        oof[te] = predict(reg_ensemble(), X[tr], y[tr], X[te])
    r2 = r2_score(y, oof); rmse = mean_squared_error(y, oof) ** 0.5
    rho = spearmanr(y, oof).correlation

    # temporal
    temporal = None
    if d["year"].notna().sum() >= 0.6 * len(d) and d["year"].nunique() > 3:
        cut = int(d["year"].quantile(0.75))
        te = (d["year"] > cut).values; tr = (~te) & d["year"].notna().values
        if te.sum() >= 40:
            p = predict(reg_ensemble(), X[tr], y[tr], X[te])
            temporal = {"cutoff_year": cut, "n_test": int(te.sum()),
                        "r2": round(float(r2_score(y[te], p)), 3),
                        "rmse": round(float(mean_squared_error(y[te], p) ** 0.5), 3)}

    # final model on all data
    finals = reg_ensemble(); [m.fit(X, y) for m in finals]
    joblib.dump({"models": finals, "n_bits": 1024}, f"{OUT}/{name}.joblib")
    meta = {"endpoint": name, "meaning": TARGETS.get(name, name), "task": "regression (pChEMBL)",
            "n": len(d), "y_min": round(float(y.min()), 2), "y_max": round(float(y.max()), 2),
            "y_std": round(float(y.std()), 2), "scaffold_cv_r2": round(float(r2), 3),
            "rmse": round(float(rmse), 3), "spearman": round(float(rho), 3), "temporal": temporal,
            "source": "ChEMBL measured pChEMBL"}
    json.dump(meta, open(f"{OUT}/{name}_meta.json", "w"), indent=2)
    print(f"  [{name:5}] n={len(d):5d} scaffold-CV R2={r2:.3f} RMSE={rmse:.2f} rho={rho:.3f} "
          f"| temporal R2={temporal['r2'] if temporal else '-'}")
    return meta


def main():
    rep = {}
    for name in TARGETS:
        f = f"data/endpoints/{name}.csv"
        if os.path.exists(f):
            m = evaluate(name, pd.read_csv(f))
            if m: rep[name] = m
    json.dump(rep, open(f"{OUT}/regression_report.json", "w"), indent=2)
    print("Saved", OUT, "regression models + report.")


if __name__ == "__main__":
    main()
