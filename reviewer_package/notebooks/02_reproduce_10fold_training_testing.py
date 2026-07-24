# %% [markdown]
# # Reproduce: Random Forest + 10-fold cross-validation (training and testing)
#
# This notebook reproduces the headline classification and regression numbers **from the raw measured
# data**, with every step of the logic written inline (no project-specific imports), so it can be
# verified independently. Fixed `random_state = 42` throughout.
#
# **Environment:** Python 3.13, RDKit 2026.03, scikit-learn 1.8, NumPy, pandas.
#
# **Inputs:** `../data/endpoints/<ENDPOINT>.csv` (columns: `smiles`, `label`, `pchembl`, `year`,
# `source`; the antioxidant file uses `y`). These are the measured, pooled ChEMBL + BindingDB (targets),
# B3DB (BBB) and ChEMBL DPPH (antioxidant) sets.

# %%
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import StratifiedKFold, KFold, GroupKFold
from sklearn.metrics import roc_auc_score, matthews_corrcoef, r2_score
from scipy.stats import spearmanr

RDLogger.DisableLog("rdApp.*")
SEED, N_SPLITS = 42, 10

ROOT = Path.cwd()
while not (ROOT / "data" / "endpoints").exists() and ROOT != ROOT.parent:
    ROOT = ROOT.parent
DATA = ROOT / "data" / "endpoints"
print("data directory:", DATA)

# %% [markdown]
# ## 1. Molecular representation: 1024-bit ECFP-4 + 12 descriptors = 1036 features
# A SMILES is reduced to its largest organic fragment (salt stripping), then encoded as a 1024-bit
# Morgan (ECFP-4, radius 2) fingerprint concatenated with twelve interpretable RDKit descriptors. The
# same 1036-feature vector is used for every endpoint.

# %%
MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
DESCRIPTORS = {
    "mw": Descriptors.MolWt, "clogp": Crippen.MolLogP, "tpsa": rdMolDescriptors.CalcTPSA,
    "hbd": rdMolDescriptors.CalcNumHBD, "hba": rdMolDescriptors.CalcNumHBA,
    "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds,
    "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings,
    "fraction_csp3": rdMolDescriptors.CalcFractionCSP3, "ring_count": rdMolDescriptors.CalcNumRings,
    "heavy_atoms": lambda m: m.GetNumHeavyAtoms(), "formal_charge": Chem.GetFormalCharge,
    "qed": QED.qed,
}


def parent(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def featurize_one(smiles):
    m = parent(smiles)
    if m is None:
        return None
    fp = MORGAN.GetFingerprintAsNumPy(m).astype(np.float32)
    desc = np.array([fn(m) for fn in DESCRIPTORS.values()], dtype=np.float32)
    return np.concatenate([fp, desc])


def featurize(smiles_list):
    rows, mask = [], np.zeros(len(smiles_list), dtype=bool)
    for i, s in enumerate(smiles_list):
        v = featurize_one(s)
        if v is not None:
            rows.append(v)
            mask[i] = True
    return np.vstack(rows), mask


def scaffold_groups(smiles_list):
    codes, mapping = [], {}
    for s in smiles_list:
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(str(s)), includeChirality=False)
        except Exception:
            scaf = ""
        codes.append(mapping.setdefault(scaf or f"_none_{len(mapping)}", len(mapping)))
    return np.asarray(codes)


print("features:", 1024 + len(DESCRIPTORS), "(1024 ECFP-4 bits +", len(DESCRIPTORS), "descriptors)")

# %% [markdown]
# ## 2. Train and test one endpoint under 10-fold cross-validation
# For classification: RandomForest (300 trees, min_samples_leaf 2, balanced weights). Random
# StratifiedKFold(10) and scaffold GroupKFold(10) on Bemis-Murcko scaffolds. Metric: AUROC and MCC.
# For regression: RandomForestRegressor; KFold / GroupKFold; R-squared and Spearman.

# %%
CLASSIFICATION = {"BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_A", "MAO_B", "hERG"}


def run_endpoint(name):
    task = "classification" if name in CLASSIFICATION else "regression"
    df = pd.read_csv(DATA / f"{name}.csv")
    target = "label" if task == "classification" else ("y" if name == "antioxidant_dpph" else "pchembl")
    df = df.dropna(subset=["smiles", target]).reset_index(drop=True)
    X, mask = featurize(df["smiles"].tolist())
    df = df.loc[mask].reset_index(drop=True)
    y = df[target].to_numpy()
    y = y.astype(int) if task == "classification" else y.astype(float)
    g = scaffold_groups(df["smiles"].tolist())
    out = {"endpoint": name, "task": task, "n_compounds": len(y), "n_scaffolds": len(set(g))}
    for split, sp in [("random", StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED)
                       if task == "classification" else KFold(N_SPLITS, shuffle=True, random_state=SEED)),
                      ("scaffold", GroupKFold(N_SPLITS))]:
        it = sp.split(X, y, g) if split == "scaffold" else sp.split(X, y)
        scores = []
        for tr, te in it:
            if task == "classification":
                mdl = RandomForestClassifier(n_estimators=300, min_samples_leaf=2,
                                             class_weight="balanced", n_jobs=-1, random_state=SEED)
                mdl.fit(X[tr], y[tr])
                scores.append(roc_auc_score(y[te], mdl.predict_proba(X[te])[:, 1]))
            else:
                mdl = RandomForestRegressor(n_estimators=300, min_samples_leaf=2, n_jobs=-1,
                                            random_state=SEED)
                mdl.fit(X[tr], y[tr])
                scores.append(r2_score(y[te], mdl.predict(X[te])))
        out[f"{split}"] = round(float(np.mean(scores)), 3)
        out[f"{split}_sd"] = round(float(np.std(scores, ddof=1)), 3)
    return out


# MAO_A is the smallest classifier: fast to run and reproduces the published scaffold AUROC 0.868.
result = run_endpoint("MAO_A")
print(result)
assert abs(result["scaffold"] - 0.868) < 0.005, "reproduction mismatch"
print("\nReproduced MAO_A scaffold AUROC =", result["scaffold"], "(published: 0.868)  -> MATCH")

# %% [markdown]
# ## 3. Cross-check against the published summary
# The full published numbers for every endpoint are in `../results/rf_cv_summary.csv`. To reproduce all
# thirteen, call `run_endpoint(name)` for each (a few minutes total on CPU).

# %%
summary = pd.read_csv(ROOT / "results" / "rf_cv_summary.csv")
print(summary[summary.split == "scaffold"][["endpoint", "task", "n", "roc_auc_mean", "r2_mean"]]
      .to_string(index=False))
