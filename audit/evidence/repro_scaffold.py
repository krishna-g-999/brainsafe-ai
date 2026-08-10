"""BBB scaffold 10-fold CV: reproduce as published, then repeat with duplicates collapsed."""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

RDLogger.DisableLog("rdApp.*")
ROOT = Path(r"D:\BRAINSAFE_AI")
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize  # noqa: E402

SEED, N_SPLITS = 42, 10
RF = dict(n_estimators=300, min_samples_leaf=2, n_jobs=-1, random_state=SEED)


def scaffold_groups(smiles):
    codes, mapping = [], {}
    for smi in smiles:
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(str(smi)),
                                                       includeChirality=False)
        except Exception:
            scaf = ""
        codes.append(mapping.setdefault(scaf or f"_none_{len(mapping)}", len(mapping)))
    return np.asarray(codes)


def cv(X, y, g):
    aucs = []
    for tr, te in GroupKFold(N_SPLITS).split(X, y, g):
        if len(set(y[te])) < 2:
            continue
        m = RandomForestClassifier(class_weight="balanced", **RF)
        m.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs)), len(aucs)


df = pd.read_csv(ROOT / "data" / "endpoints" / "BBB.csv")
smis = df["smiles"].astype(str).tolist()
X, mask = featurize(smis)
y = df.loc[mask, "label"].to_numpy().astype(int)
kept_smis = [s for s, m in zip(smis, mask) if m]
g = scaffold_groups(kept_smis)

m1, s1, k1 = cv(X, y, g)
print(f"scaffold CV as published : n={len(y)} groups={len(set(g))} "
      f"AUROC={m1:.4f} sd={s1:.4f} folds={k1}")

groups = defaultdict(list)
for i, vec in enumerate(X):
    groups[vec.tobytes()].append(i)
keep = [idxs[0] for idxs in groups.values() if len({y[i] for i in idxs}) == 1]
keep = np.array(sorted(keep))
m2, s2, k2 = cv(X[keep], y[keep], g[keep])
print(f"scaffold CV deduplicated : n={len(keep)} groups={len(set(g[keep]))} "
      f"AUROC={m2:.4f} sd={s2:.4f} folds={k2}")
print(f"difference               : {m1 - m2:+.4f}")

sizes = pd.Series(g).value_counts()
print(f"\nlargest scaffold groups: {sizes.head(5).tolist()}")
print(f"singleton scaffold groups: {(sizes == 1).sum()} of {len(sizes)}")
acyclic = sum(1 for s in kept_smis
              if MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(s),
                                                     includeChirality=False) == "")
print(f"compounds with an empty (acyclic) scaffold, each given its own group: {acyclic}")
