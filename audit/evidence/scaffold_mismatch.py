"""Do feature-identical compounds always land in the same scaffold group?

train_rf.py computes the scaffold from the RAW SMILES (Chem.MolFromSmiles at line 66), while
featurize.py strips salts to the largest fragment before featurising (line 61). If a salt and its
free base appear as separate rows, they are identical to the model but may be assigned to
different GroupKFold groups, in which case GroupKFold does not separate them and they leak.
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
ROOT = Path(r"D:\BRAINSAFE_AI")
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize  # noqa: E402


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


for ep in ["BBB", "AChE", "BACE1", "hERG", "D2"]:
    df = pd.read_csv(ROOT / "data" / "endpoints" / f"{ep}.csv")
    smis = df["smiles"].astype(str).tolist()
    X, mask = featurize(smis)
    kept = [s for s, m in zip(smis, mask) if m]
    g = scaffold_groups(kept)

    feat = defaultdict(list)
    for i, vec in enumerate(X):
        feat[vec.tobytes()].append(i)

    split_groups = 0
    split_rows = 0
    for idxs in feat.values():
        if len(idxs) > 1 and len({g[i] for i in idxs}) > 1:
            split_groups += 1
            split_rows += len(idxs)
    print(f"{ep:8s} feature-identical clusters that GroupKFold can place in DIFFERENT folds: "
          f"{split_groups:5d} clusters, {split_rows:5d} rows", flush=True)
