# %% [markdown]
# # Validation: calibration, conformal coverage, temporal, and adversarial (inversion) checks
#
# The published validation tables, plus two live checks (no model needed) the reviewer can run to
# confirm there is no leakage and no duplication.

# %%
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path.cwd()
while not (ROOT / "results").exists() and ROOT != ROOT.parent:
    ROOT = ROOT.parent
R = ROOT / "results"

# %% [markdown]
# ## 1. Probability calibration (isotonic) — expected calibration error, raw vs calibrated
# %%
print(pd.read_csv(R / "calibration.csv").to_string(index=False))

# %% [markdown]
# ## 2. Conformal prediction — empirical coverage vs the 0.90 target
# %%
print(pd.read_csv(R / "rf_conformal.csv").to_string(index=False))

# %% [markdown]
# ## 3. Temporal validation — train on the past, test on genuinely later compounds
# %%
print(pd.read_csv(R / "rf_temporal.csv").to_string(index=False))

# %% [markdown]
# ## 4. Adversarial (inversion) validation — the six checks
# %%
print(pd.read_csv(R / "inversion_validation.csv").to_string(index=False))

# %% [markdown]
# ## 5. Live re-check 1: no scaffold leakage under the scaffold split
# GroupKFold on Bemis-Murcko scaffolds; verify no scaffold is shared between any train and test fold.
# %%
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold
RDLogger.DisableLog("rdApp.*")

df = pd.read_csv(ROOT / "data" / "endpoints" / "MAO_A.csv").dropna(subset=["smiles", "label"])
scaf, mp = [], {}
for s in df.smiles:
    k = MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(str(s)), includeChirality=False) or "_"
    scaf.append(mp.setdefault(k, len(mp)))
scaf = np.array(scaf)
bad = sum(bool(set(scaf[tr]) & set(scaf[te]))
          for tr, te in GroupKFold(10).split(np.zeros((len(scaf), 1)), df.label.values, scaf))
print(f"folds with a shared scaffold: {bad} (want 0)  ->", "PASS" if bad == 0 else "FAIL")

# %% [markdown]
# ## 6. Live re-check 2: no duplicate compounds in the master table
# %%
m = pd.read_csv(ROOT / "data" / "master_compound_library.csv")
dup = len(m) - m["inchikey"].nunique()
print(f"master rows {len(m):,}, unique InChIKeys {m['inchikey'].nunique():,}, duplicates {dup}  ->",
      "PASS" if dup == 0 else "FAIL")
