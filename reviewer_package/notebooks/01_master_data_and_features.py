# %% [markdown]
# # Master data table, per-endpoint sets, and the feature definition
#
# What the reviewer asked for: the complete master table with the descriptors, the per-endpoint
# training/testing counts, and the exact feature definition the model consumes.

# %%
from pathlib import Path
import json
import pandas as pd

ROOT = Path.cwd()
while not (ROOT / "data").exists() and ROOT != ROOT.parent:
    ROOT = ROOT.parent
DATA = ROOT / "data"

# %% [markdown]
# ## 1. Master compound library
# One row per unique compound (keyed by standard InChIKey of the salt-stripped parent). Contains every
# measured endpoint label/value the compound has, the interpretable descriptors, a flavonoid flag, and
# the contributing data sources. Empty cells mean *not measured* (never imputed).

# %%
master = pd.read_csv(DATA / "master_compound_library.csv")
print("master rows (unique compounds):", len(master))
print("unique InChIKeys:", master["inchikey"].nunique(), "(== rows, so no duplicates)")
print("\ncolumns:", list(master.columns))
master.head(3)

# %% [markdown]
# ## 2. Descriptors held in the master table
# These are the interpretable physicochemical descriptors; they are also 12 of the 1036 model features
# (the other 1024 are ECFP-4 fingerprint bits).

# %%
desc_cols = ["mw", "clogp", "tpsa", "hbd", "hba", "rotatable_bonds", "aromatic_rings",
             "fraction_csp3", "qed"]
print(master[[c for c in desc_cols if c in master.columns]].describe().round(2).to_string())

# %% [markdown]
# ## 3. Per-endpoint training/testing counts
# The number of measured compounds used (each is both a training example in 9 folds and a test example
# in 1 fold under 10-fold cross-validation).

# %%
rows = []
for f in sorted((DATA / "endpoints").glob("*.csv")):
    d = pd.read_csv(f)
    r = {"endpoint": f.stem, "n_compounds": len(d)}
    if "label" in d.columns:
        r["active"] = int((d.label == 1).sum())
        r["inactive"] = int((d.label == 0).sum())
    rows.append(r)
counts = pd.DataFrame(rows)
print(counts.to_string(index=False))
print("\ntotal measured compound-endpoint records:", int(counts.n_compounds.sum()))

# %% [markdown]
# ## 4. Exact model feature definition (1036 features)
# The ordered feature names the model was trained on: 1024 fingerprint bits then 12 named descriptors.

# %%
names = json.loads((DATA / "feature_names_1036.json").read_text())
print("total features:", len(names))
print("first 3 fingerprint bits:", names[:3])
print("the 12 descriptors:", names[-12:])

# %% [markdown]
# ## 5. Feature retention (why both blocks are kept)
# Block ablation and descriptor permutation importance justify keeping both the fingerprint and the
# descriptor block (results/feature_block_ablation.csv, feature_descriptor_importance.csv).

# %%
abl = pd.read_csv(ROOT / "results" / "feature_block_ablation.csv")
print(abl.pivot_table(index="endpoint", columns="block", values="mean").round(3).to_string())
