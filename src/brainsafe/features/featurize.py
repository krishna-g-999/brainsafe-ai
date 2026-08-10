"""Turn a molecular structure into the numeric feature vector the models train on.

Every model input is numeric. A SMILES string is not fed to the model directly; it is converted
into two blocks of numbers:

  1. A 1024-bit ECFP-4 fingerprint (Morgan, radius 2). Each bit is 0/1 and records the presence of
     a particular local atomic environment. This is a collision-free-by-construction numeric
     encoding of substructure: bit k always means the same environment for every compound.
  2. Twelve interpretable physicochemical descriptors (molecular weight, logP, polar surface area,
     H-bond donors/acceptors, rotatable bonds, aromatic rings, sp3 fraction, ring count, heavy-atom
     count, formal charge, drug-likeness QED).

The result is a fixed 1036-length numeric vector with stable, named columns, identical for every
endpoint, so nothing categorical or free-text ever reaches the estimator.

Any categorical metadata that is *not* a model feature (endpoint name, assay source) is encoded
separately and reversibly in encodings.py; it is never mixed into this feature matrix.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

MORGAN_RADIUS = 2
MORGAN_BITS = 1024
_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=MORGAN_RADIUS, fpSize=MORGAN_BITS)

# Ordered so the column layout is stable and reproducible across runs.
_DESCRIPTORS = {
    "mw": Descriptors.MolWt,
    "clogp": Crippen.MolLogP,
    "tpsa": rdMolDescriptors.CalcTPSA,
    "hbd": rdMolDescriptors.CalcNumHBD,
    "hba": rdMolDescriptors.CalcNumHBA,
    "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds,
    "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings,
    "fraction_csp3": rdMolDescriptors.CalcFractionCSP3,
    "ring_count": rdMolDescriptors.CalcNumRings,
    "heavy_atoms": lambda m: m.GetNumHeavyAtoms(),
    "formal_charge": Chem.GetFormalCharge,
    "qed": QED.qed,
}

N_FEATURES = MORGAN_BITS + len(_DESCRIPTORS)


def feature_names() -> list[str]:
    """Stable column names for the full feature vector (fingerprint bits then descriptors)."""
    return [f"ecfp4_{i}" for i in range(MORGAN_BITS)] + list(_DESCRIPTORS)


def parent_mol(smiles: str):
    """The largest fragment of a structure, sanitised: the molecule the model actually sees.

    Public because the scaffold grouping in training must desalt exactly as this does, or a salt
    and its free base become identical inputs assigned to different cross-validation folds.
    """
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) > 1:  # keep the largest organic fragment (salt / counter-ion removal)
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def featurize_one(smiles: str) -> np.ndarray | None:
    """Return the numeric feature vector for one SMILES, or None if it cannot be parsed."""
    mol = parent_mol(smiles)
    if mol is None:
        return None
    fp = _MORGAN.GetFingerprintAsNumPy(mol).astype(np.float32)
    desc = np.array([fn(mol) for fn in _DESCRIPTORS.values()], dtype=np.float32)
    return np.concatenate([fp, desc])


def featurize(smiles_list) -> tuple[np.ndarray, np.ndarray]:
    """Featurize many SMILES.

    Returns (X, mask) where X has one row per *valid* structure and mask is a boolean array over the
    input marking which entries were parsed. The caller uses mask to align labels with rows.
    """
    rows, mask = [], np.zeros(len(smiles_list), dtype=bool)
    for i, smi in enumerate(smiles_list):
        v = featurize_one(smi)
        if v is not None:
            rows.append(v)
            mask[i] = True
    X = np.vstack(rows) if rows else np.empty((0, N_FEATURES), dtype=np.float32)
    return X, mask
