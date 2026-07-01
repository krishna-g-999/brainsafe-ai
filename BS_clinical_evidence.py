"""
BS_clinical_evidence.py — translational-evidence layer addressing 'engagement != efficacy'.
For a query SMILES, return the nearest clinically-advanced CNS compounds (ChEMBL ATC-N,
max_phase >= 1) with their phase and disease class. This is measured CLINICAL PRECEDENT,
not an efficacy prediction: it tells the user whether a molecule resembles compounds that
actually reached human CNS trials/approval.
"""
import os
import numpy as np
_DIR = os.path.dirname(os.path.abspath(__file__))
_REF, _FPS = None, None


def _load():
    global _REF, _FPS
    if _REF is not None:
        return
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import AllChem
    path = os.path.join(_DIR, "data", "clinical_cns_reference.csv")
    if not os.path.exists(path):
        _REF, _FPS = [], []
        return
    df = pd.read_csv(path)
    ref, fps = [], []
    for _, r in df.iterrows():
        m = Chem.MolFromSmiles(str(r["smiles"]))
        if m is None:
            continue
        fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024))
        ref.append({"name": r.get("name"), "max_phase": r.get("max_phase"),
                    "disease": r.get("disease")})
    _REF, _FPS = ref, fps


def clinical_analogs(smiles: str | None, k: int = 3) -> list:
    _load()
    if not smiles or not _FPS:
        return []
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return []
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024)
    sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, _FPS))
    out = []
    for i in np.argsort(sims)[::-1][:k]:
        r = _REF[i]
        out.append({"name": r["name"], "max_phase": r["max_phase"],
                    "disease": r["disease"], "similarity": round(float(sims[i]), 2)})
    return out
