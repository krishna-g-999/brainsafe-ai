"""Independent check: is the 306-compound 'external' BBB set really disjoint from BBB training?

Recomputes the overlap three ways:
  (a) full InChIKey        - the criterion the repository itself uses
  (b) InChIKey skeleton    - first 14 chars, i.e. constitution ignoring stereo/protonation/isotope
  (c) Tanimoto == 1.0      - identical ECFP-4 1024-bit fingerprint (what the model actually sees)
"""
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
ROOT = Path(r"D:\BRAINSAFE_AI")
sys.path.insert(0, str(ROOT / "src" / "brainsafe" / "data"))
from build_compound_library import standardise  # noqa: E402

train = pd.read_csv(ROOT / "data" / "endpoints" / "BBB.csv")
ext = pd.read_csv(ROOT / "data" / "external" / "processed" / "external_bbb_test.csv")
novel = ext[~ext["in_b3db_training"]].reset_index(drop=True)
print(f"BBB training rows: {len(train)}")
print(f"external rows: {len(ext)}  flagged novel: {len(novel)}")

train_full, train_skel, train_mols = set(), set(), []
for smi in train["smiles"]:
    _, ik = standardise(smi)
    if ik:
        train_full.add(ik)
        train_skel.add(ik.split("-")[0])
    m = Chem.MolFromSmiles(smi)
    if m is not None:
        train_mols.append(m)

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
train_fps = [gen.GetFingerprint(m) for m in train_mols]
print(f"unique full InChIKeys in training: {len(train_full)}")
print(f"unique InChIKey skeletons in training: {len(train_skel)}")
print(f"-> collapsed by stereo/salt/protonation: {len(train_full) - len(train_skel)} training rows")

hits_full, hits_skel, hits_tan = [], [], []
for r in novel.itertuples(index=False):
    _, ik = standardise(r.canonical_smiles)
    if ik is None:
        continue
    if ik in train_full:
        hits_full.append((r.name, ik))
    if ik.split("-")[0] in train_skel:
        hits_skel.append((r.name, ik))
    m = Chem.MolFromSmiles(r.canonical_smiles)
    if m is not None and train_fps:
        sims = DataStructs.BulkTanimotoSimilarity(gen.GetFingerprint(m), train_fps)
        mx = max(sims)
        if mx >= 0.9999:
            hits_tan.append((r.name, round(mx, 4)))

print()
print(f"(a) full-InChIKey overlap with training : {len(hits_full)} / {len(novel)}")
print(f"(b) skeleton-InChIKey overlap           : {len(hits_skel)} / {len(novel)}")
print(f"(c) Tanimoto==1.0 to a training compound: {len(hits_tan)} / {len(novel)}")
print()
if hits_skel:
    print("Compounds in the 'external' set whose constitution IS in training (skeleton match):")
    for n, ik in hits_skel:
        print(f"   {n}   {ik}")
print()
if hits_tan:
    print("Compounds the MODEL cannot distinguish from a training compound (identical ECFP-4):")
    for n, s in hits_tan:
        print(f"   {n}   maxT={s}")
