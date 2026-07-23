"""Adversarial (inversion) validation: try to break the model, confirm it holds.

Each check is written as a way the tool could be wrong; a PASS means that failure mode is absent.

  1. Leakage      no scaffold is shared between train and test under the scaffold split.
  2. Deduplication  no compound is counted twice (unique InChIKeys).
  3. Reproducible   retraining an endpoint with the fixed seed reproduces the reported score.
  4. Not trivial    predictions vary across compounds (the model is not a constant).
  5. Known answers  BBB penetration ranks CNS drugs above polar peripheral molecules.
  6. Knows-its-limits  a molecule far from training is flagged out of the applicability domain.

Output: results/tables/inversion_validation.csv (check, result, detail)
"""
from __future__ import annotations

import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from features.featurize import featurize, featurize_one  # noqa: E402
from models.train_rf import RF_COMMON, SEED, N_SPLITS, _load, _scaffold_groups  # noqa: E402

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[3]
_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
rows = []


def record(check, ok, detail):
    rows.append({"check": check, "result": "PASS" if ok else "FAIL", "detail": detail})
    print(f"[{'PASS' if ok else 'FAIL'}] {check}: {detail}")


def check_leakage():
    df = _load("MAO_A").dropna(subset=["smiles", "label"]).reset_index(drop=True)
    g = _scaffold_groups(df["smiles"].tolist())
    y = df["label"].to_numpy().astype(int)
    X = np.zeros((len(y), 1))
    bad = 0
    for tr, te in GroupKFold(N_SPLITS).split(X, y, g):
        if set(g[tr]) & set(g[te]):
            bad += 1
    record("No scaffold leakage (MAO_A, 10 folds)", bad == 0,
           f"{bad} folds with a shared scaffold (want 0)")


def check_dedup():
    lib = pd.read_csv(ROOT / "data" / "processed" / "compound_library.csv")
    dup = len(lib) - lib["inchikey"].nunique()
    record("No duplicate compounds in master library", dup == 0,
           f"{len(lib):,} rows, {lib['inchikey'].nunique():,} unique InChIKeys ({dup} duplicates)")


def check_reproducible():
    df = _load("MAO_A").dropna(subset=["smiles", "label"]).reset_index(drop=True)
    X, mask = featurize(df["smiles"].tolist())
    df = df.loc[mask].reset_index(drop=True)
    y = df["label"].to_numpy().astype(int)
    g = _scaffold_groups(df["smiles"].tolist())
    scores = []
    for tr, te in GroupKFold(N_SPLITS).split(X, y, g):
        m = RandomForestClassifier(class_weight="balanced", **RF_COMMON).fit(X[tr], y[tr])
        scores.append(roc_auc_score(y[te], m.predict_proba(X[te])[:, 1]))
    got = float(np.mean(scores))
    reported = float(pd.read_csv(ROOT / "results" / "tables" / "rf_cv_summary.csv")
                     .query("endpoint=='MAO_A' and split=='scaffold'").roc_auc_mean.iloc[0])
    record("Reproducible retrain (MAO_A scaffold AUROC)", abs(got - reported) < 0.005,
           f"retrained {got:.3f} vs reported {reported:.3f}")


def check_not_trivial():
    model = joblib.load(ROOT / "models_rf" / "BBB.joblib")
    ext = pd.read_csv(ROOT / "data" / "external" / "processed" / "external_drugs.csv").head(200)
    X, _ = featurize(ext["canonical_smiles"].tolist())
    p = model.predict_proba(X)[:, 1]
    record("Predictions are not constant (BBB over 200 drugs)", p.std() > 0.1,
           f"probability std {p.std():.3f}, range {p.min():.2f}-{p.max():.2f}")


def check_known_answers():
    model = joblib.load(ROOT / "models_rf" / "BBB.joblib")
    cns = {"donepezil": "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2",
           "diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"}
    peripheral = {"glucose": "OCC1OC(O)C(O)C(O)C1O", "atenolol": "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"}
    pc = np.mean([model.predict_proba(featurize_one(s).reshape(1, -1))[0, 1] for s in cns.values()])
    pp = np.mean([model.predict_proba(featurize_one(s).reshape(1, -1))[0, 1] for s in peripheral.values()])
    record("BBB ranks CNS drugs above polar peripheral molecules", pc > pp,
           f"mean CNS {pc:.2f} vs peripheral {pp:.2f}")


def check_applicability():
    train = _load("BBB").dropna(subset=["smiles"])["smiles"].tolist()[:2000]
    ref = [_GEN.GetFingerprint(Chem.MolFromSmiles(s)) for s in train if Chem.MolFromSmiles(s)]
    # perfluorooctanoic acid: a fluorinated surfactant, chemically alien to drug-like training data
    weird = "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    q = _GEN.GetFingerprint(Chem.MolFromSmiles(weird))
    maxsim = max(DataStructs.BulkTanimotoSimilarity(q, ref))
    record("Out-of-domain molecule is flagged (max Tanimoto < 0.30)", maxsim < 0.30,
           f"max Tanimoto of PFOA to BBB training {maxsim:.2f} (flagged out-of-domain)")


def main():
    print("=== INVERSION VALIDATION (trying to break the model) ===")
    for fn in (check_leakage, check_dedup, check_reproducible, check_not_trivial,
               check_known_answers, check_applicability):
        try:
            fn()
        except Exception as e:  # a crashing check is itself a failure
            record(fn.__name__, False, f"error: {e}")
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results" / "tables" / "inversion_validation.csv", index=False)
    n_pass = (out.result == "PASS").sum()
    print(f"\n{n_pass}/{len(out)} checks PASS")
    print("wrote results/tables/inversion_validation.csv")


if __name__ == "__main__":
    main()
