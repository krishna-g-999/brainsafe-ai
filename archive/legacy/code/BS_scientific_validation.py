"""
BS_scientific_validation.py — BrainSafe AI (BS)
Honest, leak-free re-validation of the neuroprotection model.

WHY THIS EXISTS
  The shipped validation_report.json reports a RANDOM-split hold-out R2 ~0.55,
  but (a) the applicability-domain reference contained all compounds incl.
  hold-out, (b) an addendum showed NEGATIVE 5-fold R2, and (c) ~35% of labels
  are themselves ML-generated. This script measures the TRUE generalization of
  the model architecture under conditions that cannot leak.

SCIENTIFIC GUARANTEES (no leakage):
  1. SCAFFOLD split — GroupKFold on Bemis-Murcko generic scaffolds, so no
     scaffold appears in both train and test (tests novel-chemotype transfer,
     not analog interpolation).
  2. Every fit-able transform (ECFP-PCA, ChemBERTa-PCA, StandardScaler) is fit
     ONLY on the training fold and applied to the test fold.
  3. Raw per-molecule featurization (ECFP bits, ChemBERTa embedding, RDKit
     descriptors) is a fixed function of structure — precomputed once; this is
     NOT leakage (no label or split information is used).
  4. ML-pseudo-labelled rows (data_source == 'chembl_ml_predicted') are
     QUARANTINED: the model is evaluated on the human/literature subset AND on
     the full set, reported separately.
  5. Honest BASELINES under the identical CV: mean predictor, descriptor-Ridge,
     and a Tanimoto k-NN (a strong, interpretable chemistry baseline). The
     93-feature ensemble must beat these to justify its complexity.

OUTPUT: BS_validation_report.json  +  a printed table.

Run:  python BS_scientific_validation.py
"""
from __future__ import annotations
import os, sys, json, warnings
os.environ.setdefault("HF_HUB_OFFLINE", "1")
os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")
os.chdir(os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.ensemble import (RandomForestRegressor, GradientBoostingRegressor,
                              ExtraTreesRegressor)
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors, QED
from rdkit.Chem.Lipinski import NumHDonors, NumHAcceptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from scorer import neuro_score

DIMS = ["antioxidant", "anti_inflammatory", "mitochondrial_support",
        "aggregation_modulation", "cognitive_enhancement", "neurogenesis",
        "synaptic_plasticity"]
DATA = "data/brainsafe_SCIENTIFIC_FIXED.csv"
N_SPLITS = 5
RANDOM_STATE = 42
HUMAN_SOURCES = {"curated", "curated_v2", "curated_smiles", "literature_derived"}
_CB_CACHE = "BS_chemberta_cache.npz"


# ---------------------------------------------------------------------------
# Featurization (fixed per-molecule; computed once)
# ---------------------------------------------------------------------------
def ecfp_bits(smiles, n_bits=1024):
    out = np.zeros((len(smiles), n_bits), dtype=np.float32)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s) if s else None
        if m:
            fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, n_bits)
            out[i] = np.array(fp, dtype=np.float32)
    return out


def ecfp_bitvects(smiles, n_bits=1024):
    bvs = []
    for s in smiles:
        m = Chem.MolFromSmiles(s) if s else None
        bvs.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, n_bits) if m else None)
    return bvs


def chemberta(smiles, batch=32):
    """Mean-pooled ChemBERTa embeddings (cached). Loads model ONCE."""
    if os.path.exists(_CB_CACHE):
        z = np.load(_CB_CACHE, allow_pickle=True)
        if list(z["smiles"]) == list(smiles):
            print("  [chemberta] using cache")
            return z["emb"]
    print("  [chemberta] computing embeddings (one-time)...")
    from transformers import AutoModel, AutoTokenizer
    import torch
    tok = AutoTokenizer.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
    mdl = AutoModel.from_pretrained("seyonec/ChemBERTa-zinc-base-v1"); mdl.eval()
    embs = []
    for i in range(0, len(smiles), batch):
        chunk = [s if s else "C" for s in smiles[i:i+batch]]
        enc = tok(chunk, padding=True, truncation=True, max_length=128, return_tensors="pt")
        with torch.no_grad():
            out = mdl(**enc)
        embs.append(out.last_hidden_state.mean(dim=1).numpy().astype(np.float32))
    emb = np.vstack(embs)
    np.savez(_CB_CACHE, smiles=np.array(smiles, dtype=object), emb=emb)
    return emb


def rdkit_descriptors(smiles):
    """7 interpretable descriptors for the Ridge baseline."""
    out = np.zeros((len(smiles), 7), dtype=np.float32)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s) if s else None
        if m:
            out[i] = [Descriptors.MolWt(m), Crippen.MolLogP(m), Descriptors.TPSA(m),
                      NumHDonors(m), NumHAcceptors(m),
                      rdMolDescriptors.CalcNumRotatableBonds(m), QED.qed(m)]
    return out


def generic_scaffold(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None:
        return "INVALID"
    try:
        core = MurckoScaffold.GetScaffoldForMol(m)
        gen = MurckoScaffold.MakeScaffoldGeneric(core)
        s = Chem.MolToSmiles(gen)
        return s if s else "ACYCLIC"
    except Exception:
        return "ACYCLIC"


# ---------------------------------------------------------------------------
# Models — production architecture (identical hyperparameters) + baselines
# ---------------------------------------------------------------------------
def make_ensemble():
    return {
        "rf":  RandomForestRegressor(n_estimators=500, max_depth=8, min_samples_leaf=2,
                                     random_state=RANDOM_STATE, n_jobs=-1),
        "gbt": GradientBoostingRegressor(n_estimators=250, learning_rate=0.04, max_depth=4,
                                         random_state=RANDOM_STATE),
        "et":  ExtraTreesRegressor(n_estimators=500, max_depth=8, min_samples_leaf=2,
                                   random_state=RANDOM_STATE, n_jobs=-1),
        "ridge": Ridge(alpha=8),
    }


def fold_features_ensemble(raw, tr, te):
    """Fit ECFP-PCA(50)+ChemBERTa-PCA(32) on TRAIN only; build 93-D matrices."""
    ep = PCA(n_components=50, random_state=RANDOM_STATE).fit(raw["ecfp"][tr])
    cp = PCA(n_components=32, random_state=RANDOM_STATE).fit(raw["cb"][tr])
    def build(idx):
        return np.hstack([ep.transform(raw["ecfp"][idx]), cp.transform(raw["cb"][idx]),
                          raw["disease"][idx], raw["bbb"][idx], raw["struct"][idx]])
    return build(tr), build(te)


def predict_ensemble(Xtr, ytr, Xte):
    """Per-dimension 4-model ensemble (matches production; GBT is single-output)."""
    sc = StandardScaler().fit(Xtr)
    Xtr_s, Xte_s = sc.transform(Xtr), sc.transform(Xte)
    out = np.zeros((Xte.shape[0], ytr.shape[1]), dtype=float)
    for di in range(ytr.shape[1]):
        preds = []
        for name, mdl in make_ensemble().items():
            if name == "ridge":
                mdl.fit(Xtr_s, ytr[:, di]); preds.append(mdl.predict(Xte_s))
            else:
                mdl.fit(Xtr, ytr[:, di]); preds.append(mdl.predict(Xte))
        out[:, di] = np.mean(preds, axis=0)
    return np.clip(out, 0, 100)


def predict_ridge_desc(desc_tr, ytr, desc_te):
    sc = StandardScaler().fit(desc_tr)
    m = Ridge(alpha=8).fit(sc.transform(desc_tr), ytr)
    return np.clip(m.predict(sc.transform(desc_te)), 0, 100)


def predict_knn(bv_tr, ytr, bv_te, k=5):
    """Tanimoto-weighted k-NN regressor on ECFP bit vectors."""
    out = np.zeros((len(bv_te), ytr.shape[1]), dtype=np.float32)
    for i, q in enumerate(bv_te):
        if q is None:
            out[i] = ytr.mean(axis=0); continue
        sims = np.array(DataStructs.BulkTanimotoSimilarity(q, bv_tr))
        idx = np.argsort(sims)[::-1][:k]
        w = sims[idx]
        out[i] = ytr[idx].mean(axis=0) if w.sum() == 0 else (w[:, None] * ytr[idx]).sum(0) / w.sum()
    return np.clip(out, 0, 100)


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------
def metrics(y_true, y_pred):
    per = {}
    for di, d in enumerate(DIMS):
        m = ~np.isnan(y_true[:, di])
        per[d] = round(r2_score(y_true[m, di], y_pred[m, di]), 3) if m.sum() > 2 else None
    # NPS from predicted dims
    nps_t = np.array([neuro_score({d: y_true[i, di] for di, d in enumerate(DIMS)}) for i in range(len(y_true))])
    nps_p = np.array([neuro_score({d: y_pred[i, di] for di, d in enumerate(DIMS)}) for i in range(len(y_pred))])
    return {
        "nps_r2": round(r2_score(nps_t, nps_p), 3),
        "nps_spearman": round(spearmanr(nps_t, nps_p).correlation, 3),
        "nps_mae": round(mean_absolute_error(nps_t, nps_p), 2),
        "per_dim_r2": per,
        "mean_dim_r2": round(np.mean([v for v in per.values() if v is not None]), 3),
    }


def evaluate(df, raw, label):
    print(f"\n========== EVALUATION: {label} (n={len(df)}) ==========")
    groups = df["scaffold"].values
    n_scaff = len(set(groups))
    print(f"  unique generic scaffolds: {n_scaff}  (GroupKFold k={N_SPLITS})")
    y = df[DIMS].values.astype(float)
    gkf = GroupKFold(n_splits=N_SPLITS)

    oof = {m: np.full_like(y, np.nan, dtype=float) for m in ["mean", "ridge_desc", "knn", "ensemble93"]}
    max_cross_sim = []
    for tr, te in gkf.split(df, groups=groups):
        # leakage audit: max Tanimoto between any test and any train compound
        bvtr = [raw["bv"][i] for i in tr if raw["bv"][i] is not None]
        for i in te:
            if raw["bv"][i] is not None and bvtr:
                max_cross_sim.append(max(DataStructs.BulkTanimotoSimilarity(raw["bv"][i], bvtr)))
        ytr = y[tr]
        # mean baseline
        oof["mean"][te] = np.nanmean(ytr, axis=0)
        # ridge on descriptors
        oof["ridge_desc"][te] = predict_ridge_desc(raw["desc"][tr], np.nan_to_num(ytr, nan=np.nanmean(ytr)), raw["desc"][te])
        # tanimoto knn
        oof["knn"][te] = predict_knn([raw["bv"][i] for i in tr], np.nan_to_num(ytr, nan=np.nanmean(ytr)),
                                     [raw["bv"][i] for i in te])
        # ensemble 93
        Xtr, Xte = fold_features_ensemble(raw, tr, te)
        oof["ensemble93"][te] = predict_ensemble(Xtr, np.nan_to_num(ytr, nan=np.nanmean(ytr)), Xte)

    cross = {"mean_test_train_max_tanimoto": round(float(np.mean(max_cross_sim)), 3),
             "median": round(float(np.median(max_cross_sim)), 3),
             "frac_test_with_T>=0.7": round(float(np.mean(np.array(max_cross_sim) >= 0.7)), 3)}
    print(f"  scaffold-split leakage audit (test->train max Tanimoto): "
          f"mean={cross['mean_test_train_max_tanimoto']}, median={cross['median']}, "
          f"frac>=0.7={cross['frac_test_with_T>=0.7']}")

    results = {}
    print(f"\n  {'model':14}  NPS_R2  NPS_rho  NPS_MAE  meanDimR2")
    for m in ["mean", "ridge_desc", "knn", "ensemble93"]:
        results[m] = metrics(y, oof[m])
        r = results[m]
        print(f"  {m:14}  {r['nps_r2']:6}  {r['nps_spearman']:7}  {r['nps_mae']:7}  {r['mean_dim_r2']:8}")
    print("\n  per-dimension R2 (scaffold-CV):")
    print(f"    {'dim':26} {'mean':>6} {'ridgeD':>7} {'kNN':>6} {'ens93':>6}")
    for d in DIMS:
        row = "    {:26}".format(d)
        for m in ["mean", "ridge_desc", "knn", "ensemble93"]:
            v = results[m]["per_dim_r2"][d]
            row += f" {('' if v is None else v):>6}" if m == "mean" else f" {('' if v is None else v):>7}" if m == "ridge_desc" else f" {('' if v is None else v):>6}"
        print(row)
    return {"n": len(df), "n_scaffolds": n_scaff, "leakage_audit": cross, "models": results}


def main():
    df = pd.read_csv(DATA)
    df = df[df["smiles"].notna() & (df["smiles"].astype(str).str.lower() != "nan")].reset_index(drop=True)
    df = df[df[DIMS].notna().all(axis=1)].reset_index(drop=True)   # require all 7 labels
    smiles = df["smiles"].astype(str).tolist()
    print(f"Loaded {len(df)} compounds with SMILES + complete labels.")

    print("Featurizing (once)...")
    raw = {
        "ecfp": ecfp_bits(smiles),
        "cb":   chemberta(smiles),
        "desc": rdkit_descriptors(smiles),
        "bv":   ecfp_bitvects(smiles),
        "disease": np.stack([pd.to_numeric(df.get(c, 0), errors="coerce").fillna(0).values
                             for c in ["ad_target_count", "pd_target_count",
                                       "als_target_count", "hd_target_count"]], axis=1).astype(np.float32),
        "struct": np.stack([pd.to_numeric(df.get(c, 0), errors="coerce").fillna(0).values /
                            n for c, n in zip(["phenolic_oh_count", "catechol_flag", "aromatic_ring_count",
                                               "hbd_count", "tpsa", "rotatable_bonds"],
                                              [10., 1., 10., 15., 200., 20.])], axis=1).astype(np.float32),
    }
    from model_config import BBB_MAP
    raw["bbb"] = df["bbb"].map(BBB_MAP).fillna(0).values.reshape(-1, 1).astype(np.float32)
    df["scaffold"] = [generic_scaffold(s) for s in smiles]

    report = {"data": DATA, "n_total_labelled": len(df),
              "cv": f"scaffold GroupKFold(k={N_SPLITS}) on Bemis-Murcko generic scaffolds",
              "evaluations": {}}

    # Eval A: human/literature labels only (pseudo-labels quarantined)
    human = df[df["data_source"].isin(HUMAN_SOURCES)].reset_index(drop=True)
    if len(human) >= 50:
        idx = df.index[df["data_source"].isin(HUMAN_SOURCES)].tolist()
        raw_h = {k: (v[idx] if isinstance(v, np.ndarray) else [v[i] for i in idx]) for k, v in raw.items()}
        report["evaluations"]["human_labels_only"] = evaluate(human, raw_h, "HUMAN/LITERATURE LABELS ONLY")

    # Eval B: full set (incl. ML pseudo-labels)
    report["evaluations"]["full_incl_ml_labels"] = evaluate(df, raw, "FULL SET (incl. ML pseudo-labels)")

    with open("BS_validation_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print("\nSaved BS_validation_report.json")


if __name__ == "__main__":
    main()
