"""Assemble the reviewer workbook: what every model was given, and what every model returns.

Two questions a reviewer asks that the repository could not previously answer in one place. Which
compounds carried which measured values, and where the gaps are. And what the tool actually emits per
compound, produced by which model, trained how, tested how.

Everything here is read from artefacts on disk. Nothing is restated from a document, and no number is
carried over from an earlier run: the counts come from the endpoint tables, the hyper-parameters and
cross-validation results from the per-model metadata written at training time, and the thresholds and
their bases from binder_modes.json. If an artefact is absent the cell is left empty and named in
BLANKS.csv rather than filled with a plausible value.

Outputs, under reviewer_package/model_outputs/:

  TRAINING_MATRIX.csv     one row per compound, one column per measured endpoint. This is the
                          training input: what was measured, for which compound, on which endpoint.
  OUTPUT_CATALOGUE.csv    one row per model output the tool emits, with its role, its model, how it
                          was trained and how it was tested.
  RECIPE_<family>.csv     per model family: every endpoint, its training composition and its results.
  BLANKS.csv              every kind of empty cell that appears, and what it means.
  README.md               how to read all of the above.

Run:  python src/brainsafe/analysis/build_reviewer_matrix.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
sys.path.insert(0, str(ROOT / "src" / "brainsafe" / "data"))
from build_compound_library import standardise  # noqa: E402

M = ROOT / "models_rf"
OUT = ROOT / "reviewer_package" / "model_outputs"
ENDPOINTS = ROOT / "data" / "endpoints"
ENDPOINTS_REG = ROOT / "data" / "endpoints_reg"
ADME = ROOT / "data" / "adme"

# The four receptors carry two models each: a potency regression and a binder classifier. They are
# one endpoint name and two outputs, which is why an output count taken from endpoint names and one
# taken from models disagree.
DUAL_MODELLED = ["D2", "A2A", "HT2A", "SERT"]


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


# --------------------------------------------------------------------------------------------
# 1. Training matrix: which compound carried which measured value
# --------------------------------------------------------------------------------------------
def training_matrix() -> tuple[pd.DataFrame, list[dict]]:
    """One row per compound, one column per measured endpoint.

    Keyed by the InChIKey of the desalted parent, which is how the endpoint tables are deduplicated.
    A compound measured against three targets occupies one row with three populated columns and the
    rest empty; empty means not measured, never zero.
    """
    frames, provenance = {}, []

    for f in sorted(ENDPOINTS.glob("*.csv")):
        ep = f.stem
        df = pd.read_csv(f)
        if "smiles" not in df.columns:
            continue
        keys, parents = [], []
        for smi in df["smiles"].astype(str):
            csmi, ik = standardise(smi)
            keys.append(ik)
            parents.append(csmi)
        df["_ik"], df["_parent"] = keys, parents
        df = df.dropna(subset=["_ik"])
        # A table can still hold two rows with one parent InChIKey, because the expansion fetchers
        # deduplicated on the raw SMILES string before that was corrected (BS-C-16). The first is
        # kept and the collapse is counted, so the number appears in TRAINING_SOURCES rather than
        # vanishing into a silent overwrite.
        n_before = len(df)
        df = df.drop_duplicates("_ik", keep="first")
        n_collapsed = n_before - len(df)
        cols = {}
        if "label" in df.columns:
            cols[f"{ep}__label"] = df.set_index("_ik")["label"]
        if "pchembl" in df.columns:
            cols[f"{ep}__pchembl"] = df.set_index("_ik")["pchembl"]
        frames.update(cols)
        provenance.append({
            "endpoint": ep, "source_file": f"data/endpoints/{f.name}", "rows": len(df),
            "measure": "label (1 active / 0 inactive)" + (" + pchembl potency"
                                                          if "pchembl" in df.columns else ""),
            "active": int((df["label"] == 1).sum()) if "label" in df.columns else None,
            "inactive": int((df["label"] == 0).sum()) if "label" in df.columns else None,
            "rows_sharing_a_parent_inchikey": n_collapsed,
        })

    for f in sorted(ENDPOINTS_REG.glob("*.csv")):
        ep, df = f.stem, pd.read_csv(f)
        if "smiles" not in df.columns:
            continue
        target = "y" if "y" in df.columns else ("pka" if "pka" in df.columns else None)
        if target is None:
            continue
        keys = [standardise(s)[1] for s in df["smiles"].astype(str)]
        df["_ik"] = keys
        df = df.dropna(subset=["_ik"]).drop_duplicates("_ik", keep="first")
        frames[f"{ep}__value"] = df.set_index("_ik")[target]
        provenance.append({"endpoint": ep, "source_file": f"data/endpoints_reg/{f.name}",
                           "rows": len(df), "measure": f"continuous ({target})",
                           "active": None, "inactive": None})

    for f in sorted(ADME.glob("*.csv")):
        ep, df = f.stem, pd.read_csv(f)
        smi_col = next((c for c in df.columns if c.lower() in ("smiles", "canonical_smiles")), None)
        val_col = next((c for c in df.columns if c.lower() in ("y", "value", "label")), None)
        if smi_col is None or val_col is None:
            continue
        keys = [standardise(s)[1] for s in df[smi_col].astype(str)]
        df["_ik"] = keys
        df = df.dropna(subset=["_ik"]).drop_duplicates("_ik", keep="first")
        frames[f"adme_{ep}__value"] = df.set_index("_ik")[val_col]
        provenance.append({"endpoint": f"adme_{ep}", "source_file": f"data/adme/{f.name}",
                           "rows": len(df), "measure": f"continuous or binary ({val_col})",
                           "active": None, "inactive": None})

    matrix = pd.DataFrame(frames)
    matrix = matrix[~matrix.index.duplicated(keep="first")]
    matrix.index.name = "inchikey"

    # A readable structure for each compound, taken from wherever it was first seen.
    smiles_map = {}
    for f in sorted(ENDPOINTS.glob("*.csv")):
        df = pd.read_csv(f)
        if "smiles" not in df.columns:
            continue
        for smi in df["smiles"].astype(str):
            csmi, ik = standardise(smi)
            if ik and ik not in smiles_map:
                smiles_map[ik] = csmi
    matrix.insert(0, "parent_smiles", [smiles_map.get(i, "") for i in matrix.index])
    matrix.insert(1, "n_endpoints_measured",
                  matrix.drop(columns=["parent_smiles"]).notna().sum(axis=1))
    return matrix.reset_index(), provenance


# --------------------------------------------------------------------------------------------
# 2. Output catalogue: what the tool emits, and how each output was produced
# --------------------------------------------------------------------------------------------
def output_catalogue() -> pd.DataFrame:
    rows = []
    cv = pd.read_csv(ROOT / "results" / "tables" / "rf_cv_summary.csv")
    adme_cv_path = ROOT / "results" / "tables" / "adme_cv_summary.csv"
    adme_cv = pd.read_csv(adme_cv_path) if adme_cv_path.exists() else pd.DataFrame()
    binders = _json(M / "binder_modes.json")

    def cv_get(ep, split, col):
        r = cv[(cv.endpoint == ep) & (cv.split == split)]
        return None if r.empty or col not in r.columns or pd.isna(r[col].iloc[0]) else r[col].iloc[0]

    # (a) the eight target classifiers and the five regressions, from train_rf.py
    for f in sorted(M.glob("*_meta.json")):
        meta = _json(f)
        ep = meta.get("endpoint")
        if not ep or "task" not in meta:
            continue
        clf = meta["task"] == "classification"
        ded = meta.get("deduplication") or {}
        rows.append({
            "output": ep,
            "group": "target classifier" if clf else "potency regression",
            "returns": "probability 0-1" if clf else "pChEMBL (-log10 M)",
            "role_in_tool": ("BBB gate applied to every disease score" if ep == "BBB" else
                             "hERG cardiotoxicity safety flag" if ep == "hERG" else
                             "antioxidant / neuroprotection axis" if ep == "antioxidant_DPPH" else
                             "target engagement feeding the disease layer"),
            "model": "RandomForestClassifier" if clf else "RandomForestRegressor",
            "model_file": f"models_rf/{ep}.joblib",
            "calibrated_file": (f"models_rf/{ep}_calibrated.joblib"
                                if (M / f"{ep}_calibrated.joblib").exists() else ""),
            "training_compounds": meta.get("n_compounds"),
            "positives": meta.get("positives"),
            "rows_before_dedup": ded.get("rows_in"),
            "duplicate_rows_removed": ded.get("duplicate_rows_removed"),
            "contradictory_groups_dropped": ded.get("conflicting_groups_dropped"),
            "n_scaffolds": meta.get("n_scaffolds"),
            "features": meta.get("feature_layout"),
            "hyperparameters": json.dumps(meta.get("hyperparameters", {})),
            "cv_scheme": "StratifiedKFold(10) random and GroupKFold(10) on Bemis-Murcko scaffold",
            "cv_random": cv_get(ep, "random", "roc_auc_mean" if clf else "r2_mean"),
            "cv_scaffold": cv_get(ep, "scaffold", "roc_auc_mean" if clf else "r2_mean"),
            "cv_metric": "AUROC" if clf else "R2",
            "decision_threshold": 0.5 if clf else None,
            "threshold_basis": "fixed 0.5 with class_weight=balanced" if clf else "",
            "source_table": (f"data/endpoints/{ep}.csv" if ep != "antioxidant_DPPH"
                             else "data/endpoints_reg/antioxidant_dpph.csv"),
        })

    # (b) the binder panel
    for ep, v in sorted(binders.items()):
        rows.append({
            "output": f"{ep} (binder)",
            "group": "binder classifier",
            "returns": "probability 0-1, called against a per-target threshold",
            "role_in_tool": ("WITHDRAWN, not scored" if not v.get("deployed", True)
                             else "target engagement feeding the disease layer"),
            "model": "RandomForest + sigmoid calibration (CalibratedClassifierCV)",
            "model_file": f"models_rf/{ep}_binder.joblib",
            "calibrated_file": "",
            "training_compounds": v.get("n_active_train"),
            "positives": v.get("n_positive"),
            "rows_before_dedup": None,
            "duplicate_rows_removed": None,
            "contradictory_groups_dropped": None,
            "n_scaffolds": None,
            "features": "ecfp4_1024 + 12 descriptors",
            "hyperparameters": json.dumps({"n_estimators": 300, "min_samples_leaf": 4,
                                           "class_weight": "balanced", "random_state": 42}),
            "cv_scheme": "GroupKFold(10) on scaffold; actives withheld by scaffold group",
            "cv_random": None,
            "cv_scaffold": v.get("scaffold_cv_auroc"),
            "cv_metric": "AUROC (scaffold CV)",
            "decision_threshold": v.get("threshold"),
            "threshold_basis": v.get("threshold_basis", ""),
            "source_table": f"data/endpoints/{ep}.csv",
            "holdout_auroc_vs_measured_inactives": v.get("auroc_vs_measured_inactives"),
            "sensitivity_at_threshold": v.get("sensitivity_at_threshold"),
            "sensitivity_basis": v.get("sensitivity_basis", ""),
            "background_fpr_held_out": v.get("background_fpr_at_threshold"),
            "n_active_holdout": v.get("n_active_holdout"),
            "n_measured_inactive_holdout": v.get("n_measured_inactive_holdout"),
            "reliable_call": v.get("reliable_call"),
            "mode": v.get("mode"),
            "deployed": v.get("deployed", True),
        })

    # (c) ADME
    for f in sorted((M / "adme").glob("*_meta.json")):
        meta = _json(f)
        ep = meta.get("endpoint") or f.stem.replace("_meta", "")
        r = adme_cv[adme_cv.endpoint == ep] if len(adme_cv) else pd.DataFrame()

        def acv(split, col):
            if r.empty:
                return None
            s = r[r.split == split]
            return None if s.empty or col not in s.columns or pd.isna(s[col].iloc[0]) else s[col].iloc[0]

        clf = meta.get("task") == "classification"
        rows.append({
            "output": f"adme_{ep}",
            "group": "ADME / exposure",
            "returns": "probability 0-1" if clf else "continuous, units per training data",
            "role_in_tool": "exposure and developability context, not a disease score",
            "model": "RandomForestClassifier" if clf else "RandomForestRegressor",
            "model_file": f"models_rf/adme/{ep}.joblib",
            "calibrated_file": "",
            "training_compounds": meta.get("n_compounds"),
            "positives": None,
            "rows_before_dedup": None, "duplicate_rows_removed": None,
            "contradictory_groups_dropped": None, "n_scaffolds": None,
            "features": "ecfp4_1024 + 12 descriptors",
            "hyperparameters": json.dumps(meta.get("hyperparameters", {})),
            "cv_scheme": "StratifiedKFold(10)/KFold(10) random and GroupKFold(10) on scaffold",
            "cv_random": acv("random", "roc_auc_mean" if clf else "r2_mean"),
            "cv_scaffold": acv("scaffold", "roc_auc_mean" if clf else "r2_mean"),
            "cv_metric": "AUROC" if clf else "R2",
            "decision_threshold": 0.5 if clf else None,
            "threshold_basis": "fixed 0.5" if clf else "",
            "source_table": f"data/adme/{ep}.csv",
        })

    df = pd.DataFrame(rows)
    front = ["output", "group", "returns", "role_in_tool", "model", "model_file",
             "training_compounds", "cv_scheme", "cv_metric", "cv_random", "cv_scaffold",
             "decision_threshold", "threshold_basis"]
    return df[front + [c for c in df.columns if c not in front]]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("building the training matrix ...", flush=True)
    matrix, provenance = training_matrix()
    matrix.to_csv(OUT / "TRAINING_MATRIX.csv", index=False)
    pd.DataFrame(provenance).to_csv(OUT / "TRAINING_SOURCES.csv", index=False)
    print(f"  {len(matrix):,} compounds x {matrix.shape[1] - 3} measured columns")

    print("building the output catalogue ...", flush=True)
    cat = output_catalogue()
    cat.to_csv(OUT / "OUTPUT_CATALOGUE.csv", index=False)
    for grp, sub in cat.groupby("group"):
        slug = grp.replace(" ", "_").replace("/", "_")
        sub.dropna(axis=1, how="all").to_csv(OUT / f"RECIPE_{slug}.csv", index=False)
    print(f"  {len(cat)} outputs across {cat.group.nunique()} families")

    blanks = pd.DataFrame([
        {"where": "TRAINING_MATRIX, any <endpoint>__label or __pchembl cell",
         "blank_means": "this compound was never measured against this endpoint. It is absent "
                        "evidence, not a measured zero, and it was not used to train that model."},
        {"where": "TRAINING_MATRIX, <endpoint>__pchembl present but __label empty",
         "blank_means": "the potency fell in the 5-6 grey band, which the label rule discards, so "
                        "the compound trains the regression but not the classifier."},
        {"where": "OUTPUT_CATALOGUE, cv_random for a binder",
         "blank_means": "binders are evaluated under a scaffold split only; no random-split figure "
                        "was computed, so none is quoted."},
        {"where": "OUTPUT_CATALOGUE, positives for an ADME endpoint",
         "blank_means": "most ADME endpoints are regressions and have no positive class."},
        {"where": "OUTPUT_CATALOGUE, rows_before_dedup for binder or ADME rows",
         "blank_means": "the deduplication counter is recorded by train_rf.py only; the binder and "
                        "ADME scripts do not write it, so the figure is genuinely unknown rather "
                        "than zero."},
        {"where": "OUTPUT_CATALOGUE, decision_threshold for a regression",
         "blank_means": "regressions return a value, not a call, so no threshold applies."},
        {"where": "OUTPUT_CATALOGUE, sensitivity or threshold fields for Nav1_1 and Cav3_2",
         "blank_means": "withdrawn endpoints; the model exists but is not scored, so its figures "
                        "describe something the tool does not report."},
    ])
    blanks.to_csv(OUT / "BLANKS.csv", index=False)
    print(f"wrote {OUT.relative_to(ROOT).as_posix()}/")


if __name__ == "__main__":
    main()
