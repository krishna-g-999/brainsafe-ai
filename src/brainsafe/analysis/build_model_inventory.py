"""One dated row per deployed estimator: what it predicts, what it was trained on, how it was tested.

The panel is described in several places and each description was assembled for a different purpose,
so the counts have drifted apart before. This builds the single table everything else should quote,
directly from the estimators and their metadata on disk, and stamps it with the date and the commit
so a reader can tell which state it describes.

One row per deployed estimator, not per endpoint. Four receptors carry both a potency regression and
a binder classifier, so an endpoint count and a model count legitimately differ, and conflating them
is how "69 models" and "72 estimators" both came to be written about the same panel.

Columns are chosen so that a reviewer asking "what is behind this number" can answer it without
opening the code: the prediction type, the training-set size after deduplication, the class balance,
the validation scheme, the headline metric with the split it came from, and the file's modification
time, which is the only honest statement of when the estimator was last fitted.

Outputs:
  results/tables/MODEL_INVENTORY.csv      one row per deployed estimator
  results/tables/MODEL_INVENTORY.md       the same, grouped, for reading

Run:  python src/brainsafe/analysis/build_model_inventory.py
"""
from __future__ import annotations

import datetime as dt
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
M = ROOT / "models_rf"
TAB = ROOT / "results" / "tables"

RF_CORE = "RandomForest, 300 trees, min_samples_leaf=2, seed 42"
RF_BINDER = "RandomForest, 300 trees, min_samples_leaf=4, seed 42"


def _json(p: Path):
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _mtime(p: Path) -> str:
    return dt.datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M") \
        if p.exists() else ""


def rows() -> list[dict]:
    cv = pd.read_csv(TAB / "rf_cv_summary.csv") if (TAB / "rf_cv_summary.csv").exists() \
        else pd.DataFrame()
    adme_cv = pd.read_csv(TAB / "adme_cv_summary.csv") if (TAB / "adme_cv_summary.csv").exists() \
        else pd.DataFrame()
    out = []

    def cv_lookup(frame, ep, split, col):
        if not len(frame):
            return None
        g = frame[(frame.endpoint == ep) & (frame.split == split)]
        return round(float(g[col].iloc[0]), 4) if len(g) and pd.notna(g[col].iloc[0]) else None

    # ---- core: classifiers and potency regressions -------------------------------------------
    for p in sorted(M.glob("*.joblib")):
        if p.name.endswith(("_binder.joblib", "_calibrated.joblib")):
            continue
        ep = p.stem
        meta = _json(M / f"{ep}_meta.json")
        task = meta.get("task")
        if task is None:
            g = cv[cv.endpoint == ep] if len(cv) else pd.DataFrame()
            task = g.task.iloc[0] if len(g) else "regression"
        is_clf = task == "classification"
        metric = "AUROC" if is_clf else "R2"
        col = "roc_auc_mean" if is_clf else "r2_mean"
        out.append({
            "model": ep,
            "family": "exposure" if ep == "BBB" else ("safety" if ep == "hERG" else "target"),
            "predicts": ("probability of activity" if is_clf else
                         ("pKa" if ep == "pka_basic" else "potency, pChEMBL")),
            "task": task,
            "algorithm": RF_CORE,
            "n_train": meta.get("n_compounds"),
            "n_positive": meta.get("positives"),
            "n_scaffolds": meta.get("n_scaffolds"),
            "validation": "10-fold random and 10-fold scaffold-grouped",
            "metric": metric,
            "random_split": cv_lookup(cv, ep, "random", col),
            "scaffold_split": cv_lookup(cv, ep, "scaffold", col),
            "calibration": "isotonic, out-of-fold" if (M / f"{ep}_calibrated.joblib").exists()
                           else "none",
            "threshold_basis": "",
            "deployed": True,
            "fitted": _mtime(p),
        })

    # ---- binder panel --------------------------------------------------------------------------
    for ep, v in sorted(_json(M / "binder_modes.json").items()):
        out.append({
            "model": f"{ep}_binder",
            "family": "binder",
            "predicts": "probability this compound binds this target",
            "task": "classification",
            "algorithm": RF_BINDER,
            "n_train": (v.get("n_active_train") or 0) + (v.get("n_decoy") or 0)
                       + (v.get("n_measured_inactive_train") or 0),
            "n_positive": v.get("n_positive"),
            "n_scaffolds": None,
            "validation": ("actives withheld by scaffold; thresholds on a disjoint background pool"),
            "metric": "AUROC vs measured non-binders",
            "random_split": None,
            "scaffold_split": v.get("auroc_vs_measured_inactives"),
            "calibration": "sigmoid, prefit",
            "threshold_basis": v.get("threshold_basis", ""),
            "deployed": bool(v.get("deployed", True)),
            "fitted": _mtime(M / f"{ep}_binder.joblib"),
        })

    # ---- ADME ----------------------------------------------------------------------------------
    for p in sorted((M / "adme").glob("*.joblib")):
        ep = p.stem
        meta = _json(M / "adme" / f"{ep}_meta.json")
        task = meta.get("task", "regression")
        is_clf = task == "classification"
        col = "roc_auc_mean" if is_clf else "r2_mean"
        out.append({
            "model": f"adme_{ep}",
            "family": "exposure",
            "predicts": meta.get("units") or ("probability" if is_clf else "measured value"),
            "task": task,
            "algorithm": RF_CORE,
            "n_train": meta.get("n_compounds"),
            "n_positive": None,
            "n_scaffolds": None,
            "validation": "10-fold random and 10-fold scaffold-grouped",
            "metric": "AUROC" if is_clf else "R2",
            "random_split": cv_lookup(adme_cv, ep, "random", col),
            "scaffold_split": cv_lookup(adme_cv, ep, "scaffold", col),
            "calibration": "none",
            "threshold_basis": "",
            "deployed": True,
            "fitted": _mtime(p),
        })
    return out


def main() -> None:
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True).stdout.strip()
    today = dt.date.today().isoformat()
    df = pd.DataFrame(rows())
    df.insert(0, "inventory_date", today)
    df.insert(1, "commit", commit)
    TAB.mkdir(parents=True, exist_ok=True)
    df.to_csv(TAB / "MODEL_INVENTORY.csv", index=False)

    fam = df.groupby("family").agg(models=("model", "size"),
                                   deployed=("deployed", "sum")).reset_index()
    lines = [f"# Model inventory, {today}", "",
             f"Commit `{commit}`. One row per deployed estimator, taken from the estimators and "
             "their metadata on disk. Four receptors carry both a potency regression and a binder "
             "classifier, so the model count exceeds the endpoint count.", "",
             f"**{len(df)} estimators, {int(df.deployed.sum())} deployed.** "
             f"Fitted between {df.fitted.min()} and {df.fitted.max()}.", "",
             "| family | estimators | deployed |", "|---|---|---|"]
    for _, r in fam.iterrows():
        lines.append(f"| {r['family']} | {r['models']} | {int(r['deployed'])} |")
    lines += ["", "## Every estimator", "",
              "| model | predicts | task | n train | metric | random | scaffold | calibration | "
              "deployed | fitted |", "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        lines.append(
            f"| {r['model']} | {r['predicts']} | {r['task']} | "
            f"{'' if pd.isna(r['n_train']) else int(r['n_train'])} | {r['metric']} | "
            f"{'' if r['random_split'] is None or pd.isna(r['random_split']) else r['random_split']} | "
            f"{'' if r['scaffold_split'] is None or pd.isna(r['scaffold_split']) else r['scaffold_split']} | "
            f"{r['calibration']} | {'yes' if r['deployed'] else 'WITHDRAWN'} | {r['fitted']} |")
    (TAB / "MODEL_INVENTORY.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"{len(df)} estimators, {int(df.deployed.sum())} deployed, inventory dated {today}")
    print(fam.to_string(index=False))
    print(f"fitted between {df.fitted.min()} and {df.fitted.max()}")
    withdrawn = df[~df.deployed].model.tolist()
    if withdrawn:
        print(f"withdrawn: {', '.join(withdrawn)}")
    print("wrote results/tables/MODEL_INVENTORY.csv and .md")


if __name__ == "__main__":
    main()
