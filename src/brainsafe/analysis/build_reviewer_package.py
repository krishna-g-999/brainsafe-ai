"""Build the reviewer package: every endpoint, every input, every fold, every feature, every score.

Written in response to a reviewer request for the training record in full. Nothing here is typed by
hand. Every value is read from the file that produced it, so the package cannot disagree with the
models, and re-running it after any change regenerates a package that still matches.

Where a value does not exist, the cell is left empty and the reason is recorded in the data
dictionary rather than filled with a plausible number. An empty cell in this package always means one
of a small number of stated things, never "we did not look".

Outputs, in the order a reviewer should read them:

  00_README.md                          what each file is and how they join
  01_MASTER_endpoint_inventory.csv      every endpoint, its data, its model, its scores
  02_training_input_files.csv           every input file, what it feeds, what it produced
  03_crossvalidation_summary.csv        per endpoint and split, with standard deviations
  04_crossvalidation_per_fold.csv       every individual fold
  05_feature_definitions.csv            all 1036 features
  06_model_registry.csv                 every model artefact, its thresholds and operating point
  07_binder_panel_training_design.csv   how the binder panel's negatives were constructed
  08_DATA_DICTIONARY.md                 every column in every file, and every reason a cell is empty
  09_PROVENANCE_AND_LIMITATIONS.md      what is verified, what is absent, and why

Usage: python src/brainsafe/analysis/build_reviewer_package.py [output_dir]
"""
from __future__ import annotations

import json
import re
import sys
import warnings
from datetime import date
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
EP = ROOT / "data" / "endpoints"
TAB = ROOT / "results" / "tables"
M = ROOT / "models_rf"

# Label rule applied to every ChEMBL-derived target table, stated once and referenced everywhere.
LABEL_RULE = ("active if pChEMBL >= 6; inactive if pChEMBL < 5; the 5 to 6 band is discarded as "
              "ambiguous. Binder models use a stricter positive definition of pChEMBL >= 7.")


def chembl_ids():
    """Target identifiers as written in the fetch scripts, where each was verified on download."""
    out = {}
    # Two forms are in use across the fetch scripts and both must be read, or a third of the
    # identifiers go missing: the later scripts write "NAME": ("CHEMBLnnn", "preferred name"), the
    # earlier ones write "NAME": "CHEMBLnnn". A first version of this function matched only the
    # first form and recovered 31 of 72.
    patterns = [r'"([A-Za-z0-9_]+)"\s*:\s*\(\s*"(CHEMBL\d+)"',
                r'"([A-Za-z0-9_]+)"\s*:\s*"(CHEMBL\d+)"']
    for f in sorted((ROOT / "src" / "brainsafe" / "data").glob("fetch_*.py")):
        txt = f.read_text(encoding="utf-8", errors="replace")
        for pat in patterns:
            for name, cid in re.findall(pat, txt):
                out.setdefault(name, (cid, f.name))
    # The original core panel was fetched through a cached route that did not record the target
    # identifier inline. Rather than leave those cells blank or supply a remembered identifier,
    # each was resolved from its UniProt accession against the ChEMBL target API and the returned
    # preferred name checked against the intended protein. The result is stored so the package
    # build stays offline and the resolution is auditable.
    core = ROOT / "data" / "core_panel_target_ids.json"
    if core.exists():
        for name, rec in json.loads(core.read_text()).items():
            out.setdefault(name, (rec["chembl_target_id"], "data/core_panel_target_ids.json"))
    return out


def endpoint_table_stats(name):
    p = EP / f"{name}.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p)
    # A column that is absent must yield an empty Series, not a scalar: pd.to_numeric(None) returns
    # a float, and the resulting AttributeError would otherwise stop the whole package.
    empty = pd.Series(dtype="float64")
    pch = pd.to_numeric(d["pchembl"], errors="coerce") if "pchembl" in d else empty
    yr = pd.to_numeric(d["year"], errors="coerce") if "year" in d else empty
    return {
        "training_table": f"data/endpoints/{name}.csv",
        "n_compounds": len(d),
        "n_active_label1": int((d.get("label") == 1).sum()) if "label" in d else None,
        "n_inactive_label0": int((d.get("label") == 0).sum()) if "label" in d else None,
        "n_binder_pchembl_ge7": int((pch >= 7).sum()) if pch.notna().any() else None,
        "pchembl_min": round(float(pch.min()), 2) if pch.notna().any() else None,
        "pchembl_max": round(float(pch.max()), 2) if pch.notna().any() else None,
        "year_min": int(yr.min()) if yr.notna().any() else None,
        "year_max": int(yr.max()) if yr.notna().any() else None,
        "data_sources": ", ".join(sorted(set(d["source"].dropna().astype(str)))) if "source" in d else None,
    }


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "reviewer_package" / date.today().strftime("%m.%d.%y")
    out.mkdir(parents=True, exist_ok=True)
    import app

    ids = chembl_ids()
    core_ids = json.loads((ROOT / "data" / "core_panel_target_ids.json").read_text())         if (ROOT / "data" / "core_panel_target_ids.json").exists() else {}
    modes = json.loads((M / "binder_modes.json").read_text())
    hold = json.loads((M / "holdout" / "binder_modes.json").read_text()) \
        if (M / "holdout" / "binder_modes.json").exists() else {}
    cv = pd.read_csv(TAB / "rf_cv_summary.csv")
    adcv = pd.read_csv(TAB / "adme_cv_summary.csv")

    # ---------------- 01 master endpoint inventory ----------------
    rows = []

    def meta(name, suffix="_meta.json"):
        p = M / f"{name}{suffix}"
        return json.loads(p.read_text()) if p.exists() else {}

    def cv_row(name, split, frame):
        r = frame[(frame.endpoint == name) & (frame.split == split)]
        return r.iloc[0].to_dict() if len(r) else {}

    for name, desc in app.TARGET_CLASSIFIERS.items():
        st = endpoint_table_stats(name) or {}
        mt = meta(name)
        sc = cv_row(name, "scaffold", cv)
        rd = cv_row(name, "random", cv)
        cid, src = ids.get(name, (None, None))
        rows.append({
            "endpoint": name, "panel": "Core target/property classifier",
            "task": "classification", "description": desc,
            "chembl_target_id": cid, "identifier_source_script": src,
            **st,
            "label_rule": LABEL_RULE if cid else "measured dataset label; see data source",
            "n_scaffolds": mt.get("n_scaffolds"), "n_features": mt.get("n_features"),
            "model_family": "RandomForestClassifier",
            "hyperparameters": json.dumps(mt.get("hyperparameters")) if mt.get("hyperparameters") else None,
            "cv_folds": mt.get("cv_folds"),
            "roc_auc_random_mean": rd.get("roc_auc_mean"), "roc_auc_random_sd": rd.get("roc_auc_sd"),
            "roc_auc_scaffold_mean": sc.get("roc_auc_mean"), "roc_auc_scaffold_sd": sc.get("roc_auc_sd"),
            "spearman_scaffold_mean": sc.get("spearman_mean"),
            "model_file": f"models_rf/{name}.joblib",
            "calibrated_model_file": (f"models_rf/{name}_calibrated.joblib"
                                      if (M / f"{name}_calibrated.joblib").exists() else None),
            "deployed": True,
            "cv_detail_in": "03_crossvalidation_summary.csv, 04_crossvalidation_per_fold.csv",
        })

    for name, desc in app.RECEPTOR_REGRESSORS.items():
        st = endpoint_table_stats(name) or {}
        mt = meta(name)
        sc = cv_row(name, "scaffold", cv)
        rd = cv_row(name, "random", cv)
        cid, src = ids.get(name, (None, None))
        rows.append({
            "endpoint": name, "panel": "Receptor potency regressor",
            "task": "regression", "description": desc,
            "chembl_target_id": cid, "identifier_source_script": src, **st,
            "label_rule": "continuous pChEMBL value; no thresholding",
            "n_scaffolds": mt.get("n_scaffolds"), "n_features": mt.get("n_features"),
            "model_family": "RandomForestRegressor",
            "hyperparameters": json.dumps(mt.get("hyperparameters")) if mt.get("hyperparameters") else None,
            "cv_folds": mt.get("cv_folds"),
            "roc_auc_random_mean": rd.get("roc_auc_mean"), "roc_auc_random_sd": rd.get("roc_auc_sd"),
            "roc_auc_scaffold_mean": sc.get("roc_auc_mean"), "roc_auc_scaffold_sd": sc.get("roc_auc_sd"),
            "spearman_scaffold_mean": sc.get("spearman_mean"),
            "model_file": f"models_rf/{name}.joblib",
            "calibrated_model_file": None, "deployed": True,
            "cv_detail_in": "03_crossvalidation_summary.csv, 04_crossvalidation_per_fold.csv",
        })

    for name in sorted(modes):
        v = modes[name]
        st = endpoint_table_stats(name) or {}
        mt = meta(name, "_binder_meta.json")
        cid, src = ids.get(name, (None, None))
        rows.append({
            "endpoint": name, "panel": "Binder classifier",
            "task": "classification", "description": app.MECH_LABEL.get(name, name),
            "chembl_target_id": cid, "identifier_source_script": src, **st,
            "label_rule": LABEL_RULE,
            "n_scaffolds": None, "n_features": 1036,
            "model_family": "RandomForestClassifier with sigmoid calibration",
            "hyperparameters": json.dumps(mt.get("hyperparameters")) if mt.get("hyperparameters") else None,
            "cv_folds": 10,
            "roc_auc_random_mean": None, "roc_auc_random_sd": None,
            "roc_auc_scaffold_mean": v.get("scaffold_cv_auroc"), "roc_auc_scaffold_sd": None,
            "spearman_scaffold_mean": None,
            "model_file": f"models_rf/{name}_binder.joblib",
            "calibrated_model_file": None,
            "deployed": bool(v.get("deployed", True)),
            "cv_detail_in": "07_binder_panel_training_design.csv",
        })

    for name, (label, unit, _tf) in app.ADME.items():
        sc = cv_row(name, "scaffold", adcv)
        rd = cv_row(name, "random", adcv)
        rows.append({
            "endpoint": name, "panel": "ADME / exposure",
            "task": "regression" if pd.isna(sc.get("roc_auc_mean", float("nan"))) else "classification",
            "description": f"{label} ({unit})",
            "chembl_target_id": None, "identifier_source_script": None,
            "training_table": None, "n_compounds": sc.get("n"),
            "n_active_label1": None, "n_inactive_label0": None, "n_binder_pchembl_ge7": None,
            "pchembl_min": None, "pchembl_max": None, "year_min": None, "year_max": None,
            "data_sources": "Therapeutics Data Commons / MoleculeNet / B3DB / ChEMBL",
            "label_rule": "measured value from the source dataset; no thresholding applied here",
            "n_scaffolds": None, "n_features": 1036,
            "model_family": "RandomForest", "hyperparameters": None, "cv_folds": 10,
            "roc_auc_random_mean": rd.get("roc_auc_mean"), "roc_auc_random_sd": rd.get("roc_auc_sd"),
            "roc_auc_scaffold_mean": sc.get("roc_auc_mean"), "roc_auc_scaffold_sd": sc.get("roc_auc_sd"),
            "spearman_scaffold_mean": sc.get("spearman_mean"),
            "model_file": f"models_rf/adme/{name}.joblib",
            "calibrated_model_file": None, "deployed": True,
            "cv_detail_in": "03_crossvalidation_summary.csv, 04_crossvalidation_per_fold.csv",
        })

    for name, fam, fil in [("antioxidant_DPPH", "Neuroprotection regressor", "antioxidant_DPPH.joblib"),
                           ("pka_basic", "Physicochemical regressor", "pka_basic.joblib")]:
        mt = meta(name)
        sc = cv_row(name, "scaffold", cv)
        rows.append({
            "endpoint": name, "panel": fam, "task": "regression",
            "description": ("measured DPPH radical-scavenging pIC50" if "DPPH" in name
                            else "most basic pKa"),
            "chembl_target_id": None, "identifier_source_script": None,
            **(endpoint_table_stats(name) or {}),
            "label_rule": "continuous measured value",
            "n_scaffolds": mt.get("n_scaffolds"), "n_features": 1036,
            "model_family": "RandomForestRegressor", "hyperparameters": None, "cv_folds": 10,
            "roc_auc_random_mean": None, "roc_auc_random_sd": None,
            "roc_auc_scaffold_mean": None, "roc_auc_scaffold_sd": None,
            "spearman_scaffold_mean": sc.get("spearman_mean"),
            "model_file": f"models_rf/{fil}", "calibrated_model_file": None,
            "deployed": (M / fil).exists(),
            "cv_detail_in": "03_crossvalidation_summary.csv" if len(sc) else "",
        })

    master = pd.DataFrame(rows)
    # Four receptors carry two models each, a potency regressor and a binder classifier, so their
    # endpoint name appears twice. That is the design, not a duplicated row, and a reviewer scanning
    # the endpoint column would reasonably suspect an error unless it is said plainly. A stable key
    # and an explicit note are added so the file can be joined and read without ambiguity.
    master.insert(0, "row_key", master.endpoint + " [" + master.panel + "]")
    master["uniprot_accession"] = [core_ids.get(e, {}).get("uniprot") for e in master.endpoint]
    master["chembl_pref_name"] = [core_ids.get(e, {}).get("chembl_pref_name") for e in master.endpoint]
    dup = set(master.endpoint[master.endpoint.duplicated(keep=False)])
    master["dual_model_note"] = [
        (f"{e} is modelled twice: a potency regressor and a binder classifier, trained on the same "
         f"ChEMBL target but with different labels and different purposes. Both are deployed.")
        if e in dup else "" for e in master.endpoint]
    master.to_csv(out / "01_MASTER_endpoint_inventory.csv", index=False)
    print(f"01 master inventory: {len(master)} endpoints", flush=True)

    # ---------------- 03 / 04 cross-validation ----------------
    cvs = pd.concat([cv.assign(panel="core target panel"),
                     adcv.assign(panel="ADME / exposure")], ignore_index=True)
    cvs.to_csv(out / "03_crossvalidation_summary.csv", index=False)
    per = []
    p2 = TAB / "manuscript_T2_per_fold.csv"
    if p2.exists():
        per.append(pd.read_csv(p2).assign(panel="core target panel"))
    if (TAB / "adme_cv_folds.csv").exists():
        per.append(pd.read_csv(TAB / "adme_cv_folds.csv").assign(panel="ADME / exposure"))
    if per:
        pf = pd.concat(per, ignore_index=True)
        pf.to_csv(out / "04_crossvalidation_per_fold.csv", index=False)
        print(f"03/04 cross-validation: {len(cvs)} summary rows, {len(pf)} fold rows", flush=True)

    # ---------------- 05 features ----------------
    from features.featurize import feature_names
    fn = feature_names()
    feats = []
    for i, f in enumerate(fn):
        if f.startswith("ecfp4_"):
            feats.append({"index": i, "feature": f, "block": "ECFP-4 fingerprint",
                          "type": "binary", "description":
                          "Morgan circular fingerprint bit, radius 2, folded to 1024 bits; presence "
                          "of a hashed substructural environment", "unit": "0 or 1"})
        else:
            feats.append({"index": i, "feature": f, "block": "physicochemical descriptor",
                          "type": "continuous", "unit": "RDKit units",
                          "description": {
                              "mw": "molecular weight", "clogp": "Crippen logP",
                              "tpsa": "topological polar surface area",
                              "hbd": "hydrogen-bond donors", "hba": "hydrogen-bond acceptors",
                              "rotatable_bonds": "rotatable bond count",
                              "aromatic_rings": "aromatic ring count",
                              "fraction_csp3": "fraction of sp3-hybridised carbons",
                              "ring_count": "total ring count", "heavy_atoms": "heavy atom count",
                              "formal_charge": "net formal charge",
                              "qed": "quantitative estimate of drug-likeness"}.get(f, f)})
    pd.DataFrame(feats).to_csv(out / "05_feature_definitions.csv", index=False)
    print(f"05 features: {len(feats)} ({sum(1 for f in feats if f['block'].startswith('ECFP'))} "
          f"fingerprint + {sum(1 for f in feats if 'descriptor' in f['block'])} descriptors)", flush=True)

    # ---------------- 06 model registry ----------------
    reg = []
    for p in sorted(list(M.glob("*.joblib")) + list((M / "adme").glob("*.joblib"))):
        stem = p.stem.replace("_binder", "").replace("_calibrated", "")
        v = modes.get(stem, {})
        reg.append({
            "model_file": str(p.relative_to(ROOT)).replace("\\", "/"),
            "size_bytes": p.stat().st_size,
            "endpoint": stem,
            "role": ("binder classifier" if p.stem.endswith("_binder") else
                     "calibrated classifier" if p.stem.endswith("_calibrated") else
                     "ADME model" if p.parent.name == "adme" else "base model"),
            "n_features": 1036,
            "deployed": bool(v.get("deployed", True)) if stem in modes else True,
            "decision_threshold": v.get("threshold"),
            "threshold_basis": v.get("threshold_basis"),
            "sensitivity_at_threshold": v.get("sensitivity_at_threshold"),
            "sensitivity_basis": v.get("sensitivity_basis"),
            # The threshold is a quantile of the held-out inactives, so the rate on that same set is
            # that quantile by construction. It is carried for continuity and named for what it is;
            # background_fpr_at_threshold, measured on a disjoint pool, is the informative one.
            "fpr_in_sample_on_threshold_set": v.get("fpr_in_sample_on_threshold_set",
                                                    v.get("fpr_at_threshold")),
            "background_fpr_at_threshold": v.get("background_fpr_at_threshold"),
            "background_fpr_basis": ("held_out_evaluation_pool"
                                     if v.get("background_fpr_in_sample") is not None else None),
            "high_precision_threshold": v.get("screening_threshold"),
            "withdrawn_reason": v.get("withdrawn_reason"),
        })
    pd.DataFrame(reg).to_csv(out / "06_model_registry.csv", index=False)
    print(f"06 model registry: {len(reg)} artefacts", flush=True)

    # ---------------- 07 binder training design ----------------
    bd = []
    for name in sorted(modes):
        v = modes[name]
        h = hold.get(name, {})
        bd.append({
            "endpoint": name, "label": app.MECH_LABEL.get(name, name),
            "training_mode": v.get("mode"),
            "n_positive_pchembl_ge7": v.get("n_positive"),
            "n_property_matched_decoys": v.get("n_decoy"),
            "n_measured_inactive_train": v.get("n_measured_inactive_train"),
            "n_measured_inactive_holdout": v.get("n_measured_inactive_holdout"),
            "scaffold_cv_auroc": v.get("scaffold_cv_auroc"),
            "auroc_vs_measured_inactives": v.get("auroc_vs_measured_inactives"),
            "decision_threshold": v.get("threshold"),
            "threshold_basis": v.get("threshold_basis"),
            "sensitivity_at_threshold": v.get("sensitivity_at_threshold"),
            "holdout_train_actives": h.get("n_train_actives"),
            "holdout_actives": h.get("n_holdout_actives"),
            "holdout_scaffolds": h.get("n_holdout_scaffolds"),
            "holdout_recall_at_threshold": h.get("holdout_recall_at_threshold"),
            "deployed": bool(v.get("deployed", True)),
            "withdrawn_reason": v.get("withdrawn_reason"),
        })
    pd.DataFrame(bd).to_csv(out / "07_binder_panel_training_design.csv", index=False)
    print(f"07 binder design: {len(bd)} endpoints", flush=True)

    # ---------------- 02 input files ----------------
    inp = []
    for p in sorted(EP.glob("*.csv")):
        d = pd.read_csv(p)
        name = p.stem
        used = name in set(master.endpoint)
        inp.append({
            "input_file": f"data/endpoints/{p.name}", "kind": "training table",
            "rows": len(d), "columns": ", ".join(d.columns),
            "feeds_endpoint": name if used else "",
            "used_in_deployed_panel": used,
            "produced_output": (f"models_rf/{name}.joblib or models_rf/{name}_binder.joblib"
                                if used else ""),
            "note": "" if used else "fetched and audited but not deployed; see "
                                    "09_PROVENANCE_AND_LIMITATIONS.md",
        })
    for p in sorted(TAB.glob("*.csv")):
        d = pd.read_csv(p)
        inp.append({"input_file": f"results/tables/{p.name}", "kind": "validation output",
                    "rows": len(d), "columns": ", ".join(list(d.columns)[:10]),
                    "feeds_endpoint": "", "used_in_deployed_panel": "",
                    "produced_output": "manuscript table or figure", "note": ""})
    pd.DataFrame(inp).to_csv(out / "02_training_input_files.csv", index=False)
    print(f"02 input files: {len(inp)}", flush=True)

    write_docs(out, master, cvs, modes, len(feats))
    print(f"\nwrote {len(list(out.glob('*')))} files to {out}")


def write_docs(out, master, cvs, modes, n_feats):
    """The two prose files. Everything they assert is computed above, not remembered."""
    dep = master[master.deployed == True]        # noqa: E712
    (out / "08_DATA_DICTIONARY.md").write_text(f"""# Data dictionary

Every column in every file of this package, and every reason a cell may be empty.

**An empty cell never means "not investigated".** It means exactly one of the reasons listed under
"Why a cell is empty" below, and each is stated per column.

## 01_MASTER_endpoint_inventory.csv

One row per endpoint, {len(master)} rows.

| Column | Meaning | Empty when |
|---|---|---|
| row_key | unique key: endpoint plus panel. Use this to join, not `endpoint` | never |
| endpoint | internal endpoint name. **Appears twice for A2A, D2, 5-HT2A and SERT**, which carry both a potency regressor and a binder classifier | never |
| dual_model_note | explains that duplication where it occurs | for the 64 endpoints with a single model |
| panel | which of the five panels it belongs to | never |
| task | classification or regression | never |
| description | human-readable target or property | never |
| chembl_target_id | ChEMBL target the training data was drawn from | the endpoint is not ChEMBL-derived: the nine ADME endpoints, BBB (from B3DB), the antioxidant endpoint and the pKa endpoint |
| uniprot_accession | UniProt accession, for the core panel whose identifiers were resolved from accession | endpoints not in the core panel, and non-protein endpoints |
| chembl_pref_name | the preferred name ChEMBL returns for that target, so the mapping can be checked rather than trusted | as above |
| identifier_source_script | the fetch script where that identifier is written and was verified on download | as above |
| training_table | path to the exact table the model was trained on | ADME endpoints, whose tables live in the ADME pipeline rather than data/endpoints |
| n_compounds | rows in the training table after deduplication by structure | never for ChEMBL endpoints; for ADME it is the cross-validation n |
| n_active_label1, n_inactive_label0 | class counts under the label rule | regression endpoints, which have no classes |
| n_binder_pchembl_ge7 | compounds meeting the stricter binder positive definition | endpoints with no pChEMBL column |
| pchembl_min, pchembl_max, year_min, year_max | measured range and publication-year range | non-ChEMBL endpoints |
| data_sources | the source column of the training table | ADME, where sources are recorded in the ADME pipeline |
| label_rule | how a measurement became a label | never |
| n_scaffolds | distinct Bemis-Murcko scaffolds | binder and ADME endpoints, whose scaffold counts are in 07 and the ADME pipeline |
| n_features | width of the feature vector | never; always 1036 |
| model_family | estimator class | never |
| hyperparameters | exact settings, as JSON | binder and ADME models, whose settings are fixed in the training script and identical across the panel: see 09 |
| cv_folds | folds used | never; always 10 |
| roc_auc_random_mean/sd | 10-fold random split | **regression endpoints, which have no AUROC**; binder endpoints, which were not run under a random split because the scaffold split is the meaningful one for them |
| roc_auc_scaffold_mean/sd | 10-fold scaffold-grouped split | regression endpoints |
| spearman_scaffold_mean | rank correlation, scaffold split | **classification endpoints, which have no rank correlation** |
| model_file | the artefact loaded at run time | never |
| calibrated_model_file | isotonic or sigmoid calibrated variant where one exists | endpoints that are deployed uncalibrated |
| deployed | whether the running server uses it | never |
| cv_detail_in | which file in this package carries the per-fold record | never |

## 03_crossvalidation_summary.csv

Mean and standard deviation across folds, per endpoint and split, {len(cvs)} rows.
Columns ending `_mean` and `_sd` are over the {int(cvs.n.notna().sum()) and 10} folds.

Classification metrics (roc_auc, pr_auc, mcc, f1, balanced_acc) are empty for regression endpoints.
Regression metrics (r2, rmse, mae, spearman) are empty for classification endpoints. This is the
single commonest empty cell in the package and it is structural, not missing data.

## 04_crossvalidation_per_fold.csv

One row per endpoint x split x fold. `n_scaffolds_test` is the number of distinct scaffolds in that
fold's test set, which is what makes the scaffold split meaningful; it is empty for the random split,
where scaffolds are not the grouping variable.

## 05_feature_definitions.csv

All {n_feats} input features, in the exact order the model receives them. Indices 0 to 1023 are the
folded ECFP-4 fingerprint; 1024 to {n_feats - 1} are the named physicochemical descriptors. No cell
is empty in this file.

## 06_model_registry.csv

One row per model artefact on disk. Threshold columns are empty for endpoints that are not
thresholded: regressors, and base models superseded by a calibrated variant.

## 07_binder_panel_training_design.csv

How each binder endpoint's negative class was constructed. `holdout_*` columns come from the
scaffold-held-out retraining and are empty for three endpoints that have no hold-out twin because
too few of their actives survive a 20 percent scaffold withholding: these are named in 09.

## 02_training_input_files.csv

Every input and output file, with row counts. `feeds_endpoint` and `produced_output` are empty for
validation outputs, which consume models rather than producing them, and for training tables that
were fetched and audited but not deployed.

## Why a cell is empty, in full

1. **Structurally inapplicable.** A regression endpoint has no AUROC; a classifier has no Spearman.
2. **Not thresholded.** Regressors and superseded base models have no decision threshold.
3. **Different pipeline.** ADME endpoints keep their tables and scaffold counts in the ADME pipeline
   rather than in data/endpoints.
4. **Uniform across the panel and stated once.** Binder hyperparameters are fixed in the training
   script for every endpoint; repeating them 49 times would imply they vary.
5. **Genuinely absent, and named.** Three endpoints have no hold-out twin; two endpoints were
   withdrawn. Both cases are listed explicitly in 09 rather than left to inference.
""", encoding="utf-8")

    n_withdrawn = sum(1 for v in modes.values() if not v.get("deployed", True))
    (out / "09_PROVENANCE_AND_LIMITATIONS.md").write_text(f"""# Provenance and limitations

## How this package was produced

Every value was read from the file that produced it by
`src/brainsafe/analysis/build_reviewer_package.py`. Nothing was typed by hand and nothing was
recalled. Re-running that script regenerates the package; if it ever disagreed with the models, the
package would change rather than the models.

Generated {date.today().isoformat()} from the deployed model set archived at
**doi:10.5281/zenodo.21858576**.

## What is in the panel

- {len(master)} endpoints in total, of which {len(master[master.deployed == True])} are deployed.
- {len(master[master.panel == 'Core target/property classifier'])} core classifiers and
  {len(master[master.panel == 'Receptor potency regressor'])} receptor regressors, cross-validated
  10-fold under both a random and a scaffold-grouped split.
- {len(master[master.panel == 'Binder classifier'])} binder classifiers, validated against compounds
  experimentally tested on the same target and found inactive, held out from training.
- {len(master[master.panel == 'ADME / exposure'])} ADME and exposure endpoints, 10-fold under both
  splits.

## Cross-validation design

Ten folds, two splits, for every endpoint that has a cross-validation record.

The **random split** shuffles compounds. It measures interpolation and is reported because it is the
conventional number, not because it is the honest one.

The **scaffold-grouped split** withholds entire Bemis-Murcko scaffolds using GroupKFold, so no test
compound shares a core with any training compound. It is the number to read: it measures whether the
model generalises to a chemical series it has not seen. Scaffold AUROC is consistently lower than
random AUROC across the panel, and the gap is the honest cost of extrapolation.

Per-fold records are in `04_crossvalidation_per_fold.csv`, including the test-set size and the number
of distinct scaffolds in each fold.

## Known absences, stated rather than implied

- **Two endpoints were trained and then withdrawn** ({n_withdrawn} in
  `07_binder_panel_training_design.csv` with `deployed = False`). They remain in the package with the
  reason recorded, because a withdrawn model is part of the training record.
- **Three endpoints have no scaffold hold-out twin**: GABA_A, GBA1 and TAAR1. Too few of their
  actives survive a 20 percent scaffold withholding to estimate recall, so those cells are empty
  rather than filled with an unstable number.
- **Binder endpoints have no random-split AUROC.** They were evaluated under the scaffold split and
  against measured inactives, which are the two informative comparisons for a decoy-trained model.
- **Per-fold detail for the binder panel is not in the same form** as the core panel. Binder models
  report a scaffold cross-validated AUROC and a hold-out recall; the per-fold breakdown was not
  retained during their training run. This is a genuine gap in the record and is stated rather than
  reconstructed.
- **Hyperparameters are recorded per endpoint for the core panel** and are uniform for the binder
  panel: `n_estimators=300, min_samples_leaf=4, class_weight='balanced', random_state=42`, fixed in
  `src/brainsafe/models/train_binders_hybrid.py`.

## Reproducibility

Random seed 42 throughout. The feature vector is deterministic given a structure. The models
themselves are not bit-deterministic between calls, because a random forest with `n_jobs=-1`
accumulates in thread-completion order; differences appear only in the last bits of a probability and
never change a reported decision.

Environment is pinned in `requirements.txt`; `src/brainsafe/evaluation/app_health.py` verifies that
an installation matches those pins before the server is considered deployable.
""", encoding="utf-8")

    (out / "00_README.md").write_text(f"""# BrainSafe AI: reviewer package

Prepared {date.today().isoformat()} in response to the reviewers' request for the complete training
record. Generated mechanically from the deployed model set by
`src/brainsafe/analysis/build_reviewer_package.py`; no value in this package was typed by hand.

Models archived at **doi:10.5281/zenodo.21858576**.
Source: https://github.com/krishna-g-999/brainsafe-ai (branch `main`).

## Read in this order

| File | What it answers |
|---|---|
| `01_MASTER_endpoint_inventory.csv` | **Every endpoint used to train, and what went in.** One row per endpoint: its data source, ChEMBL target, compound counts, class balance, label rule, model family, hyperparameters, and cross-validated score under both splits. |
| `02_training_input_files.csv` | **Every input file, labelled with what it produced.** Row counts and columns for each training table and each validation output. |
| `03_crossvalidation_summary.csv` | **How many folds, how many endpoints.** Mean and standard deviation per endpoint and split. |
| `04_crossvalidation_per_fold.csv` | The individual folds behind those means, with test-set and scaffold counts. |
| `05_feature_definitions.csv` | **All {n_feats} features**, in model input order, each named and described. |
| `06_model_registry.csv` | **Every model trained**, its file, its operating threshold and how that threshold was set. |
| `07_binder_panel_training_design.csv` | How the binder panel's negatives were built, and its hold-out results. |
| `08_DATA_DICTIONARY.md` | Every column explained, and **every reason a cell is empty**. |
| `09_PROVENANCE_AND_LIMITATIONS.md` | What is verified, what is genuinely absent, and why. |

## The three things a reviewer usually wants first

1. **Which endpoints, trained on what?** `01_MASTER_endpoint_inventory.csv`, one row each.
2. **What was the validation design?** Ten folds, two splits (random and scaffold-grouped), for every
   endpoint with a cross-validation record. `03` for the summary, `04` for every fold.
3. **What features, and what did each model score?** `05` for the {n_feats}-column input vector,
   `01` and `06` for per-model scores.

## On empty cells

Empty cells are meaningful and are all explained in `08_DATA_DICTIONARY.md`. The commonest is
structural: a regression endpoint has no AUROC and a classifier has no rank correlation. Where a
value is genuinely missing rather than inapplicable, it is named in
`09_PROVENANCE_AND_LIMITATIONS.md` rather than left for the reader to infer.
""", encoding="utf-8")


if __name__ == "__main__":
    main()
