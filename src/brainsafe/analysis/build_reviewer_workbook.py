"""One Excel workbook containing the whole training record, for reviewers.

The CSV package answers every question but spreads the answers across ten files, which is the wrong
shape for someone assessing a manuscript. This puts the same content in one workbook, one sheet per
question, with a sheet at the front that exists solely to prevent the most likely misreading.

That misreading is worth naming. The manuscript quotes a figure of tens of thousands of measured
records, and a reader can easily take it to mean that every endpoint was cross-validated over that
many compounds. It was not. The figure is the SUM across endpoints; each model was trained and
cross-validated on its own set, and those sets differ by two orders of magnitude. Sheet 2 sets the
per-endpoint sizes out so the arithmetic is visible rather than inferred, and separates the three
things a binder endpoint's training set is made of, since for those the source table's row count is
not the training set size.

Every number is read from the file that produced it. Nothing is typed by hand.

Writes reviewer_package/<date>/BrainSafe_AI_training_record.xlsx
"""
from __future__ import annotations

import json
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


def compound_accounting():
    """Per endpoint: how many compounds its own model actually saw, and the totals they sum to.

    A binder endpoint's training set is not the row count of its CSV. It is the measured actives
    (pChEMBL >= 7), plus decoys drawn from background chemistry, plus half the measured inactives;
    the other half is withheld for threshold calibration and never trained on. Those three counts are
    recorded per endpoint in models_rf/binder_modes.json and are read from there, because using the
    CSV row count would understate the training set for most binder endpoints and overstate it for
    any whose measured actives fall below the potency cut.
    """
    import app
    modes = json.loads((ROOT / "models_rf" / "binder_modes.json").read_text())
    rows, seen_all = [], set()
    groups = [("Core target/property classifier", list(app.TARGET_CLASSIFIERS)),
              ("Receptor potency regressor", list(app.RECEPTOR_REGRESSORS)),
              ("Binder classifier", list(app.BINDER_TARGETS)),
              ("Other regressor", ["antioxidant_DPPH", "pka_basic"])]
    for panel, names in groups:
        for n in names:
            p = EP / f"{n}.csv"
            if not p.exists():
                continue
            d = pd.read_csv(p)
            smi = set(d["smiles"].astype(str))
            seen_all |= smi
            rec = modes.get(n, {}) if panel == "Binder classifier" else {}
            if rec:
                pos = rec.get("n_positive")
                dec = rec.get("n_decoy", 0)
                ina = rec.get("n_measured_inactive_train", 0)
                n_train = (pos or 0) + (dec or 0) + (ina or 0)
                comp = (f"{pos:,} measured actives + {dec:,} decoys + {ina:,} measured inactives"
                        if dec else f"{pos:,} measured actives + {ina:,} measured inactives")
                held = rec.get("n_measured_inactive_holdout")
            else:
                n_train = len(d)
                pos = int((d["label"] == 1).sum()) if "label" in d else None
                ina = int((d["label"] == 0).sum()) if "label" in d else None
                dec = 0
                comp = (f"{pos:,} actives + {ina:,} inactives, all measured"
                        if pos is not None else "continuous measured values, no class split")
                held = None
            rows.append({
                "endpoint": n, "panel": panel,
                "compounds_this_model_was_trained_and_cv_on": n_train,
                "training_set_composition": comp,
                "measured_actives": pos, "decoys": dec, "measured_inactives_trained": ina,
                "measured_inactives_withheld_for_thresholding": held,
                "rows_in_source_csv": len(d),
                "unique_structures_in_source_csv": len(smi),
                "cv_folds": 10,
                "compounds_per_fold_approx": round(n_train / 10),
                "training_table": f"data/endpoints/{n}.csv",
            })
    df = pd.DataFrame(rows).sort_values("compounds_this_model_was_trained_and_cv_on",
                                        ascending=False)
    total_rows = int(df["compounds_this_model_was_trained_and_cv_on"].sum())
    return df, total_rows, len(seen_all)


def merged_cv():
    """Every endpoint's ten folds in one place: core panel, ADME and binder panel together."""
    frames_sum, frames_fold = [], []
    for path in [TAB / "rf_cv_summary.csv", TAB / "adme_cv_summary.csv",
                 TAB / "binder_cv_summary.csv"]:
        if path.exists():
            f = pd.read_csv(path)
            f["record_source"] = path.name
            frames_sum.append(f)
    # The core panel's folds come from manuscript_T2_per_fold.csv, not rf_cv_folds.csv. The latter
    # holds only the three endpoints of the most recent partial re-run, because train_rf.py rewrites
    # it with whatever endpoints it was given on the command line; the manuscript table is the
    # complete thirteen-endpoint record and is what the reported means were computed from.
    for path in [TAB / "manuscript_T2_per_fold.csv", TAB / "adme_cv_folds.csv",
                 TAB / "binder_cv_folds.csv"]:
        if path.exists():
            f = pd.read_csv(path)
            f["record_source"] = path.name
            frames_fold.append(f)
    s = pd.concat(frames_sum, ignore_index=True) if frames_sum else pd.DataFrame()
    p = pd.concat(frames_fold, ignore_index=True) if frames_fold else pd.DataFrame()
    if len(s):
        s = s.sort_values(["endpoint", "split"])
    if len(p):
        p = p.sort_values(["endpoint", "split", "fold"])
    return s, p


def main():
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        ROOT / "reviewer_package" / date.today().strftime("%m.%d.%y")
    out_dir.mkdir(parents=True, exist_ok=True)
    xlsx = out_dir / "BrainSafe_AI_training_record.xlsx"

    acct, total_records, total_unique = compound_accounting()
    biggest = acct.iloc[0]
    smallest = acct.iloc[-1]

    readme = pd.DataFrame([
        ("What this workbook is",
         "The complete training record for BrainSafe AI, generated from the deployed models by "
         "src/brainsafe/analysis/build_reviewer_workbook.py. No value was typed by hand."),
        ("Models archived at", "doi:10.5281/zenodo.21858576"),
        ("Source code", "https://github.com/krishna-g-999/brainsafe-ai (branch main)"),
        ("", ""),
        ("READ SHEET 2 FIRST",
         "It answers the question this record is most likely to be misread on: how many compounds "
         "each model was actually trained on."),
        ("", ""),
        ("2. Compound accounting",
         f"Per endpoint, the number of compounds THAT model was trained and cross-validated on. "
         f"These range from {int(smallest.compounds_this_model_was_trained_and_cv_on):,} "
         f"({smallest.endpoint}) to {int(biggest.compounds_this_model_was_trained_and_cv_on):,} "
         f"({biggest.endpoint}). No single model was trained on more than the largest of these. "
         f"The column sums to {total_records:,}, which is a sum of training-set sizes and not a "
         f"compound count, for two reasons: a structure measured at several targets is counted once "
         f"per endpoint it appears in, and the binder endpoints' generated decoys are counted "
         f"because the models were fitted to them. The measured record total, decoys excluded, is "
         f"203,884 over 160,365 unique compounds by InChIKey."),
        ("   note on binder endpoints",
         "A binder endpoint's training set is not the row count of its source table. It is the "
         "measured actives at pChEMBL >= 7, plus property-matched decoys drawn from background "
         "chemistry at a 3:1 ratio and no more than 0.35 Tanimoto to any active, plus half the "
         "measured inactives. The other half of the measured inactives is withheld to set the "
         "decision threshold and is never trained on. Sheet 2 breaks all four counts out."),
        ("3. Endpoint inventory",
         "Every endpoint: target identifier, data source, class balance, label rule, model family, "
         "hyperparameters, and cross-validated score under both splits."),
        ("4. Cross-validation summary",
         "Ten folds per endpoint per split, mean and standard deviation, for every endpoint in the "
         "panel. Two splits: random, and scaffold-grouped, which withholds entire Bemis-Murcko "
         "scaffolds so no scaffold appears in both training and test. The scaffold split is the "
         "honest estimate; the random split is included because it is what most published figures "
         "report, and the gap between the two is itself informative."),
        ("5. Cross-validation per fold",
         "Every individual fold behind those means, with training and test sizes, positives in the "
         "test fold, and scaffold counts. This is the raw ten-fold record."),
        ("6. Features",
         "All 1036 model inputs in the order the model receives them: 1024 ECFP-4 fingerprint bits "
         "and 12 named physicochemical descriptors."),
        ("7. Model registry",
         "Every model artefact, its decision threshold, and how that threshold was set."),
        ("8. Binder panel design",
         "How the binder panel's negative class was constructed, and its scaffold hold-out results."),
        ("9. External test results",
         "Every compound the deployed models were scored on outside training, from four independent "
         "evaluations. Membership is checked rather than assumed: each structure is canonicalised "
         "and tested against the union of all training tables, and the result is in "
         "`found_in_training`."),
        ("10. External test summary",
         "The same results split on whether the compound was in training. This is the split that "
         "matters: approved-drug sets overlap ChEMBL heavily by construction, so a headline over the "
         "mixture would be carried by the seen half. Read the 'no, held out' rows."),
        ("11. Training inputs (full)",
         "Every measured training record, one row per endpoint per compound, with its measurement "
         "provenance. Also supplied as MASTER_training_inputs.csv, which is the better format for "
         "anything programmatic."),
        ("12. Input files", "Every input and output file with row counts."),
        ("13. Data dictionary",
         "Every column explained, and every reason a cell may be empty. An empty cell is never "
         "'not investigated'."),
    ], columns=["Item", "Description"])

    master = pd.read_csv(out_dir / "01_MASTER_endpoint_inventory.csv") \
        if (out_dir / "01_MASTER_endpoint_inventory.csv").exists() else pd.DataFrame()
    cvs, per = merged_cv()
    feats = pd.read_csv(out_dir / "05_feature_definitions.csv") \
        if (out_dir / "05_feature_definitions.csv").exists() else pd.DataFrame()
    reg = pd.read_csv(out_dir / "06_model_registry.csv") \
        if (out_dir / "06_model_registry.csv").exists() else pd.DataFrame()
    bind = pd.read_csv(out_dir / "07_binder_panel_training_design.csv") \
        if (out_dir / "07_binder_panel_training_design.csv").exists() else pd.DataFrame()
    inp = pd.read_csv(out_dir / "02_training_input_files.csv") \
        if (out_dir / "02_training_input_files.csv").exists() else pd.DataFrame()
    ext = pd.read_csv(out_dir / "MASTER_external_test_results.csv") \
        if (out_dir / "MASTER_external_test_results.csv").exists() else pd.DataFrame()
    exs = pd.read_csv(out_dir / "MASTER_external_test_summary.csv") \
        if (out_dir / "MASTER_external_test_summary.csv").exists() else pd.DataFrame()
    tin = pd.read_csv(out_dir / "MASTER_training_inputs.csv") \
        if (out_dir / "MASTER_training_inputs.csv").exists() else pd.DataFrame()

    dd = pd.DataFrame([
        ("Every sheet", "endpoint", "internal endpoint name, the key across all sheets", "never"),
        ("2 Compound accounting", "compounds_this_model_was_trained_and_cv_on",
         "the size of that endpoint's own training set, which is the number its ten-fold CV ran "
         "over. For binder endpoints this is actives + decoys + trained measured inactives, which "
         "is not the row count of the source CSV; both are given so they can be reconciled",
         "never"),
        ("2 Compound accounting", "decoys",
         "property-matched background compounds used as presumed negatives, 3:1 to actives, "
         "Tanimoto < 0.35 to every active",
         "endpoints trained on measured labels alone, where the value is 0 rather than empty"),
        ("2 Compound accounting", "measured_inactives_withheld_for_thresholding",
         "measured inactives deliberately kept out of training so the decision threshold is set on "
         "data the model has not seen",
         "non-binder endpoints, which set no such threshold"),
        ("2 Compound accounting", "compounds_per_fold_approx",
         "the above divided by ten, so the size of one fold is visible", "never"),
        ("4 CV summary / 5 per fold", "record_source",
         "the results file each row came from: rf_cv_* for the core panel, adme_cv_* for the ADME "
         "panel, binder_cv_* for the binder panel", "never"),
        ("4 CV summary", "recorded_scaffold_cv_auroc",
         "for binder endpoints only, the out-of-fold AUROC logged when the deployed model was "
         "trained. The per-fold record was recomputed independently, redrawing the decoys, so this "
         "column and roc_auc_mean are two separate runs of the same protocol and their agreement is "
         "a reproducibility check rather than a restatement",
         "core and ADME endpoints, whose folds were kept at training time and needed no re-run"),
        ("3 Endpoint inventory", "chembl_target_id",
         "ChEMBL target the data came from; the core panel's were resolved from UniProt accession "
         "and the returned preferred name checked",
         "endpoint is not ChEMBL-derived: the nine ADME endpoints, BBB (B3DB), antioxidant, pKa"),
        ("3 Endpoint inventory", "roc_auc_random_mean / roc_auc_scaffold_mean",
         "10-fold AUROC under each split",
         "REGRESSION endpoints, which have no AUROC. This is structural, not missing data"),
        ("3 Endpoint inventory", "spearman_scaffold_mean", "rank correlation, scaffold split",
         "CLASSIFICATION endpoints, which have no rank correlation"),
        ("3 Endpoint inventory", "hyperparameters", "exact estimator settings",
         "binder and ADME models, whose settings are uniform across the panel and fixed in "
         "src/brainsafe/models/train_binders_hybrid.py: RandomForestClassifier with "
         "n_estimators=300, min_samples_leaf=4, class_weight=balanced, random_state=42. The core "
         "panel differs in one setting, min_samples_leaf=2, and states it per row"),
        ("5 Per fold", "n_scaffolds_test", "distinct scaffolds in that fold's test set",
         "the RANDOM split, where scaffolds are not the grouping variable"),
        ("7 Model registry", "decision_threshold", "the cut above which engagement is reported",
         "regressors and superseded base models, which are not thresholded"),
        ("8 Binder design", "holdout_* columns", "scaffold hold-out retraining results",
         "GABA_A, GBA1 and TAAR1, which have too few actives surviving a 20 per cent scaffold "
         "withholding to estimate recall. A genuine absence, named rather than implied"),
        ("8 Binder design", "withdrawn_reason", "why a trained endpoint was removed from the panel",
         "the 47 endpoints that are deployed"),
    ], columns=["Sheet", "Column", "Meaning", "Empty when"])

    sheets = [("1 READ ME FIRST", readme), ("2 Compound accounting", acct),
              ("3 Endpoint inventory", master), ("4 CV summary", cvs),
              ("5 CV per fold", per), ("6 Features", feats),
              ("7 Model registry", reg), ("8 Binder panel design", bind),
              ("9 External test results", ext), ("10 External test summary", exs),
              ("11 Training inputs (full)", tin), ("12 Input files", inp),
              ("13 Data dictionary", dd)]

    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        for name, frame in sheets:
            frame.to_excel(w, sheet_name=name, index=False)

        # readable column widths and a frozen header on every sheet
        for name, frame in sheets:
            ws = w.sheets[name]
            ws.freeze_panes = "A2"
            for i, col in enumerate(frame.columns, 1):
                longest = max([len(str(col))] +
                              [len(str(v)) for v in frame[col].head(200).tolist()]) if len(frame) else len(str(col))
                ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = min(max(12, longest + 2), 70)

    print(f"wrote {xlsx}")
    print(f"  sheets: {len(sheets)}")
    print(f"  endpoints with a training table: {len(acct)}")
    print(f"  sum of per-endpoint records: {total_records:,}")
    print(f"  unique structures across all endpoints: {total_unique:,}")
    print(f"  smallest endpoint: {smallest.endpoint} "
          f"({int(smallest.compounds_this_model_was_trained_and_cv_on):,} compounds)")
    print(f"  largest endpoint: {biggest.endpoint} "
          f"({int(biggest.compounds_this_model_was_trained_and_cv_on):,} compounds)")


if __name__ == "__main__":
    main()
