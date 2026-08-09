# Data dictionary

Every column in every file of this package, and every reason a cell may be empty.

**An empty cell never means "not investigated".** It means exactly one of the reasons listed under
"Why a cell is empty" below, and each is stated per column.

## 01_MASTER_endpoint_inventory.csv

One row per endpoint, 72 rows.

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

Mean and standard deviation across folds, per endpoint and split, 44 rows.
Columns ending `_mean` and `_sd` are over the 10 folds.

Classification metrics (roc_auc, pr_auc, mcc, f1, balanced_acc) are empty for regression endpoints.
Regression metrics (r2, rmse, mae, spearman) are empty for classification endpoints. This is the
single commonest empty cell in the package and it is structural, not missing data.

## 04_crossvalidation_per_fold.csv

One row per endpoint x split x fold. `n_scaffolds_test` is the number of distinct scaffolds in that
fold's test set, which is what makes the scaffold split meaningful; it is empty for the random split,
where scaffolds are not the grouping variable.

## 05_feature_definitions.csv

All 1036 input features, in the exact order the model receives them. Indices 0 to 1023 are the
folded ECFP-4 fingerprint; 1024 to 1035 are the named physicochemical descriptors. No cell
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
