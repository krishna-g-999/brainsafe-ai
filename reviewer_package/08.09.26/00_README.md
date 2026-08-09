# BrainSafe AI: reviewer package

Prepared 2026-08-09 in response to the reviewers' request for the complete training
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
| `05_feature_definitions.csv` | **All 1036 features**, in model input order, each named and described. |
| `06_model_registry.csv` | **Every model trained**, its file, its operating threshold and how that threshold was set. |
| `07_binder_panel_training_design.csv` | How the binder panel's negatives were built, and its hold-out results. |
| `08_DATA_DICTIONARY.md` | Every column explained, and **every reason a cell is empty**. |
| `09_PROVENANCE_AND_LIMITATIONS.md` | What is verified, what is genuinely absent, and why. |

## The three things a reviewer usually wants first

1. **Which endpoints, trained on what?** `01_MASTER_endpoint_inventory.csv`, one row each.
2. **What was the validation design?** Ten folds, two splits (random and scaffold-grouped), for every
   endpoint with a cross-validation record. `03` for the summary, `04` for every fold.
3. **What features, and what did each model score?** `05` for the 1036-column input vector,
   `01` and `06` for per-model scores.

## On empty cells

Empty cells are meaningful and are all explained in `08_DATA_DICTIONARY.md`. The commonest is
structural: a regression endpoint has no AUROC and a classifier has no rank correlation. Where a
value is genuinely missing rather than inapplicable, it is named in
`09_PROVENANCE_AND_LIMITATIONS.md` rather than left for the reader to infer.
