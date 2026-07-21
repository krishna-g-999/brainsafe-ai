# Data manifest

Every data file behind the models, what it contains, and how it was made. All values are measured
(ChEMBL 37, BindingDB, B3DB, ChEMBL DPPH); nothing is imputed. Regenerate any file by running the
script named beside it.

## Master table

- **`data/processed/compound_library.csv`** - the master. One row per unique compound (61,317),
  keyed by InChIKey, with every measured endpoint label and value it has, interpretable descriptors,
  a flavonoid-core flag, and the contributing data sources. Built by
  `src/brainsafe/data/build_compound_library.py`. Column definitions: `docs/DATA_DICTIONARY.md`.
- **`data/processed/endpoint_labels_long.csv`** - the same measurements in long form (one row per
  compound-endpoint measurement) for auditing.

## Per-endpoint tables (the training/testing values for each model)

- **`data/endpoints/<TARGET>.csv`** - one file per target (AChE, BChE, BACE1, GSK3B, MAO_A, MAO_B,
  D2, A2A, HT2A, SERT, hERG) and `BBB.csv`. Columns: `smiles, label, pchembl, year, source`, where
  `source` records ChEMBL, BindingDB, or both. These are exactly what each model is trained and
  tested on.
- **`data/endpoints_reg/antioxidant_dpph.csv`** - the antioxidant regression values
  (`smiles, y, year`).
- Provenance per endpoint (how many compounds from each source, actives/inactives):
  **`results/tables/endpoint_rebuild_provenance.csv`**.

## Ten-fold cross-validation - all data saved

- **`data/processed/cv_predictions/<TARGET>_random_oof.csv`** and **`_scaffold_oof.csv`** - for every
  compound, its out-of-fold prediction, the fold it was tested in, its true value, and its scaffold
  group. This is the complete record of the 10-fold training and testing, per endpoint, per split.
- **`results/tables/rf_cv_folds.csv`** - metrics for every one of the 10 folds x 2 splits x 13
  endpoints.
- **`results/tables/rf_cv_summary.csv`** - mean +/- standard deviation per endpoint and split.
- **`results/tables/expansion_audit.csv`** - pooled ChEMBL+BindingDB vs the ChEMBL-only baseline.

## Feature analysis

- **`results/tables/feature_block_ablation.csv`** - fingerprint-only vs descriptors-only vs combined.
- **`results/tables/feature_descriptor_importance.csv`** - permutation importance of the twelve
  descriptors per endpoint.

## External / evaluation sets (held out of training)

- **`data/external/processed/external_bbb_test.csv`** - 1,683 FDA-curated BBB compounds (306 novel).
- **`results/tables/external_bbb_validation.csv`** - BBB metrics on the 306 novel approved drugs.
- **`data/external/processed/external_drugs.csv`** - 11,723 DrugBank small molecules.
- **`results/tables/drugbank_screen_summary.csv`** - predicted-active fraction per endpoint.
- **`data/external/processed/flavonoid_panel.csv`** - 37,647 COCONUT flavonoids.
- **`data/external/processed/natural_products_coverage.csv`** - 214,740 CNS natural products.

## Models

- **`models_rf/<TARGET>.joblib`** - the deployed random forest (refit on all data).
- **`models_rf/<TARGET>_meta.json`** - n, positives, features, hyper-parameters, CV summary.
- **`models_rf/feature_names.json`** - the 1,036 feature names.

## Figures

- **`results/figures/fig_*.png`** - performance, compound counts, feature ablation, descriptor
  importance. Built by `src/brainsafe/viz/make_figures.py` from the tables above.

## Provenance and decisions

- **`data/raw/measured_endpoints_SOURCE.md`** - every source, URL, filter, licence.
- **`data/external/SOURCE.md`** - the external library provenance.
- **`docs/decisions_log.md`** - dated methodological decisions and their evidence.
- **`docs/ENDPOINT_JUSTIFICATION.md`** - why these twelve endpoints and where the values come from.
- **`docs/RF_CV_RESULTS.md`** - the results narrative.
