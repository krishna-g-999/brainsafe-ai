# Data manifest

Every data file behind the deployed panel, what it contains, and how it was made. All values are
measured (ChEMBL 37, BindingDB, B3DB, Therapeutics Data Commons, MoleculeNet, NPASS 3.0); nothing is
imputed. Regenerate any file by running the script named beside it. Counts below are read live by
`tools/build_data_manifest.py` and reflect the repository as of 2026-09-17.

This replaces an earlier version of this document built around `data/processed/compound_library.csv`,
a 61,317-compound, thirteen-endpoint snapshot from before the panel grew to its current scope. That
file and the scripts that built it (`src/brainsafe/data/build_compound_library.py`,
`src/brainsafe/viz/make_figures.py`) still exist on disk but are no longer part of the deployed
pipeline; nothing below points to them.

## The training data

- **`data/endpoints/<TARGET>.csv`** (63 files) - one file per target: the core
  target-potency and activity endpoints, the binder panel (including withdrawn endpoints, kept
  rather than deleted), and blood-brain-barrier labels. Columns: `smiles, label, pchembl, year,
  source`, where `source` records ChEMBL, BindingDB, or `ChEMBL_inactive` for a censored bound
  recovered as a measured non-binder. These are exactly what each model is trained and tested on.
- **`data/adme/<ENDPOINT>.csv`** (9 files) - the nine ADME and exposure endpoints
  beyond BBB (Caco-2 permeability, hepatocyte clearance, unbound brain-to-plasma ratio,
  lipophilicity, logBB, P-glycoprotein inhibition and substrate status, plasma-protein binding,
  aqueous solubility).
- Together, the **55 tables behind currently deployed models** hold
  **228,200 measured compound-endpoint records**, a median of
  3,789 rows per table (range 387 to 10,276).
- **Unique compounds.** 170,619
  distinct SMILES strings appear across the endpoint tables, of which
  2 could not be parsed. The remaining
  strings collapse to **169,341
  distinct InChIKeys** of the desalted, neutralised parent
  (1,276 SMILES strings were salts,
  tautomers or charge states of a compound already counted) — this is the compound count the
  manuscript's training-data paragraph states.
  Computed by `src/brainsafe/evaluation/unique_compound_count.py`;
  full breakdown in `results/tables/unique_compound_count.csv`.
- **Provenance per endpoint** (how many compounds from each source, actives/inactives):
  `results/tables/endpoint_rebuild_provenance.csv`.
- **Recovery of the negative class from censored bounds**, and its effect on class balance per
  endpoint: `results/tables/expansion_inactives.csv`, built by
  `src/brainsafe/data/recover_expansion_inactives.py`.
- **The actual numeric input every model trains on.** A SMILES string and a label are what the
  endpoint tables above store; the 1,036-column vector a model actually receives
  (`src/brainsafe/features/featurize.py`) is computed on the fly and was not previously written to
  any file. `results/tables/master_feature_vectors.csv` (one row per distinct compound, keyed by
  InChIKey) and `results/tables/master_training_usage.csv` (one row per endpoint-compound pair
  actually used, joinable to the vectors on `inchikey`) make it inspectable directly, built by
  `tools/build_master_training_vectors.py`. Every formula that turns a vector into a reported score
  is walked through with a live worked example in `docs/ML_METHODS_AND_FORMULAS.md`.

## Cross-validation - all data saved

- **`results/tables/rf_cv_folds.csv`** / **`rf_cv_summary.csv`** - per-fold and summary metrics for
  the core target-potency, activity and exposure endpoints, ten folds under both a random split and
  a scaffold-grouped split (Bemis-Murcko).
- **`results/tables/binder_cv_folds.csv`** / **`binder_cv_summary.csv`** - the same for the
  binder panel, scaffold-grouped only, against measured inactives.
- **`results/tables/model_comparison.csv`** - the like-for-like comparison against XGBoost,
  histogram gradient boosting, logistic regression and nearest-neighbour read-across.

## Background specificity and applicability domain

- **`models_rf/ad_reference.pkl`** and **`src/brainsafe/models/pools.py`** - the background library,
  partitioned into three disjoint pools by a stable hash of the canonical structure (one supplies
  decoys, one sets thresholds, one measures the false-positive rate), so a compound's pool is a
  property of the molecule rather than of run order.
- **`results/tables/background_specificity_disjoint.csv`** - the background false-positive rate at
  the currently deployed threshold, measured on the evaluation pool, for every deployed binder
  endpoint. Built by `src/brainsafe/evaluation/background_specificity_disjoint.py`; this is the
  figure the manuscript's Methods section quotes.
- **`results/tables/inversion_validation.csv`** - the applicability-domain adversarial check.

## External and prospective validation (held out of training)

- **`results/tables/external_bbb_validation.csv`** - the barrier model against FDA-curated approved
  drugs absent from B3DB by InChIKey.
- **`results/tables/external_prospective.csv`** - every qualifying endpoint refitted on its
  pre-cutoff rows, decision threshold frozen before the cutoff, tested on compounds first published
  afterwards.
- **`results/tables/external_novelty_strata.csv`** - recall against maximum Tanimoto similarity to
  the training actives, across three independently built test sets (by date, at random, by curator).
- **`results/tables/external_natural_products_summary.csv`** - the three targets tested for
  natural-product coverage, and the endpoints where genuinely external NPASS data existed to score.

## Models

- **`models_rf/<TARGET>.joblib` / `<TARGET>_binder.joblib`** (74 files) - the deployed
  estimators, refit on all available data after cross-validation.
- **`models_rf/<TARGET>_meta.json`** (51 files) - training size, positives, features,
  hyperparameters, calibration method.
- **`models_rf/binder_modes.json`** - the binder panel registry: threshold, AUROC and sensitivity
  against measured inactives, deployment status, and the recorded reason for every withdrawal.
- **`results/tables/MODEL_INVENTORY.csv`** - one row per trained estimator, family, task, training
  size, validation scheme, calibration, deployment status.
- **`models_manifest.json`** - SHA-256 checksum of every shipped model file, so a download is
  verified rather than trusted.

## Figures

- **`manuscript/figures/*.png`** (and `.pdf`) - built by the scripts in
  `src/brainsafe/figures/` (13 scripts:
  fig01_architecture, fig02_feature_vector, fig03_cv_design, fig04_pools_and_thresholds, fig05_negative_class, fig06_validation, fig07_binder_panel, fig08_use_case, fig09_model_atlas, fig10_endpoint_selection, fig11_external_validation, fig_decision_tree_example, fig_graphical_abstract), each reading the results tables named in
  its own docstring rather than a value typed into the script. `style.py` in the same directory is
  the single shared colour palette, font and layout convention every figure uses.

## Falsification suite

- **`inversion/results/*.csv`** - the per-hypothesis evidence for the ten-hypothesis falsification
  suite (H1-H10).
- **`inversion/results/VERDICTS.csv`** - one verdict per hypothesis, computed from the CSVs above by
  `inversion/summarise.py` rather than asserted.
- **`inversion/REPORT.md`** - the falsification narrative.

## Provenance and decisions

- **`data/raw/measured_endpoints_SOURCE.md`** - every source, URL, filter, licence.
- **`docs/decisions_log.md`** - dated methodological decisions and the evidence behind each.
- **`docs/ENDPOINT_JUSTIFICATION.md`** - why these targets, and exactly how each endpoint's training
  and test values are obtained.
- **`docs/TECHNICAL_REPORT.md`** - the full methods and results narrative, of which this manifest
  and the endpoint justification are companion reference documents.

---
*Generated by `tools/build_data_manifest.py`. 90 tables currently sit in
`results/tables/`; not all are named individually above, and the full list is browsable directly in
that directory or in `submission_package/08_VALIDATION_RESULTS/` for a reviewer working from the
assembled package.*
