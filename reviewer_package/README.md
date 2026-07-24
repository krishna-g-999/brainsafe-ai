# BrainSafe AI, reviewer reproducibility package

A self-contained bundle for independent verification: the master data, the per-endpoint training/test
sets, the complete 10-fold dataset, the model scripts in both `.py` and `.ipynb`, and every result
table. Everything here is measured public data; no value is imputed or hand-annotated. Fixed
`random_state = 42` throughout.

Repository: https://github.com/krishna-g-999/brainsafe-ai

## What a reviewer most likely wants, and where it is

| Requested item | Location |
|---|---|
| Complete master table (compounds + descriptors + labels + sources) | `data/master_compound_library.csv` (61,317 compounds) |
| Every individual measurement (long form) | `data/endpoint_measurements_long.csv` |
| Per-endpoint training/testing sets | `data/endpoints/*.csv` (13 files: `smiles,label,pchembl,year,source`) |
| The 10-fold dataset (per-compound out-of-fold prediction + fold id + scaffold group) | `data/cv_predictions/<ENDPOINT>_<random\|scaffold>_oof.csv` (26 files) |
| Exact model features (1024 ECFP-4 bits + 12 descriptors = 1036) | `data/feature_names_1036.json` |
| Model training/testing scripts (`.py`) | `scripts/` |
| Same, as runnable notebooks (`.ipynb`) | `notebooks/` |
| All result tables (CV, calibration, conformal, temporal, model comparison, validation) | `results/*.csv` |
| Figures | `figures/*.png` |

## Notebooks (both `.py` and `.ipynb` provided)

1. `01_master_data_and_features` — the master table, descriptors, per-endpoint counts, the exact
   1036-feature definition, and the feature-retention evidence.
2. `02_reproduce_10fold_training_testing` — **the core reproduction.** All logic is inline (no project
   imports): raw endpoint CSV → 1036 features → random forest → random and scaffold 10-fold →
   AUROC / R². It retrains MAO-A and asserts the scaffold AUROC equals the published **0.868**, then
   prints the full published summary for all thirteen endpoints.
3. `03_model_comparison` — random forest vs XGBoost vs gradient boosting vs graph neural network, with
   an independent live RF-vs-XGBoost re-run on one endpoint.
4. `04_validation_calibration_conformal_temporal` — the calibration, conformal-coverage, temporal and
   adversarial (inversion) tables, plus two live checks (no-leakage, no-duplication).

Each `.py` is a jupytext "percent" file: it is valid Python (run it directly) and is the source of the
matching `.ipynb`.

## Environment and how to run

```
# Python 3.13; install the pinned dependencies
pip install -r ../requirements.txt      # RDKit 2026.03, scikit-learn 1.8, XGBoost 3.3, numpy, pandas, scipy

# reproduce the headline numbers from raw data (a couple of minutes on CPU):
python notebooks/02_reproduce_10fold_training_testing.py
# or open the .ipynb in Jupyter and Run All.
```

## Verified

`02_reproduce_10fold_training_testing` was executed on assembly and reproduced **MAO-A scaffold-split
AUROC = 0.868**, matching the published value exactly. The scaffold split is verified leak-free (no
scaffold shared between train and test) and the master table verified duplicate-free (61,317 rows =
61,317 unique InChIKeys) in notebook 4.

## Method in one paragraph

Every SMILES is reduced to its largest organic fragment and encoded as a 1024-bit ECFP-4 (Morgan,
radius 2) fingerprint plus twelve interpretable descriptors. Each endpoint is a random forest (300
CART trees, minimum leaf size two, balanced class weights for classification), evaluated by random
10-fold and scaffold-grouped 10-fold cross-validation; the scaffold split is the headline estimate of
generalisation to new chemistry. Classification labels are measured potency thresholded at
pChEMBL ≥ 6 (active) / < 5 (inactive), grey zone dropped; receptors and the antioxidant assay are
regression on the measured value. Full detail: `../docs/METHODS.md`, `../docs/VALIDATION.md`.
