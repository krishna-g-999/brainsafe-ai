# Random-forest models with ten-fold cross-validation: results

Primary model requested in review: a random forest per endpoint, evaluated by ten-fold
cross-validation under a random and a scaffold-grouped split. All numbers come from
`results/tables/rf_cv_summary.csv` (per-fold values in `rf_cv_folds.csv`; per-compound out-of-fold
predictions in `data/processed/cv_predictions/`), produced by `src/brainsafe/models/train_rf.py`.
Features: 1024-bit ECFP-4 fingerprint plus twelve physicochemical descriptors (1,036 numeric
features). Hyper-parameters: 300 trees, min_samples_leaf 2, balanced class weights for
classification, random seed 42.

## Data (measured, published)

Every label is a measured experimental value. Target activity is pooled from two independent public
sources - ChEMBL 37 and BindingDB - standardised to one structure per InChIKey and to the shared
-log10(molar) potency scale; blood-brain-barrier permeability is from B3DB; the antioxidant assay is
ChEMBL DPPH. After pooling: **61,317 unique compounds** across **67,982 measured compound-endpoint
records**. Of the target compounds, **18,573 are measured by both ChEMBL and BindingDB** (independent
cross-confirmation) and **4,246 are contributed by BindingDB alone**. Provenance per endpoint is in
`results/tables/endpoint_rebuild_provenance.csv`.

## Why two splits

- **Random 10-fold** - the conventional estimate; close analogues of a test compound can sit in the
  training folds, so it flatters the model.
- **Scaffold 10-fold** - GroupKFold on the Bemis-Murcko scaffold; whole scaffolds are held out, so it
  measures generalisation to genuinely new chemistry. This is the headline number.

## Classification endpoints (area under ROC curve)

| Endpoint | n | Random 10-fold | Scaffold 10-fold | Scaffold MCC |
|---|---|---|---|---|
| BBB | 7,805 | 0.960 +/- 0.007 | 0.920 +/- 0.037 | 0.658 |
| AChE | 4,387 | 0.963 +/- 0.008 | 0.921 +/- 0.021 | 0.656 |
| BChE | 2,621 | 0.968 +/- 0.012 | 0.937 +/- 0.021 | 0.691 |
| BACE1 | 8,501 | 0.967 +/- 0.010 | 0.956 +/- 0.021 | 0.672 |
| GSK-3-beta | 4,958 | 0.969 +/- 0.013 | 0.937 +/- 0.030 | 0.559 |
| MAO-A | 2,228 | 0.947 +/- 0.017 | 0.868 +/- 0.046 | 0.564 |
| MAO-B | 3,665 | 0.954 +/- 0.008 | 0.890 +/- 0.033 | 0.626 |
| hERG | 5,875 | 0.954 +/- 0.007 | 0.921 +/- 0.035 | 0.676 |
| **Mean** | | **0.960** | **0.919** | |

## Regression endpoints (coefficient of determination, R-squared)

| Endpoint | n | Random 10-fold | Scaffold 10-fold |
|---|---|---|---|
| D2 | 7,734 | 0.601 +/- 0.018 | 0.483 +/- 0.052 |
| A2A | 6,785 | 0.682 +/- 0.023 | 0.576 +/- 0.066 |
| 5-HT2A | 5,989 | 0.636 +/- 0.028 | 0.490 +/- 0.054 |
| SERT | 4,572 | 0.602 +/- 0.038 | 0.388 +/- 0.124 |
| Antioxidant (DPPH) | 2,862 | 0.669 +/- 0.065 | 0.434 +/- 0.100 |

## Effect of the BindingDB expansion (audit)

The BindingDB affinity export contributes high-affinity binders (actives), so the data addition was
audited against the earlier ChEMBL-only models (`results/tables/expansion_audit.csv`). The
scaffold-split headline metric changed by a mean of **-0.0002** across all endpoints (range -0.015 to
+0.013): the added measured data neither inflates nor degrades cross-validated performance. The
largest single gain is A2A regression (+0.013 scaffold R-squared on +1,238 compounds). This confirms
the expansion is legitimate additional signal, not analogue inflation.

One consequence to interpret carefully: for targets whose measured data are strongly active-skewed
(for example GSK-3-beta, now ~93% active), the model's base rate of predicting "active" is
correspondingly high. Discrimination (AUROC 0.937 under scaffold split) is unaffected, but absolute
probabilities should be read with the calibration layer and this base rate in mind.

## External validation (approved drugs, never in training)

The deployed BBB model was applied to 306 FDA-curated compounds not present in B3DB
(`results/tables/external_bbb_validation.csv`): **AUROC 0.774, accuracy 0.739, sensitivity 0.798,
specificity 0.621**. The drop from the 0.92 cross-validated figure is expected and honest - these are
genuinely external approved drugs. A descriptive screen of the eight classifiers over the 11,723
DrugBank small molecules is in `results/tables/drugbank_screen_summary.csv`.

## Reading the numbers

- The eight classifiers hold up under the scaffold split (mean AUROC 0.919), so performance is real
  generalisation, not memorised analogues. MAO-A and MAO-B fall most from random to scaffold, as
  expected for the smallest, least scaffold-diverse sets.
- Receptor regressors explain 60-68% of potency variance on a random split and 39-58% on
  scaffold-new chemistry - the realistic ceiling for structure-only potency prediction, reported as
  such.

## Files

- `results/tables/rf_cv_summary.csv` - mean and standard deviation per endpoint and split.
- `results/tables/rf_cv_folds.csv` - every fold's metrics (10 folds x 2 splits x 13 endpoints).
- `data/processed/cv_predictions/<endpoint>_<split>_oof.csv` - every compound's out-of-fold
  prediction, fold id and scaffold group.
- `results/tables/expansion_audit.csv` - pooled vs ChEMBL-only comparison.
- `models_rf/<endpoint>.joblib` and `_meta.json` - deployed model and metadata.
