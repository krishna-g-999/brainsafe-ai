# Estimator choice and the value of more data (evidence)

Two questions answered with experiments rather than assertion: is a different estimator better than the
random forest, and would more data help? Date: 2026-07-21.

## 1. Random forest vs XGBoost vs histogram gradient boosting

Same features (ECFP-4 + 12 descriptors), same folds, scaffold-grouped 5-fold. Full table:
`results/tables/model_comparison.csv`.

| | Random forest | XGBoost | HistGradientBoosting |
|---|---|---|---|
| Mean scaffold AUROC (8 classifiers) | **0.914** | 0.905 | 0.901 |
| Mean scaffold R-squared (5 regressors) | 0.461 | **0.481** | 0.478 |

- Random forest is best or tied on all eight **classifiers** (it wins 7 of 13 endpoints overall).
- Gradient boosting is marginally better on the **receptor regressions** (D2, 5-HT2A, SERT gain
  ~0.02-0.04 R-squared).
- Every difference is within about 0.02. **The estimator is not the bottleneck.** A defensible small
  refinement is to keep the random forest for classification and use gradient boosting for the
  receptor regressors; the gain is real but minor.

## 2. Learning curves: does more data still help?

Each endpoint has a fixed scaffold-held-out test set; the model is trained on an increasing fraction
of the remaining scaffolds and scored on that same test set. Full table and figure:
`results/tables/learning_curve.csv`, `results/figures/fig_learning_curve.png`.

| Endpoint | score at 50% -> 100% of data | verdict |
|---|---|---|
| BBB | 0.917 -> 0.921 | **plateaued** - more data will not help; the limit is representation/task |
| BACE1 | 0.925 -> 0.942 | still improving - more data would help |
| MAO-A | 0.785 -> 0.800 | still improving (small, noisy set) - more data would help |
| A2A | 0.452 -> 0.489 | still improving - more data would help |

**Conclusion.** More data is worth collecting for BACE1, MAO-A and the receptor regressors, which are
still climbing; it is not worth collecting more BBB data, which has saturated. This is endpoint-
specific and now measured rather than guessed.

## 3. What this implies for improving the model

- Swapping estimators (XGBoost) yields at most ~0.02 and is not the lever.
- The real levers are a better molecular representation (graph neural network or a pretrained
  molecular transformer), adding measured **inactives** to fix the active-skewed targets, multi-task
  learning so low-data endpoints borrow strength, and reliability tools (conformal prediction plus the
  applicability-domain gate) so the model abstains rather than guesses out of domain.
- No method removes the task ceiling: structure alone cannot predict clinical brain effect, which
  additionally needs pharmacokinetics, dose and systems-level response.
