# Provenance and limitations

## How this package was produced

Every value was read from the file that produced it by
`src/brainsafe/analysis/build_reviewer_package.py`. Nothing was typed by hand and nothing was
recalled. Re-running that script regenerates the package; if it ever disagreed with the models, the
package would change rather than the models.

Generated 2026-08-09 from the deployed model set archived at
**doi:10.5281/zenodo.21858576**.

## What is in the panel

- 72 endpoints in total, of which 70 are deployed.
- 8 core classifiers and
  4 receptor regressors, cross-validated
  10-fold under both a random and a scaffold-grouped split.
- 49 binder classifiers, validated against compounds
  experimentally tested on the same target and found inactive, held out from training.
- 9 ADME and exposure endpoints, 10-fold under both
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

- **Two endpoints were trained and then withdrawn** (2 in
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
