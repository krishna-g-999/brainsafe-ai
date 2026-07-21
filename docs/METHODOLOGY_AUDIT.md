# Methodology audit

An independent audit of the BrainSafe random-forest pipeline: what is sound, which gaps were found,
how they were rectified, and which limitations remain and are disclosed rather than hidden. Every
claim here is backed by a file under `results/tables/`. Date: 2026-07-21.

## 1. What is methodologically sound (verified)

- **Measured labels only.** Every label is an experimental value (ChEMBL 37, BindingDB, B3DB, ChEMBL
  DPPH). No qualitative annotation is used as a label, so the earlier circularity is absent
  (`docs/decisions_log.md`).
- **Two independent measured sources.** Target activity is pooled from ChEMBL and BindingDB; 18,573
  compounds are measured by both (independent cross-confirmation), 4,246 come from BindingDB alone.
- **Honest validation.** Every endpoint is scored under random 10-fold and scaffold-grouped 10-fold
  (Bemis-Murcko); the scaffold number, which holds whole scaffolds out, is treated as the headline.
- **All fold data retained.** Per-compound out-of-fold predictions, fold ids and scaffold groups are
  saved (`data/processed/cv_predictions/`), so every number is reproducible and inspectable.
- **Expansion was audited, not assumed.** Adding BindingDB changed the scaffold headline by a mean of
  -0.0002 (`results/tables/expansion_audit.csv`): no inflation.
- **Feature retention is evidence-based.** Block ablation shows the combined feature set beats either
  block alone for all 13 endpoints; permutation importance shows physically sensible drivers (TPSA
  for BBB, logP for hERG) (`results/tables/feature_block_ablation.csv`,
  `feature_descriptor_importance.csv`).

## 2. Gaps found and rectified

### Gap 1 - probability calibration (fixed)
Random-forest probabilities are biased toward the centre, so a raw score is not a true probability.
The pipeline reported raw scores. **Fix:** isotonic recalibration, measured honestly on out-of-fold
predictions. Expected calibration error fell for seven of eight classifiers (mean ~0.072 to ~0.012;
hERG 0.107 to 0.006, AChE 0.092 to 0.008). BBB was already well calibrated and is left as is.
Calibrated deployment models are saved as `models_rf/<endpoint>_calibrated.joblib`
(`results/tables/calibration.csv`).

### Gap 2 - applicability domain (fixed)
The pipeline had no domain-of-reliability check, so it would return a confident-looking number for a
compound unlike anything it was trained on. **Fix:** nearest-neighbour Tanimoto similarity (ECFP-4,
2048-bit); a compound is in domain at maximum Tanimoto >= 0.30. This surfaced an important, honest
limitation: DrugBank is 72% in domain for BBB but only 20-37% in domain for the target-specific
enzyme models, so predictions on arbitrary libraries are largely extrapolation and are now flagged.
On the external BBB drugs the in- and out-of-domain performance is comparable (AUROC 0.77 vs 0.82,
the latter on only 48 compounds), consistent with BBB being governed by global physicochemical
properties rather than scaffold similarity (`results/tables/applicability_bbb_validation.csv`,
`applicability_coverage.csv`).

## 3. Limitations retained and disclosed

These are inherent to the measured data or the task and are stated plainly rather than hidden.

1. **Class imbalance and base rate.** Several targets are strongly active-skewed in the public data
   (GSK-3-beta ~93% active); the BindingDB affinity export contributes actives, which sharpens this.
   Discrimination is unaffected (scaffold AUROC 0.94), but the model's base rate of predicting
   "active" is high, so a raw predicted-active fraction on a random library over-states prevalence.
   Mitigations in place: balanced class weights, isotonic calibration, applicability-domain flagging,
   and reporting the balance-sensitive Matthews correlation and PR-AUC alongside AUROC.
2. **BindingDB contributes actives only.** Its affinity export returns high-affinity binders, so it
   adds actives, not inactives. The effect on performance was audited (no inflation). hERG could not
   be retrieved (server rate-limiting) and remains ChEMBL-only.
3. **Censored affinities.** BindingDB values reported as ">" or "<" are treated by their numeric
   bound. For labelling this is directionally correct; for regression it adds minor noise. Disclosed.
4. **Regression ceiling.** Structure-only potency prediction explains 39-58% of variance under the
   scaffold split; SERT (0.39) and the antioxidant model (0.43) are the weakest and are flagged as
   lower-confidence rather than presented as equal to the classifiers.
5. **Narrow target-model domains.** The enzyme/receptor models are trained on target-focused
   chemistry and therefore have narrow applicability domains; this is why the domain flag is
   mandatory for any library-wide prediction.
6. **Fingerprint size.** Features use a 1024-bit ECFP-4 (possible bit collisions); the applicability
   domain uses 2048-bit to reduce collision in the similarity calculation. A move to 2048-bit
   features is a low-risk future change.

## 4. Deployment readiness per endpoint

All eight classifiers exceed the deployment gate of scaffold Matthews correlation >= 0.45 (range
0.56-0.69). Among regressors, D2, A2A and 5-HT2A are usable (scaffold R-squared 0.48-0.58); SERT and
the antioxidant model are reported as lower-confidence. Every prediction should carry its calibrated
probability and its applicability-domain flag.

## 5. Evidence index

| Claim | File |
|---|---|
| 10-fold metrics, all folds | `results/tables/rf_cv_folds.csv`, `rf_cv_summary.csv` |
| Per-compound fold predictions | `data/processed/cv_predictions/*.csv` |
| Expansion did not inflate | `results/tables/expansion_audit.csv` |
| Feature retention | `results/tables/feature_block_ablation.csv`, `feature_descriptor_importance.csv` |
| Calibration | `results/tables/calibration.csv` |
| Applicability domain | `results/tables/applicability_bbb_validation.csv`, `applicability_coverage.csv` |
| External approved-drug test | `results/tables/external_bbb_validation.csv` |
| Data provenance | `results/tables/endpoint_rebuild_provenance.csv`, `data/raw/measured_endpoints_SOURCE.md` |
