# Decisions log

A dated record of methodological decisions and their rationale. Newest first. Every entry states the
decision, the reason, and where the supporting evidence lives.

## 2026-07-21, Measured inactives (PubChem) tested and reverted from the primary model
Decision: do not add bulk PubChem high-throughput inactives to the classifiers. Reason: on GSK-3-beta
(the worst-skewed endpoint, 93% active) adding 4,276 measured inactives corrected the DrugBank base
rate (71.6% -> 16.4% predicted active, genuine) but inflated scaffold AUROC 0.937 -> 0.989 because the
added negatives are chemically unlike inhibitors (median Tanimoto to actives 0.29), an easy-negative /
decoy artefact, not real discrimination. The honest discrimination estimate stays 0.937; the base-rate
skew is handled by calibration and the applicability-domain flag and disclosed as a data limitation.
Correct future approach: similarity-matched hard negatives only. Evidence:
`docs/INACTIVES_EXPERIMENT.md`, `results/tables/inactives_audit.csv`.

## 2026-07-21, Measured data expanded with BindingDB (audited, no inflation)
Decision: pool a second independent measured source (BindingDB) with ChEMBL for the eleven protein
targets, keeping labels measured-only and retraining the random forests. Reason: the review asked for
more data; ChEMBL was verified near-complete for these targets (live totals within ~3% of what we
held), so genuine growth required an independent source. Result: 57,088 to 61,317 unique compounds
(67,982 measured records); 4,246 compounds added by BindingDB alone, 18,573 measured by both sources.
Because the BindingDB affinity export returns actives, the change was audited against the ChEMBL-only
baseline: scaffold-split headline metric moved by a mean of -0.0002 (range -0.015 to +0.013), i.e. the
addition neither inflates nor degrades performance; A2A regression gained the most (+0.013 R-squared on
+1,238 compounds). hERG remained ChEMBL-only (BindingDB rate-limited). Evidence:
`results/tables/expansion_audit.csv`, `endpoint_rebuild_provenance.csv`, `docs/RF_CV_RESULTS.md`.

## 2026-07-21, External validation on approved drugs held out of training
Decision: hold the 306 FDA-curated BBB compounds that are absent from B3DB out as an external test
rather than adding them to training. Reason: an independent approved-drug test is more informative than
a small training augmentation. Result: BBB model AUROC 0.7645 (this entry read 0.774 until an audit compared it with results/tables/external_bbb_validation.csv) (accuracy 0.739, sensitivity 0.798,
specificity 0.621) on these unseen approved drugs. Evidence:
`results/tables/external_bbb_validation.csv`, `src/brainsafe/evaluation/external_validation.py`.

## 2026-07-20, Random-forest models with ten-fold cross-validation (primary model)
Decision: adopt a random forest per endpoint (300 trees, ECFP-4 plus twelve descriptors) evaluated
by ten-fold cross-validation under both a random and a scaffold-grouped split, as the primary
predictive model. Reason: this was the specific model and validation asked for in review; it
reproduces the strong, honestly-degrading behaviour seen with the earlier ensemble. Mean
classification AUROC is 0.960 random and 0.918 scaffold; receptor regression R-squared is 0.61-0.66
random and 0.39-0.56 scaffold. MAO-A and MAO-B degrade most under the scaffold split, consistent
with their being the smallest, least scaffold-diverse sets. Evidence:
`docs/RF_CV_RESULTS.md`, `results/tables/rf_cv_summary.csv`, `results/tables/rf_cv_folds.csv`,
`models_rf/*_meta.json`.

## 2026-07-20, External libraries integrated as evaluation and coverage, not training labels
Decision: fold the HPC compound libraries in as held-out evaluation (DrugBank 11,723 drugs; an
FDA-curated BBB set of 1,683, of which 306 are novel versus B3DB) and coverage (37,647 COCONUT
flavonoids; 214,740 CNS natural products; the ChEMBL 34 pool), never as training labels. Reason:
they carry structures, not measured endpoint values, so using them as labels would reintroduce the
circularity removed earlier. Evidence: `data/external/SOURCE.md`,
`data/external/processed/external_summary.csv`, `src/brainsafe/data/integrate_external.py`.

## 2026-07-20, Repository reorganised into a standard research layout
Decision: adopt the `data/{raw,interim,processed,external}`, `src/brainsafe`, `results`, `docs`,
`manuscript`, `presentation`, `archive/legacy` structure (see `PROJECT_STRUCTURE.md`). Superseded
scripts and old-version folders moved to `archive/legacy/`; the working pipeline and the application
left intact and smoke-tested after the move. Reason: the flat layout (49 root scripts, loose reports,
backups, duplicate app and figure versions) was not traceable for review. Evidence: git history of the
reorganisation commit; smoke test confirming core imports and data paths intact.

## 2026-07-20, Random forest and fold-count comparison
Decision: report scaffold-grouped 5-fold cross-validation as the primary estimate, with random 5-fold,
random 10-fold and random-forest-only reported alongside. Reason: mean classification AUROC is 0.912
(scaffold-5) versus 0.958 (random-5) and 0.964 (random-10); random-forest-only is 0.960. The fold
count barely changes the result while the split type changes it by about 0.05, showing the random-split
gain is analogue leakage between folds, not a stronger model. Evidence:
`results/tables/STable15_cv_comparison.csv`, `BS_cv_comparison.json`.

## 2026-07-20, Ensemble-versus-baseline significance
Decision: report DeLong tests and paired bootstrap on the ensemble-minus-kNN AUROC difference. Reason:
point deltas alone (for example MAO-A +0.014) do not establish significance. The paired test, which
accounts for the correlation between the two ROC curves on identical compounds, finds all eight
endpoints significant at p<0.05, decisively for seven (p<0.001) and marginally for MAO-A (p=0.033).
Evidence: `results/tables/STable14_significance.csv`, `BS_significance_report.json`.

## Earlier, Measured labels only (circularity lesson)
Decision: train every supervised endpoint on measured bioactivity (ChEMBL pChEMBL, B3DB, DPPH), never
on qualitative annotation derived from the same knowledge as the features. Reason: an earlier prototype
that used curated annotation scores was shown by feature ablation to be reading the answer back out of
disease-association features rather than learning from structure; structure-only performance collapsed
to near zero. Evidence: `docs/BS_MODEL_CARD.md`, sections 5 and 10.

## Earlier, Receptors served as regression, not classification
Decision: D2, A2A, 5-HT2A and SERT are modelled as potency regression rather than binary
classification. Reason: their ChEMBL sets are 96-98% active (only binders are reported), so a binary
task is ill-posed (Matthews correlation 0.21-0.44) and fails the deployment quality gate; regression on
pChEMBL is the appropriate task. Evidence: `docs/BS_MODEL_CARD.md` section 14;
`results/tables/STable2_receptor_regression.csv`.
