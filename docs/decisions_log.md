# Decisions log

A dated record of methodological decisions and their rationale. Newest first. Every entry states the
decision, the reason, and where the supporting evidence lives.

## 2026-07-20 — Repository reorganised into a standard research layout
Decision: adopt the `data/{raw,interim,processed,external}`, `src/brainsafe`, `results`, `docs`,
`manuscript`, `presentation`, `archive/legacy` structure (see `PROJECT_STRUCTURE.md`). Superseded
scripts and old-version folders moved to `archive/legacy/`; the working pipeline and the application
left intact and smoke-tested after the move. Reason: the flat layout (49 root scripts, loose reports,
backups, duplicate app and figure versions) was not traceable for review. Evidence: git history of the
reorganisation commit; smoke test confirming core imports and data paths intact.

## 2026-07-20 — Random forest and fold-count comparison
Decision: report scaffold-grouped 5-fold cross-validation as the primary estimate, with random 5-fold,
random 10-fold and random-forest-only reported alongside. Reason: mean classification AUROC is 0.912
(scaffold-5) versus 0.958 (random-5) and 0.964 (random-10); random-forest-only is 0.960. The fold
count barely changes the result while the split type changes it by about 0.05, showing the random-split
gain is analogue leakage between folds, not a stronger model. Evidence:
`results/tables/STable15_cv_comparison.csv`, `BS_cv_comparison.json`.

## 2026-07-20 — Ensemble-versus-baseline significance
Decision: report DeLong tests and paired bootstrap on the ensemble-minus-kNN AUROC difference. Reason:
point deltas alone (for example MAO-A +0.014) do not establish significance. The paired test, which
accounts for the correlation between the two ROC curves on identical compounds, finds all eight
endpoints significant at p<0.05, decisively for seven (p<0.001) and marginally for MAO-A (p=0.033).
Evidence: `results/tables/STable14_significance.csv`, `BS_significance_report.json`.

## Earlier — Measured labels only (circularity lesson)
Decision: train every supervised endpoint on measured bioactivity (ChEMBL pChEMBL, B3DB, DPPH), never
on qualitative annotation derived from the same knowledge as the features. Reason: an earlier prototype
that used curated annotation scores was shown by feature ablation to be reading the answer back out of
disease-association features rather than learning from structure; structure-only performance collapsed
to near zero. Evidence: `docs/BS_MODEL_CARD.md`, sections 5 and 10.

## Earlier — Receptors served as regression, not classification
Decision: D2, A2A, 5-HT2A and SERT are modelled as potency regression rather than binary
classification. Reason: their ChEMBL sets are 96-98% active (only binders are reported), so a binary
task is ill-posed (Matthews correlation 0.21-0.44) and fails the deployment quality gate; regression on
pChEMBL is the appropriate task. Evidence: `docs/BS_MODEL_CARD.md` section 14;
`results/tables/STable2_receptor_regression.csv`.
