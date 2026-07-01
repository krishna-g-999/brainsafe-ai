# BrainSafe AI — Publication Package

Single submission folder for the manuscript and all supporting materials. Every reported value is
computed by the released scripts and stored in Supplementary/Reports (no manual transcription).

## Manuscript/
- BrainSafe_AI_Manuscript.docx  — submission-standard Word manuscript (line-numbered): Highlights,
  structured Abstract, Abbreviations, Methods (equations, hyperparameter & software tables),
  Results (AUROC with 95% bootstrap CIs), Discussion, Conclusion, Ethics, Author contributions,
  40 Harvard references. Embeds the graphical abstract and Figures 1–7; Tables 1–5.
- BS_MANUSCRIPT_FINAL.md — markdown source.

## Figures/  (300 dpi, Okabe–Ito colorblind-safe, consistent fonts, panel labels)
- graphical_abstract.png
- fig1_workflow.png  — pipeline schematic
- fig2_dataset.png   — dataset size / class balance
- fig3_validation.png — AUROC across random/scaffold/cluster/temporal (scaffold 95% CI error bars)
- fig4_roc_calibration.png — (A) ROC, (B) calibration reliability
- fig5_conformal_comparison.png — (A) conformal coverage, (B) ensemble vs baselines
- fig6_benchmark.png — vs published random-split ranges
- fig7_regression.png — predicted vs measured (antioxidant + 4 receptors, panels A–E)

## Supplementary/
- Tables/  STable0 provenance · S1 classification metrics · S2 receptor regression ·
  S3 antioxidant · S4 threshold sensitivity · S5 similarity-binned AUROC · S6 clinical composition ·
  S7 benchmark vs literature
- Reports/ raw validation JSON (endpoints, external/cluster+conformal, temporal+PR, regression,
  random-split benchmark, AUROC 95% CIs, baselines, antioxidant)
- Datasets/ measured training labels (ChEMBL_37 per-endpoint, B3DB BBB, DPPH, clinical reference)
- Documentation/ model card; benchmark/competitive analysis

## Provenance
Data: ChEMBL_37 (release 2026-05-01) + B3DB. Toolchain: Python 3.13, RDKit 2026.03.2,
scikit-learn 1.8.0. Random seed 42. Research use; pending peer review.
