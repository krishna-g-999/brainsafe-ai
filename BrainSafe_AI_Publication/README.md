# BrainSafe AI — Submission Package

Single folder containing the manuscript and **all** supporting materials for submission.
Every reported value is computed by the released scripts and stored under `Supplementary/`
(no manual transcription); figures are regenerated from out-of-fold predictions.

**Integrity statement.** All numbers in the manuscript, figures, tables and presentations were
cross-checked against the deployed model artifacts (`models_brain/*_meta.json`, validation JSON
reports) on the date of this package. No value is fabricated, estimated, or assumed.

---

## Manuscript/
- **BrainSafe_AI_Manuscript.docx** — submission-standard Word manuscript (continuous line numbering):
  Highlights, structured Abstract, Abbreviations, Methods (equations, hyperparameter & software
  tables), Results (non-comparative + comparative, with 95% bootstrap CIs), the executed LLM
  head-to-head (Table 8), Discussion (including a threats-to-validity self-audit), Conclusion, Ethics,
  Author contributions, and **44 Harvard references**. Embeds the graphical abstract, Figures 1–7, and
  Tables 1–8.
- **BS_MANUSCRIPT_FINAL.md** — Markdown source of the manuscript.

## Figures/  (300 dpi · Okabe–Ito colorblind-safe · consistent fonts · panel labels)
graphical_abstract · fig1 workflow · fig2 dataset size/balance · fig3 AUROC across the four
validation regimes (scaffold 95% CI error bars) · fig4 ROC + calibration · fig5 conformal coverage +
ensemble-vs-baselines · fig6 vs published ranges · fig7 predicted-vs-measured (antioxidant + 4 receptors).

## Presentation/
- **BrainSafe_AI_Presentation.pptx** — 16-slide dual-audience talk (expert + commercial). Every
  technical term is glossed in plain language on-slide; **full speaker notes** on every slide with
  narration and anticipated expert and general/commercial Q&A.
- **BrainSafe_AI_Supplementary.pptx** — supplementary deck: Part A all 8 figures; Part B every
  supplementary table (S0–S13) rendered from the CSVs.

## Supplementary/
- **Tables/** — S0 provenance · S1 classification metrics · S2 receptor regression · S3 antioxidant
  (DPPH) · S4 threshold sensitivity · S5 similarity-binned AUROC · S6 clinical composition ·
  S7 benchmark vs literature · S8 BrainSafe-vs-LLM capability matrix · S9 ablation vs baselines ·
  S10 label-definition robustness · S11 assay-type composition · S12 IC50-only-vs-pooled sensitivity ·
  S13 LLM head-to-head scoreboard.
- **Reports/** — raw validation JSON (endpoints; external/cluster + conformal; temporal + PR;
  regression; random-split benchmark; AUROC 95% CIs; baselines; antioxidant; flaw-fix analyses;
  LLM benchmark ground truth, responses, scoreboard, capability bundle).
- **Datasets/** — measured training labels (ChEMBL_37 per-endpoint, B3DB BBB, DPPH, clinical reference).
- **Documentation/** — model card; benchmark/competitive analysis; pre-registered LLM benchmark protocol.

---

## Headline verified results
- **Classification AUROC** — random 0.94–0.98 · scaffold/cluster 0.87–0.95 · temporal 0.61–0.92;
  Brier 0.04–0.14; conformal coverage 0.885–0.905 (target 0.90).
- **Comparative** — deployed ensemble mean scaffold AUROC 0.912 vs kNN-Tanimoto 0.867 vs logistic
  0.808; best on all 8 endpoints.
- **Robustness** — label-cut spread ≤ 0.109; IC50-only-vs-pooled ΔAUROC ≤ 0.006; AUROC 0.958→0.770
  as nearest-neighbour similarity falls (justifies the applicability-domain flag).
- **LLM head-to-head (executed, 4 models)** — LLMs match/beat on famous-drug classification, but 45%
  (14/31) of the ChEMBL identifiers they cited were fabricated or resolved to the wrong molecule, and
  all four confabulated on an unpublished compound; BrainSafe: 0 fabricated, grounded, honestly uncertain.

## Provenance
Data: ChEMBL_37 (release 2026-05-01) + B3DB + ChEMBL DPPH + ChEMBL ATC-N; 64,474 measured records
(53,301 ChEMBL targets + 7,807 B3DB + 2,862 DPPH + 504 clinical). Toolchain: Python 3.13,
RDKit 2026.03.2, scikit-learn 1.8.0; fixed random seed 42. Research use; pending peer review.
Full source and regeneration scripts: https://github.com/krishna-g-999/brainsafe-ai
