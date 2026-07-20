# Project structure

This document is the map of the repository. It is kept current whenever files move.

```
BrainSafe-AI/
├── README.md                 Project overview and quick start
├── LICENSE                   MIT
├── CITATION.cff              How to cite the software and dataset
├── requirements.txt          Python dependencies (pinned)
│
├── data/
│   ├── raw/                  Data exactly as downloaded, never edited. Every source folder
│   │                         carries a SOURCE.md giving the database, version, URL and access date.
│   ├── interim/              Standardised intermediates (canonical SMILES, InChIKey dedup).
│   ├── processed/            Final modelling tables, one per endpoint (features + label), CSV.
│   ├── external/             Reference and evaluation sets (clinical compounds, flavonoid panels).
│   ├── endpoints/            [current pipeline] per-target measured training labels (SMILES, label, pChEMBL, year).
│   └── endpoints_reg/        [current pipeline] regression labels (receptors, antioxidant).
│       (endpoints/ and endpoints_reg/ migrate into processed/ during the rebuild; see decisions_log.md.)
│
├── src/brainsafe/            Go-forward code as an importable package.
│   ├── data/                 Download and curation.
│   ├── features/             Featurisation and numeric encoding of non-numeric fields.
│   ├── models/               Random forest and ensemble training with k-fold cross-validation.
│   ├── evaluation/           Cross-validation, significance testing, external comparison.
│   └── viz/                  Figures.
│
├── models/                   Trained model artifacts (binaries are git-ignored; regenerable from src).
│   models_brain/, models_brain_reg/, models_genuine/   [current pipeline] deployed model metadata.
│
├── results/
│   ├── metrics/              Validation reports (JSON).
│   ├── tables/               Supplementary tables (Sx), CSV.
│   └── figures/              Publication figures (PNG).
│
├── docs/                     Model card, data dictionary, decisions log, protocols, this file.
├── manuscript/               Manuscript source and compiled document.
├── presentation/             Slide decks.
├── BrainSafe_AI_Publication/ Assembled submission package (manuscript + figures + supplementary).
│
└── archive/legacy/           Superseded scripts and outputs, retained for provenance.
    ├── code/                 Retired scripts (earlier app and figure versions).
    └── _bulky/               Large old-version folders, backups and archives (git-ignored).
```

## Two code lineages currently coexist

1. **Measured-data endpoint pipeline (current).** `BS_predictive_model.py`, `BS_train_endpoints.py`,
   `BS_train_regression.py`, the `BS_fetch_*`, validation (`BS_external_validation.py`,
   `BS_temporal_pr.py`, `BS_significance.py`, `BS_cv_comparison.py`), figures (`BS_figures_v3.py`,
   `BS_fig_generalisation.py`) and reporting (`BS_build_docx.py`, `BS_supplementary.py`). This lineage
   produces the manuscript.

2. **Application (`app_v6_final.py`) and its engine** (`ml_v5_engine.py`, `BS_predict.py`,
   `dose_response.py`, `model_config.py`, `brain_region_mapper.py`). This lineage serves the
   interactive tool and is kept intact.

During the data expansion (see `EXPANSION_PLAN.md`) the measured-data pipeline is rewritten as clean
modules under `src/brainsafe/`, with paths and imports updated and each step tested before the old
scripts are retired to `archive/legacy/`.
```
