# BrainSafe AI v5

**Multi-Dimensional Neuroprotective Scoring Platform for
Neurodegenerative Disease Drug Discovery**

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://YOUR-URL.streamlit.app)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## About
BrainSafe AI v5 computes a Neuroprotective Score (NPS) and 7 mechanistic
dimension scores for any compound against all four major NDDs
(Alzheimer's, Parkinson's, ALS, Huntington's).

**Hold-out R² = 0.782 | Spearman ρ = 0.880 | n=325 gold compounds**

## Web Interface
→ https://YOUR-URL.streamlit.app
No login required. Enter compound name, SMILES, or InChI.

## Repository Structure
brainsafe_ai/
├── app.py # Streamlit web application
├── brainsafe_v5_training_set.csv # Gold-standard database (via Zenodo)
├── manuscript_final/ # Manuscript figures and tables
│ ├── figures/ # Main + supplementary figures
│ ├── validation_figures/ # Validation figures
│ ├── data/ # Tables S1-S2 + validation reports
│ └── manuscript/ # Manuscript draft (NAR Web Server)
├── models_v5/ # Trained models (via Zenodo)
└── requirements.txt

## Citation
> Gunanathan K et al. (2026) BrainSafe AI: A Multi-Dimensional
> Neuroprotective Scoring Platform. *Nucleic Acids Research* (submitted).

## Zenodo DOI
[![DOI](https://zenodo.org/badge/DOI/PLACEHOLDER.svg)](https://doi.org/PLACEHOLDER)
