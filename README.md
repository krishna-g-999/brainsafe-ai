# BrainSafe AI

A computational platform for neuroprotective compound scoring, BBB permeability prediction, and ML-based compound expansion for neurodegenerative disease research.

## Live App
[**https://brainsafe-ai.streamlit.app**](https://brainsafe-ai.streamlit.app)

## Features
- Neuroprotection Score (NPS) for 134 curated + 191 ML-predicted compounds
- Blood-brain barrier permeability prediction
- Morgan fingerprint-based ML feature vectors (98.4% score uniqueness)
- Supports 325 compounds with live PubChem + ChEMBL lookup
- Interactive Streamlit web application

## Run Locally
Install dependencies and launch:

    pip install -r requirements.txt
    streamlit run app.py

## Tech Stack
- Python 3.11, RDKit 2025.09.5, scikit-learn, Streamlit
- Random Forest classifier with Morgan circular fingerprints (radius=2, 64-bit)

## Citation
Gunanathan K. et al. BrainSafe AI: Brain Health Compound Explorer — a SAI-Net Translational Module (2026)

## Source Code
https://github.com/krishna-g-999/brainsafe-ai
