# BrainSafe AI

**Repository:** https://github.com/krishna-g-999/brainsafe-ai

**An evidence-grounded, calibrated, blood–brain-barrier-gated multi-endpoint predictor of
small-molecule effects on the human brain.**

BrainSafe AI predicts, from chemical structure alone, a compound's profile of brain-relevant
effects: blood–brain-barrier (BBB) penetration; engagement of disease-relevant CNS targets
(AChE, BChE, BACE1, GSK-3β, MAO-A, MAO-B; receptor potencies for D2, A2A, 5-HT2A, SERT);
an hERG cardiotoxicity safety flag; a measured-data antioxidant model; and a deterministic
druggability/CNS-MPO layer. Every machine-learning endpoint is trained on **measured public
bioactivity data** (ChEMBL_37 and B3DB, 64,474 measured records). Predictions are
probability-calibrated, carry conformal prediction sets, are grounded in nearest measured
analogues, and are integrated into BBB-gated per-disease scores.

> **Research use, pending peer review.** This tool predicts molecular target engagement and
> physicochemical properties, not clinical efficacy, and has not undergone wet-lab or clinical
> validation. It is not for medical, diagnostic, or treatment decisions.

## Model and validation (honest, multi-regime)
Deployed model: a **random forest** per endpoint (chosen after comparison with XGBoost, gradient
boosting and a graph neural network), features are a 1024-bit ECFP-4 fingerprint plus 12 descriptors.
Target panel (13 endpoints) plus a 9-endpoint ADME/exposure layer; 61,317 unique measured compounds.

| Regime | Classifier AUROC |
|---|---|
| Random 10-fold | 0.95–0.97 (mean 0.96) |
| Scaffold 10-fold (GroupKFold) | 0.87–0.96 (mean 0.92) |
| Temporal (future compounds) | 0.61–0.91 |
| External (306 FDA drugs, BBB) | 0.774 |

Isotonic-calibrated probabilities (mean ECE 0.072→0.012); conformal coverage 0.89–0.92 (target 0.90);
applicability-domain flag on every prediction. **Validated by inversion (six adversarial checks, all
pass): [`docs/VALIDATION.md`](docs/VALIDATION.md).** Methods and results:
[`docs/METHODS.md`](docs/METHODS.md), [`docs/RF_CV_RESULTS.md`](docs/RF_CV_RESULTS.md),
[`docs/ADME_RESULTS.md`](docs/ADME_RESULTS.md).

## Repository structure
Standard research layout (`data/`, `src/brainsafe/{data,features,models,evaluation,viz,gnn,adme}`,
`results/`, `docs/`, `manuscript/`, `archive/legacy/`). Superseded code is in `archive/legacy/`.

Key entry points: **`app.py`** (interactive application over the current models),
`src/brainsafe/models/train_rf.py` (training), `src/brainsafe/evaluation/` (validation, comparison,
conformal/temporal, inversion checks), `src/brainsafe/adme/` (ADME/exposure + K_p,uu),
`src/brainsafe/viz/make_figures.py` (figures), `src/brainsafe/build_manuscript_docx.py` (manuscript).

## Reproducibility
Trained-model binaries and large structure libraries are **not** stored in git; they are regenerated
by the released scripts. Full recipe in [`docs/VALIDATION.md`](docs/VALIDATION.md) §5; in brief:
```bash
python src/brainsafe/data/rebuild_endpoints.py          # pool ChEMBL + BindingDB
python src/brainsafe/models/train_rf.py                 # RF + 10-fold (all endpoints)
python src/brainsafe/models/calibrate.py                # isotonic calibration
python src/brainsafe/evaluation/rf_conformal_temporal.py
python src/brainsafe/adme/fetch_adme.py && python src/brainsafe/adme/train_adme.py
python src/brainsafe/evaluation/validate_inversion.py   # adversarial checks
streamlit run app.py                                    # interactive tool
```

## Using the web server

Three ways in, all from the same models:

| Mode | What it is for |
|---|---|
| **Compound Search** | one compound, full report: exposure, mechanism network, target engagement profile, disease relevance, ADME, read-across, CNS MPO and applicability domain |
| **Batch Screening** | up to 300 compounds pasted or uploaded as CSV, TSV or plain text, returned as one ranked row each with a CSV download |
| **Export** | every result downloadable as a tidy CSV, a self-contained HTML report, a structured JSON object, or the network as vector SVG |

The HTML report inlines its own figures, structure image and styling, so it opens offline in any
browser and prints to PDF without a network connection. The JSON carries thresholds, the screening
mode and the caveats, so a stored answer stays interpretable away from the interface.

## Running it

Locally:

```bash
pip install -r requirements.txt
streamlit run app.py
```

As a container, which is how it should be deployed:

```bash
docker build -t brainsafe-ai . && docker run -p 8501:8501 brainsafe-ai
```

The image runs as an unprivileged user, reads only from disk, and exposes Streamlit's health
endpoint at `/_stcore/health` so an orchestrator can restart a wedged container. Server settings and
the theme live in `.streamlit/config.toml`.

Before deploying, verify the environment and the application in one step:

```bash
python src/brainsafe/evaluation/app_health.py
```

It checks that every declared dependency resolves and matches its pin, that all 63 model artefacts
load, that the knowledge graph is internally consistent, that chemically unrelated compounds produce
distinct and directionally correct profiles, that every export format is well formed and
self-contained, and that no red hue has entered the palette. It exits non-zero on any failure, so it
can gate a release.

## Environment
Python 3.13, RDKit 2026.03.2, scikit-learn 1.8.0, NumPy 2.4.6, pandas 3.0.3, SciPy 1.17.1,
matplotlib 3.10.9, Plotly 6.7.0. See `requirements.txt`. Random seed 42 throughout.

## Data sources
ChEMBL 37 and BindingDB (target activity), B3DB (BBB; Meng et al.,
2021), and ChEMBL DPPH radical-scavenging assays. Structures for user-entered compounds are
resolved via PubChem.

## Citation
Manuscript in preparation (Sri Sathya Sai Institute of Higher Learning). See
`BrainSafe_AI_Publication/Manuscript/`.

## License
See `LICENSE`.
