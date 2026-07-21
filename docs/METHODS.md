# Methods

A consolidated, manuscript-ready account of the data, model and validation as they stand after the
2026-07 revision. Every number traces to a file under `results/tables/` or `data/`; scripts are named
so each step is reproducible.

## Data sources (measured only)

All supervised labels are measured experimental values; no qualitative annotation is used as a label.

- **Target activity** — ChEMBL 37 (release 2026-05-01) and BindingDB, two independent public
  databases, pooled at the compound level for eleven human targets (AChE, BChE, BACE1, GSK-3-beta,
  MAO-A, MAO-B, hERG; and D2, A2A, 5-HT2A, SERT as regressors). Activities with a defined potency
  (ChEMBL pChEMBL, or a BindingDB IC50/Ki/Kd/EC50 in nM converted to -log10 molar) are retained and
  the per-compound median across both sources is taken. `src/brainsafe/data/{fetch_bindingdb,
  rebuild_endpoints}.py`.
- **Measured inactives** — for the protein-target classifiers, genuine negatives from PubChem
  BioAssay high-throughput screens are added, excluding any compound ever active for that target and
  never overriding a dose-response active. `src/brainsafe/data/{fetch_pubchem_inactives,
  rebuild_with_inactives}.py`.
- **Blood-brain barrier** — B3DB measured permeability labels (Meng et al., Sci Data 2021). A further
  306 FDA-curated approved-drug compounds absent from B3DB are held out for external validation.
- **Antioxidant** — ChEMBL DPPH radical-scavenging assays; IC50/EC50 converted to pIC50.

## Curation and standardisation

Every structure is passed through one pipeline: parse with RDKit, keep the largest organic fragment
(salt/counter-ion removal), sanitise, generate the canonical SMILES and standard InChIKey, and
deduplicate by InChIKey. The original source record and the standardised training variant are both
retained (`data/README.md`). Classification labels: active at pChEMBL >= 6, inactive at < 5; the 5-6
grey zone is dropped so the two classes are unambiguous.

## Feature representation

Each compound is encoded as 1,036 numeric features: a 1024-bit ECFP-4 Morgan fingerprint (radius 2)
and twelve interpretable physicochemical descriptors (molecular weight, cLogP, TPSA, H-bond
donors/acceptors, rotatable bonds, aromatic rings, sp3 fraction, ring count, heavy-atom count, formal
charge, QED). The fingerprint encoding is collision-free by construction (bit k always denotes the
same substructure). Categorical metadata is encoded reversibly and separately and never enters the
feature matrix (`src/brainsafe/features/{featurize,encodings}.py`).

## Models

- **Primary: random forest** per endpoint (300 trees, min_samples_leaf 2, balanced class weights for
  classification), one classifier for each of the eight classification endpoints and one regressor
  for each of the five regression endpoints (`src/brainsafe/models/train_rf.py`).
- **Estimator comparison** — random forest was compared against XGBoost and histogram gradient
  boosting under identical features and folds. Random forest is best or tied on all eight classifiers
  (mean scaffold AUROC 0.914 vs 0.905 vs 0.901); gradient boosting is marginally better on the
  receptor regressions. All differences are within about 0.02, so the estimator is not the limiting
  factor (`src/brainsafe/evaluation/model_comparison.py`).
- **Graph neural network** — a Graph Isomorphism Network (pure PyTorch) was trained on the raw
  molecular graph and compared to the random forest on an identical scaffold hold-out. At this data
  scale (2,000-8,500 compounds) the random forest wins on every endpoint tested, consistent with the
  literature that graph networks require far larger data or pretraining to overtake fingerprints
  (`src/brainsafe/gnn/`).

## Validation

- **Cross-validation** — every endpoint is evaluated by random 10-fold and scaffold-grouped 10-fold
  (Bemis-Murcko GroupKFold). The scaffold split, which holds whole scaffolds out, is the headline
  estimate of generalisation to new chemistry. Every compound's out-of-fold prediction and fold id
  are saved (`data/processed/cv_predictions/`).
- **Probability calibration** — isotonic calibration per classifier, measured honestly on out-of-fold
  predictions; mean expected calibration error 0.072 -> 0.012
  (`src/brainsafe/models/calibrate.py`).
- **Applicability domain** — nearest-neighbour ECFP-4 Tanimoto; a prediction is in domain at maximum
  Tanimoto >= 0.30. Library-wide predictions outside the domain are flagged as extrapolation
  (`src/brainsafe/evaluation/applicability_domain.py`).
- **External validation** — the BBB model applied to 306 FDA-approved drugs absent from training gives
  AUROC 0.774 (`src/brainsafe/evaluation/external_validation.py`).
- **Data-addition audits** — each data change is compared against the prior baseline. Adding BindingDB
  moved the scaffold headline by a mean of -0.0002 (no inflation); the PubChem-inactive addition is
  audited both for the metric change and for easy-negative bias (similarity of added inactives to
  actives) (`src/brainsafe/evaluation/{audit_expansion,audit_inactives}.py`).
- **Feature retention** — block ablation (fingerprint vs descriptors vs both) and descriptor
  permutation importance; both blocks are retained on measured evidence and the important descriptors
  are physically sensible (TPSA for BBB, cLogP for hERG)
  (`src/brainsafe/evaluation/feature_analysis.py`).
- **Learning curves** — scaffold-honest curves show BBB is data-saturated whereas BACE1, MAO-A and the
  receptors are still improving, quantifying where more data would help
  (`src/brainsafe/evaluation/learning_curve.py`).

## Reproducibility

Random seed 42 throughout. Environment: Python 3.13, RDKit 2026.03.2, scikit-learn 1.8, XGBoost 3.3,
PyTorch 2.12 (CPU). Trained-model binaries and large structure libraries are regenerable and are not
committed; every table and figure is produced by the scripts above.
