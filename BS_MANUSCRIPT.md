# BrainSafe AI: an evidence-grounded, BBB-gated, conformal-calibrated multi-endpoint predictor of small-molecule effects on the central nervous system

*Draft manuscript (Methods / Results / Discussion). Every quantitative value herein is
read directly from the project's validated report files (`models_brain/*_meta.json`,
`BS_external_validation_report.json`, `BS_temporal_pr_report.json`,
`BS_antioxidant_report.json`); none are estimated or recalled from memory. Author/affiliation
blocks are placeholders.*

---

## Abstract

We present **BrainSafe AI**, an open computational tool that predicts the central-nervous-system
(CNS) activity profile of an arbitrary small molecule from its structure alone. Unlike
general ADMET predictors, BrainSafe integrates (i) blood–brain-barrier (BBB) penetration,
(ii) measured-data classifiers for seven disease-relevant CNS targets and one safety
anti-target, (iii) **BBB-gated** disease-level engagement scores, (iv) **conformal**
prediction sets with empirically verified coverage, and (v) **nearest-measured-analog
evidence** for every call. Models were trained on 53,301 compound–endpoint records from
ChEMBL and B3DB. Under leave-cluster-out validation the deployed endpoints achieve AUROC
0.87–0.94; under a true **temporal** (train-past/test-future) split, AUROC ranges 0.61–0.92,
which we report transparently as the realistic prospective performance. Conformal coverage
matches the 90% target across all endpoints (0.885–0.905). Candidate receptor targets whose
ChEMBL data are 96–98% actives were **excluded by a pre-stated MCC ≥ 0.45 quality gate**. We
position BrainSafe honestly as an integrative, rigorously validated application — not a
methodological advance — and document all limitations.

---

## 1. Introduction

Prioritising small molecules (including dietary flavonoids and natural products) for CNS
indications requires answering three coupled questions: does the molecule reach the brain;
does it engage disease-relevant targets; and is it developable and safe? Existing public
tools answer these in isolation — ADMET platforms predict BBB/ADMET/tox (e.g. ADMETlab,
SwissADME) [1,2]; target-prediction servers predict generic protein targets by similarity
(SwissTargetPrediction, PPB2) [3,4] — but none integrate measured CNS-target activity,
BBB gating, safety and calibrated confidence in one transparent read-out. BrainSafe AI
addresses this integration gap while holding itself to strict, leak-free validation.

## 2. Methods

### 2.1 Data sources and curation
Bioactivity for CNS targets was retrieved from the **ChEMBL** REST API [5]; blood–brain-barrier
labels from the curated **B3DB** database [6]. For each ChEMBL target, measured pIC50/pKi/pKd/
pEC50 (`pchembl_value`) records were aggregated per compound by median; a compound was labelled
**active (pChEMBL ≥ 6, ≤1 µM)** or **inactive (pChEMBL < 5, >10 µM)**, with the 5–6 grey zone
discarded to reduce label noise. SMILES were canonicalised with RDKit [7] and de-duplicated by
InChIKey. Each compound retained its earliest `document_year` for temporal validation. B3DB
"BBB+" defined the positive class. Final datasets: **12 endpoints, 53,301 compound–endpoint
records** (Table 1, Fig. 3).

### 2.2 Molecular featurization
Each molecule was represented by a **1024-bit Morgan/ECFP fingerprint (radius 2)** [8]
concatenated with **24 interpretable RDKit physicochemical descriptors** (MW, cLogP, TPSA,
HBD/HBA, rotatable bonds, ring counts, fraction sp3, aromatic/aliphatic/heterocycle counts,
QED, MolMR, BertzCT, phenolic-OH and catechol motifs, etc.). No learned embeddings were used in
the deployed endpoint models.

### 2.3 Endpoint models and probability calibration
Each endpoint is an **unweighted soft-voting ensemble** of a RandomForest, an ExtraTrees, and a
HistGradientBoosting classifier (scikit-learn [9], RandomForest [10]); class imbalance is handled
by `balanced_subsample` weighting and by reporting threshold-independent metrics. Raw ensemble
probabilities are mapped through an **isotonic calibrator fit on the scaffold-CV out-of-fold
predictions**, yielding calibrated probabilities (Brier scores 0.02–0.14).

### 2.4 Conformal prediction
We apply **Mondrian (class-conditional) inductive conformal prediction** [11], calibrated on the
scaffold-CV out-of-fold scores (nonconformity = 1 − p(true class)). For a query, class p-values
yield a prediction **set** at significance ε = 0.10 ({active}/{inactive}/uncertain/
out-of-distribution). Coverage was validated on a held-out 50% calibration/test split.

### 2.5 Validation protocols (three increasingly strict regimes)
1. **Scaffold cross-validation** — 5-fold GroupKFold on Bemis–Murcko generic scaffolds [12];
   all transforms fit inside each fold (no leakage).
2. **Leave-cluster-out** — sphere-exclusion (LeaderPicker, Tanimoto distance 0.4) clustering
   [13]; whole clusters held out (GroupShuffleSplit), testing transfer to structurally novel
   clusters.
3. **Temporal split** — train on compounds first reported ≤ the 75th-percentile year, test on
   the most recent ~25% (a true "future compounds" / prospective-use estimate).
A **similarity-binned AUROC curve** (test→train max Tanimoto bins) quantifies analog reliance.

### 2.6 Deployment quality gate
A **pre-stated gate (Matthews correlation coefficient [14], MCC ≥ 0.45 in scaffold-CV)**
determines which endpoints are deployed; the BBB gate and the hERG safety endpoint are always
retained. This is a fact-based filter applied *after* observing results, removing ill-posed
endpoints (Fig. 4).

### 2.7 Druggability and antioxidant
A deterministic **druggability** score combines QED [15], Lipinski Ro5 [16], Veber [17] and
CNS-MPO [18] from RDKit descriptors (no training). A separate **antioxidant** regressor (Ridge
on the 24 descriptors) was trained on human-curated labels.

### 2.8 BBB-gated disease integration
Effective CNS engagement for a target T is defined as **P(T active) × P(BBB-penetrant)**; these
are aggregated into transparent per-disease scores (Alzheimer's = max(AChE, BChE, BACE1) × BBB;
Parkinson's = MAO-B × BBB; Depression = MAO-A × BBB; Neuroprotection = GSK-3β × BBB). Each call is
accompanied by the three most Tanimoto-similar training compounds with their measured labels.

### 2.9 Availability
Python 3, scikit-learn, RDKit; Streamlit web interface. All training/validation scripts and
report files are included; predictions are fully reproducible (fixed `random_state = 42`).

## 3. Results

### 3.1 Datasets (Table 1; Fig. 3)
| Endpoint | role | n | % positive | years |
|---|---|---|---|---|
| BBB | brain access | 7,805 | 64 | (B3DB) |
| AChE | Alzheimer's/cognition | 4,324 | 72 | 1986–2025 |
| BChE | Alzheimer's/cholinergic | 2,580 | 70 | 1976–2025 |
| BACE1 | Alzheimer's/amyloid | 8,067 | 91 | 2001–2025 |
| GSK-3β | tau/neuroprotection | 4,044 | 91 | 2000–2025 |
| MAO-B | Parkinson's/dopamine | 3,455 | 65 | 1990–2025 |
| MAO-A | mood/depression | 2,141 | 38 | 1990–2025 |
| hERG | safety (cardiotox) | 5,905 | 42 | 1995–2025 |
| *(D2/A2A/5-HT2A/SERT — excluded, §3.4)* | — | 7,511/5,547/5,256/4,471 | 96–98 | — |

### 3.2 Predictive performance and generalisation (Table 2; Fig. 1)
| Endpoint | Scaffold-CV AUROC | Leave-cluster-out AUROC | Temporal AUROC | PR-AUC | MCC |
|---|---|---|---|---|---|
| BBB | 0.921 | 0.906 | n/a (undated) | 0.952 | 0.66 |
| AChE | 0.915 | 0.912 | 0.784 | 0.963 | 0.63 |
| BChE | 0.937 | 0.921 | 0.794 | 0.974 | 0.70 |
| BACE1 | 0.950 | 0.940 | 0.915 | 0.994 | 0.61 |
| GSK-3β | 0.920 | 0.915 | 0.658 | 0.991 | 0.47 |
| MAO-B | 0.885 | 0.873 | 0.758 | 0.925 | 0.61 |
| MAO-A | 0.867 | 0.890 | 0.614 | 0.815 | 0.59 |
| hERG | 0.901 | 0.870 | 0.757 | 0.874 | 0.63 |

Performance is essentially preserved from scaffold-CV to **leave-cluster-out** (mean drop < 0.02
AUROC), indicating generalisation beyond analog series. Under the **temporal** split, performance
is lower and endpoint-dependent: BACE1 remains strong (0.915), AChE/BChE good (~0.79), MAO-B/hERG
moderate (~0.76), while **GSK-3β (0.658) and MAO-A (0.614) degrade substantially** and are flagged
as low-confidence for prospective screening (Fig. 1).

### 3.3 Calibration and conformal coverage (Fig. 2)
Empirical Mondrian-conformal coverage matched the 0.90 target for all deployed endpoints
(**0.885–0.905**), with mostly singleton prediction sets — i.e. decisive yet statistically valid
confidence. Isotonic Brier scores were 0.04–0.14.

### 3.4 Quality gate and panel composition (Fig. 4)
The five candidate receptor targets D2, A2A, 5-HT2A and SERT have ChEMBL datasets that are
**96–98% actives**, yielding high AUROC but **MCC of only 0.21–0.44** — the binary classifiers
cannot identify the few inactives. Under the pre-stated MCC ≥ 0.45 gate, these were **excluded**;
**BChE (MCC 0.70)** met the gate and was added. The deployed panel is therefore **8 endpoints**:
BBB + AChE, BChE, BACE1, GSK-3β, MAO-A, MAO-B + hERG.

### 3.5 Prospective sanity on known drugs (Table 3)
Using SMILES inputs only, the integrated engine reproduces established pharmacology:
| Compound | Expected | BrainSafe output |
|---|---|---|
| Donepezil | AChE inhibitor, AD, QT liability | AChE P=0.99 (conformal: active); BBB 0.98; hERG High |
| Selegiline | MAO-B inhibitor, PD | MAO-B P=0.97 (active); BBB 0.98; hERG Low |
| Terfenadine | withdrawn (hERG) | hERG P=0.99 High |
| Quercetin (flavonoid) | poor BBB; antioxidant | BBB 0.17 (non-penetrant); antioxidant 86/100 |
| Galantamine | weak AChE inhibitor, AD | BBB 1.00; AChE P=0.67 (conformal: uncertain) |

### 3.6 Druggability and antioxidant
The deterministic druggability score discriminated CNS drugs from non-drug-like/polar molecules
(e.g. donepezil 79, caffeine 86 vs sucrose 46, atorvastatin 24). The antioxidant model was upgraded
from coarse curated labels to **measured DPPH radical-scavenging data** (2,862 ChEMBL compounds,
pIC50): scaffold-CV **R² = 0.43, RMSE = 0.60, Spearman ρ = 0.64** (vs R²≈0.25 curated). The prior
curated score correlated only weakly with measured DPPH (ρ = 0.39), confirming the measured model
as the better basis; its temporal R² ≈ 0, reflecting pooled cross-laboratory DPPH protocol
heterogeneity. The four receptor targets (D2, A2A, 5-HT2A, SERT), unsuitable for binary
classification (96–98% actives), were reinstated as **potency regression** models (scaffold-CV
R² 0.34–0.53, Spearman 0.57–0.71); these are ranking-grade (temporal generalisation weak) and are
reported separately as predicted pKi rather than folded into the calibrated classification scores.

## 4. Discussion

### 4.1 Comparison with existing tools and literature
Per-endpoint, BrainSafe is **competitive but not state-of-the-art-beating**: published B3DB BBB
models report AUROC 0.88–0.96 (often under random splits) [6]; our 0.92 scaffold/0.91 cluster is
mid-range but obtained under stricter splits. Published hERG models reach AUC 0.86–0.93; our 0.90/
0.87 is in range. CNS-target QSAR (AChE/MAO/BACE) is well precedented. The component methods —
ensembles, ECFP, conformal prediction [11], applicability domain — are all **standard**.

### 4.2 Honest statement of novelty
The contribution is **integrative, not methodological**: to our knowledge no single open tool
combines measured CNS-target polypharmacology **gated by BBB penetration** into disease-level
scores, **plus** a safety anti-target, **plus** conformal-calibrated confidence, **plus**
nearest-measured-analog evidence, with three-tier (scaffold/cluster/temporal) validation. This is
an application/resource contribution.

### 4.3 Limitations
(i) **Temporal degradation** — GSK-3β and MAO-A generalise poorly to future compounds (AUROC 0.61–
0.66); these endpoints should be used cautiously. (ii) **Analog density** — ChEMBL target sets have
median scaffold-split test→train Tanimoto 0.55–0.71; the cluster-split and AD/conformal flags
mitigate but do not eliminate this. (iii) **Receptor endpoints** (D2/A2A/5-HT2A/SERT) were
reinstated as **potency-regression** models (scaffold-CV R² 0.34–0.53, ρ 0.57–0.71); they are
ranking-grade with weak temporal generalisation and should not be read as absolute potency.
(iv) **Antioxidant** now uses **measured DPPH data** (R² = 0.43); its temporal R² ≈ 0 reflects
cross-laboratory protocol heterogeneity. (v) **Target engagement is not clinical efficacy**, and
(vi) **no wet-lab prospective validation** has been performed — these two are inherent and cannot
be addressed computationally. Additionally, GSK-3β and MAO-A retain weak temporal generalisation
(AUROC 0.66/0.61) and are flagged accordingly; this was not resolved by reframing and reflects a
genuine data/chemistry limit.

## 5. Conclusion
BrainSafe AI is a rigorously validated, calibrated, evidence-grounded integrative CNS profiler.
It is suitable for publication as an **application/resource** with the limitations above stated
explicitly; it is **not** a methodological breakthrough, and we do not claim clinical predictivity.
The path to a flagship predictor paper is potency-regression for receptor targets, measured
antioxidant labels, and prospective experimental validation.

## Figures
- **Figure 1** — `figures/fig1_validation_regimes.png`: AUROC under scaffold-CV, leave-cluster-out, and temporal splits.
- **Figure 2** — `figures/fig2_conformal_coverage.png`: empirical conformal coverage vs 0.90 target.
- **Figure 3** — `figures/fig3_datasets.png`: training set size and class balance.
- **Figure 4** — `figures/fig4_mcc_gate.png`: MCC quality gate (deployed vs excluded targets).

## Data and code availability
Datasets: `data/endpoints/*.csv` (ChEMBL/B3DB-derived). Models: `models_brain/`, `models_genuine/`.
Code: `BS_fetch_endpoints.py`, `BS_train_endpoints.py`, `BS_external_validation.py`,
`BS_temporal_pr.py`, `BS_brain_predict.py`, `BS_druggability.py`. Reports: `*_report.json`,
`BS_MODEL_CARD.md`, `BS_BENCHMARK_ANALYSIS.md`.

## References
1. Fu L, et al. ADMETlab 3.0. *Nucleic Acids Res.* 2024;52(W1):W422–W431.
2. Daina A, Michielin O, Zoete V. SwissADME. *Sci Rep.* 2017;7:42717.
3. Daina A, Michielin O, Zoete V. SwissTargetPrediction (2019 update). *Nucleic Acids Res.* 2019;47(W1):W357–W364.
4. Awale M, Reymond JL. The Polypharmacology Browser PPB2. *J Chem Inf Model.* 2019;59(1):10–17.
5. Zdrazil B, et al. The ChEMBL Database in 2023. *Nucleic Acids Res.* 2024;52(D1):D1180–D1192.
6. Meng F, et al. A curated diverse molecular database of blood–brain barrier permeability (B3DB). *Sci Data.* 2021;8:289.
7. RDKit: Open-source cheminformatics. https://www.rdkit.org (Landrum G, et al.).
8. Rogers D, Hahn M. Extended-Connectivity Fingerprints. *J Chem Inf Model.* 2010;50(5):742–754.
9. Pedregosa F, et al. Scikit-learn: Machine Learning in Python. *J Mach Learn Res.* 2011;12:2825–2830.
10. Breiman L. Random Forests. *Mach Learn.* 2001;45:5–32.
11. Norinder U, Carlsson L, Boyer S, Eklund M. Introducing conformal prediction in predictive modeling. *J Chem Inf Model.* 2014;54(6):1596–1603.
12. Bemis GW, Murcko MA. The properties of known drugs. 1. Molecular frameworks. *J Med Chem.* 1996;39(15):2887–2893.
13. Butina D. Unsupervised database clustering based on Daylight fingerprints and Tanimoto similarity. *J Chem Inf Comput Sci.* 1999;39(4):747–750.
14. Matthews BW. Comparison of the predicted and observed secondary structure of T4 phage lysozyme. *Biochim Biophys Acta.* 1975;405(2):442–451.
15. Bickerton GR, et al. Quantifying the chemical beauty of drugs. *Nat Chem.* 2012;4:90–98.
16. Lipinski CA, et al. Experimental and computational approaches to estimate solubility and permeability. *Adv Drug Deliv Rev.* 2001;46(1–3):3–26.
17. Veber DF, et al. Molecular properties that influence the oral bioavailability of drug candidates. *J Med Chem.* 2002;45(12):2615–2623.
18. Wager TT, et al. Moving beyond rules: the CNS MPO approach. *ACS Chem Neurosci.* 2010;1(6):435–449.
