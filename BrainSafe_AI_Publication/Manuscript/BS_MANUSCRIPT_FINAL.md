# BrainSafe AI: an evidence-grounded, calibrated, BBB-gated multi-endpoint predictor of small-molecule effects on the human brain

*Complete methodology, results, and discussion. All figures/values are produced by the
reproducible scripts in this repository; numbers are read from the saved validation
artifacts, not estimated.*

---

## Abstract
We present BrainSafe AI, an open computational tool that, from a chemical structure alone,
estimates a compound's profile of effects relevant to the human brain. It integrates eight
machine-learning endpoints trained on **measured** public bioactivity data (ChEMBL
pChEMBL and the B3DB blood–brain-barrier database): blood–brain-barrier (BBB) penetration;
inhibition of AChE, BChE, BACE1, GSK-3β, MAO-A and MAO-B; and the hERG cardiotoxicity
liability. Four additional receptor targets (D2, A2A, 5-HT2A, SERT) are modelled as
potency regressions. A deterministic druggability/CNS-MPO layer, a measured-data
antioxidant (DPPH) model, and a clinical-precedent layer (504 nervous-system compounds with
real clinical-phase data) complete the system. Predictions are **isotonic-calibrated**, carry
**Mondrian conformal prediction sets** with empirically verified ~90% coverage, are
**grounded in the nearest real measured analogues**, and are integrated into **BBB-gated
per-disease scores**. Validation is reported across a full rigour hierarchy — random,
scaffold, leave-cluster-out, and temporal splits. On like-for-like random splits the
classifiers reach AUROC 0.94–0.98 (at/above published state of the art); under strict
scaffold and cluster splits 0.87–0.95; and under true temporal (future-compound) splits
0.61–0.92, transparently exposing where generalisation is limited. The methodological
contribution is the **integration** — calibrated, evidence-grounded, BBB-gated,
safety-aware CNS profiling from measured data — rather than any single new algorithm.

## 1. Introduction
Assessing whether a small molecule (drug or natural product/flavonoid) is likely to affect
the brain requires several questions answered together: can it cross the BBB; does it engage
disease-relevant CNS targets; is it developable; is it safe; and is there clinical precedent?
Existing web tools answer subsets — general ADMET (SwissADME, ADMETlab, admetSAR, pkCSM) or
generic target prediction (SwissTargetPrediction, PPB2) — but none unify measured-data CNS
target activity, BBB gating, calibrated uncertainty, safety, and clinical precedent in one
transparent, evidence-grounded tool. BrainSafe AI provides this integration.

## 2. Methods

### 2.1 Data sources (all measured, public)
- **CNS target bioactivity** — ChEMBL REST API. For each target we retrieved activities with
  a defined pChEMBL value (standard types IC50/Ki/Kd/EC50/Potency): AChE (CHEMBL220), BChE
  (CHEMBL1914), BACE1 (CHEMBL4822), GSK-3β (CHEMBL262), MAO-A (CHEMBL1951), MAO-B (CHEMBL2039),
  the hERG safety anti-target (CHEMBL240), and receptors D2 (CHEMBL217), A2A (CHEMBL251),
  5-HT2A (CHEMBL224), SERT (CHEMBL228).
- **Blood–brain barrier** — B3DB classification dataset (7,807 measured compounds).
- **Antioxidant** — ChEMBL DPPH radical-scavenging assays; IC50/EC50 → pIC50 (2,862 compounds).
- **Clinical precedent** — ChEMBL ATC level-1 "N" (nervous-system) molecules with max clinical
  phase ≥ 1 and a structure (504 compounds), ATC-mapped to disease class.

### 2.2 Curation and labelling
SMILES were canonicalised and deduplicated by InChIKey (RDKit). For classification, activities
were aggregated per compound by median pChEMBL and labelled **active = pChEMBL ≥ 6 (≤1 µM)**,
**inactive = pChEMBL < 5 (>10 µM)**, with the 5–6 grey zone discarded to reduce label noise.
Receptor targets (96–98 % active) were unsuitable for binary classification and were instead
modelled as **potency regression** on median pChEMBL. Antioxidant DPPH IC50/EC50 were converted
to pIC50 (−log10 M) and aggregated per compound.

### 2.3 Molecular representation
Each molecule was represented by a **1,024-bit Morgan (ECFP, radius 2) fingerprint** concatenated
with **24 interpretable RDKit physicochemical descriptors** (MW, cLogP, TPSA, HBD, HBA,
rotatable bonds, ring counts, FractionCSP3, QED, phenolic-OH and catechol flags, etc.).

### 2.4 Models
Classification endpoints use an unweighted-mean ensemble of **RandomForest (300 trees,
class-balanced), ExtraTrees (300, class-balanced), and HistGradientBoosting**. Regression
endpoints (receptors, antioxidant) use the RandomForest/ExtraTrees/HistGradientBoosting
**regressor** ensemble. A fact-based **quality gate (Matthews correlation coefficient ≥ 0.45
under scaffold CV)** governs deployment; endpoints failing the gate (D2/A2A/5-HT2A/SERT as
binary) are excluded from classification and served as regression instead.

### 2.5 Calibration and conformal prediction
Probabilities were **isotonic-calibrated** on scaffold-CV out-of-fold predictions. **Mondrian
(class-conditional) inductive conformal prediction** produces per-compound prediction sets at
the 90 % level; empirical coverage was verified on held-out calibration splits.

### 2.6 Validation hierarchy (no single split relied upon)
1. **Random** stratified split (like-for-like with most literature).
2. **Scaffold** GroupKFold(5) on Bemis–Murcko generic scaffolds (all transforms fit in-fold).
3. **Leave-cluster-out** (LeaderPicker sphere-exclusion clusters held out whole).
4. **Temporal** — train on compounds reported ≤ 75th-percentile ChEMBL document year, test on
   the most recent ~25 % (a true future-compound test).
Metrics: AUROC, PR-AUC, balanced accuracy, MCC, Brier; for regression R²/RMSE/Spearman.

### 2.7 Integration layers
- **BBB-gated disease synthesis:** effective engagement = P(target) × P(BBB); per-disease scores
  take the max over relevant targets (classifier probability, or regression pKi mapped to
  engagement), each driver tagged by provenance.
- **Evidence grounding:** every call returns the nearest real measured analogues (Tanimoto +
  measured outcome + pChEMBL).
- **Druggability:** deterministic composite of QED, Lipinski Ro5, Veber, CNS-MPO, PAINS (RDKit).
- **Clinical precedent:** nearest ATC-N clinical compounds with phase + disease class.

### 2.8 Software
Python 3, RDKit 2026.03, scikit-learn 1.8, NumPy/Pandas/SciPy; Streamlit front end. All code and
data-fetch scripts are in the repository; fixed random_state = 42 throughout.

## 3. Results

### 3.1 Classification endpoints — full validation hierarchy
| Endpoint (n) | Random | Scaffold | Cluster | **Temporal** | Brier | Conformal cov. |
|---|---|---|---|---|---|---|
| BBB (7,805) | 0.963 | 0.921 | 0.906 | – | 0.105 | 0.897 |
| AChE (4,324) | 0.975 | 0.915 | 0.912 | 0.784 | 0.099 | 0.899 |
| BChE (2,580) | 0.976 | 0.937 | 0.921 | 0.794 | 0.090 | 0.899 |
| BACE1 (8,067) | 0.956 | 0.950 | 0.940 | 0.915 | 0.039 | 0.902 |
| GSK-3β (4,044) | 0.943 | 0.920 | 0.915 | 0.658 | 0.044 | 0.885 |
| MAO-B (3,455) | 0.960 | 0.885 | 0.873 | 0.758 | 0.122 | 0.895 |
| MAO-A (2,141) | 0.950 | 0.867 | 0.890 | 0.614 | 0.136 | 0.905 |
| hERG (5,905) | 0.950 | 0.901 | 0.870 | 0.757 | 0.123 | 0.896 |

AUROC values. Conformal coverage targets 0.90; empirical 0.885–0.905 across endpoints.
*(Figure 1: validation hierarchy; Figure 2: scaffold-CV ROC curves; Figure 3: calibration
reliability; Figure 4: conformal coverage; Figure 7: dataset size/balance. Full metrics:
Supplementary Table S1; per-threshold precision/recall/F1: S4; similarity-binned generalisation: S5.)*

### 3.2 Receptor potency regression (scaffold-CV)
| Receptor (n) | R² | RMSE | Spearman | Temporal R² |
|---|---|---|---|---|
| A2A (5,547) | 0.526 | 0.74 | 0.706 | 0.326 |
| 5-HT2A (5,256) | 0.460 | 0.76 | 0.684 | 0.085 |
| D2 (7,511) | 0.425 | 0.71 | 0.652 | −0.007 |
| SERT (4,471) | 0.338 | 0.84 | 0.573 | 0.171 |

*(Predicted-vs-measured scatter: Figure 6; Supplementary Table S2.)*

### 3.3 Antioxidant (measured DPPH) and druggability
Measured DPPH regression (n=2,862): scaffold-CV **R² = 0.43, RMSE = 0.60, Spearman = 0.636**
(vs curated R²≈0.25). The prior curated score correlated only weakly with measured DPPH
(Spearman 0.39), confirming the measured model as the superior basis. Druggability/CNS-MPO is
deterministic and discriminates CNS drugs from polar non-drugs (donepezil 79, caffeine 86 vs
sucrose 46, atorvastatin 24).

### 3.4 Prospective sanity (chemistry-only inputs)
Donepezil → Alzheimer's via AChE (P=0.99) + BBB-penetrant + **hERG high** (matches its QT
liability); Selegiline → Parkinson's via MAO-B; Terfenadine → **hERG 0.99** (withdrawn for
cardiotoxicity, correctly flagged); Fluoxetine → Depression via SERT with clinical precedent
(Fluoxetine, Phase 4, Depression); Quercetin → BBB non-penetrant + high antioxidant. The system
reproduces known pharmacology and safety.

### 3.5 Benchmark vs literature
On random splits, BrainSafe AUROC (0.94–0.98) is at/above published ranges (BBB 0.88–0.96;
hERG 0.86–0.93). The same models additionally report the stricter scaffold/cluster/temporal
numbers most studies omit. *(Figure 5; Supplementary Table S7.)*

## 4. Discussion

**Contribution.** Individual components (ECFP/RF ensembles, BBB/hERG/target QSAR, conformal
prediction, QED/CNS-MPO) are standard and not novel. The contribution is their **integration**
into a single, transparent, measured-data CNS profiler that is calibrated, conformal,
evidence-grounded, BBB-gated, safety-aware, and clinically contextualised — a configuration no
existing single tool provides.

**Validation honesty.** We deliberately report a four-level split hierarchy. The collapse from
random (0.94–0.98) to temporal (0.61–0.92) quantifies real prospective difficulty: 71–91 % of
recent compounds carry scaffolds unseen in training. Where a temporal AUROC is high (BACE1
0.92) it is partly because the recent test set is 93 % active; where it is balanced (MAO-A,
45 % active) the honest number is 0.61. We surface, rather than hide, this.

**Limitations (explicit and, where inherent, unfixable computationally).**
(i) Models predict **engagement/binding, not direction** (agonist vs antagonist).
(ii) **Engagement is not efficacy**; the clinical layer provides *precedent from real trial
data*, not an efficacy prediction.
(iii) **No wet-lab prospective validation** — requires experiments.
(iv) **Temporal generalisation to novel scaffolds is bounded** by covariate shift; receptor
and pooled-assay (DPPH) endpoints generalise across time only weakly and are flagged.
(v) GSK-3β and MAO-A degrade temporally and are marked lower-confidence.

**Intended use.** Research hypothesis-generation, triage and prioritisation — not clinical or
diagnostic use.

## 5. Conclusion
BrainSafe AI is a scientifically validated, calibrated, evidence-grounded multi-endpoint CNS
profiler built entirely on measured public data, with state-of-the-art-grade per-endpoint
performance on like-for-like splits and fully transparent harder-split and prospective numbers.
Everything computationally validatable has been done; the residual gaps (efficacy prediction,
agonist/antagonist direction, wet-lab confirmation) are inherent and stated plainly. It is
suitable for an application/resource publication and as a usable research tool.

## Data and code availability
All models (`models_brain/`, `models_brain_reg/`, `models_genuine/`), datasets
(`data/endpoints/`, `data/endpoints_reg/`, `data/clinical_cns_reference.csv`), validation
reports (`*_report.json`, `BS_randomsplit_benchmark.json`), the model card (`BS_MODEL_CARD.md`),
and all fetch/train/validation scripts are in the repository. App: `app_v6_final.py`.

## Supplementary materials
**Figures** (`figures/`, 300 dpi, regenerated from out-of-fold predictions):
Fig 1 validation hierarchy · Fig 2 ROC curves · Fig 3 calibration reliability ·
Fig 4 conformal coverage · Fig 5 benchmark vs literature · Fig 6 regression scatter ·
Fig 7 dataset overview.
**Supplementary tables** (`supplementary/`, exact values from the validation artifacts):
S1 classification metrics (random/scaffold/cluster/temporal, PR-AUC, BA, MCC, Brier, conformal) ·
S2 receptor regression · S3 antioxidant (measured DPPH) · S4 threshold sensitivity ·
S5 similarity-binned generalisation · S6 clinical-reference composition · S7 benchmark vs literature.
**Model card:** `BS_MODEL_CARD.md` (full provenance, diagnosis, and limitations).

## Key references
1. Mendez D, et al. ChEMBL: towards direct deposition of bioassay data. *Nucleic Acids Res* 2019;47:D930.
2. Meng F, et al. A curated diverse molecular database of blood–brain barrier permeability (B3DB). *Sci Data* 2021;8:289.
3. Rogers D, Hahn M. Extended-connectivity fingerprints. *J Chem Inf Model* 2010;50:742.
4. Bemis GW, Murcko MA. The properties of known drugs. 1. Molecular frameworks. *J Med Chem* 1996;39:2887.
5. Pedregosa F, et al. Scikit-learn: machine learning in Python. *JMLR* 2011;12:2825.
6. Norinder U, Carlsson L, Boyer S, Eklund M. Introducing conformal prediction in predictive modeling. *J Chem Inf Model* 2014;54:1596.
7. Bickerton GR, et al. Quantifying the chemical beauty of drugs (QED). *Nat Chem* 2012;4:90.
8. Lipinski CA, et al. Rule of five. *Adv Drug Deliv Rev* 2001;46:3.
9. Veber DF, et al. Molecular properties affecting oral bioavailability. *J Med Chem* 2002;45:2615.
10. Wager TT, et al. CNS multiparameter optimization (CNS-MPO). *ACS Chem Neurosci* 2010;1:435.
11. RDKit: Open-source cheminformatics. https://www.rdkit.org
