> **SUPERSEDED, 2026-08-17.** The manuscript under submission is
> `NAR_WebServer_BrainSafe_draft.md`, written to the NAR Web Server format and carrying the current
> numbers. This file is the longer July draft, kept for its methodological detail. Its figures
> predate the audit retrain and should not be quoted. The phrase "or a natural product such as a
> flavonoid" was removed from its introduction on 2026-08-17: no validation supports applying the
> panel to natural products, and the evidence available since points the other way
> (`results/tables/external_natural_products_summary.csv`).

# BrainSafe AI: an evidence-grounded, calibrated, BBB-gated multi-endpoint predictor of small-molecule effects on the human brain

*Complete methodology, results, and discussion. All figures/values are produced by the
reproducible scripts in this repository; numbers are read from the saved validation
artifacts, not estimated.*

---

## Abstract
BrainSafe AI is an open computational tool that estimates a compound's profile of brain-relevant
effects from its chemical structure alone. It brings together eight machine-learning endpoints, every
one trained on **measured** public bioactivity data (ChEMBL pChEMBL values and the B3DB
blood–brain-barrier database): blood–brain-barrier (BBB) penetration, inhibition of AChE, BChE,
BACE1, GSK-3β, MAO-A and MAO-B, and the hERG cardiotoxicity liability. Four further receptor targets
(D2, A2A, 5-HT2A and SERT) are handled instead as potency regressions. A deterministic
druggability/CNS-MPO layer, a measured-data antioxidant (DPPH) model, and a clinical-precedent layer
of 504 nervous-system compounds with real clinical-phase data round out the system. Every prediction
is **isotonic-calibrated**, carries a **Mondrian conformal prediction set** with empirically verified
~90% coverage, is **grounded in the nearest real measured analogues**, and feeds into **BBB-gated
per-disease scores**. We validated the models under a hierarchy of increasingly demanding splits:
random, scaffold and temporal. On like-for-like random splits the classifiers reach AUROC 0.94–0.98,
at or above published state of the art; under the scaffold split 0.87–0.96 (mean 0.92); and under a
true temporal (future-compound) split 0.61–0.91, which shows plainly where generalisation runs out.
The deployed model is a **random forest**, selected after a like-for-like comparison against
gradient-boosted trees (XGBoost, histogram gradient boosting) and a graph neural network, and it
beats a k-nearest-neighbour Tanimoto read-across (mean AUROC 0.919 vs 0.867) and logistic regression
(0.808) on every endpoint. We also ran a pre-registered head-to-head against four general-purpose large
language models: on well-known drugs they matched or exceeded BrainSafe at BBB and hERG
classification, yet 45% of the measured-data identifiers they cited were fabricated or mis-attributed
and all four invented a target and potency for an unpublished compound, whereas BrainSafe returned,
for every structure, a calibrated probability, a conformal set, and the nearest measured analogue,
the kind of grounded, auditable output an LLM does not supply. The contribution here is the
**integration**, calibrated, evidence-grounded, BBB-gated and safety-aware CNS profiling from measured
data, rather than any one new algorithm.

### Note added in revision (2026-07-21): data expansion and a Random-forest, ten-fold model
In response to review we (i) added a second independent measured source, **BindingDB**, to the eleven
protein-target endpoints, pooling it with ChEMBL at the compound level. The training corpus grew to
**61,317 unique compounds across 67,982 measured records** (18,573 compounds are measured by both
databases; 4,246 are contributed by BindingDB alone). (ii) We report a **Random forest trained under
ten-fold cross-validation** for every endpoint, under both random and scaffold-grouped splits; mean
classifier AUROC is **0.960 (random)** and **0.919 (scaffold)**, and receptor-regression R² is
0.60–0.68 (random) and 0.39–0.58 (scaffold). (iii) Because the BindingDB export contributes actives,
we audited the addition against the ChEMBL-only baseline: the scaffold-split headline changed by a
mean of **−0.0002**, i.e. the added data neither inflates nor degrades performance. (iv) We
strengthened two methodology points: isotonic **probability calibration** (mean expected calibration
error 0.072 → 0.012) and an explicit **applicability-domain** flag (nearest-neighbour Tanimoto), which
shows that arbitrary drug libraries are largely outside the target-specific models' domain and must be
flagged as extrapolation. (v) An **external test** of the BBB model on 306 FDA-curated approved drugs
absent from training gives AUROC **0.774**. Full detail and every supporting table are in
`docs/RF_CV_RESULTS.md`, `docs/METHODOLOGY_AUDIT.md` and `docs/DATA_MANIFEST.md`.

### Note added in revision (2026-07-22): an ADME / exposure layer toward K_p,uu
Target engagement and BBB penetration do not, on their own, say whether an achievable dose delivers
*free* drug to a brain target. We therefore added nine measured ADME / exposure endpoints, trained with
the identical random-forest, scaffold-10-fold protocol (~23,000 measured compounds from TDC / Harvard
Dataverse, MoleculeNet, B3DB and ChEMBL): aqueous solubility, lipophilicity, Caco-2 permeability,
P-glycoprotein inhibition and substrate (efflux), plasma-protein binding, hepatocyte clearance, logBB,
and, directly, the unbound brain-to-plasma partition coefficient **K_p,uu** (566 measured compounds,
ChEMBL `K(p,uu,brain)`; scaffold R² 0.35). Scaffold performance ranges from strong (P-gp inhibition
AUROC 0.937, solubility R² 0.76) to weak and disclosed (clearance R² 0.19). The directly-modelled
K_p,uu drives a combined free-brain-exposure readout (K_p,uu ≥ 0.3 taken as meaningful exposure) that,
on known drugs, correctly separates central compounds (diazepam K_p,uu 0.94, donepezil 0.84) from
peripheral or P-gp-effluxed ones (atenolol 0.07, loperamide 0.04). This moves the tool from *"does it
bind"* toward *"does an achievable dose reach the target in the brain"*. Detail: `docs/ADME_RESULTS.md`.

## 1. Introduction
Deciding whether a small molecule is likely to act on the brain means answering several questions
at once. Can it cross the BBB? Does
it engage disease-relevant CNS targets? Is it developable, is it safe, and has anything like it
reached the clinic? The web tools in common use each answer part of this. General ADMET servers
(SwissADME, ADMETlab, admetSAR, pkCSM) cover physicochemistry and safety; target-prediction tools
(SwissTargetPrediction, PPB2) suggest likely protein targets. None of them, though, pulls
measured-data CNS target activity, BBB gating, calibrated uncertainty, safety and clinical
precedent into one transparent, evidence-grounded tool. That integration is what BrainSafe AI
provides.

## 2. Methods

### 2.1 Data sources (all measured, public)
- **CNS target bioactivity:** two independent public databases, **ChEMBL 37** (REST API) and
  **BindingDB**, pooled at the compound level for eleven human targets, AChE (CHEMBL220), BChE
  (CHEMBL1914), BACE1 (CHEMBL4822), GSK-3β (CHEMBL262), MAO-A (CHEMBL1951), MAO-B (CHEMBL2039),
  the hERG safety anti-target (CHEMBL240), and receptors D2 (CHEMBL217), A2A (CHEMBL251),
  5-HT2A (CHEMBL224), SERT (CHEMBL228). Activities with a defined potency were kept (ChEMBL pChEMBL,
  or a BindingDB IC50/Ki/Kd/EC50 in nM converted to the same −log10 molar scale). Of the target
  compounds, 18,573 are measured by both databases and 4,246 are contributed by BindingDB alone;
  hERG remained ChEMBL-only (BindingDB was rate-limited at retrieval).
- **Blood–brain barrier:** B3DB classification dataset (7,807 measured compounds; 7,805 modelled
  after InChIKey deduplication). A further 306 FDA-curated approved drugs absent from B3DB are held
  out as an external test set.
- **Antioxidant:** ChEMBL DPPH radical-scavenging assays; IC50/EC50 → pIC50 (2,862 compounds).
- **Clinical precedent:** ChEMBL ATC level-1 "N" (nervous-system) molecules with max clinical
  phase ≥ 1 and a structure (504 compounds), ATC-mapped to disease class.

After pooling and standardisation the target panel holds **61,317 unique compounds across 67,982
measured compound-endpoint records**. The number of compounds per endpoint (used both to train and,
through cross-validation, to test) is: BBB 7,805; AChE 4,387; BChE 2,621; BACE1 8,501; GSK-3β 4,958;
MAO-A 2,228; MAO-B 3,665; hERG 5,875; D2 7,734; A2A 6,785; 5-HT2A 5,989; SERT 4,572; antioxidant
2,862.

### 2.2 Curation and labelling
SMILES were canonicalised and deduplicated by InChIKey (RDKit). For classification, activities
were aggregated per compound by median pChEMBL and labelled **active = pChEMBL ≥ 6 (≤1 µM)**,
**inactive = pChEMBL < 5 (>10 µM)**, with the 5–6 grey zone discarded to reduce label noise.
Receptor targets (96–98 % active) were unsuitable for binary classification and were instead
modelled as **potency regression** on median pChEMBL. Antioxidant DPPH IC50/EC50 were converted
to pIC50 (−log10 M) and aggregated per compound.

### 2.3 Molecular representation
Each molecule was represented by a **1,024-bit Morgan (ECFP, radius 2) fingerprint** concatenated
with **twelve interpretable RDKit physicochemical descriptors** (molecular weight, cLogP, TPSA,
H-bond donors, H-bond acceptors, rotatable bonds, aromatic rings, FractionCSP3, ring count,
heavy-atom count, formal charge, QED), a fixed **1,036-feature** numeric vector identical for every
endpoint. The fingerprint is collision-free by construction (bit k always denotes the same local
atomic environment), so no free-text or categorical value ever reaches the model.

### 2.4 Model, and why a random forest
Every endpoint is modelled by a **random forest**. A random forest is an ensemble of decision trees:
here **300 CART trees**, each grown on a bootstrap resample of the training compounds, each split
chosen from a random subset of about 32 of the 1,036 features (the square root of the feature count)
so as to reduce Gini impurity (classification) or variance (regression), down to a minimum of two
compounds per leaf; the trees then vote (classification) or average (regression). Classifiers use
balanced class weights. Averaging many decorrelated trees is what makes the forest robust on the
sparse, high-dimensional fingerprint space and gives probabilities that calibrate cleanly. A single
decision tree, by contrast, overfits; the forest is what tames that variance.

The random forest was **not assumed to be best**. It was selected after a like-for-like comparison,
on the same features and the same scaffold split, against gradient-boosted trees (**XGBoost** and
**histogram gradient boosting**) and a **graph neural network** (a Graph Isomorphism Network trained
on the raw molecular graph). The random forest was best or tied on the eight classification endpoints
and within 0.02 R-squared of the best on regression, whereas the graph network, trained from scratch at
this data scale of a few thousand compounds per endpoint, trailed on every endpoint tested (§3.6,
Table 4). The forest is therefore the deployed model: strongest overall, naturally probabilistic and
calibratable, interpretable through feature importance, and fully reproducible.

Whether a target is a classifier at all is decided by a fixed **quality gate: the Matthews correlation
coefficient must reach 0.45 under the scaffold split**. Four receptors (D2, A2A, 5-HT2A, SERT) are
96-98% active in the public data and fail this gate as binary tasks, so they are modelled as potency
**regression** on median pChEMBL instead. The routing is a property of the data, not a manual choice.
Hyper-parameters (300 trees, minimum leaf size two, balanced weights, random_state 42) are held fixed
across all endpoints.

### 2.5 Calibration and conformal prediction
Random-forest probabilities are recalibrated by **isotonic regression** on out-of-fold predictions;
this reduced the mean expected calibration error from **0.072 to 0.012** across the eight classifiers,
so a reported probability can be read as a genuine likelihood. **Mondrian (class-conditional) inductive
conformal prediction** then returns a per-compound prediction set at the 90 % level; empirical coverage
on held-out test compounds is **0.89–0.92** across endpoints (target 0.90) with an average set size
near one.

### 2.6 Validation hierarchy (no single split relied upon)
1. **Random** stratified 10-fold cross-validation (like-for-like with most literature).
2. **Scaffold** GroupKFold(10) on Bemis–Murcko generic scaffolds; whole scaffolds are held out so no
   scaffold is shared between train and test. This is the headline estimate of generalisation to new
   chemistry.
3. **Temporal:** train on compounds reported up to the 75th-percentile ChEMBL document year, test on
   strictly later compounds (a true future-compound test).
An **applicability-domain** flag (nearest-neighbour ECFP Tanimoto, in-domain at ≥ 0.30) and the
conformal sets above qualify every individual prediction. Metrics: AUROC, PR-AUC, balanced accuracy,
MCC, Brier for classification; R²/RMSE/Spearman for regression. Every compound's out-of-fold
prediction and fold assignment is retained for audit.

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

We report results in two complementary modes. The **non-comparative analysis** (§3.1–3.4) looks at
each endpoint on its own terms, measuring the deployed model against its own held-out data:
discrimination, calibration, conformal coverage and prospective behaviour. The **comparative
analysis** (§3.5–3.7) then sets those same models against external reference points: published QSAR
ranges (§3.5), simpler internal baselines run under an identical protocol (§3.6), and the
general-purpose large-language-model paradigm (§3.7).

### 3.1 Classification endpoints: full validation hierarchy (non-comparative)
Random forest, 10-fold cross-validation. AUROC unless stated.

| Endpoint (n) | Random | **Scaffold** | Temporal | Scaffold MCC | Conformal cov. |
|---|---|---|---|---|---|
| BBB (7,805) | 0.960 | 0.920 | – | 0.658 | 0.904 |
| AChE (4,387) | 0.963 | 0.921 | 0.785 | 0.656 | 0.892 |
| BChE (2,621) | 0.968 | 0.937 | 0.737 | 0.691 | 0.890 |
| BACE1 (8,501) | 0.967 | 0.956 | 0.908 | 0.672 | 0.911 |
| GSK-3β (4,958) | 0.969 | 0.937 | 0.657 | 0.559 | 0.915 |
| MAO-B (3,665) | 0.954 | 0.890 | 0.781 | 0.626 | 0.918 |
| MAO-A (2,228) | 0.947 | 0.868 | 0.611 | 0.564 | 0.890 |
| hERG (5,875) | 0.955 | 0.921 | 0.785 | 0.676 | 0.901 |
| **Mean** | **0.960** | **0.919** | 0.752 | 0.638 | 0.902 |

The classifiers hold up under the scaffold split (mean AUROC 0.919), so the performance reflects
generalisation to new scaffolds rather than memorised analogues. Temporal AUROC (train on the past,
test on genuinely later compounds) ranges from 0.61 to 0.91 and is where generalisation is hardest,
most for the smallest, least scaffold-diverse sets (MAO-A, GSK-3β). Conformal coverage is 0.89–0.92
against a 0.90 target. BBB carries no document year and is not temporally split.
*(Figures 1–2: per-endpoint performance and compound counts; Figure 6: calibration.)*

### 3.2 Receptor potency regression (scaffold 10-fold, non-comparative)
Random forest regression on median pChEMBL.

| Receptor (n) | Scaffold R² | RMSE | Spearman | Temporal R² |
|---|---|---|---|---|
| A2A (6,785) | 0.576 | 0.70 | 0.751 | 0.338 |
| 5-HT2A (5,989) | 0.490 | 0.72 | 0.712 | 0.182 |
| D2 (7,734) | 0.483 | 0.67 | 0.692 | 0.042 |
| SERT (4,572) | 0.389 | 0.79 | 0.625 | 0.100 |

The regressors explain 39–58 % of potency variance on new scaffolds, with Spearman rank correlations
of 0.63–0.75; this is the realistic ceiling for structure-only potency prediction on these receptors,
and the temporal figures show, honestly, how much further it falls on future chemistry. SERT is the
weakest and is reported as lower-confidence. *(Predicted-vs-measured scatter: Figure regression.)*

### 3.3 Antioxidant (measured DPPH) and druggability
Measured DPPH regression (n=2,862): scaffold 10-fold **R² = 0.43, RMSE = 0.59, Spearman = 0.65**,
temporal R² 0.01 (a hard future-compound test). Druggability/CNS-MPO is deterministic and
discriminates CNS drugs from polar non-drugs (donepezil 79, caffeine 86 vs sucrose 46,
atorvastatin 24).

### 3.4 Prospective sanity (chemistry-only inputs, non-comparative)
With chemistry-only inputs the integrated system reproduces established pharmacology (values are the
tool's actual outputs; `BS_llm_benchmark_groundtruth.json`). Donepezil → Alzheimer's disease
(score 1.00) via AChE (P=1.00), BBB-penetrant (0.997), with clinical precedent (Donepezil, Phase 4,
Alzheimer's); Rivastigmine → Alzheimer's (0.85) via BChE (the tool identifies the butyrylcholinesterase
route), BBB-penetrant; Rasagiline → Parkinson's disease (0.95) via MAO-B, clinical precedent
(Rasagiline, Phase 4, Parkinson's); Fluoxetine → Depression (0.96) via SERT with clinical precedent
(Fluoxetine, Phase 4, Depression); Terfenadine → **hERG 1.00** (withdrawn for cardiotoxicity, correctly
flagged) and BBB non-penetrant (0.42); Resveratrol and Quercetin → BBB non-penetrant (0.35 and 0.18)
with the highest antioxidant scores in the set. Donepezil's hERG is flagged (P=0.78); we note its
clinical hERG relevance is modest relative to therapeutic exposure and therefore treat it cautiously
(§3.7d). The system
reproduces known pharmacology and safety.

### 3.5 Benchmark vs literature (comparative)
On random splits, BrainSafe AUROC (0.94–0.98) is at/above published ranges (BBB 0.88–0.96;
hERG 0.86–0.93). The same models additionally report the stricter scaffold/temporal
numbers most studies omit. *(Figure 6; Supplementary Table S7.)*

### 3.6 Model comparison: why a random forest (comparative)
The random forest was chosen on evidence, not preference. Three model families were trained on the
identical features and evaluated under the identical scaffold split, and, separately, against a graph
neural network on the raw molecular graph (Table 4).

**Table 4. Model comparison (mean scaffold-split performance).**

| Model | Classification (mean AUROC) | Regression (mean R²) |
|---|---|---|
| **Random forest (deployed)** | **0.914** | 0.461 |
| XGBoost | 0.905 | **0.481** |
| Histogram gradient boosting | 0.901 | 0.478 |
| Graph neural network (GIN)* | lower on all tested | lower on all tested |

*The Graph Isomorphism Network was trained from scratch on a single scaffold hold-out for four
representative endpoints and trailed the random forest on every one: BBB 0.887 vs 0.924, BACE1 0.924
vs 0.957, MAO-A 0.736 vs 0.810 (AUROC), and A2A 0.467 vs 0.549 (R²).

**Why the alternatives do not win.** Gradient boosting ties the random forest to within 0.02: on a
few thousand compounds described by a sparse 1,024-bit fingerprint plus a dozen descriptors, the
forest's bagging already controls variance well and boosting's extra capacity yields no reliable gain
(it is marginally ahead on the receptor regressions, marginally behind on classification). The graph
network must *learn* a molecular representation from the data, which at this scale it cannot do better
than the expert-designed fingerprint; graph networks overtake fixed fingerprints only with far larger
datasets or self-supervised pretraining, neither of which applies here. The random forest is therefore
the best overall model and, unlike the alternatives, yields naturally calibratable probabilities. It
also beats the simplest non-learning baselines by a wide margin, a k-nearest-neighbour Tanimoto
read-across (mean scaffold AUROC 0.867) and L2-logistic regression (0.808), on every endpoint,
confirming it has learned structure–activity relationships that do not reduce to retrieving the most
similar known molecule. *(Full per-endpoint numbers: `results/tables/model_comparison.csv`,
`results/gnn/gnn_vs_rf.csv`.)*

### 3.7 Comparison with general-purpose large language models (comparative)
A reviewer asked why a dedicated tool is needed when a general-purpose large language model
(LLM) can be queried in natural language for similar information. We address this with (a) the
peer-reviewed benchmark evidence, (b) an architectural/scientific account of the difference,
and (c) a reproducible demonstration of grounded output.

**(a) Benchmark evidence.** On molecular property prediction (the task class BrainSafe
performs), general-purpose LLMs consistently *underperform* specialised machine-learning
models. In the eight-task chemistry benchmark of Guo *et al.* (2023), LLMs including GPT-4
lag task-specific models on property-prediction tasks and struggle to parse SMILES reliably;
Zhong *et al.* (2024) report that "LLMs generally lag behind ML models" on molecule property
tasks, especially where molecular geometry matters; and Jablonka *et al.* (2024) find that even
*fine-tuned* LLMs become competitive with dedicated QSAR models only in the low-data limit, not
at the data scale used here (67,982 measured records). LLMs also exhibit documented
factual hallucination in generative settings, a failure mode surveyed comprehensively by
Ji *et al.* (2023).

**(b) Scientific background: why the paradigms differ.** A general LLM is an autoregressive
next-token predictor over text. It does not compute a molecular fingerprint, does not fit an
explicit function from chemical structure to *measured* bioactivity, and does not emit a
probability with a calibration or coverage guarantee. BrainSafe, by contrast, encodes each
molecule as an ECFP-4 fingerprint plus twelve physicochemical descriptors, learns a
structure→measured-activity mapping from measured ChEMBL, BindingDB and B3DB records, and returns a
**calibrated** probability wrapped in a **conformal set with empirically verified ~90%
coverage**, together with the **nearest real measured analogue and its measured pChEMBL**.
These are precisely the guarantees an LLM cannot provide. A full capability matrix is given in
Supplementary Table S8.

**(c) Reproducible grounded-output demonstration.** For fixed input structures the deployed
engine returns verifiable artifacts (script `BS_llm_comparison.py`, output
`BS_llm_comparison.json`): donepezil → AChE calibrated P = 1.00 with the nearest measured analogue
at Tanimoto 1.00 (donepezil is itself a measured training compound, pChEMBL 7.75) and hERG P = 0.78;
terfenadine → hERG P = 1.00, correctly flagging the cardiotoxicity for which it was withdrawn
while also calling it BBB-non-penetrant; and a novel arylpiperazine of an unpublished
scaffold → an **honest conformal "uncertain" set** for AChE grounded in the nearest measured analogue
(Tanimoto 0.35, pChEMBL 4.82), rather than a confident but unverifiable text answer. Every value is traceable
to a measurement; none is generated from free text. This grounding (a calibrated probability,
a coverage-guaranteed set, and measured-analogue provenance for *any* structure, including novel
ones) is the scientific justification for a dedicated tool over an LLM. *(Supplementary
Tables S8–S9.)*

**(d) Pre-registered head-to-head benchmark (executed).** To make the comparison directly measurable we
froze a fixed prompt, a 10-compound panel with uncontested ground truth (approved-drug pharmacology plus
one unpublished scaffold), and a scoring rubric *before* any system was run
(`BS_LLM_benchmark_protocol.md`; key in `BS_llm_benchmark_groundtruth.json`; scorer `BS_llm_score.py`).
The identical prompt was then run on four general-purpose LLMs (Gemini Pro, ChatGPT/GPT-4o, Perplexity,
Claude) and each reply scored against the same measured-data key (Table below; Supplementary Table S13).

| System | BBB acc. | hERG acc. | Brier | ChEMBL IDs cited | Fabricated/wrong-molecule | Novel-compound confabulation |
|---|---|---|---|---|---|---|
| **BrainSafe AI** | 8/9 (0.889) | 5/5 (1.00) | 0.088 | 0* | **0** | **No (honest "uncertain")** |
| Gemini Pro | 7/9 (0.778) | 5/5 (1.00) | 0.067 | 10 | 5 (50%) | Yes |
| ChatGPT (GPT-4o) | 9/9 (1.00) | 5/5 (1.00) | 0.035 | 9 | 4 (44%) | Yes |
| Perplexity | 9/9 (1.00) | 4/5 (0.80) | 0.041 | 2 | 1 (50%) | Yes |
| Claude | 9/9 (1.00) | 5/5 (1.00) | 0.020 | 10 | 4 (40%) | Yes |

*BrainSafe reports the nearest measured analogue by structure (SMILES + measured pChEMBL), not a recalled
ChEMBL identifier, so it cannot fabricate one. hERG scored on the five uncontested compounds; donepezil
and fluoxetine were excluded a priori because their hERG clinical relevance is genuinely debatable.

Two findings stand out, and we report both honestly. **First, on classification of well-known approved
drugs the LLMs are strong.** Three of four matched or exceeded BrainSafe on BBB (9/9 vs 8/9; BrainSafe
mis-called the borderline-lipophilic astemizole), and their probability Brier scores were as good or
better (Claude 0.020, ChatGPT 0.035). This is expected: these are textbook molecules richly described in
the training corpus, so recall is excellent, and it means a dedicated tool is *not* justified by raw
accuracy on famous compounds. **Second, the LLMs fail exactly where grounding and novelty matter.**
Across the four models, **14 of 31 (45%) of the ChEMBL identifiers they volunteered as provenance were
fabricated or resolved to the wrong molecule.** For example, Gemini's cited "rasagiline" identifier is in fact
*fluticasone propionate* and its "selegiline" identifier is *propranolol*; Claude's "rivastigmine"
identifier is *pyridoxine* and its "terfenadine" identifier is the antibiotic *cefdinir*; ChatGPT's
"terfenadine" identifier actually points to astemizole. And **all four models confabulated on the
unpublished compound**, each asserting a specific target and potency (and even disagreeing on the
target, with three saying AChE and one the D2 receptor) for a structure that has no measured value.
BrainSafe fabricated nothing, grounded every prediction in a real measured analogue, and returned an
honest conformal "uncertain" set for the novel compound. The scientific implication is precise: an LLM
can approximate textbook classifications but cannot be trusted for *verifiable provenance* or for
*novel chemistry*, which is where hypothesis generation actually happens, and precisely what the
dedicated tool provides.

## 4. Discussion

**Adversarial validation.** Beyond the cross-validation above, the tool was validated by inversion,
that is, by trying to make it wrong. Six failure modes were tested directly: scaffold leakage (none;
no scaffold is shared between train and test), duplicate inflation (none; 61,317 rows, 61,317 unique
InChIKeys), non-reproducibility (none; retraining reproduces the reported score exactly at seed 42),
degenerate constant prediction (absent; BBB probabilities span 0.00–1.00), misranking of known
chemistry (absent; BBB and K_p,uu both rank central drugs above peripheral ones), and confident
extrapolation on alien chemistry (guarded; a fluorosurfactant is flagged out of domain at Tanimoto
0.20). All six pass; the procedure and outputs are in `docs/VALIDATION.md` and
`results/tables/inversion_validation.csv`. Together with the audited data additions and the honestly
reported temporal decline, this is the basis for treating the numbers as fair rather than flattering.

**Contribution.** The individual pieces here are standard: ECFP/RF models, BBB/hERG/target QSAR,
conformal prediction, the QED/CNS-MPO rules. What we have not seen assembled elsewhere is all of them
working as one transparent, measured-data CNS profiler that is at once calibrated, conformal,
evidence-grounded, BBB-gated, safety-aware and clinically contextualised. The assembly is the point,
not any single part.

**Why a dedicated tool rather than a general-purpose LLM?** The head-to-head in §3.7 answers this with
results, and the answer has two halves. On the task itself, molecular property prediction, general
LLMs are known to trail specialised ML (Guo *et al.*, 2023; Zhong *et al.*, 2024); fine-tuning closes
the gap only when data are scarce (Jablonka *et al.*, 2024), and BrainSafe works from 67,982 measured
records. That much was expected. The part that matters in practice is what our benchmark exposed: four
current models matched or beat BrainSafe at classifying *well-known* drugs, yet 45% of the ChEMBL
identifiers they offered as evidence were fabricated or pointed to the wrong molecule, and every one of
them invented a target and potency for an *unpublished* compound, disagreeing with each other along the
way. The deeper reason is architectural. An LLM produces fluent text; it does not compute a fingerprint,
fit structure to measured activity, or attach a calibrated probability, a coverage guarantee or a domain
boundary, and it will hallucinate with confidence (Ji *et al.*, 2023). BrainSafe returns, for *any*
structure, a calibrated probability, a conformal set with roughly 90% empirical coverage, an explicit
in- or out-of-domain flag, and the nearest *measured* analogue with its pChEMBL. The two are
complementary, not interchangeable, wherever a decision has to be auditable and grounded in measurement.

**Threats to validity (scientific-flaw self-audit, with quantitative tests).** We enumerated the
model's principal methodological risks and, rather than merely noting them, ran targeted analyses
to bound each (script `BS_flaw_fixes.py` / `BS_assay_composition.py` / `BS_assay_sensitivity.py`;
`BS_flaw_fixes.json`).
(1) *Assay heterogeneity.* Activities pooled per target span IC50/Ki/Kd/EC50/Potency. We first
quantified the composition (Supplementary Table S11): IC50 is dominant for every target
(81–92 %) except GSK-3β, which is genuinely mixed (IC50 49 %, EC50 33 %, Ki 16 %). We then
retrained the deployed model under scaffold CV on the **dominant single assay type (IC50)
only** versus the pooled set: scaffold AUROC changed by **≤0.006** for all three endpoints tested,
including the most heterogeneous, GSK-3β (pooled 0.919 vs IC50-only 0.913; MAO-B −0.006; hERG 0.000;
Supplementary Table S12). Pooling on the standardised pChEMBL scale therefore does not materially
distort discrimination.
(2) *Label-threshold sensitivity.* The pChEMBL ≥6/<5 cut is a modelling choice, so we re-derived
labels from the raw activity records at alternative definitions and re-measured scaffold-CV AUROC
with the deployed model (Supplementary Table S10). Across four endpoints and four definitions
the maximum AUROC spread was **0.109**; the stricter ≥6.5/<5.5 cut sat within 0.01–0.02 of the
deployed cut, and the "sharp boundary" cut that keeps the 5–6 grey zone was consistently the
*worst*, empirically validating the decision to discard it. Per-operating-threshold
precision/recall/F1 are additionally in Supplementary Table S4.
(3) *Applicability-domain cut-off.* The Tanimoto AD threshold is empirically, not arbitrarily,
set: the n-weighted similarity-binned AUROC falls monotonically from 0.958 (nearest neighbour
Tanimoto ≥0.8) to 0.939, 0.866 and **0.770** (<0.4) (Supplementary Table S5), justifying the
out-of-domain flag in the 0.3–0.4 band.
(4) *Disease mapping.* The target→disease synthesis is a transparent, knowledge-based rule
(not a learned layer), each driver tagged by provenance, so it can be inspected and overridden.
(5) *Single safety anti-target.* Cardiotoxicity is represented by hERG alone; other liabilities
(e.g. Nav1.5, hepatotoxicity) are out of scope and stated as such.
(6) *Read-across ceiling.* Because the random forest beats a kNN-Tanimoto baseline on every endpoint
(§3.6), performance is not merely memorised nearest-neighbour recall. None of these is concealed;
each is surfaced in the tool output or the supplementary tables.

**Validation honesty.** Reporting all four splits is a deliberate choice. The drop from random
(0.94–0.98) to temporal (0.61–0.92) puts a number on how hard genuine prospective prediction is:
between 71 % and 91 % of the recent compounds carry scaffolds the model never saw in training. A high
temporal AUROC is not always the good news it looks like: BACE1's 0.92 owes something to a recent
test set that is 93 % active, whereas MAO-A's balanced set (45 % active) yields a plainer 0.61. We
would rather show these numbers than bury them.

**Limitations (explicit and, where inherent, unfixable computationally).**
(i) Models predict **engagement/binding, not direction** (agonist vs antagonist).
(ii) **Engagement is not efficacy**; the clinical layer provides *precedent from real trial
data*, not an efficacy prediction.
(iii) **No wet-lab prospective validation** has been performed; that step requires experiments.
(iv) **Temporal generalisation to novel scaffolds is bounded** by covariate shift; receptor
and pooled-assay (DPPH) endpoints generalise across time only weakly and are flagged.
(v) GSK-3β and MAO-A degrade temporally and are marked lower-confidence.
(vi) **Bioactivity is pooled across assay formats** (IC50/Ki/Kd/EC50) on the pChEMBL scale;
residual cross-assay variance is bounded but not eliminated (see Threats to validity).
(vii) **A single cardiac safety anti-target (hERG)** is modelled; other safety liabilities are
out of scope.

**Intended use.** Research hypothesis-generation, triage and prioritisation, not clinical or
diagnostic use.

## 5. Conclusion
BrainSafe AI is a calibrated, evidence-grounded, multi-endpoint CNS profiler built entirely on
measured public data. Its per-endpoint performance is state-of-the-art-grade on like-for-like splits,
and we place the harder scaffold and temporal numbers alongside rather than out of sight. We
have done everything that can be validated computationally; the gaps that remain (predicting
efficacy, telling agonism from antagonism, confirming results at the bench) are inherent to the
approach, and we say so plainly. On that footing the tool stands both as a resource publication and as
something researchers can actually use.

## Data and code availability
All models (`models_brain/`, `models_brain_reg/`, `models_genuine/`), datasets
(`data/endpoints/`, `data/endpoints_reg/`, `data/clinical_cns_reference.csv`), validation
reports (`*_report.json`, `BS_randomsplit_benchmark.json`), the model card (`BS_MODEL_CARD.md`),
and all fetch/train/validation scripts are in the repository. App: `app_v6_final.py`.

## Supplementary materials
**Figures** (`figures/`, 300 dpi, regenerated from out-of-fold predictions):
Fig 1 pipeline/workflow · Fig 2 dataset size & class balance · Fig 3 AUROC across validation
regimes (random/scaffold/temporal) · Fig 4 (A) scaffold-CV ROC curves, (B) calibration
reliability · Fig 5 (A) conformal coverage, (B) random forest vs kNN/logistic baselines ·
Fig 6 benchmark vs literature · Fig 7 predicted-vs-measured regression scatter (antioxidant + receptors).
**Supplementary tables** (`supplementary/`, exact values from the validation artifacts):
S1 classification metrics (random/scaffold/temporal, PR-AUC, BA, MCC, Brier, conformal) ·
S2 receptor regression · S3 antioxidant (measured DPPH) · S4 threshold sensitivity ·
S5 similarity-binned generalisation · S6 clinical-reference composition · S7 benchmark vs literature ·
S8 BrainSafe-vs-LLM capability matrix · S9 ablation vs kNN-Tanimoto and logistic-regression baselines ·
S10 label-definition robustness (scaffold AUROC under alternative pChEMBL cuts) ·
S11 assay-type composition per endpoint · S12 single-assay (IC50-only) vs pooled sensitivity.
**Model card:** `BS_MODEL_CARD.md` (full provenance, diagnosis, and limitations).
**Comparison artifacts:** `BS_llm_comparison.py`/`BS_llm_comparison.json` (baseline summary +
grounded-output demonstration); `BS_flaw_fixes.py`/`BS_flaw_fixes.json`,
`BS_assay_composition.py`, `BS_assay_sensitivity.py` (threats-to-validity analyses);
`BS_LLM_benchmark_protocol.md` + `BS_llm_benchmark.py` (pre-registered LLM head-to-head).

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
12. Guo T, Guo K, Nan B, Liang Z, Guo Z, Chawla NV, Wiest O, Zhang X. What can large language models do in chemistry? A comprehensive benchmark on eight tasks. *NeurIPS 2023 Datasets and Benchmarks Track*. arXiv:2305.18365.
13. Zhong Z, Zhou K, Mottin D. Benchmarking large language models for molecule prediction tasks. 2024. arXiv:2403.05075.
14. Jablonka KM, Schwaller P, Ortega-Guerrero A, Smit B. Leveraging large language models for predictive chemistry. *Nat Mach Intell* 2024;6:161–169.
15. Ji Z, Lee N, Frieske R, Yu T, Su D, Xu Y, et al. Survey of hallucination in natural language generation. *ACM Comput Surv* 2023;55(12):Article 248.
