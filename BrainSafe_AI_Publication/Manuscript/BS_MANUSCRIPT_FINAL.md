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
0.61–0.92, transparently exposing where generalisation is limited. Under an identical protocol the deployed ensemble outperforms both a k-nearest-neighbour
Tanimoto read-across (mean AUROC 0.912 vs 0.867) and logistic regression (0.808) on every
endpoint. In a pre-registered head-to-head, four general-purpose large language models matched or
exceeded BrainSafe on blood–brain-barrier/hERG classification of well-known drugs but fabricated or
mis-attributed 45% of the measured-data identifiers they cited and confabulated a specific target and
potency for an unpublished compound, whereas BrainSafe returns, for any structure, a calibrated
probability, a conformal set, and the nearest measured analogue — grounded, auditable output an LLM
does not provide. The methodological
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
- **Blood–brain barrier** — B3DB classification dataset (7,807 measured compounds; 7,805 modelled
  after InChIKey deduplication).
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

Results are reported in two complementary modes requested for a rigorous evaluation.
**Non-comparative analysis** (§3.1–3.4) reports each endpoint's *standalone* performance
against its own held-out measured data — the absolute discrimination, calibration,
conformal coverage and prospective behaviour of the deployed models. **Comparative
analysis** (§3.5–3.7) benchmarks the same models *against external references*: published
QSAR performance ranges (§3.5), simpler internal baselines under an identical protocol
(§3.6), and the general-purpose large-language-model paradigm (§3.7).

### 3.1 Classification endpoints — full validation hierarchy (non-comparative)
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
*(Figure 3: validation hierarchy; Figure 4A: scaffold-CV ROC curves; Figure 4B: calibration
reliability; Figure 5A: conformal coverage; Figure 2: dataset size/balance. Full metrics:
Supplementary Table S1; per-threshold precision/recall/F1: S4; similarity-binned generalisation: S5.)*

### 3.2 Receptor potency regression (scaffold-CV, non-comparative)
| Receptor (n) | R² | RMSE | Spearman | Temporal R² |
|---|---|---|---|---|
| A2A (5,547) | 0.526 | 0.74 | 0.706 | 0.326 |
| 5-HT2A (5,256) | 0.460 | 0.76 | 0.684 | 0.085 |
| D2 (7,511) | 0.425 | 0.71 | 0.652 | −0.007 |
| SERT (4,471) | 0.338 | 0.84 | 0.573 | 0.171 |

*(Predicted-vs-measured scatter: Figure 7; Supplementary Table S2.)*

### 3.3 Antioxidant (measured DPPH) and druggability
Measured DPPH regression (n=2,862): scaffold-CV **R² = 0.43, RMSE = 0.60, Spearman = 0.636**
(vs curated R²≈0.25). The prior curated score correlated only weakly with measured DPPH
(Spearman 0.39), confirming the measured model as the superior basis. Druggability/CNS-MPO is
deterministic and discriminates CNS drugs from polar non-drugs (donepezil 79, caffeine 86 vs
sucrose 46, atorvastatin 24).

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
hERG 0.86–0.93). The same models additionally report the stricter scaffold/cluster/temporal
numbers most studies omit. *(Figure 6; Supplementary Table S7.)*

### 3.6 Ablation against simpler baselines (comparative)
To confirm the ensemble contributes signal beyond structural look-up, we compared it, under an
identical scaffold-split protocol and feature set, against (i) a **k-nearest-neighbour Tanimoto
"read-across" baseline** — the closest algorithmic analogue to associative structure recall —
and (ii) **L2-regularised logistic regression**. Across the eight deployed classification
endpoints the deployed ensemble attains a mean scaffold-split AUROC of **0.912**, versus
**0.867** for kNN-Tanimoto (mean Δ = **+0.045**) and **0.808** for logistic regression
(mean Δ = **+0.104**); the ensemble is best on **every one of the eight endpoints**
(Table below; Supplementary Table S9). That the ensemble consistently exceeds a pure
nearest-neighbour read-across shows it captures structure–activity relationships that are
not reducible to "find the most similar known molecule."

| Endpoint | Ensemble | kNN-Tanimoto | Logistic reg. | Δ vs kNN |
|---|---|---|---|---|
| BBB | 0.921 | 0.882 | 0.792 | +0.039 |
| AChE | 0.915 | 0.874 | 0.828 | +0.041 |
| BChE | 0.937 | 0.899 | 0.866 | +0.038 |
| BACE1 | 0.950 | 0.885 | 0.883 | +0.065 |
| GSK-3β | 0.920 | 0.842 | 0.831 | +0.078 |
| MAO-B | 0.885 | 0.833 | 0.742 | +0.052 |
| MAO-A | 0.867 | 0.853 | 0.731 | +0.014 |
| hERG | 0.901 | 0.871 | 0.792 | +0.030 |
| **Mean** | **0.912** | **0.867** | **0.808** | **+0.045** |

### 3.7 Comparison with general-purpose large language models (comparative)
A reviewer asked why a dedicated tool is needed when a general-purpose large language model
(LLM) can be queried in natural language for similar information. We address this with (a) the
peer-reviewed benchmark evidence, (b) an architectural/scientific account of the difference,
and (c) a reproducible demonstration of grounded output.

**(a) Benchmark evidence.** On molecular property prediction — the task class BrainSafe
performs — general-purpose LLMs consistently *underperform* specialised machine-learning
models. In the eight-task chemistry benchmark of Guo *et al.* (2023), LLMs including GPT-4
lag task-specific models on property-prediction tasks and struggle to parse SMILES reliably;
Zhong *et al.* (2024) report that "LLMs generally lag behind ML models" on molecule property
tasks, especially where molecular geometry matters; and Jablonka *et al.* (2024) find that even
*fine-tuned* LLMs become competitive with dedicated QSAR models only in the low-data limit, not
at the data scale used here (64,474 measured records). LLMs also exhibit documented
factual hallucination in generative settings, a failure mode surveyed comprehensively by
Ji *et al.* (2023).

**(b) Scientific background — why the paradigms differ.** A general LLM is an autoregressive
next-token predictor over text. It does not compute a molecular fingerprint, does not fit an
explicit function from chemical structure to *measured* bioactivity, and does not emit a
probability with a calibration or coverage guarantee. BrainSafe, by contrast, encodes each
molecule as an ECFP-4 fingerprint plus 24 physicochemical descriptors, learns a
structure→measured-activity mapping from 64,474 ChEMBL/B3DB records, and returns a
**calibrated** probability wrapped in a **conformal set with empirically verified ~90%
coverage**, together with the **nearest real measured analogue and its measured pChEMBL**.
These are precisely the guarantees an LLM cannot provide. A full capability matrix is given in
Supplementary Table S8.

**(c) Reproducible grounded-output demonstration.** For fixed input structures the deployed
engine returns verifiable artifacts (script `BS_llm_comparison.py`, output
`BS_llm_comparison.json`): donepezil → AChE calibrated P = 1.00 with the nearest measured analogue
at Tanimoto 1.00 (donepezil is itself a measured training compound, pChEMBL 7.75) and hERG P = 0.78;
terfenadine → hERG P = 1.00, correctly flagging the cardiotoxicity for which it was withdrawn,
while correctly calling it BBB-non-penetrant; and a novel arylpiperazine of an unpublished
scaffold → an **honest conformal "uncertain" set** for AChE grounded in the nearest measured analogue
(Tanimoto 0.35, pChEMBL 4.82), rather than a confident but unverifiable text answer. Every value is traceable
to a measurement; none is generated from free text. This grounding — calibrated probability,
coverage-guaranteed set, and measured-analogue provenance for *any* structure, including novel
ones — is the scientific justification for a dedicated tool over an LLM. *(Supplementary
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
drugs the LLMs are strong** — three of four matched or exceeded BrainSafe on BBB (9/9 vs 8/9; BrainSafe
mis-called the borderline-lipophilic astemizole), and their probability Brier scores were as good or
better (Claude 0.020, ChatGPT 0.035). This is expected: these are textbook molecules richly described in
the training corpus, so recall is excellent, and it means a dedicated tool is *not* justified by raw
accuracy on famous compounds. **Second, the LLMs fail exactly where grounding and novelty matter.**
Across the four models, **14 of 31 (45%) of the ChEMBL identifiers they volunteered as provenance were
fabricated or resolved to the wrong molecule** — e.g. Gemini's cited "rasagiline" identifier is in fact
*fluticasone propionate* and its "selegiline" identifier is *propranolol*; Claude's "rivastigmine"
identifier is *pyridoxine* and its "terfenadine" identifier is the antibiotic *cefdinir*; ChatGPT's
"terfenadine" identifier actually points to astemizole. And **all four models confabulated on the
unpublished compound**, each asserting a specific target and potency (and even disagreeing on the
target — three said AChE, one said the D2 receptor) for a structure that has no measured value.
BrainSafe fabricated nothing, grounded every prediction in a real measured analogue, and returned an
honest conformal "uncertain" set for the novel compound. The scientific implication is precise: an LLM
can approximate textbook classifications but cannot be trusted for *verifiable provenance* or for
*novel chemistry* — which is where hypothesis generation actually happens, and precisely what the
dedicated tool provides.

## 4. Discussion

**Contribution.** Individual components (ECFP/RF ensembles, BBB/hERG/target QSAR, conformal
prediction, QED/CNS-MPO) are standard and not novel. The contribution is their **integration**
into a single, transparent, measured-data CNS profiler that is calibrated, conformal,
evidence-grounded, BBB-gated, safety-aware, and clinically contextualised — a configuration no
existing single tool provides.

**Why a dedicated tool rather than a general-purpose LLM?** The comparative evidence (§3.7)
gives a results-backed answer. First, on the exact task class in question — molecular property
prediction — general LLMs are documented to *underperform* specialised ML (Guo *et al.*, 2023;
Zhong *et al.*, 2024), with fine-tuned LLMs matching QSAR only in the low-data limit
(Jablonka *et al.*, 2024), whereas BrainSafe operates at a 64,474-record data scale. Second, and
more fundamentally, the two paradigms differ in kind: an LLM emits fluent text without a
calibrated probability, a coverage guarantee, an applicability-domain boundary, or a link to a
specific measurement, and is prone to confident hallucination (Ji *et al.*, 2023). Our executed
head-to-head (§3.7d) makes this concrete: four current LLMs matched or beat BrainSafe on BBB/hERG
classification of *well-known* drugs, yet 45% of the ChEMBL identifiers they cited as evidence were
fabricated or resolved to the wrong molecule, and all four confabulated a specific target and potency
for an *unpublished* compound (disagreeing with one another). The lesson is not that LLMs are
inaccurate in general — on memorised compounds they are strong — but that they cannot be trusted for
verifiable provenance or for novel chemistry. BrainSafe
returns, for *any* structure, a calibrated probability, a conformal set with empirically verified
~90% coverage, an explicit in/out-of-domain flag, and the nearest *measured* analogue with its
pChEMBL. It is therefore complementary to, not replaced by, general LLMs where a decision must be
auditable and grounded in measured data.

**Threats to validity (scientific-flaw self-audit, with quantitative tests).** We enumerated the
model's principal methodological risks and, rather than merely noting them, ran targeted analyses
to bound each (script `BS_flaw_fixes.py` / `BS_assay_composition.py` / `BS_assay_sensitivity.py`;
`BS_flaw_fixes.json`).
(1) *Assay heterogeneity.* Activities pooled per target span IC50/Ki/Kd/EC50/Potency. We first
quantified the composition (Supplementary Table S11): IC50 is dominant for every target
(81–92 %) except GSK-3β, which is genuinely mixed (IC50 49 %, EC50 33 %, Ki 16 %). We then
retrained the deployed ensemble under scaffold CV on the **dominant single assay type (IC50)
only** versus the pooled set: scaffold AUROC changed by **≤0.006** for all three endpoints tested,
including the most heterogeneous, GSK-3β (pooled 0.919 vs IC50-only 0.913; MAO-B −0.006; hERG 0.000;
Supplementary Table S12). Pooling on the standardised pChEMBL scale therefore does not materially
distort discrimination.
(2) *Label-threshold sensitivity.* The pChEMBL ≥6/<5 cut is a modelling choice, so we re-derived
labels from the raw activity records at alternative definitions and re-measured scaffold-CV AUROC
with the deployed ensemble (Supplementary Table S10). Across four endpoints and four definitions
the maximum AUROC spread was **0.109**; the stricter ≥6.5/<5.5 cut sat within 0.01–0.02 of the
deployed cut, and the "sharp boundary" cut that keeps the 5–6 grey zone was consistently the
*worst*, empirically validating the decision to discard it. Per-operating-threshold
precision/recall/F1 are additionally in Supplementary Table S4.
(3) *Applicability-domain cut-off.* The Tanimoto AD threshold is empirically, not arbitrarily,
set: the n-weighted similarity-binned AUROC falls monotonically from 0.958 (nearest neighbour
Tanimoto ≥0.8) to 0.939, 0.866 and **0.770** (<0.4) (Supplementary Table S5), justifying the
out-of-domain flag in the 0.3–0.4 band.
(4) *Disease mapping* — the target→disease synthesis is a transparent, knowledge-based rule
(not a learned layer), each driver tagged by provenance, so it can be inspected and overridden.
(5) *Single safety anti-target* — cardiotoxicity is represented by hERG alone; other liabilities
(e.g. Nav1.5, hepatotoxicity) are out of scope and stated as such.
(6) *Read-across ceiling* — because the ensemble beats a kNN-Tanimoto baseline on every endpoint
(§3.6), performance is not merely memorised nearest-neighbour recall. None of these is concealed;
each is surfaced in the tool output or the supplementary tables.

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
(vi) **Bioactivity is pooled across assay formats** (IC50/Ki/Kd/EC50) on the pChEMBL scale;
residual cross-assay variance is bounded but not eliminated (see Threats to validity).
(vii) **A single cardiac safety anti-target (hERG)** is modelled; other safety liabilities are
out of scope.

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
Fig 1 pipeline/workflow · Fig 2 dataset size & class balance · Fig 3 AUROC across validation
regimes (random/scaffold/cluster/temporal) · Fig 4 (A) scaffold-CV ROC curves, (B) calibration
reliability · Fig 5 (A) conformal coverage, (B) ensemble vs kNN/logistic baselines ·
Fig 6 benchmark vs literature · Fig 7 predicted-vs-measured regression scatter (antioxidant + receptors).
**Supplementary tables** (`supplementary/`, exact values from the validation artifacts):
S1 classification metrics (random/scaffold/cluster/temporal, PR-AUC, BA, MCC, Brier, conformal) ·
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
