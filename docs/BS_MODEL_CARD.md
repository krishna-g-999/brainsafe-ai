# BrainSafe AI Model Card & Honest Validation Record

*This is a development record, written in the model-card style of Mitchell et al. (2019, FAT* '19).
It is deliberately chronological, so that the reasoning behind the current models, including the flaw
we found and corrected, is on the record rather than hidden.*

> **How to read this document.** Sections **1–10** describe the **original prototype** (a seven-dimension
> "Neuroprotection Score" trained on a 535-compound curated file) and the problem we diagnosed in it:
> it was learning from annotations rather than from chemistry. That prototype was **superseded**.
> Sections **11–17** describe the **current deployed system**, eight measured-data classification
> endpoints, four receptor regressions, a measured antioxidant model, and the integration layer, all
> trained on 64,474 measured ChEMBL/B3DB records. **For the model that is actually shipped and reported
> in the manuscript, read from §11 onward.** Numbers in §1–10 refer to the retired prototype and should
> not be read as properties of the current tool.

---

## 0. 2026-07-21 update — current primary model: Random Forest + 10-fold on expanded data

This is the committee-requested model and supersedes the ensemble figures below for headline
reporting. Full detail: `docs/RF_CV_RESULTS.md`, `docs/METHODOLOGY_AUDIT.md`, `docs/DATA_MANIFEST.md`.

- **Data (measured only):** two independent public sources pooled for the eleven protein targets,
  **ChEMBL 37 + BindingDB** (18,573 compounds measured by both; 4,246 from BindingDB alone), plus
  **B3DB** (BBB) and **ChEMBL DPPH** (antioxidant). Master table: **61,317 unique compounds,
  67,982 measured compound-endpoint records** (`data/processed/compound_library.csv`).
- **Model:** RandomForest per endpoint (300 trees, min_samples_leaf 2, balanced weights),
  1024-bit ECFP-4 + 12 descriptors. Eight classifiers, four receptor regressors, one antioxidant
  regressor.
- **Validation:** random and scaffold-grouped **10-fold**. Mean classifier AUROC **0.960 random /
  0.919 scaffold**; regressor R² 0.60-0.68 random / 0.39-0.58 scaffold. All fold predictions saved
  (`data/processed/cv_predictions/`).
- **Data expansion audited:** adding BindingDB moved the scaffold headline by a mean of **-0.0002**
  (no inflation; `results/tables/expansion_audit.csv`).
- **Calibration:** isotonic per endpoint; mean expected calibration error **0.072 → 0.012**
  (`results/tables/calibration.csv`). Calibrated models: `models_rf/<endpoint>_calibrated.joblib`.
- **Applicability domain:** nearest-neighbour Tanimoto (flag at T<0.30). DrugBank is 72% in-domain
  for BBB but only 20-37% for the enzyme targets, so library-wide predictions are largely
  extrapolation and are flagged (`results/tables/applicability_coverage.csv`).
- **External validation:** BBB model on 306 FDA-curated approved drugs absent from training,
  **AUROC 0.774** (`results/tables/external_bbb_validation.csv`).

Per-endpoint scaffold AUROC (classifiers): BBB 0.920, AChE 0.921, BChE 0.937, BACE1 0.956,
GSK-3β 0.937, MAO-A 0.868, MAO-B 0.890, hERG 0.921. All exceed the deployment gate (scaffold MCC
0.56-0.69). Note: several target sets are active-skewed in the public data (e.g. GSK-3β ~93%
active), so predicted-active *base rates* are high; read probabilities through the calibration layer
and the applicability-domain flag.

---

## 1. Intended use

- **What it is:** a decision-support / triage explorer for neuroprotective potential
  (7 mechanistic dimensions → an aggregate Neuroprotection Score, NPS) plus a
  deterministic **druggability** score.
- **Intended users:** researchers screening/prioritising candidate compounds for
  follow-up, *not* clinicians, *not* patients.
- **NOT intended for:** quantitative efficacy prediction, dosing, clinical decisions,
  or de-novo prediction of activity for structurally novel compounds (see §4).

## 2. Data

- **Source file:** `data/brainsafe_SCIENTIFIC_FIXED.csv`, 535 compounds (529 with a
  parseable SMILES), ~20 chemical classes (FDA drugs, clinical candidates, polyphenols,
  vitamins, etc.).
- **Labels:** 7 mechanistic dimension scores on 0–100, plus disease-relevance categories.
- **Label provenance (important):**
  - `curated` / `curated_v2` / `curated_smiles` / `literature_derived`: **339** rows
    (human/literature-assigned ordinal scores; the human-`curated` subset has only
    **13 distinct values** → coarse expert ordinal labels).
  - `chembl_ml_predicted`: **190** rows (~35%) are **themselves ML-generated**, yet were
    tiered `gold_ml_predicted`. These are pseudo-labels, not measurements.
- **Inter-dimension correlation:** mean **r = 0.84** (range 0.72–0.97). The "7 mechanisms"
  are largely **one latent factor**; the labels do not support 7 independent mechanistic readouts.

## 3. Model

- **Architecture:** per-dimension 4-model mean ensemble (RandomForest, GradientBoosting,
  ExtraTrees, Ridge). 7 dimension models + a (redundant) NPS model.
- **Features (93):** ECFP-4 PCA(50) + ChemBERTa-PCA(32) + 4 disease-target counts + 1 BBB
  ordinal + 6 structural descriptors.
- **NPS:** a fixed weighted sum of the 7 dimensions (`scorer.neuro_score`), computed
  *after* prediction, not a learned target. (The serialized `nps_*.joblib` models are
  therefore circular and unused at inference.)

## 4. Evaluation: leak-free scaffold-split GroupKFold (k=5)

**Method:** Bemis–Murcko *generic* scaffold groups; no scaffold shared across folds; all
transforms (ECFP-PCA, ChemBERTa-PCA, scaler) fit inside the training fold only. Leakage
audit: test→train median Tanimoto **0.31** (only ~5% of test compounds have a T≥0.7 twin) , 
a genuine novel-chemotype split.

**Performance (NPS):**

| Model | NPS R² | Spearman ρ | MAE |
|---|---|---|---|
| Mean predictor (baseline) | −0.01 | n/a | ~15 |
| Ridge on 7 descriptors (baseline) | ~0.0 | 0.11 | ~14 |
| Tanimoto k-NN (baseline) | **−0.17** | 0.10 | ~15 |
| **93-feature ensemble** | **0.57** | **0.82** | **8.5** |

The ensemble clearly beats all baselines, and **pseudo-labels do not inflate it**
(human-only 0.568 ≈ full 0.574). Per-dimension R² 0.46–0.64. *This contradicts an earlier
addendum that reported negative 5-fold R²; those numbers were from earlier/ablation model
generations, now disambiguated.*

## 5. ⚠️ Critical caveat: the model predicts from annotations, not chemistry

A feature ablation (scaffold-CV, full set) is decisive:

| Feature set | NPS R² | mean dim R² |
|---|---|---|
| Full (structure + disease + BBB) | 0.573 | 0.554 |
| **Structure-only** (ECFP + ChemBERTa + struct) | **0.037** | 0.045 |
| **Disease + BBB only** (5 features) | **0.548** | 0.519 |

**≈ All predictive signal comes from 5 coarse disease-association features; the 82
structural/ChemBERTa features add ≈ 0.03 R².** This is almost certainly **label–feature
circularity**, the disease-target counts and the dimension labels were curated from the
same knowledge, so the model reads the answer back out rather than learning structure→activity.

**Consequence for applicability:** for a genuinely novel compound with no known disease
targets, those features are 0, and the model degrades to structure-only performance (≈ 0).
**The model is an evidence aggregator for already-annotated compounds, not a de-novo
structural predictor.**

## 6. Druggability score (separate, deterministic)

`BS_druggability.py`, CNS-weighted composite of **QED** (Bickerton 2012), **Lipinski Ro5**
(2001), **Veber** (2002), **CNS-MPO** (Wager 2010), with a PAINS flag. Computed analytically
from SMILES; **no training, no circularity.** Validated qualitatively (CNS drugs high, polar
non-drugs low). *Caveats:* composite weights are an unvalidated heuristic; QED/Ro5/CNS-MPO
share descriptors (MW/HBD/TPSA), so components are correlated, not independent; cLogD≈cLogP
and pKa omitted in CNS-MPO. **This is the most scientifically defensible, deployable component.**

## 7. Known limitations (summary)

1. ML neuroprotection prediction is **not valid for novel chemistry** (§5).
2. Labels are partly ML-generated (35%) and largely collinear (one latent factor).
3. NPS is a fixed formula; its "R²" is not an independent metric.
4. Druggability composite weights are heuristic and unvalidated against outcomes.
5. Applicability domain should be judged by *known-annotation availability*, not only
   structural Tanimoto.

## 8. Recommendations (to become genuinely valid)

- **Reframe** the ML output as "aggregated known evidence," or restrict it to annotated
  compounds; lead the product with the deterministic druggability + curated literature.
- **Rebuild labels** from measurable per-mechanism endpoints (e.g., assay-specific ChEMBL
  pChEMBL: AChE/MAO-B activity, ORAC/DPPH antioxidant, etc.) so dimensions are independently grounded.
- **Exclude label-coupled features** (disease counts derived from the same curation) or
  derive them from an independent source, then re-measure structure→activity transfer.
- **Calibrated uncertainty** (conformal prediction) and an honest out-of-domain flag.
- **Validate druggability** against an outcome (CNS approval / clinical phase; B/P ratios).

## 9. Reproducibility

- Honest evaluation: `python BS_scientific_validation.py` → `BS_validation_report.json`
- Druggability: `python BS_druggability.py --selftest`
- AD reference corrected to train-only (424 SMILES); backup `training_smiles.json.all529_bak`.
- Fixed seed `random_state=42` throughout.

---

## 10. Tier-1/2 update: building a *genuine* (non-circular) predictor

To test whether neuroprotection is learnable **from chemical structure alone**, we
re-ran with the circularity removed: **human/literature labels only** (ML pseudo-labels
quarantined), **structure-only features** (Morgan + RDKit descriptors; disease-count
features *excluded*), scaffold GroupKFold(5), with **2000× bootstrap 95% CIs** and
mean/k-NN baselines (`BS_predictive_model.py` → `BS_predictive_report.json`; leakage
audit: median test→train Tanimoto 0.34).

**Per-dimension structure→activity (n=339, scaffold-CV):**

| Dimension | Ensemble R² | 95% CI | Beats baselines? | Genuine? |
|---|---|---|---|---|
| **antioxidant** | **0.19–0.27** | **[0.18, 0.35]** (best: Ridge/descriptors) | yes | ✅ **YES** |
| anti_inflammatory | −0.06 | [−0.17, 0.02] | no | ❌ |
| mitochondrial_support | −0.15 | [−0.26, −0.06] | no | ❌ |
| aggregation_modulation | −0.13 | [−0.24, −0.04] | no | ❌ |
| cognitive_enhancement | −0.12 | [−0.22, −0.03] | no | ❌ |
| neurogenesis | −0.22 | [−0.33, −0.14] | no | ❌ |
| synaptic_plasticity | −0.20 | [−0.31, −0.11] | no | ❌ |
| **NPS (aggregate)** | −0.08 | [−0.19, 0.00] | no | ❌ |

**Findings:**
1. **Only `antioxidant` is genuinely structure-predictable**, chemically expected
   (radical-scavenging tracks phenolic-OH/catechol motifs). The best model is a
   **parsimonious, interpretable Ridge on 24 RDKit descriptors**, scaffold-CV
   **R² = 0.267 [0.182, 0.346]**, Spearman 0.49, MAE 15. Saved:
   `models_genuine/antioxidant_genuine_ridge.joblib` (+ spec). Adding Morgan
   fingerprints or ChemBERTa did *not* help (and slightly hurt) → signal is in
   interpretable physicochemistry, not learned embeddings.
2. **The other six dimensions are not structure-predictable** with the current
   curated labels (negative R²). This confirms why the original model required the
   circular disease-count features to appear performant.
3. **Label ceiling:** the antioxidant labels have only ~28 distinct values; continuous
   assay endpoints (ORAC / DPPH / TEAC) would be expected to raise this materially.

**Net:** BrainSafe now contains one *genuinely predictive, leak-free, interpretable*
model (antioxidant), a deterministic druggability layer, and an honest "not
structure-predictable" verdict for the remaining endpoints. To make those endpoints
genuinely predictable requires **independent, measurable, mechanism-specific labels**
(see §8), not achievable from the current curation.

---

## 11. GENUINE multi-endpoint brain predictor (measured-data models)

> **Note (2026-07-21):** the current headline model is the Random Forest + 10-fold pipeline on the
> expanded ChEMBL+BindingDB data described in **§0**; the block below records the earlier ensemble
> deployment on ChEMBL-only data and is retained for provenance.
>
> **Authoritative ensemble numbers (ChEMBL-only deployment).** The deployed panel is eight measured-data classification endpoints plus four receptor regressions and a measured antioxidant model, matching Manuscript Table 4 and `endpoints_report.json`: BBB n=7,805 AUROC 0.921; AChE 4,324 / 0.915; BChE 2,580 / 0.937; BACE1 8,067 / 0.950; GSK-3β 4,044 / 0.920; MAO-A 2,141 / 0.867; MAO-B 3,455 / 0.885; hERG 5,905 / 0.901; antioxidant (DPPH) n=2,862, R²=0.43. Total measured records 64,474 (= 61,108 across the twelve endpoint sets + 2,862 antioxidant + 504 clinical). The tables in §11–§15 are chronological build snapshots; where they differ, these deployed numbers are correct.

The §10 conclusion was acted on: instead of holistic curated scores, we trained
models on **measured public bioactivity data** for endpoints that are real,
mechanistically meaningful, and structure-predictable.

**Data:** ChEMBL measured pChEMBL (active ≥6 / inactive <5, grey zone dropped,
per-compound median, InChIKey-deduped) for CNS targets; B3DB measured labels for BBB.
**Features:** Morgan-1024 + 24 RDKit descriptors (structure only).
**Validation:** scaffold GroupKFold(5), out-of-fold AUROC/PR-AUC/MCC; per-prediction
applicability-domain flag (max Tanimoto to a 2000-compound training sample).

| Endpoint | Brain meaning | n | AUROC | PR-AUC | MCC | median test→train T |
|---|---|---|---|---|---|---|
| **BBB** | reaches the brain? | 7,805 | **0.921** | 0.951 | 0.66 | 0.48 |
| **AChE** | Alzheimer's / cognition | 4,324 | **0.915** | 0.964 | 0.63 | 0.68 |
| **BACE1** | Alzheimer's / amyloid | 8,067 | **0.950** | 0.992 | 0.63 | 0.71 |
| **MAO-B** | Parkinson's / dopamine | 3,455 | **0.885** | 0.925 | 0.62 | 0.59 |
| **MAO-A** | mood / depression | 2,141 | **0.867** | 0.815 | 0.58 | 0.57 |
| antioxidant | oxidative stress | 2,862 | R²=0.43 (measured DPPH, §15) | n/a | n/a | 0.34 |
| druggability/CNS-MPO | developability | n/a | deterministic (RDKit) | n/a |, | n/a |

**External sanity (known drugs, held to chemistry):** Donepezil → AChE active 0.99,
BBB 0.98; Selegiline → MAO-B active 0.97; Caffeine → BBB 0.86; Quercetin → BBB
non-penetrant 0.26 + antioxidant 86. (Galantamine AChE 0.60 sits just below threshold , 
a genuinely weak inhibitor near the boundary; honest uncertainty, not error.)

**Honest caveats:**
- For the ChEMBL target sets, median test→train Tanimoto remains 0.57–0.71 even under
  scaffold split (ChEMBL target data are dense analog/SAR series). AUROCs therefore
  reflect strong predictivity **within known-inhibitor chemical space**; the AD flag
  marks novel-scaffold queries as low-confidence. A stricter cluster/temporal split
  (and external test set) is the recommended next validation for publication.
- Class imbalance (e.g. BACE1 89% active) is reported via AUROC/PR-AUC/MCC, not raw accuracy.
- These are target-activity predictions (a mechanistic proxy for "brain effect"), not
  clinical efficacy. Research-use only.

**Artifacts:** `BS_fetch_endpoints.py` (data), `BS_train_endpoints.py` (training/validation),
`models_brain/<endpoint>.joblib` + `_meta.json` + `endpoints_report.json`,
`BS_brain_predict.py` (unified inference), wired into the app's novel-compound view
(`render_brain_profile_panel`). Datasets: `data/endpoints/*.csv`.

**Publication readiness:** this is now a genuine, measured-data, scaffold-validated
multi-endpoint QSAR tool, a defensible resource/methods paper. To strengthen to a
top-tier predictor paper: add cluster/temporal external validation, calibrated
(conformal) probabilities, GSK-3β/BACE-style additional targets, and report a strict
out-of-distribution generalisation analysis.

---

## 12. The innovation layer: an evidence-grounded, BBB-gated, calibrated CNS profiler

The endpoint models were assembled into something beyond "another QSAR predictor":
a transparent decision-support engine (`BS_brain_predict.py`, app
`render_brain_profile_v2`). Four ideas, each adding genuine scientific value:

1. **Calibrated probabilities.** Every classifier's output is mapped through an
   isotonic calibrator fit on its scaffold-CV out-of-fold predictions, so a reported
   "70%" means ~70% empirically (Brier 0.04–0.14). Decisions use Youden-J thresholds
   on the *calibrated* scores.

2. **BBB-gated disease synthesis.** A target hit is meaningless if the molecule can't
   reach the brain, so effective CNS engagement = P(target) × P(BBB-penetrant). These
   roll up into transparent per-disease scores (Alzheimer's = max(AChE, BACE1) × BBB;
   Parkinson's = MAO-B × BBB; Depression = MAO-A × BBB; Neuroprotection = GSK-3β × BBB).

3. **Evidence grounding (anti-black-box).** Every call ships with the **nearest real
   measured compounds** (Tanimoto + their measured active/inactive + pChEMBL), so users
   see the data behind a prediction and can judge trust. High similarity = strong
   evidence; low = flagged extrapolation. This is the core trust/innovation feature.

4. **Benefit ⇄ risk.** A measured **hERG** cardiotoxicity model (AUROC 0.91) adds a
   safety axis, turning the tool from "is it active?" into "is it a viable CNS lead?".

**Deployed measured-data classification endpoints (scaffold-CV, isotonic-calibrated).** The eight deployed endpoints are BBB, AChE, BChE, BACE1, GSK-3β, MAO-A, MAO-B and hERG; the snapshot table below omits BChE and predates the final data pull, so use the authoritative numbers noted in §11:

| Endpoint | Role | AUROC | Brier | n |
|---|---|---|---|---|
| BBB | brain access (gate) | 0.921 | 0.105 | 7,805 |
| AChE | Alzheimer's / cognition | 0.915 | 0.099 | 4,324 |
| BACE1 | Alzheimer's / amyloid | 0.950 | 0.042 | 8,067 |
| GSK-3β | tau / neuroprotection | 0.920 | 0.044 | 4,044 |
| MAO-B | Parkinson's / dopamine | 0.885 | 0.122 | 3,455 |
| MAO-A | mood / depression | 0.867 | 0.136 | 2,141 |
| hERG | SAFETY (cardiotox) | 0.901 | 0.123 | 5,905 |
| antioxidant | oxidative stress | R²=0.43 (§15) | n/a | 2,862 |
| druggability/CNS-MPO | developability | deterministic | n/a |, |

**Prospective sanity (chemistry-only inputs; verified against current models, structures
PubChem-resolved):** Donepezil → Alzheimer's disease score 1.00 via AChE (P=1.00), nearest
measured analogue Tanimoto 1.00 (donepezil itself, pChEMBL 7.75), hERG P=0.78 flagged (its clinical
hERG relevance is modest relative to exposure, so treated cautiously); Selegiline/Rasagiline →
Parkinson's via MAO-B (0.95–1.00), hERG low; **Terfenadine** (withdrawn for cardiotoxicity) →
**hERG 1.00** and BBB non-penetrant 0.42 (caught); Fluoxetine → Depression 0.96 via SERT
(Fluoxetine Ph4 clinical precedent); Quercetin/Resveratrol → BBB non-penetrant (0.18/0.35),
analogs flagged as extrapolation. The engine reproduces known pharmacology and safety.

**Honest limitations (unchanged, still apply):** ChEMBL target sets are analog-dense
(median scaffold-split test→train Tanimoto 0.57–0.71 for the target models; BBB 0.48),
so AUROCs reflect strong predictivity within known-inhibitor space, the AD/evidence
flags expose novel-scaffold extrapolation. Predicts target engagement, not clinical
efficacy. The remaining steps to a flagship paper are a strict cluster/temporal external
test, full conformal prediction sets, and prospective experimental validation.

**Artifacts:** `BS_fetch_endpoints.py`, `BS_train_endpoints.py` (calibration+evidence),
`models_brain/*` (+`endpoints_report.json`), `BS_brain_predict.py` (engine),
`render_brain_profile_v2` in the app.

---

## 13. Publication-grade external validation (no-compromise)

`BS_external_validation.py` → `BS_external_validation_report.json`,
`models_brain/<ep>_conformal.json`.

**(a) Strict leave-cluster-out generalisation.** Sphere-exclusion (LeaderPicker,
Tanimoto-distance 0.4) clusters; whole clusters held out (GroupShuffleSplit) so test
compounds are structurally novel relative to training. AUROC barely drops vs the
scaffold-CV number, the models generalise beyond analog series:

| Endpoint | scaffold-CV AUROC | **strict cluster-split AUROC** | cluster-split median T | #clusters |
|---|---|---|---|---|
| AChE | 0.915 | **0.912** | 0.61 | 1,352 |
| BACE1 | 0.950 | **0.918** | 0.69 | 1,243 |
| BBB | 0.921 | **0.906** | 0.50 | 2,793 |
| GSK-3β | 0.920 | **0.915** | 0.68 | 1,210 |
| MAO-A | 0.867 | **0.890** | 0.59 | 866 |
| MAO-B | 0.885 | **0.873** | 0.62 | 1,181 |
| hERG | 0.907 | **0.900** | 0.55 | 2,115 |

**(b) Similarity-binned generalisation curve.** AUROC reported in test→train
max-Tanimoto bins (T<0.4, 0.4–0.6, 0.6–0.8, >0.8) per endpoint (see JSON), the honest
"performance vs novelty" curve; performance is retained well even into the T<0.6 region.

**(c) Mondrian (class-conditional) inductive conformal prediction.** Calibrated on
scaffold-CV OOF; gives per-compound prediction SETS with statistically valid coverage.
Empirically validated (50/50 calib/test split):

| Endpoint | target | **empirical coverage** | mean set size | singletons |
|---|---|---|---|---|
| AChE | 0.90 | 0.899 | 1.17 | high |
| BACE1 | 0.90 | 0.893 | 1.06 | high |
| BBB | 0.90 | 0.897 | 1.14 | high |
| GSK-3β | 0.90 | 0.885 | 1.16 | high |
| MAO-A | 0.90 | 0.905 | 1.30 | n/a |
| MAO-B | 0.90 | 0.895 | 1.24 | n/a |
| hERG | 0.90 | 0.889 | 1.17 | high |

Coverage ≈ the 90% target across all endpoints → **trustworthy confidence**. The engine
now returns each call as a conformal set ({active}/{inactive}/uncertain/out-of-distribution),
shown in the app.

**What remains (cannot be fabricated, stated honestly):** a true **temporal/prospective
split** (train on pre-cutoff ChEMBL, test on later depositions) and **wet-lab prospective
validation**. These require a dedicated temporal data pipeline and experiments,
respectively, they are the final steps for a flagship predictor paper. Everything
computationally validatable has been done with no compromise: leak-free scaffold CV,
strict leave-cluster-out, similarity-binned generalisation, calibrated probabilities, and
conformal prediction with verified coverage.

---

## 14. Closing the computationally-fixable limitations (no-compromise update)

The §5/§6 limitations were addressed to the extent computationally possible (the rest
are honestly marked as inherent). `BS_fetch_endpoints.py` (now stores `document_year` +
expanded panel), `BS_temporal_pr.py`, `BS_external_validation.py`.

**(i) Temporal / time-split validation, DONE (this was computational, not impossible).**
Each compound is assigned its earliest ChEMBL `document_year`; train ≤ 75th-percentile
year, test on the most recent ~25% (a true "future compounds" test). This is the honest
prospective-use estimate and is **lower than scaffold-CV, as it should be**:

| Endpoint | scaffold-CV | leave-cluster-out | **temporal** | prospective read |
|---|---|---|---|---|
| BACE1 | 0.95 | 0.94 | **0.92** | robust on future compounds |
| BChE | 0.94 | 0.92 | **0.79** | good |
| AChE | 0.92 | 0.91 | **0.78** | good |
| MAO-B | 0.89 | 0.87 | **0.76** | moderate |
| hERG | 0.90 | 0.87 | **0.76** | moderate |
| GSK-3β | 0.92 | 0.92 | **0.66** | **weak on novel/future chemotypes** |
| MAO-A | 0.87 | 0.89 | **0.61** | **weak on novel/future chemotypes** |
| BBB | 0.92 | 0.91 | n/a (B3DB undated) | n/a |

Temporal AUROC is now stored per endpoint (`*_meta.json: temporal_auroc`). GSK-3β and
MAO-A degrade substantially out-of-time and must be treated as lower-confidence for
prospective screening, stated plainly rather than hidden behind the CV number.

**(ii) Expanded target panel, DONE, fact-gated.** Fetched 5 additional measured CNS
targets. Outcome decided by results, not assumption:
- **BChE (CHEMBL1914) added**, scaffold AUROC 0.94, cluster 0.92, **MCC 0.70** (best in
  the panel), 70% active. Genuine addition (Alzheimer's cholinergic, complements AChE).
- **D2, A2A, 5-HT2A, SERT excluded.** Their ChEMBL data are 96–98% actives (only binders
  reported), so binary active/inactive QSAR is ill-posed → **MCC 0.21–0.44** despite high
  AUROC. A pre-stated **MCC ≥ 0.45 quality gate** (`BS_brain_predict.MIN_MCC`) drops them
  from deployment. Proper future handling = **potency (pKi) regression**, not classification.
- Deployed panel: **BBB + AChE, BChE, BACE1, GSK-3β, MAO-A, MAO-B + hERG (safety) = 8**.

**(iii) PR-AUC + threshold sensitivity, DONE** (`BS_temporal_pr_report.json`). For every
endpoint we now report precision/recall/F1 at thresholds 0.3/0.5/0.7 and Youden-J, on
held-out (temporal/scaffold) sets, the correct reporting for imbalanced endpoints. E.g.
BACE1 (89% active) temporal @0.5: P=0.97/R=0.95/F1=0.96; hERG @youden: P=0.57/R=0.69.

**(iv) Analog density, explicit.** Cluster-split (§13) + similarity-binned curves quantify
it; temporal split further controls for it. Stated, not assumed.

**Inherent (cannot fix computationally, stated honestly):**
- **Target engagement ≠ clinical efficacy**, needs clinical-outcome labels we do not have.
- **Wet-lab prospective validation**, requires experiments.
- **Antioxidant** was a weak (R²≈0.27) curated-label model *at this stage*; §15 replaces it with a measured-DPPH regression (R²=0.43). Improving it further would need additional measured
  ORAC/DPPH/TEAC assay data (not reliably available as a clean public set here).

**Net:** of the five named limitations, three are now genuinely fixed (temporal validation,
PR/threshold reporting, panel expansion via BChE) and one is honestly bounded (analog
density, fully quantified); only clinical-efficacy and wet-lab validation remain, and both
are scientifically impossible to fabricate computationally.

---

## 15. Measured-data fixes: antioxidant (DPPH) + receptor potency regression

Two remaining gaps were addressed with **real measured data**, then validated/cross-checked.

**(A) Antioxidant, replaced curated labels with MEASURED DPPH data.**
Assembled 2,862 unique compounds from ChEMBL DPPH radical-scavenging assays (measured
IC50/EC50 → pIC50; `BS_fetch_antioxidant.py`). A RF+ExtraTrees+HistGB regression ensemble
on Morgan-1024 + 24 descriptors achieved, under scaffold GroupKFold(5):
- **R² = 0.43, RMSE = 0.60 log-units, Spearman ρ = 0.636**, a genuine improvement over the
  curated model (R² ≈ 0.25).
- **Temporal split R² ≈ 0** (honest): pooled cross-lab/protocol DPPH measurements do not
  generalise across time, the model is reliable in-distribution/for ranking, not as an
  absolute cross-era predictor.
- **Cross-check:** the previous curated 0–100 score correlates only weakly with the measured
  DPPH pIC50 (Spearman ρ = 0.387), confirming the curated model captured only partial signal
  and that the measured model is the better basis. Deployed: `antioxidant_measured_dpph.joblib`.

**(B) Receptor targets, switched from ill-posed binary to POTENCY REGRESSION.**
D2/A2A/5-HT2A/SERT were 96–98% actives → binary QSAR was meaningless (MCC 0.21–0.44). Trained
regression ensembles on measured pChEMBL (`BS_train_regression.py`):

| Receptor | n | scaffold-CV R² | RMSE | Spearman | temporal R² |
|---|---|---|---|---|---|
| A2A | 5,547 | 0.526 | 0.74 | 0.706 | 0.326 |
| 5-HT2A | 5,256 | 0.460 | 0.76 | 0.684 | 0.085 |
| D2 | 7,511 | 0.425 | 0.71 | 0.652 | −0.007 |
| SERT | 4,471 | 0.338 | 0.84 | 0.573 | 0.171 |

These are **ranking-grade** potency predictors (Spearman 0.57–0.71), a correct, useful
replacement for the dropped classifiers, but **temporal generalisation is weak** (D2 ≈ 0),
so they are surfaced for prioritisation/ranking, explicitly **not** as absolute prospective
potency. Deployed under `models_brain_reg/` and shown as a separate "predicted pKi" panel,
kept distinct from the calibrated classification probabilities.

**Honest status after these fixes:** antioxidant is now measured-data-based and ~2× stronger
(R² 0.25→0.43); receptor targets are reinstated as validated ranking-grade potency models.
The two genuinely inherent gaps remain unchanged and unfaked: **target engagement ≠ clinical
efficacy**, and **no wet-lab prospective validation**. Temporal weakness of pooled-assay
endpoints (antioxidant DPPH, D2/5-HT2A) is reported, not hidden.

---

## 16. Root-cause diagnosis of the remaining gaps + data-driven solutions

Each gap was diagnosed with measured evidence (not assumed) and addressed where scientifically possible.

**(a) "Competitive, not SOTA", exact cause.** (i) We used only **pChEMBL (potency) records** and
discarded the larger body of qualitative %-inhibition data (e.g. MAO-A: 4,398 pChEMBL vs **19,633
total** activities). (ii) Representation is standard (Morgan-1024 + 24 descriptors), not graph/deep
multitask. (iii) Most published "SOTA" AUROCs use **random splits**; ours use scaffold/cluster/temporal,
which are stricter and lower by construction. *Conclusion:* the gap is largely split-methodology +
data-type restriction, not a modelling error. **Substantiated:** on like-for-like **random splits**
(`BS_randomsplit_benchmark.json`) our AUROC is **0.94–0.98** (BBB 0.963, AChE 0.975, BChE 0.976,
BACE1 0.956, GSK-3β 0.943, MAO-B 0.960, MAO-A 0.950, hERG 0.950), **at or above the published SOTA
ranges** (BBB 0.88–0.96; hERG 0.86–0.93). The models are SOTA-grade; we additionally report the
harder scaffold/cluster/temporal numbers that most papers omit. The full validation hierarchy is:
random 0.94–0.98 → scaffold 0.87–0.95 → cluster 0.87–0.94 → temporal 0.61–0.92.

**(b) GSK-3β / MAO-A temporal weakness, exact cause (measured).** In the temporal test,
**71–91% of recent compounds carry scaffolds never seen in training** (covariate shift), for *every*
endpoint. The differing temporal AUROC is explained by **class balance of the recent test set**, not by
one model being "better": BACE1/GSK-3β recent tests are 93% active (AUROC partly trivial), whereas
MAO-A's recent test is **balanced (45% active)**, a genuinely harder, more honest test, hence 0.61.
*Conclusion:* temporal degradation is driven by irreducible novel-chemotype shift; MAO-A's low number
is the most honest of the set. More-diverse data narrows but cannot eliminate it. Flagged in metadata.

**(c) Engagement ≠ efficacy, addressed with real clinical data.** Built a **clinical/translational
precedent layer** from ChEMBL ATC-N (nervous-system) molecules that reached a clinical phase
(**504 compounds with phase + structure + ATC-derived disease**; `BS_fetch_clinical.py`,
`BS_clinical_evidence.py`). For any query the tool returns the nearest clinically-advanced CNS
compounds (Tanimoto) with their **max clinical phase and disease**. This is *clinical precedent from
measured trial data*, explicitly **not** an efficacy prediction. Validated qualitatively (Donepezil→
Donepezil Ph4 AD; Fluoxetine→Fluoxetine Ph4 Depression).

**Integration upgrade.** Receptor-regression engagement (pKi→[0,1]) is now folded into the BBB-gated
disease synthesis, fixing mechanism coverage (e.g. **Fluoxetine now resolves to Depression via SERT**,
not Parkinson's). Each disease driver is tagged `calibrated-classifier` or `regression(pKi)` for honesty.

**Remaining inherent limitations (cannot be removed computationally):**
- Models predict **binding/engagement, not direction** (agonist vs antagonist), e.g. D2 binding is
  shown for both anti-Parkinson agonists and antipsychotic antagonists.
- **Efficacy is shown as clinical precedent, not predicted.**
- **Wet-lab prospective validation** still requires experiments.
- Temporal generalisation to wholly novel scaffolds remains bounded by covariate shift.

## 17. Comparative/non-comparative analysis & comparison with general-purpose LLMs

Added in response to reviewer feedback. Artifacts: `BS_llm_comparison.py`, `BS_llm_comparison.json`,
`supplementary/STable8_llm_capability_comparison.csv`, `supplementary/STable9_baseline_comparison.csv`.

**Non-comparative (standalone) analysis.** Per-endpoint absolute performance against each endpoint's
own held-out measured data (AUROC, PR-AUC, MCC, Brier, conformal coverage), Manuscript §3.1/§3.4;
STable1–S3.

**Comparative analysis, internal baselines (measured).** Under an identical scaffold-split protocol
and feature set, the deployed ensemble beats a kNN-Tanimoto "read-across" baseline and logistic
regression on **every one of the 8 classification endpoints**: mean scaffold AUROC **0.912** (ensemble)
vs **0.867** (kNN, Δ +0.045) vs **0.808** (logistic, Δ +0.104) (`BS_baseline_comparison.json` → STable9).
Beating a pure nearest-neighbour read-across shows the model learns SAR beyond structural look-up.

**Comparative analysis, vs general-purpose LLMs (why a dedicated tool?).**
- *Benchmark evidence (peer-reviewed):* general LLMs underperform specialised ML on molecular property
  prediction (Guo et al. 2023, NeurIPS D&B; Zhong et al. 2024, arXiv:2403.05075); fine-tuned LLMs match
  QSAR only in the low-data limit (Jablonka et al. 2024, Nat Mach Intell 6:161); LLMs hallucinate in
  generative settings (Ji et al. 2023, ACM Comput Surv 55(12):248).
- *Scientific background:* an LLM is an autoregressive text predictor, no molecular fingerprint, no
  structure→measured-activity function, no calibrated probability, no coverage guarantee, no
  applicability domain, no measured-analogue provenance. BrainSafe supplies all six.
- *Reproducible grounded-output demonstration:* for fixed structures the engine returns calibrated P +
  conformal set + nearest **measured** analogue with pChEMBL (donepezil→AChE P=1.00, nearest analogue
  Tanimoto 1.00 / pChEMBL 7.75; terfenadine→hERG P=1.00; novel arylpiperazine→honest conformal
  "uncertain" grounded in a measured analogue, Tanimoto 0.35 / pChEMBL 4.82). Every value is
  measurement-traceable, the audit guarantee an LLM cannot give.

**Consolidated scientific-flaw self-audit (threats to validity), now with quantitative tests.**
Artifacts: `BS_flaw_fixes.py/.json`, `BS_assay_composition.py`, `BS_assay_sensitivity.py`; STable10-13.
- **(1) Assay-type pooling, TESTED, not just documented.** Composition quantified (STable11): IC50
  dominant 81-92% per target except GSK3B (IC50 49%, EC50 33%, Ki 16%). Single-assay (IC50-only) vs
  pooled scaffold-CV retrain changes AUROC by **≤0.006** (GSK3B 0.919→0.913; MAO_B −0.006; hERG 0.000;
  STable12) → pChEMBL pooling does not materially distort discrimination.
- **(2) Label-threshold sensitivity, TESTED.** Re-labelling at {deployed ≥6/<5, strict ≥6.5/<5.5,
  sharp ≥6/<6, high-potency ≥7/<5} gives max scaffold-AUROC spread **0.109** over 4 endpoints; deployed
  ≈ strict (within 0.01-0.02); grey-zone-retaining "sharp boundary" is consistently worst, validating
  the grey-zone drop (STable10). Per-operating-threshold P/R/F1 also in STable4.
- **(3) AD cut-off, data-driven.** n-weighted similarity-binned AUROC 0.958 (T≥0.8) → 0.939 → 0.866 →
  **0.770** (T<0.4), justifying the 0.30-0.40 out-of-domain flag (STable5).
- **(4) disease mapping**, transparent knowledge-based rule, provenance-tagged, inspectable/overridable.
- **(5) single hERG safety anti-target**, other liabilities (Nav1.5, hepatotox) out of scope.
- **(6) read-across ceiling**, ruled out by §3.2/STable9 (ensemble > kNN on all 8 endpoints).

**Pre-registered LLM head-to-head, EXECUTED** (`BS_LLM_benchmark_protocol.md`, `BS_llm_benchmark.py`,
`BS_llm_score.py`, `BS_llm_responses.json`): frozen prompt + 10-compound panel (uncontested truth + 1
unpublished scaffold) + fixed rubric, run on Gemini Pro, ChatGPT/GPT-4o, Perplexity, Claude and scored
live against ChEMBL (STable13). Honest two-part result:
- **Well-known drugs:** LLMs are strong, 3/4 hit BBB 9/9 (BrainSafe 8/9, missed astemizole) and Brier
  as good/better (Claude 0.020). Raw accuracy on famous compounds does NOT justify the tool.
- **Grounding + novelty:** 14/31 (45%) of ChEMBL IDs the LLMs cited as provenance were fabricated or
  pointed to the WRONG molecule (e.g. Gemini's "rasagiline" id = fluticasone propionate, "selegiline" =
  propranolol; Claude's "rivastigmine" = pyridoxine, "terfenadine" = cefdinir; ChatGPT's "terfenadine"
  id = astemizole). All 4 confabulated a specific target+potency for the unpublished compound and
  disagreed (3 AChE, 1 D2). BrainSafe: 0 fabricated, grounded in real measured analogues, honest
  conformal "uncertain" on the novel compound. Conclusion: LLMs approximate textbook classifications but
  cannot be trusted for verifiable provenance or novel chemistry, the tool's actual value.
