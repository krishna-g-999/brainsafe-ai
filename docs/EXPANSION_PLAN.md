# Data expansion and modelling plan

Prepared in response to the review committee. Each section states what will be done, the scientific
reason for it, and how it will be checked. This is a living document; results are recorded in
`decisions_log.md` as each stage completes. No step is executed before its design is agreed.

A single principle runs through the plan and was learned the hard way in an earlier version of this
project (documented in `docs/BS_MODEL_CARD.md`): **labels for a supervised model must come from
measurement, not from annotation that was itself derived from the same knowledge used to build the
features.** Adding compounds broadens chemical coverage; adding *labels* requires an experimental
measurement. The two are kept separate throughout.

---

## 1. Additional data (DrugBank, FDA, flavonoids, and measured bioactivity)

**What.** Assemble a larger, more chemically diverse compound collection and, separately, a larger set
of measured training labels.

- **Measured bioactivity labels** for every endpoint continue to come from **ChEMBL** (pChEMBL for
  target activity; B3DB for blood-brain-barrier permeability; DPPH assays for antioxidant capacity),
  because these are experimental measurements on a common, comparable scale. We will re-pull each
  target with the current ChEMBL release, record the exact query and access date, and expand where
  additional measured records exist.
- **DrugBank and the FDA approved-drug list** provide structures and regulatory status but their
  target/CNS annotations are largely qualitative. They are therefore used for three specific,
  defensible purposes, **not** as a source of training labels:
  (i) an **external evaluation set** of approved drugs, independent of the ChEMBL training splits;
  (ii) **chemical-space diversification** to widen the applicability domain and to make the
  similarity-binned generalisation analysis more informative;
  (iii) **coverage of natural products and flavonoids** (including dietary sources) as
  domain-relevant test compounds, with measured activity drawn from ChEMBL or primary literature
  where it exists.

**Why.** Mixing qualitative annotation with measured pChEMBL is what produced label-feature
circularity in the earlier prototype. Keeping measured labels (ChEMBL/B3DB/DPPH) distinct from
annotation-only sources (DrugBank/FDA) preserves the validity of the supervised task while still
letting the added compounds improve diversity and external testing.

**How checked.** Every source folder under `data/raw/` gets a `SOURCE.md` (database, version, URL,
access date, record count). A provenance table records, for each compound, its origin and whether it
contributes a *label* or only a *structure/evaluation* role.

## 2. Endpoint selection: why these, and what else is possible

**Current 12 endpoints** span the mechanistic axes relevant to neurodegeneration and CNS drug
viability: brain access (BBB); Alzheimer-related enzymes (AChE, BChE, BACE1, GSK-3beta);
Parkinson/mood monoamine oxidases (MAO-A, MAO-B); dopaminergic and serotonergic receptors (D2, A2A,
5-HT2A, SERT); and a cardiac safety anti-target (hERG).

**Selection criteria (to be stated explicitly in the manuscript):**
1. mechanistic relevance to neurodegenerative or neuronal disease, or to CNS drug viability;
2. availability of enough **measured** ChEMBL/B3DB records to train and validate a model at a
   pre-set data threshold (minimum n and minimum minority-class fraction);
3. a deployment quality gate (Matthews correlation coefficient >= 0.45 under scaffold cross-validation);
   endpoints that fail as classifiers are served as potency regression instead.

**Candidate additional endpoints** to evaluate against the same criteria (inclusion decided by data,
not assumption): NMDA receptor (excitotoxicity), nicotinic acetylcholine receptor, adenosine A1,
sigma-1 receptor, tau and alpha-synuclein aggregation, Nrf2/antioxidant-response signalling,
neuroinflammatory markers (COX-2, iNOS), P-glycoprotein efflux (complements the passive BBB model),
and additional safety anti-targets (Nav1.5, CYP inhibition). Each candidate will be data-scoped first;
only those clearing criterion 2 and the quality gate will be deployed, and the rest will be reported
as attempted-but-insufficient-data so the selection is transparent.

## 3. Feature retention: statistical justification

**What.** For each endpoint, justify which molecular descriptors are retained rather than asserting a
fixed feature set.

**Method.** The representation is a 1,024-bit ECFP-4 fingerprint (local substructure, the basis of
structure-activity relationships) concatenated with interpretable physicochemical descriptors
(molecular weight, cLogP, topological polar surface area, hydrogen-bond donors and acceptors,
rotatable bonds, aromatic ring count, fraction sp3, and related quantities that govern permeability
and drug-likeness). We will report, per endpoint:
- random-forest impurity importance and **permutation importance** (model-agnostic, less biased);
- a **feature-ablation** comparison (fingerprint only, descriptors only, combined), extending the
  analysis already run for assay-type and label-threshold robustness;
- removal of near-zero-variance and highly collinear descriptors, with the correlation threshold and
  the dropped features listed.

**Why.** This turns "we kept these features" into a measured statement about which features carry
signal for which endpoint, and documents any that are redundant.

## 4. Random forest with k-fold cross-validation, per module

Random forest is already one of three ensemble members. We will additionally report, per endpoint, a
**random-forest-only** model and both **5-fold and 10-fold** cross-validation, alongside the primary
scaffold-grouped split. An initial run of this comparison (`BS_cv_comparison.py`,
`results/tables/STable15_cv_comparison.csv`) shows the mean classification AUROC is 0.912 under
scaffold 5-fold, 0.958 under random 5-fold and 0.964 under random 10-fold, with random-forest-only at
0.960. The purpose of reporting all of these is to make explicit that the fold count changes little
whereas the split type changes a lot: the higher random-split numbers reflect analogue leakage
between random folds, not a better model, which is why the scaffold-grouped number is reported as the
honest estimate.

## 5. Numeric encoding of non-numeric fields

**What.** Any categorical or textual field used in modelling is converted to a numeric encoding that
does not collide with other fields.

**Method.**
- **Structure** is already numeric via the fingerprint and descriptors.
- **Ordered categories** (for example clinical phase 0-4) use integer ordinal encoding.
- **Unordered categories** (for example ATC class, data source, disease group) use one-hot encoding,
  with a documented column name per category so no two categories share an index.
- Encodings, their value maps and their column ranges are recorded in `docs/DATA_DICTIONARY.md`, and
  the encoder is fit on training folds only to avoid leakage.

**Why.** A single shared integer code for distinct categories would let the model read spurious order
into unordered classes; explicit one-hot or ordinal encoding, documented, prevents that.

## 6. Auditing every step

- Each stage is designed and reviewed **before** execution and recorded in `decisions_log.md` with its
  rationale and outcome.
- Every reported number is produced by a released script from released data; nothing is hand-entered.
- Data provenance is captured in `data/raw/**/SOURCE.md` and a per-compound provenance table.
- Code is written as tested, documented modules under `src/brainsafe/`, with cross-validation splits
  fixed by seed for reproducibility.

---

## Staging

1. **Structure and documentation (this stage).** Clean layout, provenance conventions, this plan.
2. **Data assembly.** Re-pull measured endpoints from the current ChEMBL release; add DrugBank/FDA/
   flavonoid structures to `data/raw/` with provenance; build `data/processed/` tables.
3. **Endpoint scoping.** Data-scope the candidate endpoints; decide inclusion by the stated criteria.
4. **Feature and encoding analysis.** Importance, permutation, ablation, collinearity; finalise the
   encoded feature tables and the data dictionary.
5. **Modelling.** Random forest and ensemble, 5- and 10-fold and scaffold splits, per module.
6. **Validation and write-up.** Significance, external DrugBank/FDA test, generalisation curve; update
   the manuscript.

Each stage produces clear CSV tables under `results/tables/` and a short entry in `decisions_log.md`.
