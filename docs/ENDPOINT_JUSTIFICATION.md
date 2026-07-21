# The twelve prediction endpoints: what they are, why these, and where the values come from

This document answers three questions raised in review: (1) what each of the twelve modelled
endpoints is, (2) why these twelve were chosen and whether more are possible, and (3) exactly how
the training and test values for each endpoint are obtained. Every count below is the number of
measured compounds actually held in `data/endpoints/` and `data/endpoints_reg/` on 2026-07-20.

## 1. The twelve endpoints

The tool models twelve structure-based endpoints, grouped by the biological question each answers.
Eight are binary classifiers (active / inactive) and four are potency regressors; a measured
antioxidant regressor and a rule-based druggability layer sit alongside them but are counted
separately (see section 4).

### Exposure gate (1)
| # | Endpoint | Task | Compounds | Why it is in the panel |
|---|----------|------|-----------|------------------------|
| 1 | BBB penetration | classification | 7,807 | A CNS drug must first cross the blood-brain barrier; this endpoint gates every downstream brain score. |

### Alzheimer's-disease target engagement (4)
| # | Endpoint | Task | Compounds | Rationale |
|---|----------|------|-----------|-----------|
| 2 | AChE (acetylcholinesterase) | classification | 4,387 | The validated symptomatic AD target (donepezil, rivastigmine, galantamine). |
| 3 | BChE (butyrylcholinesterase) | classification | 2,621 | Rises as AD progresses; a co-target of dual cholinesterase inhibitors. |
| 4 | BACE1 (beta-secretase 1) | classification | 8,501 | The rate-limiting enzyme of amyloid-beta production; the main disease-modifying AD target. |
| 5 | GSK-3-beta | classification | 4,958 | Drives tau hyperphosphorylation; the principal tau-side AD target. |

### Parkinson's / mood target engagement (2 classification + 4 regression)
| # | Endpoint | Task | Compounds | Rationale |
|---|----------|------|-----------|-----------|
| 6 | MAO-B (monoamine oxidase B) | classification | 3,665 | Established PD target (selegiline, rasagiline, safinamide). |
| 7 | MAO-A (monoamine oxidase A) | classification | 2,228 | Antidepressant target; the A/B selectivity split is clinically important. |
| 8 | D2 dopamine receptor | regression | 7,734 | Central to PD therapy and antipsychotic action; potency, not a yes/no call, is what matters. |
| 9 | A2A adenosine receptor | regression | 6,785 | Non-dopaminergic PD target (istradefylline). |
| 10 | 5-HT2A serotonin receptor | regression | 5,989 | Psychiatric and PD-psychosis target (pimavanserin). |
| 11 | SERT (serotonin transporter) | regression | 4,572 | The primary antidepressant (SSRI) target. |

### Safety liability (1)
| # | Endpoint | Task | Compounds | Rationale |
|---|----------|------|-----------|-----------|
| 12 | hERG channel | classification | 5,875 | The dominant cardiac-safety flag in CNS drug discovery; a standard early counter-screen. |

Together these twelve cover the full decision a medicinal chemist makes about a candidate CNS
molecule: can it get in (BBB), does it hit the disease target (AChE/BChE/BACE1/GSK-3-beta for AD;
MAO-A/MAO-B/D2/A2A/5-HT2A/SERT for PD and mood), and will it be safe (hERG).

## 2. Why classification for eight and regression for four

D2, A2A, 5-HT2A and SERT are modelled as potency regression rather than binary classification
because their ChEMBL activity sets are 96-98% "active" (only binders tend to be reported for these
receptors). A binary split on such data is ill-posed: it collapses to a near-constant label and
fails the deployment quality gate (Matthews correlation 0.21-0.44). Regression on the continuous
pChEMBL value is the statistically appropriate task and is what a chemist actually needs (how
potent, not merely whether it binds). This is recorded in `docs/decisions_log.md`.

## 3. Are there more? Candidate endpoints considered and deferred

The framework is deliberately extensible; twelve is the current scope, not a ceiling. Endpoints
that were considered and the reason each is deferred:

| Candidate | Category | Why deferred (for now) |
|-----------|----------|------------------------|
| NMDA / NR2B, GABA-A, alpha-7 nAChR, sigma-1 | CNS target | Measured potency data are sparser and assay formats more heterogeneous; needs careful assay harmonisation before a defensible model. |
| COMT, alpha-synuclein aggregation | PD target | COMT sets are small; anti-aggregation assays are not standardised into a single comparable readout. |
| CYP2D6 / CYP3A4 inhibition, P-glycoprotein efflux | ADME/safety | Important but outside the current "brain target + core safety" scope; planned as a second safety module. |
| Aqueous solubility, microsomal clearance, plasma-protein binding | ADME | Physicochemical/PK properties rather than target engagement; a natural next module once the target panel is locked. |

The selection rule applied throughout: an endpoint is included only when (a) it answers one of the
three questions above for AD or PD, and (b) enough measured, assay-comparable public data exist to
train and honestly validate a model. Endpoints failing (b) are listed here rather than modelled on
thin data.

## 4. Additional layers (not counted among the twelve)

- **Antioxidant (DPPH)** — regression, 2,862 measured compounds; a neuroprotection-relevant
  property with a genuine measured assay, kept separate because it is a chemical property rather
  than a target-engagement call.
- **Druggability / CNS-MPO** — a deterministic, rule-based physicochemical layer (no training); it
  is a transparent scoring formula, not a machine-learning endpoint.

## 5. How the train/test values are obtained (per endpoint)

No value is invented or hand-annotated. Every label is a measured experimental value from a public
source, standardised identically before use.

- **Target-engagement classifiers (AChE, BChE, BACE1, GSK-3-beta, MAO-A, MAO-B).**
  Sources: two independent public databases, ChEMBL version 37 (release 2026-05-01) and BindingDB,
  pooled at the compound level. Kept: activities with a defined potency (pChEMBL, or a BindingDB
  IC50/Ki/Kd/EC50 converted to the same -log10 molar scale) against the named target. Per compound
  the median across both sources is taken. Label: active if pChEMBL >= 6 (<= 1 uM), inactive if
  pChEMBL < 5 (> 10 uM); the 5-6 grey zone is dropped so the two classes are unambiguous. hERG uses
  ChEMBL only (BindingDB was rate-limited at retrieval).
- **BBB penetration.** Source: B3DB, a curated database of measured blood-brain-barrier
  permeability (Meng et al., Sci Data 2021). Label: the experimentally determined BBB+ / BBB- class.
  A further 306 FDA-curated approved-drug compounds are held out as an external validation set.
- **Potency regressors (D2, A2A, 5-HT2A, SERT).** Sources: ChEMBL 37 and BindingDB, same filtering;
  the target value is the continuous per-compound median potency across both sources.
- **Antioxidant.** Source: ChEMBL DPPH radical-scavenging assays; IC50/EC50 converted to pIC50
  (-log10 molar) as the regression target.

**Train/test discipline.** For every endpoint the same compound never appears in both train and
test. Performance is reported under four increasingly strict regimes — random split, scaffold-
grouped cross-validation (Bemis-Murcko), leave-cluster-out, and temporal (train on older ChEMBL
records, test on newer) — so the number quoted is honest about how the model behaves on genuinely
novel chemistry, not just on close analogues of the training set.

## 6. Standardisation (applied identically to all endpoints)

Every SMILES is passed through the same pipeline before use: keep the largest organic fragment
(salt/solvent stripping), sanitise, generate the canonical SMILES and InChIKey, then deduplicate on
InChIKey. The original source record and the standardised training variant are both retained so the
transformation is fully traceable (`data/README.md`).

---
*Counts are the pooled ChEMBL + BindingDB compound totals per endpoint, verified against
`data/endpoints/*.csv`, `results/tables/endpoint_rebuild_provenance.csv` and the master
`data/processed/compound_library.csv` (61,317 unique compounds, 67,982 measured records) on
2026-07-21. Sources and licences: `data/raw/measured_endpoints_SOURCE.md`.*
