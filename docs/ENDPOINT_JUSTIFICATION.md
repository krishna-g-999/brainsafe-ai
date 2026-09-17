# Endpoint justification: what is modelled, why, and where the training values come from

Every count in this document is read live from `results/tables/MODEL_INVENTORY.csv`,
`models_rf/binder_modes.json` and `results/tables/np_target_survey.csv` by
`tools/build_endpoint_justification.py`, and reflects the panel as deployed on 2026-09-17. It replaces
an earlier version of this document describing a twelve-endpoint panel that predates the binder
architecture and the target-panel expansion; nothing in this file is carried forward from that
version without being re-verified against the current artefacts.

## 1. Overview

The deployed server answers four questions about a candidate molecule, in the order a medicinal
chemist would ask them: does it reach the brain (*exposure*), does it engage a mechanism relevant to
a CNS condition (*target engagement*), is it safe (*safety*), and is it developable (*ADME*). This is
covered by **75 trained estimators**, of which **70 are deployed**
and 5 were trained, tested and withdrawn after failing a specificity or
discrimination check described in section 3. The deployed panel spans **54 molecular
targets**, a 10-endpoint exposure and ADME layer, and one cardiac-safety classifier.

Target engagement is covered two ways, and both are counted above rather than one being treated as
supplementary. **12 core endpoints** (A2A, AChE, BACE1, BChE, D2, GSK3B, HT2A, MAO_A, MAO_B, SERT, antioxidant_DPPH, pka_basic) are trained on
measured potency or activity, most reported as a decade-old, well-characterised assay. **The
52-endpoint binder panel** (47 deployed) extends this to receptors and
enzymes across the monoaminergic, opioid, cannabinoid, histaminergic, adenosine, purinergic,
glutamatergic and nicotinic systems, each validated against compounds experimentally tested at that
target and found inactive, never against property-matched decoys. The rationale for each therapeutic
axis is set out in the manuscript's Methods section ("Endpoint selection") and is not repeated in
full here to avoid the two documents drifting apart; this document's job is the table a reviewer
asked for: every endpoint, its training size, and its deployment status, in one place.

## 2. The exposure and ADME layer (10 endpoints)

| Endpoint | Task | Training compounds | Deployed |
|---|---|---|---|
| BBB | classification | 3,901 | yes |
| caco2_permeability | regression | 897 | yes |
| clearance_hepatocyte | regression | 1,020 | yes |
| kpuu | regression | 566 | yes |
| lipophilicity | regression | 4,200 | yes |
| logbb | regression | 1,058 | yes |
| pgp_inhibition | classification | 1,212 | yes |
| pgp_substrate | classification | 1,371 | yes |
| plasma_protein_binding | regression | 1,797 | yes |
| solubility | regression | 9,573 | yes |

Blood-brain barrier penetration is modelled first among these because it gates every downstream
target score: a target score is admitted only in proportion to predicted exposure, so potency at a
target the compound cannot reach contributes nothing to the server's output. The remaining nine
endpoints (Caco-2 permeability, hepatocyte clearance, unbound brain-to-plasma ratio, lipophilicity,
logBB, P-glycoprotein inhibition and substrate status, plasma-protein binding, and aqueous
solubility) cover developability, the fourth of the four questions in section 1.

## 3. Core target-potency and activity endpoints (12)

| Endpoint | Task | Training compounds | Deployed |
|---|---|---|---|
| A2A | regression | 6,743 | yes |
| AChE | classification | 5,125 | yes |
| BACE1 | classification | 8,207 | yes |
| BChE | classification | 3,278 | yes |
| D2 | regression | 7,905 | yes |
| GSK3B | classification | 5,439 | yes |
| HT2A | regression | 6,075 | yes |
| MAO_A | classification | 3,585 | yes |
| MAO_B | classification | 4,534 | yes |
| SERT | regression | 4,479 | yes |
| antioxidant_DPPH | regression | 2,782 | yes |
| pka_basic | regression | 6,384 | yes |

`antioxidant_DPPH` and `pka_basic` are carried in this family for training-pipeline reasons (a
measured chemical property and a physicochemical property, respectively) rather than as
target-engagement calls, and are described as auxiliary endpoints in the manuscript's architecture
figure.

## 4. Safety

| Endpoint | Task | Training compounds | Deployed |
|---|---|---|---|
| hERG | classification | 9,933 | yes |

hERG blockade is the dominant cardiac-safety flag in CNS drug discovery and a standard early
counter-screen; it is reported here rather than folded into the target-engagement count because a
liability, not a mechanism, is what it answers.

## 5. The binder panel (52 endpoints, 47 deployed)

Each binder classifier answers "does this compound bind this target," validated on compounds
experimentally tested at that target and found inactive rather than on property-matched decoys, at a
threshold constrained simultaneously by held-out measured inactives and by the false-positive rate on
a disjoint pool of unrelated chemistry (`results/tables/background_specificity_disjoint.csv`). AUROC
and sensitivity below are both measured against the target's own held-out measured inactives, by
scaffold.

| Target | Training compounds | AUROC vs. measured inactives | Sensitivity | Status |
|---|---|---|---|---|
| A1 | 7,352 | 0.914 | 0.736 | deployed |
| A2A | 15,636 | 0.949 | 0.824 | deployed |
| CB1 | 10,198 | 0.949 | 0.862 | deployed |
| CGRP | 2,251 | 0.985 | 0.993 | deployed |
| COX2 | 4,429 | 0.783 | 0.338 | deployed, low-power flag |
| CSF1R | 6,772 | 0.950 | 0.901 | deployed |
| Cav3_2 | 1,325 | 0.982 | 0.975 | deployed |
| D2 | 13,861 | 0.872 | 0.654 | deployed |
| D3 | 15,831 | 0.956 | 0.860 | deployed |
| DAT | 5,113 | 0.959 | 0.864 | deployed |
| DHODH | 3,619 | 0.966 | 0.948 | deployed |
| GABA_A | 771 | 0.719 | 0.303 | deployed, low-power flag |
| GBA1 | 441 | 0.932 | 0.795 | deployed |
| GluA2 | 68 | 0.696 | 0.103 | withdrawn |
| GluN2B | 3,461 | 0.761 | 0.374 | deployed, low-power flag |
| H3 | 12,295 | 0.978 | 0.953 | deployed |
| HDAC1 | 11,728 | 0.960 | 0.900 | deployed |
| HDAC6 | 14,155 | 0.975 | 0.930 | deployed |
| HT1A | 14,179 | 0.936 | 0.754 | deployed |
| HT2A | 14,369 | 0.925 | 0.670 | deployed |
| HT6 | 10,532 | 0.960 | 0.931 | deployed |
| HT7 | 5,623 | 0.947 | 0.830 | deployed |
| KEAP1 | 463 | 0.880 | 0.520 | deployed |
| LRRK2 | 4,470 | 0.968 | 0.919 | deployed |
| MT1 | 2,670 | 0.896 | 0.678 | deployed |
| NET | 5,933 | 0.920 | 0.732 | deployed |
| NFKB1 | 102 | 0.459 | 0.000 | withdrawn |
| NLRP3 | 843 | 0.930 | 0.867 | deployed |
| NR3C1 | 37 | 0.410 | 0.167 | withdrawn |
| NRF2 | 285 | 0.789 | 0.545 | withdrawn |
| Nav1_1 | 227 | 0.952 | 0.120 | withdrawn |
| Nav1_5 | 870 | 0.921 | 0.800 | deployed |
| Nav1_6 | 1,027 | 0.862 | 0.532 | deployed |
| Nav1_7 | 10,397 | 0.957 | 0.877 | deployed |
| Nav1_8 | 1,030 | 0.956 | 0.892 | deployed |
| OPRK1 | 11,642 | 0.945 | 0.801 | deployed |
| OPRM1 | 14,058 | 0.954 | 0.872 | deployed |
| OX1 | 12,724 | 0.964 | 0.922 | deployed |
| OX2 | 14,679 | 0.964 | 0.885 | deployed |
| P2X7 | 11,122 | 0.813 | 0.465 | deployed, low-power flag |
| PDE10A | 15,463 | 0.962 | 0.875 | deployed |
| PDE4B | 4,332 | 0.966 | 0.867 | deployed |
| RIPK1 | 5,619 | 0.966 | 0.938 | deployed |
| SERT | 11,429 | 0.945 | 0.835 | deployed |
| SIRT1 | 276 | 0.792 | 0.378 | deployed, low-power flag |
| Sigma1 | 7,427 | 0.881 | 0.600 | deployed |
| TAAR1 | 82 | 0.780 | 0.355 | deployed, low-power flag |
| a3b4nAChR | 760 | 0.974 | 0.906 | deployed |
| a4b2nAChR | 1,835 | 0.923 | 0.754 | deployed |
| a7nAChR | 1,285 | 0.763 | 0.571 | deployed |
| mGluR5 | 4,511 | 0.893 | 0.723 | deployed |
| mTOR | 11,510 | 0.983 | 0.941 | deployed |

6 deployed endpoints fall below the reliability gate on the corrected sensitivity figure
and carry a low-power marker on any negative call in the interface, rather than being withdrawn:
each holds its background false-positive rate at or below target, so each is weak rather than
misleading.

## 6. Endpoints withdrawn, and why

An endpoint is deployed only if a threshold exists that recovers real ligands without firing on
unrelated chemistry, and only if it discriminates better than chance against its own held-out measured
inactives. 5 endpoints were trained, tested against both conditions, and withdrawn.
They are recorded here rather than silently dropped, because a panel that reports only what survived
is a selection, not an inventory.

- **GluA2.** fires on glucose and atenolol at its calibrated threshold of 0.629, reaching 0.719 on a trivial molecule, with a random-chemistry false-positive rate of 0.072 and a sensitivity of 0.103. Withdrawn on the current fits, having previously passed.
- **NFKB1.** added to test natural-product coverage. Fires on glucose, urea, acetate, glycine and lactate at its calibrated threshold of 0.416, and scores AUROC 0.459 against its own held-out measured inactives at a sensitivity of 0.000: it recovers no active while calling five trivial metabolites binders. Of its labelled NPASS records 0.3 per cent are a direct binding constant; 344 of 362 are Potency. Same cause as NRF2.
- **NR3C1.** added to test natural-product coverage. It passes the specificity audit, firing on no trivial molecule, and fails on discrimination instead: AUROC 0.410 against its own held-out measured inactives is below chance, at a sensitivity of 0.167. The only one of the three with real binding data (14.8 per cent Ki or Kd) and the one with too few compounds to fit: 140 after deduplication.
- **NRF2.** added to test natural-product coverage. Its random-chemistry false-positive rate is 0.057, above the 5 per cent the panel holds to, and reaching that rate would need a threshold of 0.303 against a calibrated 0.301. On the refit it discriminates better than it did (AUROC 0.789 against its own held-out measured inactives, sensitivity 0.545) but it still cannot be given a cut that controls false positives. Of its labelled NPASS records 0.0 per cent are a direct binding constant (Ki or Kd); 1,019 of 1,029 are Potency, a pooled functional readout that does not define a binding class a ligand fingerprint can separate.
- **Nav1_1.** fires on glucose, urea, glycine, lactate and atenolol at its calibrated threshold of 0.571, with a random-chemistry false-positive rate of 0.080. Holding 5 per cent on random chemistry would need 0.594, and its sensitivity at the calibrated cut is already 0.120, so no threshold separates trivial metabolites from real ligands. Its AUROC of 0.952 measures ranking, which is not the quantity a deployed cut needs.

## 7. The candidate-target survey

Beyond the panel described above, 1,688 candidate targets were surveyed for whether enough
measured data exists to train and honestly validate a binder classifier (at least 60 compounds, at
least 15 per class). 66 cleared that bar
(`results/tables/np_target_survey.csv`; Figure 10, panel C, `figures/Figure10_endpoint_selection.png`).
Three of those, NRF2, NFKB1, NR3C1, were added specifically to test whether the panel's methodology
extends to targets relevant to natural-product coverage. The outcome is reported whichever way it
fell rather than only where it succeeded: NRF2 (withdrawn), NFKB1 (withdrawn), NR3C1 (withdrawn). Section 6 above gives the specific
reason each failed.

Figure 10 also shows the two conditions that decide deployment more generally: panel A shows that
more training data buys discrimination with diminishing returns, and panel B shows that
discrimination is not sufficient on its own, a usable threshold is what actually gates deployment,
which is why endpoints that rank well by AUROC can still be withheld.

## 8. Data sources and standardisation

Target-engagement and exposure labels are pooled from ChEMBL (version 37, release 2026-05-01) and
BindingDB at compound level; blood-brain-barrier labels come from B3DB augmented with FDA-curated
approved drugs; the exposure and ADME endpoints beyond BBB use measured sets from Therapeutics Data
Commons, MoleculeNet, B3DB and ChEMBL; the three natural-product-coverage endpoints (section 7) use
NPASS 3.0. No value is imputed, hand-annotated, or drawn from qualitative curator annotation; every
label is a measured experimental value, standardised identically before use: reduced to the largest
organic fragment, neutralised, sanitised, and keyed by the InChIKey of that standardised parent.
Chirality is excluded from the representation, so rows identical in feature space are collapsed
before any split is drawn. Full provenance, including which compounds came from which source per
endpoint, is in `results/tables/endpoint_rebuild_provenance.csv`, and source licences are recorded in
`data/raw/measured_endpoints_SOURCE.md`.

A compound assayed and found inactive is frequently deposited only as a censored bound
(`standard_relation` of `>` with a concentration), which the conventional pChEMBL-only query
discards. This project recovers those rows as measured non-binders wherever the whole interval they
define falls on one side of the activity cut, rather than training the negative class on
property-matched decoys; the recovery and its effect on class balance are documented in
`results/tables/expansion_inactives.csv`.

## 9. Train/test discipline

Every core and exposure endpoint is cross-validated ten-fold under both a random split and a
scaffold-grouped split (Bemis-Murcko), so the distance between the two is the honest statement of how
far a model generalises to structurally novel chemistry
(`results/tables/rf_cv_folds.csv`, `rf_cv_summary.csv`). Binder classifiers are cross-validated
scaffold-grouped only, against their measured inactives
(`results/tables/binder_cv_folds.csv`, `binder_cv_summary.csv`). For every endpoint the same compound
never appears on both sides of a fold. Beyond cross-validation, the barrier model and the wider panel
are further tested against FDA-curated approved drugs absent from the training source
(`results/tables/external_bbb_validation.csv`) and under a prospective, date-based refit
(`results/tables/external_prospective.csv`); both are reported in the manuscript's Results section
and are not repeated here.

Exactly which compounds trained which endpoint, and the numeric feature vector each one was reduced
to, is not asserted here but laid out directly: `results/tables/master_training_usage.csv` gives one
row per endpoint-compound pair actually used, and `results/tables/master_feature_vectors.csv` gives
the full 1,036-column vector for every distinct compound, keyed by InChIKey so the two join directly.
Every formula in this document that turns such a vector into a reported score, prediction set or
disease score is walked through with a live worked example in `docs/ML_METHODS_AND_FORMULAS.md`.

---
*Generated by `tools/build_endpoint_justification.py` from the artefacts named beside each figure
above. Regenerate after any change to the panel with
`python tools/build_endpoint_justification.py`.*
