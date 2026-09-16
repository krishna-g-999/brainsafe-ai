# BrainSafe AI: a calibrated, exposure-gated web server for multi-endpoint prediction of small-molecule action in the human brain

**Authors:** Krishnasalini Gunanathan¹, Raghunatha Sarma¹, Sai Shyam¹, Ramya E. M.¹,
Venketesh Sivaramakrishnan¹

¹Sri Sathya Sai Institute of Higher Learning (SSSIHL), Puttaparthi, Andhra Pradesh 515134, India

**Correspondence:** Prof. Venketesh Sivaramakrishnan, svenketesh@sssihl.edu.in

**Manuscript type:** NAR Web Server Issue.

**Graphical abstract:** submitted as a separate file, `manuscript/figures/GraphicalAbstract.png`
(and `.pdf`).

**Key Points**

- BrainSafe AI profiles small-molecule mechanism in the human brain from structure alone, gating
  every target score by predicted exposure so potency at a target a compound cannot reach
  contributes nothing to the output.
- Its 74 cross-validated estimators are each scaffold-grouped cross-validated and carry a
  calibrated probability; the 22 target-potency, exposure and safety estimators are additionally
  validated under a random split, the core classifiers reporting a conformal interval as well.
- The negative class is recovered from compounds measured and found inactive rather than simulated
  with decoys, and every validation is reported whichever way it falls, which led to withdrawing
  two endpoints that could not separate a real ligand from an unrelated metabolite.

## Abstract

BrainSafe AI predicts, from structure alone, how a small molecule may act on the human brain. For a
submitted SMILES string or compound name it returns engagement of 54 molecular targets spanning the
principal neurodegenerative, psychiatric, neuroinflammatory, analgesic and sleep-related mechanisms,
predicted blood-brain barrier penetration, nine ADME and exposure endpoints including a directly
modelled unbound brain-to-plasma ratio, and two cardiac safety liabilities. Every target score is
admitted only in proportion to predicted brain exposure, so potency at a target a compound cannot
reach contributes nothing, and engaged targets are traced through a curated pathway graph to the
conditions those mechanisms touch. The server is built on 75 estimators, 70 deployed, trained on
228,200 measured compound-endpoint records from ChEMBL [@chembl], BindingDB [@bindingdb] and B3DB
[@b3db]. Under scaffold-grouped 10-fold cross-validation the measured-label classifiers reach a mean
AUROC of 0.925 (0.958 under a random split), and expected calibration error falls from 0.0801 to
0.0147 after isotonic calibration [@calibration]. The binder panel is validated against compounds
tested at the same target and found inactive rather than against decoys, reaching a mean AUROC of
0.917. Every prediction carries a calibrated probability, a conformal interval, and an
applicability-domain distance to the nearest measured analogue, and the server reports silence
rather than a guess for compounds outside its competence: on non-CNS chemistry its specificity is
0.925 (95% CI 0.907 to 0.940). Disease-level scores are presented as a route from a mechanism to the
conditions it touches, not as an indication prediction. BrainSafe AI is freely available without
registration at https://huggingface.co/spaces/Krishnag999/brainsafe-ai, with source code, trained
models and every validation artefact at https://github.com/krishna-g-999/brainsafe-ai.

## Introduction

Central nervous system drug discovery fails more often, and later, than any other therapeutic area.
Two reasons dominate: a compound potent at its target may never reach the brain, and a compound that
reaches the brain may engage more than the target it was designed for. Existing servers address these
questions separately. ADMET platforms predict barrier penetration without saying what the compound
would do on arrival; target-prediction servers rank probable targets without asking whether the
compound can reach them.

BrainSafe AI answers both in one pass, and couples them: a target score is admitted only in
proportion to predicted exposure (Figure 1). It reports what it cannot do as prominently as what it
can. Five endpoints were trained, tested and withdrawn because no threshold separated real ligands
from trivial chemistry; every validation is written so that it could fail and is reported whichever
way it falls; and the recall achieved on chemistry distant from the training set is reported beside
every result rather than left to be discovered.

## Materials and Methods

**Training data.** Labels are measured experimental values only, never qualitative annotation.
Potency data are ChEMBL pChEMBL values augmented with BindingDB affinities pooled at compound level;
blood-brain barrier labels come from B3DB augmented with FDA-curated approved drugs; the nine ADME
endpoints use measured sets from Therapeutics Data Commons [@tdc], MoleculeNet [@moleculenet], B3DB
and ChEMBL. The panel
holds 228,200 measured compound-endpoint records over 169,341 unique compounds keyed by the InChIKey
of the desalted parent. Each endpoint is trained on its own measured set alone; across the deployed
panel those sets hold a median of 3,789 rows and span 387 compounds (KEAP1) to 10,276 (hERG).
That set is every deployed model owning a table in data/endpoints, the barrier model included;
excluding it, since it is an exposure rather than a target endpoint, gives 54 tables and a
median of 3,587.

**Recovery of the negative class.** A compound assayed and found inactive is often deposited only as
a censored bound, and the conventional query discards exactly those rows, leaving a training set of
actives against synthetic decoys. A bound settles a label when the whole interval it defines falls on
one side of the activity cut, and is discarded as undecidable when it spans both. Recovering these
returned experimentally tested non-binders to 57 endpoints.

**Representation.** Each compound is a 1,036-column vector: a 1,024-bit folded ECFP-4 fingerprint
[@ecfp] and twelve physicochemical descriptors. Structures are reduced to the largest organic fragment and
neutralised. Neutralisation is part of the representation rather than a detail of it, because a drug
and its salt must give the same answer: removing a counter-ion without it leaves the parent carrying
the salt's charge, and haloperidol hydrochloride then scored a barrier probability of 0.613 against
0.993 for the free base. A permanent charge is retained, because a quaternary ammonium's charge is
precisely what prevents it crossing. Chirality is excluded, so two enantiomers give identical rows;
rows identical in feature space are collapsed before any split. Of 228,198 training structures, 40.8
per cent carry a stereocentre, but where one skeleton appears as several stereoisomers measured at
the same endpoint the labels agree in 94.6 per cent of 8,013 cases, so the share of the panel where
chirality could change a class call is 0.19 per cent.

**Models.** A random forest [@random_forest] is fitted per endpoint. That choice was made on a
like-for-like comparison over thirteen core endpoints against XGBoost [@xgboost], histogram gradient
boosting, L2 logistic regression and a nearest-neighbour read-across, and against a graph neural
network on four of them, which the forest won on all four. Under the scaffold split the forest is
best on seven of eight classification endpoints, losing AChE to histogram gradient boosting, and on
none of the five regressions, where boosting scores higher. It was deployed for its stability under
hyperparameters, for not extrapolating beyond the training range, and because TreeSHAP [@shap_trees]
is exact for it rather than approximate. Classifiers are isotonically calibrated [@calibration] on
out-of-fold predictions, so no compound contributes to the calibrator that scores it. Binder models
use Platt scaling [@platt], the withheld set for one target often being too small to fit a step
function.

**Thresholds.** The background library is partitioned into three disjoint pools by a stable hash of
the structure: one supplies property-matched decoys [@dude], one sets thresholds, one measures the
false-positive rate. Choosing a threshold as a quantile of a sample and then measuring the rate on
that same sample restates the target rather than testing it.

**Exposure gating and the disease layer.** A target score is admitted in proportion to predicted
barrier penetration. Engaged targets are traced through a curated graph anchored to KEGG [@kegg],
Reactome [@reactome] and IUPHAR [@iuphar] to the conditions they touch. The graph's edge weights were
ablated and carry no predictive information beyond the topology (curated 0.7901, uniform 0.7897,
permuted 0.7874), so they are reported as structure rather than as tuned parameters.

## Results

**Cross-validation.** The panel and its per-endpoint performance are shown in the model atlas
(`figures/Figure9_model_atlas.png`). Under random 10-fold cross-validation the measured-label
classifiers reach a mean AUROC of 0.958 (0.899 to 0.976); under a scaffold-grouped split that
withholds entire structural classes, 0.925 (0.878 to 0.965). Expected calibration error falls from
0.0801 to 0.0147 after isotonic calibration, and conformal prediction [@conformal] on the eight core
classifiers, on the deduplicated matrix the classifiers are trained on, achieves empirical coverage
of 0.876 to 0.933 against a 0.90 target, with mean set size from 0.956 to 1.215 on a two-class
problem where 1.0 is a confident single label.

**The binder panel.** The 52 binder classifiers are validated not against the decoys used to train
them but against compounds experimentally tested at the same target and found inactive. Across the 47
deployed they reach a mean AUROC of 0.917 and a mean sensitivity of 0.764, both on actives withheld
by scaffold. Both are means over 47 endpoints and the spread is wide: AUROC ranges from 0.719 at
GABA-A to 0.985 at CGRP, sensitivity from 0.303 at GABA-A to 0.993 at CGRP with a median of 0.835.

That sensitivity figure carries a correction we report rather than leave for a reader to find. Four
scripts write it in sequence, and the last scored every active in the endpoint table, so roughly four
fifths of the scoring set were compounds the model had been fitted on, while the field's own basis
annotation went on reading that the measurement was held out. The panel mean was published as 0.898
against a held-out 0.764, an inflation of 0.134. The registry has been recomputed from the
scaffold-held-out partition, it agrees with the threshold table on all 47 deployed endpoints, and a
regression test fails if any endpoint reports sensitivity on an undeclared basis. Six endpoints fall
below the reliability gate on the corrected figure, against one before; each holds its background
false-positive rate at or below target, so each is weak rather than misleading and stays deployed
with a low-power marker on any negative call. Five further endpoints are withdrawn outright: two for
firing on trivial metabolites at every usable threshold, and three added specifically to test
natural-product coverage, reported in the limitations.

**Leakage and null models.** On the deduplicated matrix the pipeline fits, no InChIKey, no feature
vector and no scaffold appears on both sides of any fold. With labels permuted the same pipeline
returns a mean AUROC of 0.496 (random) and 0.503 (scaffold) over the eight core classifiers, every
one of the sixteen values within 0.02 of chance and a smallest margin over its own null 19-fold
larger than the largest such departure.

**External and prospective validation.** The barrier model was tested on FDA-curated approved drugs
absent from B3DB by InChIKey: AUROC 0.764 on all 306, and 0.767 on the 227 also distinguishable from
training in feature space, once compounds the featuriser cannot tell apart from a training row are
also excluded. For the target panel no external set of comparable size exists, because
for most of these targets the public measured chemistry is the training set. Two kinds of
independence were therefore constructed, with every model refitted rather than scored. By date: each
endpoint refitted on its pre-cutoff rows with its decision threshold also frozen before the cutoff,
and tested on compounds first published afterwards; 39 of 47 deployed endpoints qualify, giving
45,244 test compounds. By curator: compounds deposited in BindingDB and absent from ChEMBL, withheld
entirely (`figures/Figure11_external_validation.png`).

Read in aggregate the time split suggests prospective decay, with mean AUROC 0.823 against 0.951 for
a size-matched random control. It is not decay. The false-positive rate on background chemistry is
unchanged, and the gap closes once test compounds are stratified by maximum Tanimoto similarity to
the training actives. A random split of medicinal-chemistry data holds out 83 per cent close
analogues of its own training set, because the published record is series; a time split holds out 28
per cent chemistry below Tanimoto 0.40. Three test sets built by unrelated rules trace one recall
curve: 0.16, 0.55, 0.74 and 0.86 by date across four novelty bands, against 0.12, 0.52, 0.77 and 0.93
at random and 0.05, 0.46, 0.83 and 0.90 by curator. Recall is therefore a function of chemical
distance rather than of publication date, and the expected recall for a submitted compound is
reported at query time.

**Specificity.** One thousand compounds with no recorded activity at any modelled target were scored
through the deployed pipeline; 925 returned no actionable disease signal, a specificity of 0.925 (95%
CI 0.907 to 0.940). These compounds are presumed inactive because nothing is recorded about them, not
proven inactive, so this is a lower bound.

**Adversarial checks and falsification.** Six checks were written so that each could fail, and all
six pass. The domain-flag check initially failed and passes only after its control set was corrected:
28 of the original controls are measured compounds inside the flag's own reference library, where
calling them in domain is truthful rather than a failure. The passing criterion was not moved. A
further falsification analysis found that the curated edge weights add nothing beyond topology
(curated 0.7901, uniform 0.7897, permuted 0.7874); that silence at a target reflects the operating
point rather than a non-discriminative model, every target separating its own actives at AUROC 0.91
or better; and that engaged targets are not independent observations, 37 firing targets spanning
only 16 independent directions. It also found a deployed endpoint calling common metabolites
binders at its calibrated threshold, and that endpoint is withdrawn.

**Use case.** A worked profile and the silence behaviour are shown in Figure 2. For donepezil the
server returns AChE engagement at 1.00 against a training base rate of 0.596, barrier penetration
0.991, and Alzheimer's disease as the top condition at 0.991 with AChE named as the driver. It also
returns a hERG probability of 0.734, an enrichment of 0.652 over a 0.236 base rate, so the compound
that scores highest on the efficacy axis carries a liability a medicinal chemist would need to see.
For atenolol, a beta-blocker optimised not to enter the brain, every target probability sits below
its base rate and the top condition scores exactly 0, well below the reporting threshold: the server
is correctly silent.

## Discussion

BrainSafe AI is, to our knowledge, the first freely available server to couple predicted brain
exposure and target engagement in one calibrated pass rather than reporting them side by side. The
coupling holds up under scrutiny built to break it: recall traces one curve across three test sets
constructed by unrelated rules (by date, at random, and by curator), the binder panel separates real
ligands from measured non-binders at a mean AUROC of 0.917 without ever training against a decoy, and
all six adversarial checks pass, including one that failed on its first run and was fixed at the
control set rather than at the criterion. On approved CNS drugs it recovers the pharmacologically
correct mechanism every time it is asked, and on compounds built not to enter the brain it stays
silent rather than guessing. Its distinguishing design choices, that the negative class is recovered
from measurement rather than simulated, that thresholds are measured on a pool disjoint from the one
that set them, and that every target score is gated by predicted exposure, are what make that
combination possible rather than coincidental.

Five limitations bound its use, reported with the same rigour as the results above. First, the
applicability-domain flag is a weak discriminator of
non-drug-like chemistry: in the adversarial check it scores genuinely absent chemistry at a median
maximum similarity of 0.47 against 0.57 for unseen approved drugs (n = 25, one-sided Mann-Whitney
p = 1.8e-03), but at a threshold rejecting a tenth of genuine drugs it catches only a fifth of
distant chemistry. What it
does predict well is sensitivity, the distance it measures being the variable recall tracks. Recall
on chemistry beyond Tanimoto 0.40 is also near 0.16, so a negative result on a novel scaffold is close
to uninformative, and the server reports the expected recall beside it. The specificity estimate
rests on compounds presumed rather than proven inactive. Terpenoid and steroidal natural products,
in turn, are largely outside the training library, whose median fraction-sp3 is 0.34, and such
compounds are flagged as outside the domain. Last, the disease layer does not predict indication:
27 of the 51 targets in the pathway graph drive more than one condition, and what selects among them,
dose, regimen and patient population, is not present in a structure.

The server does not distinguish an agonist from an antagonist. The training label is a potency value
measuring affinity, which an agonist and an antagonist at the same receptor can share, and ChEMBL's
action_type field appears nowhere in this project's data. The honest description of what the panel
predicts is engagement, not modulation.

## Data availability

BrainSafe AI is freely available without registration or login at
https://huggingface.co/spaces/Krishnag999/brainsafe-ai, served over HTTPS; no submitted structure is
retained beyond the lifetime of its request. Source code, trained models, every training table and
every validation artefact are at https://github.com/krishna-g-999/brainsafe-ai under the MIT
licence; underlying data retain their own sources' licences. Trained estimators and the raw API
responses are deposited separately with a manifest recording the SHA-256 of every file. **[TO BE
SUPPLIED BEFORE SUBMISSION]** the archive deposit's own DOI.

## Author contributions

**[TO BE SUPPLIED BEFORE SUBMISSION]** NAR requires a CRediT-style contribution statement naming
what each author did. This must be supplied by the author team.

## Funding

[TO BE SUPPLIED]

## Conflict of interest

None declared.

<!-- REFERENCES -->



## Figure legends

![Figure 1](figures/Figure1_architecture.png)

**Figure 1.** How a query is answered. (**A**) A submitted structure is standardised and scored by
four model families: a 52-endpoint binder panel (47 deployed), twelve target-potency endpoints, a
ten-endpoint exposure/ADME layer, and one safety classifier. Every target score is admitted only in
proportion to predicted brain exposure and ranked by enrichment over its own base rate, and every
value carries a calibrated probability and an applicability-domain distance. (**B**) Each deployed
estimator is the last of several fits kept only to measure held-out behaviour: 960 fits stand behind
the 70-estimator deployed panel.

![Figure 2](figures/Figure8_use_case.png)

**Figure 2.** A worked profile and the silence behaviour. For four approved CNS drugs the server
recovers the pharmacologically correct driving mechanism and condition; for four peripherally acting
compounds no disease score reaches the reporting threshold. A target score is admitted only in
proportion to predicted barrier penetration, so a compound that does not arrive cannot generate a
call.

Three further validation figures, referenced in the Results above, are given in full in the primary
manuscript and reproduce from the repository named beside each: the per-estimator panel with no
claim resting on a mean a reader cannot check (`figures/Figure9_model_atlas.png`), the four
calibration, coverage, specificity and adversarial-check validations
(`figures/Figure6_validation.png`), and the external and prospective validation showing that an
apparent temporal decay is a chemical-distance effect rather than model drift
(`figures/Figure11_external_validation.png`).
