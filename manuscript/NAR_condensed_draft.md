# BrainSafe AI: a calibrated, exposure-gated web server for multi-endpoint prediction of small-molecule action in the human brain

**Authors:** Krishnasalini Gunanathan¹, Raghunatha Sarma¹, Sai Shyam¹, Ramya E. M.¹,
Venketesh Sivaramakrishnan¹

¹Sri Sathya Sai Institute of Higher Learning (SSSIHL), Puttaparthi, Andhra Pradesh 515134, India

**Correspondence:** Prof. Venketesh Sivaramakrishnan, svenketesh@sssihl.edu.in

## Abstract

BrainSafe AI predicts, from structure alone, how a small molecule may act on the human brain. For a
submitted SMILES string or compound name it returns engagement of 54 molecular targets spanning the
principal neurodegenerative, psychiatric, neuroinflammatory, analgesic and sleep-related mechanisms,
predicted blood-brain barrier penetration, nine ADME and exposure endpoints including a directly
modelled unbound brain-to-plasma ratio, and a cardiac safety liability. Every target score is
admitted only in proportion to predicted brain exposure, so potency at a target a compound cannot
reach contributes nothing, and engaged targets are traced through a curated pathway graph to the
conditions those mechanisms touch. The server is built on 75 estimators, 70 deployed, trained on
228,200 measured compound-endpoint records from ChEMBL, BindingDB and B3DB. Under scaffold-grouped
10-fold cross-validation the measured-label classifiers reach a mean AUROC of 0.925, and expected
calibration error falls from 0.080 to 0.015 after isotonic calibration. The binder panel is validated
against compounds tested at the same target and found inactive rather than against decoys, reaching a
mean AUROC of 0.917. Every prediction carries a calibrated probability, an applicability-domain
distance to the nearest measured analogue, and the recall the panel achieves at that distance, which
falls from 0.86 for close analogues to 0.16 for a genuinely novel scaffold. Disease-level scores are
presented as a route from a mechanism to the conditions it touches, not as an indication prediction.
BrainSafe AI is freely available without registration at
https://huggingface.co/spaces/Krishnag999/brainsafe-ai, with source code, trained models and every
validation artefact at https://github.com/krishna-g-999/brainsafe-ai.

## Introduction

Central nervous system drug discovery fails more often, and later, than any other therapeutic area.
Two reasons dominate: a compound potent at its target may never reach the brain, and a compound that
reaches the brain may engage more than the target it was designed for. Existing servers address these
questions separately. ADMET platforms predict barrier penetration without saying what the compound
would do on arrival; target-prediction servers rank probable targets without asking whether the
compound can reach them.

BrainSafe AI answers both in one pass, and couples them: a target score is admitted only in
proportion to predicted exposure (Figure 1). It reports what it cannot do as prominently as what it can. Five
endpoints were trained, tested and withheld because no threshold separated real ligands from trivial
metabolites; four of nine falsification hypotheses were refuted; and the recall achieved on chemistry
distant from the training set is reported beside every result rather than left to be discovered.

## Materials and Methods

**Training data.** Labels are measured experimental values only, never qualitative annotation.
Potency data are ChEMBL pChEMBL values augmented with BindingDB affinities pooled at compound level;
blood-brain barrier labels come from B3DB augmented with FDA-curated approved drugs; the nine ADME
endpoints use measured sets from Therapeutics Data Commons, MoleculeNet, B3DB and ChEMBL. The panel
holds 228,200 measured compound-endpoint records over 169,341 unique compounds keyed by the InChIKey
of the desalted parent. Each endpoint is trained on its own measured set alone; across the deployed
panel those sets hold a median of 3,789 rows and span 387 compounds (KEAP1) to 10,276 (hERG).

**Recovery of the negative class.** A compound assayed and found inactive is often deposited only as
a censored bound, and the conventional query discards exactly those rows, leaving a training set of
actives against synthetic decoys. A bound settles a label when the whole interval it defines falls on
one side of the activity cut, and is discarded as undecidable when it spans both. Recovering these
returned experimentally tested non-binders to 57 endpoints.

**Representation.** Each compound is a 1,036-column vector: a 1,024-bit folded ECFP-4 fingerprint and
twelve physicochemical descriptors. Structures are reduced to the largest organic fragment and
neutralised. Neutralisation is part of the representation rather than a detail of it, because a drug
and its salt must give the same answer: removing a counter-ion without it leaves the parent carrying
the salt's charge, and haloperidol hydrochloride then scored a barrier probability of 0.613 against
0.993 for the free base. A permanent charge is retained, because a quaternary ammonium's charge is
precisely what prevents it crossing. Chirality is excluded, so two enantiomers give identical rows;
rows identical in feature space are collapsed before any split. Of 228,198 training structures, 40.8
per cent carry a stereocentre, but where one skeleton appears as several stereoisomers measured at
the same endpoint the labels agree in 94.6 per cent of 8,013 cases, so the share of the panel where
chirality could change a class call is 0.19 per cent.

**Models.** A random forest is fitted per endpoint. That choice was made on a like-for-like
comparison over thirteen core endpoints against XGBoost, histogram gradient boosting, L2 logistic
regression and a nearest-neighbour read-across, and against a graph neural network on four of
them, which the forest won on all four. Under the scaffold split the forest is best on seven of
eight classification endpoints, losing AChE to histogram gradient boosting, and on none of the
five regressions, where boosting scores higher. It was deployed for its stability under
hyperparameters, for not extrapolating beyond the training range, and because TreeSHAP is exact
for it rather than approximate. Classifiers are isotonically calibrated on out-of-fold
predictions, so no compound contributes to the calibrator that scores it. Binder models use Platt
scaling, the withheld set for one target often being too small to fit a step function.

**Thresholds.** The background library is partitioned into three disjoint pools by a stable hash of
the structure: one supplies decoys, one sets thresholds, one measures the false-positive rate.
Choosing a threshold as a quantile of a sample and then measuring the rate on that same sample
restates the target rather than testing it.

**Exposure gating and the disease layer.** A target score is admitted in proportion to predicted
barrier penetration. Engaged targets are traced through a curated graph anchored to KEGG, Reactome
and IUPHAR to the conditions they touch. The graph's edge weights were ablated and carry no
predictive information beyond the topology (curated 0.7901, uniform 0.7897, permuted 0.7874), so they
are reported as structure rather than as tuned parameters.

## Results

**Cross-validation.** The panel and its per-endpoint performance are shown in Figure 2. Under random 10-fold cross-validation the measured-label classifiers reach a
mean AUROC of 0.958 (0.899 to 0.976); under a scaffold-grouped split that withholds entire structural
classes, 0.925 (0.878 to 0.965). Expected calibration error falls from 0.080 to 0.015 after isotonic
calibration, and conformal prediction on the eight core classifiers achieves empirical coverage of
0.89 to 0.92 against a 0.90 target (Figure 3A).

**The binder panel.** The 52 binder classifiers are validated not against the decoys used to train
them but against compounds experimentally tested at the same target and found inactive. Across the 47
deployed they reach a mean AUROC of 0.917 and a mean sensitivity of 0.898 on actives withheld by
scaffold. Both are means over 47 endpoints and the spread is wide: AUROC ranges from 0.719 at GABA-A
to 0.985, sensitivity from 0.639 at COX-2 to 0.997. Five endpoints were withdrawn for firing on
trivial metabolites at every usable threshold.

**Leakage and null models.** On the deduplicated matrix the pipeline fits, no InChIKey, no feature
vector and no scaffold appears on both sides of any fold. With labels permuted the same pipeline
returns a mean AUROC of 0.4938 (random) and 0.4921 (scaffold).

**External and prospective validation.** The barrier model was tested on FDA-curated approved drugs
absent from B3DB by InChIKey: AUROC 0.764 on all 306, and 0.793 on the 241 also distinguishable from
training in feature space. For the target panel no external set of comparable size exists, because
for most of these targets the public measured chemistry is the training set. Two kinds of
independence were therefore constructed, with every model refitted rather than scored. By date: each
endpoint refitted on its pre-cutoff rows with its decision threshold also frozen before the cutoff,
and tested on compounds first published afterwards; 39 of 47 deployed endpoints qualify, giving
45,244 test compounds. By curator: compounds deposited in BindingDB and absent from ChEMBL, withheld
entirely (Figure 5).

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
through the deployed pipeline; 949 returned no actionable disease signal, a specificity of 0.949 (95%
CI 0.934 to 0.961). These compounds are presumed inactive because nothing is recorded about them, not
proven inactive, so this is a lower bound.

**Adversarial checks and falsification.** Six checks were written so that each could fail, and all
six pass. The domain-flag check initially failed and passes only after its control set was corrected:
28 of the original controls are measured compounds inside the flag's own reference library, where
calling them in domain is truthful rather than a failure. The passing criterion was not moved. Nine
falsification hypotheses were tested and four refuted, all reported: the curated edge weights add
nothing beyond topology; exposure gating cannot discriminate between diseases, being a filter;
silence at a target reflects the operating point rather than a non-discriminative model, every target
separating its own actives at AUROC 0.91 or better; and engaged targets are not independent
observations, 36 firing targets spanning only 16 independent directions.

**Use case.** A worked profile and the silence behaviour are shown in Figure 4. For donepezil the server returns AChE engagement at 1.00 against a training base rate of
0.724, barrier penetration 0.991, and Alzheimer's disease as the top condition at 0.991 with AChE
named as the driver. It also returns a hERG probability of 0.734, an enrichment of 0.652 over a 0.236
base rate, so the compound that scores highest on the efficacy axis carries a liability a medicinal
chemist would need to see. For atenolol, a beta-blocker optimised not to enter the brain, every
target probability sits below its base rate and the top condition scores 0.0009, three hundred times
below the reporting threshold: the server is correctly silent.

## Discussion

BrainSafe AI couples exposure and engagement in one pass, and reports the confidence of each. Its
distinguishing choices are that the negative class is recovered from measurement rather than
simulated, that thresholds are measured on a pool disjoint from the one that set them, and that
target scores are gated by predicted exposure.

Five limitations bound its use. First, the applicability-domain flag is a weak discriminator of
non-drug-like chemistry: in the adversarial check it scores genuinely absent chemistry at a median maximum
similarity of 0.47 against 0.59 for unseen approved drugs (n = 25, one-sided Mann-Whitney
p = 1.1e-03), but at a
threshold rejecting a tenth of genuine drugs it catches only a fifth of distant chemistry. What it
does predict well is sensitivity, the distance it measures being the variable recall tracks. Second,
recall on chemistry beyond Tanimoto 0.40 is near 0.16, so a negative result on a novel scaffold is
close to uninformative, and the server reports the expected recall beside it. Third, the specificity
estimate rests on compounds presumed rather than proven inactive. Fourth, terpenoid and steroidal
natural products are largely outside the training library, whose median fraction-sp3 is 0.34, and
such compounds are flagged as outside the domain. Fifth, the disease layer does not predict
indication: 27 of the 51 targets in the pathway graph drive more than one condition, and what selects
among them, dose, regimen and patient population, is not present in a structure.

The server does not distinguish an agonist from an antagonist. The training label is a potency value
measuring affinity, which an agonist and an antagonist at the same receptor can share, and ChEMBL's
action_type field appears nowhere in this project's data. The honest description of what the panel
predicts is engagement, not modulation.

## Data availability

The server is freely available without registration at
https://huggingface.co/spaces/Krishnag999/brainsafe-ai. Source code, trained models, every training
table and every validation artefact are at https://github.com/krishna-g-999/brainsafe-ai under the
MIT licence; underlying data retain their sources' licences. A permanent archive is deposited at
[DOI TO BE SUPPLIED].

## Funding

[TO BE SUPPLIED]

## Conflict of interest

None declared.

## Figure legends

![Figure 1](figures/Figure1_architecture.png)

**Figure 1.** How a query is answered. (**A**) A submitted structure is standardised and represented
as one fixed 1,036-column vector, scored by four model families: nine exposure and ADME endpoints,
twelve target potency and activity endpoints, the 52-endpoint binder panel of which 47 are deployed,
and two auxiliary regressions. Every target score is admitted only in proportion to the predicted probability that the
compound reaches the brain, and surviving scores are ranked by enrichment over each endpoint's base
rate rather than by raw probability. Every reported value carries a calibrated probability, a
conformal interval and an applicability-domain distance. (**B**) The counts in (**A**) are trained
estimators, 75 in total, of which 70 are deployed. Each was preceded by twenty fits that never serve
a prediction and exist only to measure how the twenty-first behaves on withheld compounds, 1,480
across the panel.

![Figure 2](figures/Figure9_model_atlas.png)

**Figure 2.** The panel, one mark per estimator, so that no claim rests on a mean a reader cannot
check. (**A**) Every estimator, deployed or withdrawn, placed by the number of compounds it was
trained on
and by the performance claimed for it, coloured by model family. Marker shape carries the metric:
AUROC and R² both run to 1.0 and are not the same quantity, since 0.5 is chance for one and a
respectable fit for the other, so they are distinguished rather than averaged. The five estimators
withdrawn after specificity testing are drawn in outline, because a panel showing only what survived
is a selection rather than an inventory. Training sets span two orders of magnitude, from 68 to 15,723 rows;
binder training sets include property-matched decoys while the others are measured compounds only,
which the axis states. (**B**) The same population by family, with the
family median marked. The spread is the point: the binder classifiers have a median of 0.935 while
the exposure and ADME regressions have a median R² of 0.575, and a single panel average would
describe neither. The complete inventory, with training-set composition, validation scheme,
calibration and fitting date for every estimator, is Supplementary Table S1.

![Figure 3](figures/Figure6_validation.png)

**Figure 3.** Four validations that a cross-validated score cannot replace. (**A**) Expected
calibration error before and after isotonic regression fitted on out-of-fold predictions, so no
compound contributes to the calibrator that scores it. (**B**) Recall on whole scaffold classes
withheld before training, with 95 per cent Wilson intervals [@wilson_ci] and marker area proportional
to the number of withheld actives, so an interval that is wide because the evidence is thin looks
thin. (**C**) Specificity on chemistry the server should stay quiet about, and external
discrimination on approved drugs absent from the training source. (**D**) The adversarial suite, in
which each check was written so that it could fail. All six pass, one of them only after its controls were corrected; each is shown at the same
size as the rest, is the applicability-domain flag, reported rather than retuned.

![Figure 4](figures/Figure8_use_case.png)

**Figure 4.** A worked profile and the silence behaviour. For four approved CNS drugs the server
recovers the pharmacologically correct driving mechanism and the corresponding condition; for four
peripherally acting compounds no disease score reaches the reporting threshold. Bars are disease
scores after exposure gating, and the driving target is named beside each. The two behaviours are the
same design decision seen from opposite sides: a target score is admitted only in proportion to
predicted barrier penetration, so a compound that does not arrive cannot generate a disease call.

Supplementary figures, each named by the file it is generated into so that the number and the
artefact cannot come apart:

![Figure 5](figures/Figure11_external_validation.png)

**Figure 5.** External validation, and an apparent temporal decay that is a composition effect.
(**A**) Per endpoint, AUROC under a size-matched random split against AUROC under a time split that
withholds the most recent quarter of the data and freezes the decision threshold before the cutoff.
The size match matters: a time split trains on less data as well as none of the future, so without a
control at the same n a drop cannot be attributed to either. (**B**) The same comparison for
sensitivity at the frozen operating point, where the gap is larger. (**C**) Why the gap exists. A
random split of medicinal-chemistry data holds out mostly close analogues of its own training set,
because the published record is series; a time split does not. The two are not testing comparable
populations. (**D**) The resolution. Recall against maximum Tanimoto similarity to the training
actives, for three test sets built by unrelated rules: withheld by publication date, withheld at
random, and withheld by curator, the last being compounds deposited in BindingDB and absent from
ChEMBL. They trace one curve, so recall is a function of chemical distance rather than of how the set
was held out, and the expected sensitivity for a submitted compound is knowable at query time.
