# Chapter 1. Introduction

> Every quantitative statement in this chapter is read from an artefact in this repository. Where a
> claim requires a figure from the external literature that this repository does not yet hold a
> verified citation for, the gap is marked `[REF NEEDED]` rather than filled from memory. Citation
> numbers in square brackets refer to `manuscript/references.md`, whose entries were each resolved by
> exact-title query against CrossRef or Europe PMC at a normalised title similarity above 0.82.

---

## 1.1 The problem this thesis addresses

Drug discovery for the central nervous system fails more often, and later, than discovery for any
other therapeutic area `[REF NEEDED: a citable attrition series, for example Kola and Landis, or a
CNS-specific success-rate analysis]`. The failures are expensive precisely because they are late: a
compound that survives to a clinical readout has already consumed the medicinal chemistry, the
toxicology and the manufacturing that a first-in-human study requires.

Two distinct mechanisms account for a large share of that loss, and it is worth separating them
carefully, because the rest of this thesis is organised around the fact that they are usually
addressed apart.

The first is **exposure**. A compound may be potent at a well-validated target and never reach it.
The blood-brain barrier is not a passive membrane with a permeability constant but an active
interface, with efflux transporters that return compounds to the circulation and tight junctions that
exclude them in the first place [8]. A programme that optimises affinity without tracking exposure
can produce a compound that is excellent against isolated protein and inert in an animal.

The second is **engagement**, understood broadly. A compound that does arrive may engage more than
the target it was designed for. Sometimes that is fatal on the safety axis: blockade of the hERG
potassium channel prolongs the QT interval and remains a leading cause of late cardiovascular
attrition [13]. Sometimes it is fatal on the efficacy axis, when the mechanism reached is not the
mechanism that drives the disease. Either way the question is the same: given that the compound is
present in brain tissue, what does it touch there?

Neither question is new, and neither is unaddressed. What is unusual is how rarely they are answered
in the same breath, by the same instrument, on the same molecule, with a statement of confidence that
covers both.

## 1.2 Exposure and engagement are not separable questions

The argument for coupling them is not one of convenience. It is that the quantity a medicinal chemist
actually needs is a product, and neither factor is interpretable alone.

Consider what a target-prediction result means without an exposure term. A statement that a compound
binds a named central target with high probability is, for a compound that does not cross the
barrier, not a weaker claim than it appears. It is a claim about a counterfactual world in which the
compound is delivered directly to the tissue. The published record contains many such compounds:
peripherally restricted agents designed precisely so that their central pharmacology never occurs. A
tool that reports the affinity without the restriction has not made a small error, it has answered a
question nobody asked.

The converse is equally true. A barrier-penetration probability, on its own, says nothing about
whether the compound should be pursued. Caffeine and a solvent may both cross. The exposure axis
becomes a decision variable only when there is a mechanism on the other side of it to be gated.

This thesis therefore treats the disease-level output as a gated quantity. The formal statement is in
Chapter 7, but the shape of it belongs here, because it determines what the rest of the work has to
establish. Writing $s_t(x)$ for the engagement signal at target $t$ and $\gamma(x)$ for predicted
barrier penetration, the score for a condition $d$ is

$$\tilde S_d(x) \;=\; \gamma(x)\cdot \max_{(t,w)\in G(d)} w\, s_t(x),$$

with $\gamma(x)=1$ for the two conditions in the graph whose mechanism is peripheral, multiple
sclerosis and migraine, and $\gamma(x)=\hat q_{\mathrm{BBB}}(x)$ otherwise. Two properties of that
expression matter from the outset.

First, the gate is **multiplicative and shared**. Potency at an unreachable target contributes
essentially nothing, which is the pharmacological content of the design. This is what makes the
system a coupled instrument rather than two instruments printed on one page.

Second, and this is a result of the falsification programme rather than an intention of the design,
the gate **cannot discriminate between conditions**. Because $\gamma$ does not depend on $d$ across
the non-peripheral conditions, multiplying every condition's score by the same number leaves their
ranking untouched. Gating decides whether anything is reported at all; it never decides which
condition is reported. This was posed as a falsifiable hypothesis (H3) and recorded as refuted by
construction, and it is stated here, in the introduction, because the honest description of the
architecture is available from the algebra and there is no reason to defer it to a limitations
section. The system is an exposure filter over a mechanism ranking, and it should be described that
way throughout.

A third consequence follows from taking exposure seriously as a physiological quantity rather than a
convenient binary. The variable that governs whether a free concentration is available at the target
is the unbound brain-to-plasma partition coefficient, $K_{p,uu}$, not the total brain-to-plasma ratio
[29]. Total ratio counts drug bound to tissue lipid, which cannot engage a receptor. This thesis
therefore models $K_{p,uu}$ directly. The honest report of that model is that it is the weakest
member of the exposure layer: under a scaffold-grouped split it reaches $R^2 = 0.3523$ with a
fold-to-fold standard deviation of 0.1583 on 566 measured compounds, against $R^2 = 0.4131$ for
logBB on 1,058 and $R^2 = 0.7263$ for aqueous solubility on 9,573
(`results/tables/adme_cv_summary.csv`). Modelling the right quantity badly is preferable to modelling
the wrong quantity well, but only if the weakness is reported, and Chapter 7 reports it.

## 1.3 What existing servers do, and what they do not do

Publicly available prediction servers for this problem fall into three groups, and the boundary
between them is the boundary this thesis is trying to cross.

**Physicochemical and ADMET predictors** estimate barrier penetration, permeability, efflux liability
and related properties from structure. They answer the exposure question, often well, and they say
nothing about pharmacology. The best-known heuristic in this family, the CNS multiparameter
optimisation score, is explicitly a desirability function over six physicochemical properties and was
presented as a design aid rather than a model of activity [11]. Used as intended it is valuable. Used
as a proxy for central activity it is a category error, since it contains no information about any
target.

**Target-prediction servers** rank probable protein targets for a submitted structure, typically by
similarity to annotated ligands. They answer the engagement question, and they are agnostic about
whether the compound reaches the tissue in which those targets sit. Similarity-based target
prediction is a strong baseline and this thesis treats it as such rather than as a straw man: a
five-nearest-neighbour Tanimoto read-across is one of the five model families benchmarked in Chapter
3, and the falsification suite shows it recovering the correct target for 0.9726 of held-out
compounds against 0.0587 for a frequency baseline when the target family is represented in the index
(`inversion/results/H5_readacross_value.csv`). That is a genuine capability. What it does not
include is exposure.

**Single-endpoint barrier predictors** model blood-brain penetration alone, usually on one of the
curated permeability datasets [27]. They are the exposure axis in isolation.

`[REF NEEDED: verified citations for the specific servers named in the comparison, namely SwissADME,
ADMETlab, pkCSM, admetSAR and SwissTargetPrediction. None is currently in references.md, and the
comparison table in docs/BS_BENCHMARK_ANALYSIS.md carries no citations. These must be resolved
through the same CrossRef and Europe PMC pipeline that produced the existing thirty-two entries
before this section can be submitted.]`

Three capabilities are absent across all three groups, and together they define the contribution
claimed here.

The first is **coupling**. No server in these families admits a target score in proportion to
predicted exposure. Where both quantities are produced they are produced side by side, and the
integration is left to the reader.

The second is **calibrated, compound-specific uncertainty**. A ranked list without a probability is
not actionable, and a probability that is not calibrated is not a probability. This thesis reports a
calibrated probability, an empirically verified conformal coverage statement [14], and an
applicability-domain distance [7] for every value it returns.

The third, and the one that has proved most consequential in the work, is a **quantified statement of
what silence means**. When a target-prediction tool returns nothing, the user does not know whether
the compound is inactive or whether the model cannot see it. This system reports, beside every
result, the recall the panel actually achieves at that compound's own distance from the training
chemistry. That figure runs from 0.862 for a close analogue to 0.161 for a compound below Tanimoto
0.40 from anything measured (`results/tables/external_novelty_strata.csv`). A silence accompanied by
a recall of 0.16 is not a negative result, and the interface says so.

## 1.4 The approach taken in this thesis

The system described in the following chapters is a panel of independently fitted estimators combined
by a stated rule rather than a learned one.

At the time of writing, `submission_package/07_MODELS/model_inventory.csv` holds **75 fitted
estimators, 70 of them deployed**: 52 binder classifiers, 7 exposure regressions, 6 target
regressions, 6 target classifiers, 3 exposure classifiers and 1 safety classifier. Counting distinct
molecular targets rather than estimators, and excluding the four receptor regressions that duplicate
binder endpoints together with the antioxidant assay and the pKa property model, the deployed panel
covers **54 molecular targets**: the 47 deployed binders plus acetylcholinesterase,
butyrylcholinesterase, BACE1, GSK-3β, MAO-A, MAO-B and hERG.

Training data are measured experimental values only. Counting the endpoint tables directly gives
**228,200 measured compound-endpoint records over 63 tables**, drawn from ChEMBL [30], BindingDB
[19], B3DB [27], Therapeutics Data Commons [28] and MoleculeNet [23]. Each endpoint is fitted on its
own set alone; across the 55 deployed classification tables the median is 3,789 rows, from 387 for
KEAP1 to 10,276 for hERG. No label comes from curator annotation, a decision whose origin is recorded
in `docs/decisions_log.md`: an earlier prototype trained on curated annotation scores was shown by
feature ablation to be reading the answer back out of its own features.

Three methodological decisions distinguish the work, and each is the subject of a later chapter.

**The negative class is recovered from measurement rather than simulated.** A compound assayed and
found inactive is frequently deposited only as a censored bound, and the conventional pChEMBL query
discards exactly those rows, leaving a positive class drawn from measurement and a negative class
drawn from property-matched decoys [12]. A censored bound settles a label whenever the whole interval
it defines falls on one side of the activity cut, and is discarded as undecidable when it spans both.
Recovering these yields **29,751 measured non-binders across 57 endpoints**, counted from
`submission_package/06_TRAINING_DATA/endpoints/`. Chapter 2 develops this.

**Decision thresholds are measured on a sample disjoint from the one that set them.** Choosing a cut
as a quantile of a sample and then measuring the false-positive rate on that same sample cannot fail,
because the rate restates the quantile. The background library is therefore partitioned into three
disjoint pools by a stable hash of the canonical structure, one supplying decoys, one setting
thresholds and one measuring the rate achieved. That the measured rate can disagree with its target
is the evidence that it is a measurement. Chapter 5 develops this.

**Target scores are gated by predicted exposure**, in the sense set out in section 1.2. Chapter 7
develops this, and Chapter 9 reports what the falsification programme found the gate can and cannot
do.

## 1.5 What is claimed, and what is not

A thesis introduction that lists only capabilities produces a document a reader cannot check. The
principal negative results are therefore stated here rather than deferred, since several of them
constrain how the positive results should be read.

**Discrimination is strong under cross-validation and the spread is wide.** Over the eight
measured-label classifiers, mean AUROC is 0.958 under a random split, ranging 0.899 for the barrier
model to 0.976 for BACE1, and 0.925 under a scaffold-grouped split that withholds entire structural
classes, ranging 0.878 to 0.965 (`results/tables/rf_cv_summary.csv`). Over the 47 deployed binder
classifiers, validated against compounds tested at the same target and found inactive rather than
against decoys, mean AUROC is 0.9174 with a median of 0.947, ranging from 0.719 at GABA-A to 0.985 at
the CGRP receptor; mean sensitivity is 0.8983 with a median of 0.932, ranging from 0.639 at COX-2 to
0.997 (`submission_package/07_MODELS/binder_panel_registry.json`). Quoting either mean without its
range would flatter the panel, and the endpoints at the bottom of those ranges are named wherever the
means appear.

**Five endpoints were trained, tested and withdrawn**, and the panel is reported as an inventory of
52 binder endpoints rather than as a selection of 47. GluA2 and Nav1.1 fired on glucose, urea and
atenolol at their calibrated thresholds. NRF2, NFKB1 and NR3C1 were added specifically to test
natural-product coverage and all three failed, NFKB1 reaching an AUROC of 0.459 against its own
held-out measured inactives.

**Four of nine falsification hypotheses were refuted** and are reported as refuted
(`inversion/results/VERDICTS.csv`). The curated pathway-graph edge weights carry no measurable
information beyond the graph's topology, scoring 0.7901 curated against 0.7897 uniform and 0.7874
randomly permuted. Exposure gating cannot discriminate between conditions. Silence at a target
reflects the operating point rather than a non-discriminative model. And engaged targets are not
independent observations: 36 targets fire across approved drugs but span only 16 independent
directions.

**Specificity is a lower bound.** On 1,000 compounds with no recorded activity at any modelled
target, 949 returned no actionable disease signal, a specificity of 0.949 with a 95 per cent interval
of 0.9336 to 0.961 (`results/tables/noncns_specificity_summary.csv`). Those compounds are presumed
inactive because nothing is recorded about them, not proven inactive, and the artefact itself labels
the paired false-positive rate an upper bound.

**Recall on genuinely novel chemistry is poor, and the finding of this thesis is that it is
predictable rather than that it is better than it looked.** Below Tanimoto 0.40 from the training
actives, recall is 0.161. No analysis presented here improves that number.

**Three things the system does not do.** It does not resolve chirality, so two enantiomers receive
identical predictions, a limitation bounded in Chapter 10 at 0.19 per cent of the panel. It does not
distinguish an agonist from an antagonist, because the training label is an affinity that both can
share. And it does not predict clinical indication: 27 of the 51 targets in the pathway graph drive
more than one of its 16 conditions, and what selects among them is dose, regimen and patient
population, none of which is present in a structure.

The honest summary of the contribution is that it is one of integration, validation and uncertainty
reporting rather than of algorithm. The estimator is a random forest [4] over an ECFP-4 fingerprint
[10] and twelve descriptors. None of those components is new. What has not been assembled elsewhere,
so far as this work has been able to establish, is the coupling of a measured-label CNS target panel
to a predicted exposure term, with a calibrated probability, a verified coverage statement and a
distance-conditioned recall attached to every value, and with the failures reported at the same size
as the successes.

## 1.6 Structure of the thesis

Chapter 2 sets out the data: the sources, the recovery of the negative class from censored bounds,
the deduplication that must precede any split, and the per-endpoint sizes rather than the aggregate.

Chapter 3 covers the representation and the estimator, including the like-for-like comparison against
four alternative model families over thirteen endpoints, which is also where the limits of the random
forest are stated: it leads classification under the scaffold split and it does not lead regression.

Chapter 4 covers calibration, conformal prediction and the applicability domain, the three components
of the uncertainty statement.

Chapter 5 covers thresholds and the three disjoint background pools.

Chapter 6 covers the binder panel and its validation against measured inactives rather than decoys.

Chapter 7 covers exposure gating and the pathway graph.

Chapter 8 covers validation, including the external barrier test, the prospective simulation across
39 of the 47 deployed endpoints, and the composition finding that explains an apparent temporal decay
as a change in the population being tested rather than a decay in the models.

Chapter 9 covers the falsification suite in full, including the four refutations and what each cost.

Chapter 10 covers the limitations and the work that follows from them.

---

## Outstanding items for this chapter

1. `[REF NEEDED]` A citable CNS attrition series for section 1.1.
2. `[REF NEEDED]` Verified citations for the named comparison servers in section 1.3.
3. The count of unique compounds underlying the 228,200 records is quoted elsewhere in the project as
   169,341, keyed by the InChIKey of the desalted parent. Counting distinct SMILES strings across the
   endpoint tables in this session gives 170,619, which is consistent with the two being different
   quantities but does not verify the InChIKey figure. It has not been used in this chapter and
   should be re-derived from the tables before it is used anywhere.
