# BrainSafe AI: a calibrated, applicability-aware web server for multi-endpoint prediction of small-molecule effects on the human brain

**Authors:** [Author list to be finalised]

**Affiliation:** SAI-Net initiative, Sri Sathya Sai Institute of Higher Learning (SSSIHL), Prasanthi Nilayam, India.

**Correspondence:** [email]

**Manuscript type:** NAR Web Server Issue.

**Repository:** https://github.com/krishna-g-999/brainsafe-ai

---

## Abstract

Deciding whether a small molecule is likely to act on the brain requires answering several
questions at once: can it cross the blood-brain barrier, which disease-relevant targets does it
engage, does an achievable dose deliver free drug to the central nervous system, and is it safe.
BrainSafe AI is an open web server that answers these questions from chemical structure alone. It
integrates 72 deployed estimators trained on 227,146 measured compound-endpoint records drawn from
ChEMBL (1), BindingDB (2) and the B3DB blood-brain-barrier database (3): 60 molecular-target
endpoints spanning blood-brain-barrier penetration, the principal neurodegenerative, psychiatric,
neuroinflammatory, epileptic, analgesic, migraine, demyelinating and sleep-related target classes,
and two cardiac safety liabilities, together with a nine-endpoint ADME and exposure layer that
includes a directly modelled unbound brain-to-plasma partition coefficient (4). Every compound is
represented by one fixed 1,036-column vector, a 1,024-bit folded ECFP-4 fingerprint (5) with twelve
physicochemical descriptors, and every endpoint is a random forest (6) validated under 10-fold
cross-validation in two regimes: a random split, and a split grouped on Bemis-Murcko scaffolds (7)
that withholds entire chemical series. Across 71 cross-validated endpoints this is 1,420 fitted
models behind the 72 that are deployed. Predictions are probability-calibrated, carry an
endpoint-specific applicability-domain flag with the nearest measured analogue, and are combined into
per-disease relevance scores, filtered by blood-brain-barrier exposure and traced through a curated
target-to-pathway-to-disease knowledge graph. The measured-label classifier panel reaches a mean
AUROC of 0.958 on the random split and 0.925 on the scaffold split, and mean expected calibration
error falls from 0.0795 to 0.0161 after isotonic calibration (8). The 49 binder classifiers are validated
not against the decoys used to train them (9) but against compounds experimentally tested on the same
target and found inactive, giving a mean AUROC of 0.904 and a mean sensitivity of 0.866 on actives
withheld by scaffold, at thresholds constrained simultaneously by held-out measured inactives and by
the false-positive rate on a disjoint pool of unrelated chemistry.

A negative class recovered from censored measurements is central to these figures. Public databases
record what bound; a compound assayed and found inactive is often stored only as a bound such as
"IC50 > 10 uM", which the standard potency query discards. Recovering the bounds that settle a label
whatever the true value is added 21,994 measured non-binders and reduced the endpoints exceeding 90
per cent active from 35 to 13. Classification became marginally harder as a result and regression
improved, both reported here rather than only the improvement.

Beyond conventional validation, the server was subjected to a systematic falsification analysis in
which each of its central claims was paired with a null model capable of reproducing the same
apparent success by accident. This recovered results in both directions. The disease layer is
informative: top-3 accuracy 0.804 against a permutation null of 0.157 and a frequency null of 0.548.
Its curated edge weights are not: uniform and randomly permuted weights score 0.8025 and 0.8023
against 0.8043, so the predictive content lies in the graph topology and the weights are reported as
structure rather than as tuned parameters. Validated against clinical indications drawn from ChEMBL
rather than against the tool's own target-to-disease map, and restricted to the 162 approved drugs
whose structure appears nowhere in the training chemistry, top-3 accuracy is 0.352 (95% CI
0.283-0.428) against a permutation null of 0.145 (p=0.0005) and a frequency null of 0.654. The
output therefore depends on the compound and is not memorisation, but it does not beat a constant
answer naming the commonest CNS indications, and this is stated rather than omitted. The same
analysis identified a deployed endpoint, Nav1.1, that assigned binder probabilities between 0.80 and
0.82 to glucose, urea, acetate and glycine against a threshold of 0.796; it was withdrawn, as was
Cav3.2, leaving 47 of 49 binder endpoints deployed. Of six adversarial checks designed so that each
could fail, five pass; the sixth shows that the applicability-domain flag does not separate
non-drug-like chemistry from unseen drugs, and it is reported as a limitation rather than retuned.
The server returns an auditable, mechanistically interpretable brain-relevance profile rather than a
single opaque score, exportable as a tidy data table, a self-contained report, a structured record or
a vector figure, and supports batch screening of compound sets. Source code, trained models and the
raw source caches are available at https://github.com/krishna-g-999/brainsafe-ai.

---

## Introduction

Central-nervous-system drug discovery has a distinctive failure profile: a candidate can be potent
against its intended target yet never reach the brain, or reach it but carry an unacceptable safety
liability, or engage unintended targets that reshape its clinical profile. Pharmacokinetics and the
barrier itself, rather than target affinity, account for a large share of central attrition
(10), and the parameter that governs central action is not total brain concentration but
the unbound brain-to-plasma ratio (4). Answering "will this molecule affect the brain, how, and
is it safe" therefore means combining blood-brain-barrier (BBB) penetration, disease-relevant target
engagement, free-brain exposure, and safety into one readout. Property-based schemes such as CNS MPO
(11) score the first of these well but do not name a mechanism; single-endpoint QSAR models name
one mechanism but not exposure. Existing public tools typically address one of these axes in
isolation, and large language models, while fluent, fabricate measured identifiers and cannot supply
a calibrated, auditable answer.

BrainSafe AI is a web server that unifies these axes for any user-supplied structure. Its design
priorities are three. First, every machine-learning endpoint is trained only on *measured*
experimental data, never on qualitative annotation, so predictions are not circular. Second, every
prediction is reported with its uncertainty: a calibrated probability, an applicability-domain flag,
and the nearest measured analogue behind it, so the user always knows whether a number is
interpolation or extrapolation. Third, predictions are made mechanistically interpretable: individual
target engagements are traced through a curated target-to-pathway-to-disease knowledge graph into
per-disease relevance scores, filtered by predicted blood-brain-barrier exposure, so a result is an
explanation rather than a bare number.

## Materials and Methods

### Endpoint selection and scientific rationale

The endpoint panel is organised around the four sequential questions a CNS candidate must satisfy
(Figure 1). *Exposure*: a molecule that reaches no free concentration in brain tissue cannot act
centrally however potent it is, so BBB penetration, the unbound brain-to-plasma ratio (K_p,uu), total
brain distribution (logBB), P-glycoprotein efflux and passive permeability are modelled first.
*Target engagement*: the panel covers the cholinergic axis (AChE, BChE, alpha-7 nicotinic receptor),
where acetylcholinesterase inhibition remains the mainstay symptomatic treatment in Alzheimer's
disease (12), and the amyloid and tau axes (BACE1, GSK-3beta) (13); monoamine oxidase
B (14) and LRRK2 (15) for Parkinson's disease; the serotonergic, dopaminergic,
noradrenergic, opioid, cannabinoid, histaminergic, adenosine and sigma-1 systems that underlie
depression, anxiety, psychosis, addiction, attention deficit, chronic pain and sleep regulation, the
last including the orexin receptors (16). Three axes are included because they are
mechanistically implicated across several neurodegenerative conditions rather than tied to one:
NLRP3-driven neuroinflammation (17), KEAP1-NRF2 antioxidant signalling (18), and
histone deacetylase activity, whose genetic removal modifies pathology in Huntington's disease models
(19). Glutamatergic targets are included on the same basis, riluzole being the long-standing
example of an approved agent acting on that axis (20). *Safety*: hERG blockade is a
leading cause of late-stage cardiovascular attrition through QT prolongation (21) and is
modelled as an explicit liability. *Developability*: solubility, lipophilicity, plasma protein binding
and hepatocyte clearance determine whether an achievable dose sustains exposure, and a measured
antioxidant endpoint captures the oxidative-stress axis common to neurodegeneration.

In total the server comprises 62 endpoints (53 molecular-target endpoints and 9 ADME endpoints),
realised as 69 trained models because four receptor targets are represented both as potency
regressions and as binder classifiers, and because a pKa regression and an antioxidant regression
support the CNS-likeness and neuroprotection axes. Two further endpoints were trained and withdrawn
after deployment testing and are not counted here; the reasons are given in the Results.

### Training data

Protein-target activity is pooled from ChEMBL (1) (pChEMBL values) and BindingDB (2)
at the compound level; BBB penetration uses the B3DB database (3) augmented with FDA-curated
approved drugs; the antioxidant endpoint uses measured DPPH pIC50 values; and the nine ADME endpoints
use measured sets from Therapeutics Data Commons (22), MoleculeNet (23), B3DB and ChEMBL.
The full panel holds 227,146 measured compound-endpoint records over 193,536 unique compounds keyed
by the InChIKey of the desalted parent; the thirteen core target endpoints account for 76,850 of
those records. No value is imputed and no source overrides a measurement.

These totals are sums across endpoints and are not the size of any training set. Each endpoint is
trained and cross-validated on its own measured set alone, and those sets span nearly two orders of
magnitude, from 234 compounds at GluA2 to 10,276 at hERG. A compound measured at several targets
contributes one record to each and is counted once per endpoint. Per-endpoint compound counts,
scaffold counts and class balance are given in Tables 1 to 3, and the complete per-endpoint
accounting is in the Supplementary training record.

**The negative class.** A public bioactivity record describes what was found to bind. A compound
assayed and found inactive is frequently deposited only as a censored bound, `standard_relation` `>`
with no pChEMBL value, and the conventional query, which filters on pChEMBL, discards precisely those
rows. Training on what survives that filter yields a positive class drawn from measurement and a
negative class drawn from property-matched decoys, and it left 35 of the 60 endpoints above 90 per
cent active, which is a property of the query rather than of the chemistry. A censored bound settles
a label whenever the entire interval it defines falls on one side of the activity cut: `IC50 > 10 uM`
places the true potency strictly below pChEMBL 5.0 and is therefore a measured non-binder, whereas
`IC50 > 100 nM` spans both classes and is discarded as undecidable rather than guessed at (Figure 5).
Applying this recovered 21,994 measured non-binders across 57 endpoints and reduced the endpoints
above 90 per cent active from 35 to 13. Bounds are never converted to potencies: a bound is used only
to assign a class, never as a value in any regression.

### Methods compared and model selection

Molecules are represented by a 1,024-bit ECFP-4 fingerprint (5) plus twelve interpretable
physicochemical descriptors (molecular weight, cLogP, topological polar surface area, hydrogen-bond
donors and acceptors, rotatable bonds, aromatic rings, fraction sp3, ring count, heavy-atom count,
formal charge, and QED drug-likeness), giving the 1,036-column vector shown in Figure 2. Five model
families were benchmarked under identical 5-fold cross-validation on both split regimes, on the
deduplicated matrix the deployed pipeline fits, using scikit-learn (24). Two are baselines a
reader is entitled to demand: a five-nearest-neighbour read-across on Tanimoto similarity, which is
what a medicinal chemist does by eye and which any model must beat to justify itself, and
L2-regularised logistic regression. Three are ensembles: a random forest (6), XGBoost
(25), and histogram gradient boosting.

On the scaffold split the random forest gives the highest mean classifier AUROC, 0.9228, ahead of
histogram gradient boosting (0.9160), XGBoost (0.9144), the read-across (0.8844) and logistic
regression (0.8338). It exceeds the read-across on all eight classifier endpoints and all five
regression endpoints, and it is the single best classifier on seven of eight. The margin over the
boosted ensembles is small, 0.007 mean AUROC, and on regression the random forest is not the leader:
XGBoost and histogram gradient boosting reach mean scaffold R2 of 0.5453 and 0.5452 against 0.5186
for the random forest, an advantage of about 0.027 that holds on four of five endpoints.

The random forest is therefore deployed for every endpoint as a considered choice rather than as a
clean win. It leads decisively where the panel's principal claims are made, calibrates stably under
isotonic regression, supplies the out-of-bag vote distribution that the conformal layer consumes, is
interpretable through feature importance, and needs no GPU, which keeps the server lightweight. One
estimator across all endpoints also means a single well-characterised failure mode rather than five.
The cost of that uniformity is roughly 0.03 R2 on the potency regressions, and it is stated here
rather than left for a reader to discover.

A graph isomorphism network trained directly on molecular graphs was compared against the random
forest on four endpoints during an earlier build of the panel (2026-07-21) and did not exceed it on
any of them. That comparison has not been repeated on the current models, so it is reported as
indicative of the representation question at this data scale rather than as a current result.

### Cross-validation design

Every endpoint is evaluated by 10-fold cross-validation under two regimes (Figure 5a, 5b). In the
**random** regime folds are drawn at random, so a test compound will usually have close analogues in
the training folds; this measures interpolation within known chemistry. In the **scaffold-grouped**
regime folds are formed by GroupKFold over generic Bemis-Murcko scaffolds, so every compound sharing a
scaffold is assigned to the same fold and each fold holds out entire chemical series unseen during
training; this measures extrapolation to new chemotypes. Fold sizes are equal to within one compound in
both regimes, so the two are directly comparable. All reported means and standard deviations are across
the ten folds, and per-fold values for every endpoint are provided in the repository
(`results/tables/manuscript_T2_per_fold.csv`, `adme_cv_folds.csv` and `binder_cv_folds.csv`).

Cross-validation here fits models in order to measure, not in order to deploy. Each endpoint is fitted
twenty times during evaluation, ten in each regime, and each of those twenty models is scored once on
the fold withheld from it and then discarded. The model that is served is a twenty-first, refitted on
the endpoint's full set after the estimate is fixed, so no compound used to report a score was in the
training set of the model that scored it. The panel holds 71 cross-validated models over 67
endpoint names, four targets (D2, A2A, 5-HT2A and SERT) carrying both a potency regressor and a
binder classifier, so the evaluation comes to 71 x 2 x 10 = 1,420 fits standing behind the deployed
panel.

A third, more demanding regime is also reported. In the **temporal** split the model is trained only on
compounds published before a cutoff year and tested on compounds published after it (Table 5), which is
the closest available analogue of prospective use. Aggregate performance falls in this regime:
classifier AUROC averages 0.752 (range 0.611 to 0.908) against 0.919 under the scaffold split, and the
receptor potency regressions average R2 0.134. Because this figure determines how the tool should be
used, we investigated its cause rather than simply reporting it (next section).

### What limits prospective performance, and how the server handles it

We tested four candidate explanations for the temporal drop.

*Range restriction.* R2 is normalised by the variance of the test set, so a narrowing potency window
over time would depress R2 without any loss of accuracy. This is refuted: the variance of the future
test set is equal to or larger than that of the training set for every regression endpoint (variance
ratio 0.92 to 1.89). The low R2 is not a normalisation artifact.

*Reduced training data.* The temporal model is fitted on the pre-cutoff subset only. Refitting on a
random subset of the same size recovers R2 0.59 to 0.65 and rank correlation 0.76 to 0.80, so the loss
is caused by temporal distribution shift and not by sample size.

*Recency weighting.* Exponentially up-weighting recent training compounds changed classifier AUROC by a
mean of only +0.004 and did not help the regressions, so drift is not a simple recency effect.

*Applicability domain.* This is the explanation that holds. Stratifying future compounds by their
maximum Tanimoto similarity to the training set (Figure 7) separates the aggregate figure into two very
different regimes. For future compounds that fall **inside** the domain the classifiers retain a mean
AUROC of **0.827** and the potency models a mean rank correlation of **0.56**; in the intermediate band
these fall to 0.742 and 0.30; and **outside** the domain the classifiers approach chance (**0.573**) and
the potency models retain essentially no rank information (0.06). The aggregate temporal number is
therefore a mixture of a usable regime and an uninformative one, and the ratio between them is decided
by a quantity the server can compute at query time.

The server acts on this rather than merely disclosing it. Every result opens with a reliability
statement that places the query in one of the three regimes and quotes the measured prospective
performance for that regime, so a user is told, before reading any number, how much of it to believe.
This converts the applicability domain from a caveat into an operational safeguard, and it is the
reason the receptor targets are deployed as binder classifiers rather than as potency regressions.
Rank correlation is reported alongside R2 for the regression endpoints because ranking, not absolute
potency, is the decision-relevant quantity for triage.

### Why the scaffold error bars are larger

The standard deviation reported with each endpoint mixes two distinct sources. The first is *sampling
noise*: each fold metric is computed on a finite test set and carries its own standard error even if
every fold were equally difficult. The second is *fold heterogeneity*: genuine variation in how hard
each held-out fold is. For k folds with observed fold metrics m_i these combine as

    Var_between(m) = Var_heterogeneity + mean_i Var_sampling(m_i)

so the heterogeneity component is recovered by subtraction, with Var_sampling estimated by
bootstrapping within each fold's own test set (400 resamples). Applying this decomposition to every
endpoint (Figure 5c, Table 4) shows a clear separation. Under the random split a mean of 15% of the
between-fold variance is heterogeneity (median 8%, and indistinguishable from zero for six of thirteen
endpoints), confirming that folds are statistically exchangeable and the error bar is essentially pure
sampling noise. Under the scaffold split a mean of 71% of the variance is genuine heterogeneity (range
53% to 92%).

The wider scaffold error bar is therefore not model instability, and it is not a defect of the
protocol: it is a measurement of how much predictive performance depends on which chemical series is
being queried. Endpoints with the largest heterogeneity share (hERG 90%, BBB 92%, SERT 89%) are those
whose training data spans the most disconnected chemical space, and they are precisely the endpoints
for which an applicability-domain flag matters most. This result is the empirical justification for
reporting a per-endpoint domain flag with every prediction rather than a single global confidence.

### Decoy-aware binder classifiers

The receptor, transporter and kinase targets are reported in ChEMBL almost entirely as actives (82% to
99% of records), so a naive potency regressor learns to predict the training median for any input and
cannot discriminate a real binder from an arbitrary molecule. We therefore model each as a
binder-versus-decoy classifier: positives are measured binders (pChEMBL >= 7); negatives are
property-matched decoys sampled from a 158,890-compound background library with a maximum ECFP-4
Tanimoto below 0.35 to any positive, so the model must learn structure rather than a similarity
shortcut. Decoys are drawn only from the 95,515-compound decoy pool, one of three disjoint
partitions of that library (Figure 4A); the other two carry the threshold and the evaluation sets,
so no compound used to train a classifier can later be used to set or to test its threshold. Each classifier is a compact random forest with prefit sigmoid calibration (26), which is used
in place of isotonic regression here because the withheld set for a binder endpoint is often too
small to fit a step function without overfitting it. Decoy-based validation is reported here only for completeness; the figures the panel is judged on
are those against measured non-binders in the next section, because an AUROC against decoys measures
separation from chemistry chosen to be dissimilar and is inflated by that choice.

### Validation against experimentally measured inactives, and threshold calibration

Decoy-based validation is convenient but circular in an important sense: negatives are defined by
structural dissimilarity, so a high AUROC may only demonstrate that the model recognises the chemical
neighbourhood of a target's ligand series. We therefore tested every binder model against a negative
set it had never seen and that carries experimental authority: compounds measured against the same
target and found inactive (pChEMBL below 5). These compounds are present in our own endpoint tables
and are excluded from binder-model training by construction.

Two findings follow. First, the ranking is genuine. Across 34 targets with sufficient measured
inactives the mean AUROC of measured binders against held-out measured inactives is 0.956, no target falls
below 0.88, and the mean difference from the near-miss decoy estimate is only +0.028. The concern that
the models solve a merely lexical task is therefore not supported.

Second, the decision threshold was wrong. At the conventional 0.5 cut, a large fraction of
experimentally confirmed inactives is called a binder, reaching 81% for HDAC1 and exceeding 50% for
several other targets, because the property-matched decoys are easier to reject than genuine tested
non-binders and the calibration inherits that optimism. Ranking quality and threshold quality are
distinct, and only the former is captured by AUROC. We therefore set each target's operating point on
its measured inactives, choosing the threshold at which 10% of them would be called binders. Across
the panel this yields a mean false-positive rate of 0.11 and a mean sensitivity of 0.88 (Table 2).
Two targets, KEAP1 and melatonin MT1, initially retained a sensitivity below 0.6 at that operating
point. Adding the measured inactives to training as hard negatives raised KEAP1 from 0.33 to 0.83 and
MT1 from 0.52 to 0.78, and no target now falls below the reliability threshold. The server retains the
machinery to mark a low-sensitivity target as possibly under-called should future data reintroduce one.

### Calibration, applicability domain and exposure

Every classifier is isotonically calibrated (8) on out-of-fold predictions, so the
calibrator never sees a compound in its own fit; mean expected calibration error falls from 0.0795
to 0.0161 (Figure 6A). Each prediction additionally carries a Mondrian conformal interval, which
turns the applicability domain from a caveat into a coverage statement (27); empirical
coverage is 0.887 to 0.920 against a 0.90 target. For every prediction the server also computes the
maximum ECFP-4 Tanimoto similarity of the query to each endpoint's own measured chemistry, a
nearest-neighbour definition of the applicability domain (28), and reports an endpoint-specific
in-domain, near-domain or out-of-domain flag together with the nearest measured analogue and its
structure. The ADME layer
provides a free-brain-exposure verdict; on known drugs this correctly separates central compounds
(diazepam K_p,uu 0.94, donepezil 0.84) from peripheral or efflux-limited ones (atenolol 0.07,
loperamide 0.04).

### Base-rate enrichment and the knowledge graph

Because several classifier training sets are active-heavy (BACE1 91%, GSK-3beta 93% positive), a raw
calibrated probability is not evidence of engagement unless it exceeds the base rate. The server
therefore scores each target by its *enrichment* over the endpoint base rate, and each receptor by its
floored binder probability, so that a prediction near chance contributes nothing. A curated, versioned
knowledge graph maps each target through a biological pathway to the diseases it informs, anchored to
KEGG synapse and disease maps (hsa04725, hsa04726, hsa04728, hsa04723, hsa05032, hsa04080, hsa05010,
hsa05012) (29), the Reactome KEAP1-NFE2L2 oxidative-stress response (R-HSA-9755511) (30)
and IUPHAR/BPS Guide to Pharmacology associations (31). A per-disease relevance score is the strongest engaged target for that
disease, scaled by predicted BBB penetration; taking the strongest rather than an average prevents
unrelated mechanisms from diluting a real signal. Coverage spans fourteen brain conditions (Figure 1).

Two properties of this construction were established by ablation rather than assumed, and both
constrain how it should be described. First, the BBB term multiplies every disease identically and
therefore cannot alter which condition ranks highest; it is an exposure filter that determines
whether any signal is reported, not a discriminative term. Second, the curated edge weights do not
measurably affect the disease ranking: over 15,609 scaffold-held-out compounds, curated, uniform and
randomly permuted weights give top-3 accuracies of 0.7917, 0.7911 and 0.7899. The predictive content
of the graph lies in its topology, in which target is connected to which disease, and the weights are
retained as a mechanistic prior expressing directness of linkage rather than as fitted parameters.

### Web implementation

The server is a single-page application built with Streamlit (32) and RDKit (33). The
user submits either a compound name, resolved to a structure through PubChem, or a SMILES string.
Models are loaded once and cached; a complete profile across all 72 deployed estimators, including
both applicability-domain calculations against the 158,890-compound reference library, returns in a
few seconds on a single CPU core. The interface is self-contained and requires no installation or
GPU.

## Results

### Per-endpoint performance

Tables 1 to 3 report, for every endpoint, the measured data size, the number of distinct scaffolds,
class balance where applicable, the training and test set size per fold, and the mean and standard
deviation of the appropriate metric under both cross-validation regimes. Figure 6 summarises the same
results graphically.

The eight target classifiers reach a mean AUROC of 0.958 (random) and 0.925 (scaffold). BACE1 is the
most robust to chemotype change, losing 0.012 between the two splits (0.978 to 0.965), and MAO-A the
least, losing 0.065 (0.964 to 0.899). An external test of the BBB model on 306 FDA-curated approved
drugs absent from B3DB by InChIKey gives AUROC 0.761; restricted to the 241 of those that are also
distinguishable from the training set in feature space, which is the subset that supports an external
claim, it gives 0.788. The receptor potency regressions reach R2 0.64 to 0.72 (random) and 0.46 to
0.61 (scaffold). The 49 binder classifiers reach a mean AUROC of 0.904 against measured non-binders,
reported in place of any decoy-based figure. Within the ADME layer, performance ranges from strong
(solubility R2 0.73, P-gp inhibition AUROC 0.94 under the scaffold split) to weak and explicitly
disclosed (hepatocyte clearance R2 0.21, K_p,uu R2 0.35).

## Table 1. Core target panel (13 endpoints), 10-fold cross-validation

| Endpoint | Target | Task | Compounds | Scaffolds | Active fraction | Train/fold | Test/fold | Metric | Random 10-fold | Scaffold 10-fold | Why this endpoint |
|---|---|---|---|---|---|---|---|---|---|---|---|
| BBB | Blood-brain barrier penetration | classification | 7,805 | 2,439 | 0.64 | 7,024 | 780 | AUROC | 0.961 ± 0.007 | 0.920 ± 0.037 | Gate: no CNS effect without brain entry |
| AChE | Acetylcholinesterase | classification | 4,387 | 2,116 | 0.72 | 3,948 | 439 | AUROC | 0.963 ± 0.008 | 0.921 ± 0.021 | Symptomatic Alzheimer therapy (donepezil class) |
| BChE | Butyrylcholinesterase | classification | 2,621 | 1,345 | 0.70 | 2,359 | 262 | AUROC | 0.968 ± 0.012 | 0.937 ± 0.021 | Rises as Alzheimer progresses; selective-inhibitor target |
| BACE1 | Beta-secretase 1 | classification | 8,501 | 3,088 | 0.91 | 7,651 | 850 | AUROC | 0.967 ± 0.010 | 0.956 ± 0.021 | Rate-limiting step of amyloid-beta generation |
| GSK3B | Glycogen synthase kinase-3 beta | classification | 4,958 | 2,118 | 0.93 | 4,462 | 496 | AUROC | 0.969 ± 0.013 | 0.937 ± 0.030 | Tau hyperphosphorylation; neuroprotection |
| MAO_A | Monoamine oxidase A | classification | 2,228 | 827 | 0.41 | 2,005 | 223 | AUROC | 0.947 ± 0.017 | 0.868 ± 0.046 | Serotonin/noradrenaline catabolism; depression |
| MAO_B | Monoamine oxidase B | classification | 3,665 | 1,196 | 0.68 | 3,298 | 366 | AUROC | 0.955 ± 0.007 | 0.889 ± 0.033 | Dopamine catabolism; Parkinson therapy (selegiline) |
| hERG | hERG potassium channel | classification | 5,875 | 3,377 | 0.41 | 5,287 | 587 | AUROC | 0.954 ± 0.007 | 0.921 ± 0.035 | Cardiotoxicity liability; principal safety filter |
| D2 | Dopamine D2 receptor | regression | 7,734 | 3,329 | n/a | 6,961 | 773 | R2 | 0.601 ± 0.018 | 0.483 ± 0.051 | Antipsychotic efficacy; motor control |
| A2A | Adenosine A2A receptor | regression | 6,785 | 2,654 | n/a | 6,106 | 678 | R2 | 0.682 ± 0.023 | 0.576 ± 0.066 | Non-dopaminergic Parkinson target (istradefylline) |
| HT2A | Serotonin 5-HT2A receptor | regression | 5,989 | 2,254 | n/a | 5,390 | 599 | R2 | 0.636 ± 0.028 | 0.490 ± 0.054 | Atypical antipsychotics; psychedelics |
| SERT | Serotonin transporter | regression | 4,572 | 1,459 | n/a | 4,115 | 457 | R2 | 0.602 ± 0.038 | 0.389 ± 0.124 | SSRI antidepressant target |
| antioxidant_DPPH | Radical-scavenging capacity | regression | 2,862 | 919 | n/a | 2,576 | 286 | R2 | 0.669 ± 0.065 | 0.434 ± 0.100 | Oxidative stress in neurodegeneration |


## Table 2. Binder classifiers validated against held-out measured inactives

| Endpoint | Target | Measured binders | Training negatives | AUROC (held-out measured inactives) | Threshold | Sensitivity | Why this endpoint |
|---|---|---|---|---|---|---|---|
| D2 | Dopamine D2 receptor | 3,676 | 10,895 | 0.938 | 0.986 | 0.785 | Antipsychotic efficacy |
| A2A | Adenosine A2A receptor | 4,352 | 12,941 | 0.960 | 0.994 | 0.861 | Non-dopaminergic Parkinson target |
| HT2A | Serotonin 5-HT2A receptor | 3,951 | 11,813 | 0.953 | 0.983 | 0.850 | Atypical antipsychotic profile |
| SERT | Serotonin transporter | 3,028 | 9,012 | 0.983 | 0.768 | 0.970 | SSRI antidepressant target |
| HT1A | 5-HT1A receptor | 3,703 | 11,012 | 0.989 | 0.758 | 0.971 | Anxiety, depression (buspirone) |
| HT6 | 5-HT6 receptor | 2,741 | 8,195 | 0.977 | 0.947 | 0.942 | Cognition enhancement in Alzheimer |
| HT7 | 5-HT7 receptor | 1,478 | 4,377 | 0.960 | 0.921 | 0.913 | Mood, circadian rhythm, sleep |
| H3 | Histamine H3 receptor | 3,212 | 9,548 | 0.990 | 0.919 | 0.964 | Wakefulness (pitolisant), cognition |
| DAT | Dopamine transporter | 1,337 | 3,855 | 0.986 | 0.531 | 0.984 | ADHD, addiction, stimulant liability |
| NET | Noradrenaline transporter | 1,572 | 4,583 | 0.990 | 0.715 | 0.971 | Depression, ADHD (atomoxetine) |
| Sigma1 | Sigma-1 receptor | 1,928 | 5,762 | 0.940 | 0.955 | 0.856 | Neuroprotection, ER-stress chaperone |
| CB1 | Cannabinoid CB1 receptor | 2,680 | 7,936 | 0.965 | 0.988 | 0.878 | Pain, appetite, mood |
| OPRK1 | Kappa-opioid receptor | 3,059 | 9,022 | 0.968 | 0.991 | 0.903 | Analgesia, dysphoria, mood |
| OPRM1 | Mu-opioid receptor | 3,718 | 10,971 | 0.981 | 0.874 | 0.969 | Analgesia, addiction liability |
| D3 | Dopamine D3 receptor | 4,127 | 12,289 | 0.978 | 0.960 | 0.938 | Addiction, Parkinson motor complications |
| A1 | Adenosine A1 receptor | 1,943 | 5,720 | 0.954 | 0.978 | 0.847 | Neuroprotection, epilepsy, sedation |
| a7nAChR | Alpha-7 nicotinic receptor | 337 | 921 | 0.982 | 0.642 | 0.932 | Cognition, neuroinflammation in Alzheimer |
| LRRK2 | LRRK2 kinase | 1,173 | 3,510 | n/a | 0.400 | n/a | Most common genetic cause of Parkinson disease |


## Table 3. ADME / exposure layer (9 endpoints)

| Endpoint | Property | Task | Compounds | Metric | Random 10-fold | Scaffold 10-fold | Why this endpoint |
|---|---|---|---|---|---|---|---|
| kpuu | Unbound brain/plasma ratio (Kp,uu) | regression | 566 | R2 | 0.404 ± 0.099 | 0.352 ± 0.158 | Free drug available to CNS targets |
| logbb | Total brain/plasma ratio (logBB) | regression | 1,058 | R2 | 0.577 ± 0.081 | 0.423 ± 0.150 | Bulk brain distribution |
| caco2_permeability | Caco-2 permeability | regression | 897 | R2 | 0.734 ± 0.051 | 0.584 ± 0.158 | Passive membrane permeability |
| pgp_substrate | P-glycoprotein substrate | classification | 1,371 | AUROC | 0.858 ± 0.032 | 0.806 ± 0.038 | Active efflux out of the brain |
| pgp_inhibition | P-glycoprotein inhibition | classification | 1,212 | AUROC | 0.955 ± 0.018 | 0.935 ± 0.029 | Efflux-mediated drug interactions |
| solubility | Aqueous solubility (logS) | regression | 9,573 | R2 | 0.804 ± 0.017 | 0.728 ± 0.077 | Formulation and absorption |
| lipophilicity | Lipophilicity (logD) | regression | 4,200 | R2 | 0.639 ± 0.028 | 0.566 ± 0.042 | Permeability/promiscuity balance |
| plasma_protein_binding | Plasma protein binding | regression | 1,797 | R2 | 0.434 ± 0.090 | 0.364 ± 0.056 | Determines free fraction |
| clearance_hepatocyte | Hepatocyte clearance | regression | 1,020 | R2 | 0.230 ± 0.104 | 0.206 ± 0.064 | Metabolic stability, exposure duration |


## Table 4. Between-fold error-bar decomposition

Observed between-fold SD separated into sampling noise (finite test set) and genuine fold-to-fold heterogeneity, by within-fold bootstrap.


| Endpoint | Random SD | of which sampling | heterogeneity share | Scaffold SD | of which sampling | heterogeneity share |
|---|---|---|---|---|---|---|
| BBB | 0.0075 | 0.0062 | 32% | 0.0367 | 0.0102 | 92% |
| AChE | 0.0083 | 0.0094 | 0% | 0.0206 | 0.0137 | 56% |
| BChE | 0.0118 | 0.0096 | 34% | 0.0212 | 0.0145 | 53% |
| BACE1 | 0.0099 | 0.0106 | 0% | 0.0209 | 0.0118 | 68% |
| GSK3B | 0.0128 | 0.0112 | 23% | 0.0302 | 0.0204 | 54% |
| MAO_A | 0.0169 | 0.0147 | 24% | 0.0464 | 0.0251 | 71% |
| MAO_B | 0.0075 | 0.0103 | 0% | 0.0330 | 0.0185 | 69% |
| hERG | 0.0066 | 0.0082 | 0% | 0.0351 | 0.0109 | 90% |
| D2 | 0.0177 | 0.0255 | 0% | 0.0515 | 0.0299 | 66% |
| A2A | 0.0228 | 0.0251 | 0% | 0.0662 | 0.0289 | 81% |
| HT2A | 0.0276 | 0.0264 | 8% | 0.0537 | 0.0311 | 67% |
| SERT | 0.0379 | 0.0329 | 25% | 0.1239 | 0.0404 | 89% |
| antioxidant_DPPH | 0.0648 | 0.0457 | 50% | 0.0995 | 0.0592 | 65% |


## Table 5. Temporal (future-compound) validation

Models are trained only on compounds published before the cutoff year and tested on compounds published after it. This is the most demanding regime and the closest analogue of prospective use.


| Endpoint | Target | Cutoff year | Train | Test | Metric | Score |
|---|---|---|---|---|---|---|
| AChE | Acetylcholinesterase | 2020 | 4,141 | 1,028 | AUROC | 0.760 |
| BChE | Butyrylcholinesterase | 2021 | 2,696 | 617 | AUROC | 0.801 |
| BACE1 | Beta-secretase 1 | 2017 | 6,901 | 1,673 | AUROC | 0.910 |
| GSK3B | Glycogen synthase kinase-3 beta | 2021 | 4,205 | 546 | AUROC | 0.776 |
| MAO_A | Monoamine oxidase A | 2020 | 2,859 | 819 | AUROC | 0.727 |
| MAO_B | Monoamine oxidase B | 2020 | 3,719 | 923 | AUROC | 0.843 |
| hERG | hERG potassium channel | 2020 | 8,158 | 2,059 | AUROC | 0.715 |
| D2 | Dopamine D2 receptor | 2019 | 6,796 | 1,556 | R2 | 0.080 |
| A2A | Adenosine A2A receptor | 2020 | 4,600 | 1,516 | R2 | 0.226 |
| HT2A | Serotonin 5-HT2A receptor | 2020 | 4,741 | 986 | R2 | 0.157 |
| SERT | Serotonin transporter | 2015 | 3,647 | 1,147 | R2 | 0.111 |
| antioxidant_DPPH | Radical-scavenging capacity | 2016 | 2,340 | 522 | R2 | 0.009 |


Classifier endpoints: mean AUROC 0.790 (range 0.715 to 0.910). Regression endpoints: mean R2 0.117 (range 0.009 to 0.226).


## Table 6. Prospective validation under a scaffold hold-out

Twenty per cent of Bemis-Murcko scaffolds were withheld per target and every model retrained on the remainder, so no held-out compound shares a scaffold with anything its model saw. Thresholds were recalibrated on held-out negatives and an independent background sample. Targets marked excluded produced a threshold at the permitted floor, meaning no separation from background chemistry, and do not contribute to the pooled estimate.


| Target | Train actives | Held-out actives | Held-out scaffolds | Threshold | Recall | 95% CI | Note |
|---|---|---|---|---|---|---|---|
| OX2 | 3,228 | 624 | 233 | 0.050 | 0.990 | [0.979, 0.996] | excluded |
| LRRK2 | 752 | 421 | 90 | 0.050 | 0.957 | [0.933, 0.973] | excluded |
| PDE4B | 954 | 174 | 99 | 0.262 | 0.954 | [0.912, 0.977] |  |
| mTOR | 2,429 | 555 | 216 | 0.494 | 0.950 | [0.928, 0.965] |  |
| CSF1R | 1,510 | 274 | 152 | 0.190 | 0.927 | [0.890, 0.952] |  |
| H3 | 2,615 | 597 | 308 | 0.758 | 0.926 | [0.903, 0.945] |  |
| GABA_A | 166 | 35 | 16 | 0.098 | 0.914 | [0.776, 0.970] |  |
| DAT | 1,130 | 207 | 93 | 0.315 | 0.903 | [0.855, 0.937] |  |
| PDE10A | 3,141 | 943 | 326 | 0.979 | 0.902 | [0.882, 0.920] |  |
| SERT | 2,423 | 548 | 185 | 0.560 | 0.901 | [0.874, 0.924] |  |
| Nav1_5 | 194 | 48 | 27 | 0.367 | 0.875 | [0.753, 0.941] |  |
| NLRP3 | 176 | 46 | 27 | 0.953 | 0.870 | [0.743, 0.939] |  |
| OPRM1 | 2,931 | 787 | 302 | 0.842 | 0.859 | [0.833, 0.882] |  |
| HT6 | 2,023 | 718 | 167 | 0.961 | 0.858 | [0.830, 0.882] |  |
| a7nAChR | 266 | 71 | 35 | 0.494 | 0.845 | [0.743, 0.911] |  |
| Nav1_7 | 2,172 | 572 | 179 | 0.957 | 0.844 | [0.812, 0.872] |  |
| D3 | 3,283 | 844 | 351 | 0.965 | 0.842 | [0.816, 0.865] |  |
| HDAC1 | 2,421 | 634 | 292 | 0.905 | 0.841 | [0.810, 0.867] |  |
| OX1 | 2,736 | 576 | 178 | 0.983 | 0.833 | [0.801, 0.862] |  |
| KEAP1 | 105 | 17 | 15 | 0.430 | 0.824 | [0.590, 0.938] |  |
| NET | 1,307 | 265 | 90 | 0.646 | 0.819 | [0.768, 0.861] |  |
| HT1A | 2,902 | 801 | 327 | 0.980 | 0.814 | [0.786, 0.839] |  |
| A2A | 3,353 | 748 | 337 | 0.993 | 0.813 | [0.783, 0.839] |  |
| OPRK1 | 2,575 | 484 | 234 | 0.985 | 0.800 | [0.762, 0.833] |  |
| CB1 | 2,161 | 519 | 196 | 0.969 | 0.775 | [0.737, 0.808] |  |
| A1 | 1,535 | 408 | 177 | 0.980 | 0.762 | [0.719, 0.801] |  |
| HDAC6 | 2,928 | 787 | 335 | 0.995 | 0.729 | [0.697, 0.759] |  |
| HT7 | 1,166 | 312 | 133 | 0.954 | 0.721 | [0.669, 0.768] |  |
| Sigma1 | 1,571 | 357 | 186 | 0.961 | 0.675 | [0.625, 0.722] |  |
| MT1 | 571 | 126 | 48 | 0.971 | 0.667 | [0.581, 0.743] |  |
| SIRT1 | 141 | 21 | 17 | 0.766 | 0.667 | [0.454, 0.828] |  |
| HT2A | 3,102 | 715 | 306 | 0.996 | 0.639 | [0.603, 0.674] |  |
| GluN2B | 742 | 162 | 35 | 0.986 | 0.617 | [0.541, 0.689] |  |
| D2 | 2,989 | 668 | 340 | 0.992 | 0.554 | [0.516, 0.591] |  |
| TAAR1 | 78 | 15 | 9 | 0.397 | 0.533 | [0.301, 0.752] |  |
| P2X7 | 2,450 | 488 | 139 | 0.999 | 0.531 | [0.486, 0.575] |  |
| GBA1 | 101 | 19 | 10 | 0.453 | 0.526 | [0.317, 0.727] |  |
| COX2 | 945 | 237 | 67 | 0.951 | 0.371 | [0.312, 0.434] |  |
| mGluR5 | 986 | 214 | 106 | 0.999 | 0.308 | [0.250, 0.373] |  |
| GluA2 | 78 | 19 | 9 | 0.695 | 0.105 | [0.029, 0.314] |  |


Pooled recall 11,794/15,011 = 0.786; median per-target 0.817; 22 of 38 targets at or above 0.80.


## Table 7. Falsification analysis: each claim paired with a null model

Every hypothesis was stated so that it could fail. Where predictive power was at issue, scoring used scaffold hold-out models that never saw the compounds they scored. A refuted hypothesis is reported as prominently as a confirmed one; the purpose of the exercise was to find failures.


| Hypothesis | Verdict | Evidence |
|---|---|---|
| H1 the disease score is informative | **SUPPORTED** | top-3 accuracy 0.769 vs permutation null 0.152 (p=0.005) and frequency null 0.560 |
| H2 the curated edge weights add value | **REFUTED** | curated 0.7691, uniform 0.7678, permuted 0.7670 |
| H3 BBB gating discriminates between diseases | **REFUTED (by construction)** | the gate multiplies every disease equally and cannot change their order |
| H4 specificity transfers to novel chemistry | **SUPPORTED** | false-positive rate 0.033 on 61 distant compounds against 0.080 measured on library chemistry |
| H5 read-across beats a frequency baseline | **SUPPORTED** | recall 0.973 against 0.060 |
| H6 the disease scores match real clinical indications | **WEAKENED** | top-3 accuracy 0.352 on 162 drugs never seen in training, against permutation null 0.145 (p=0.001) and frequency null 0.654 |
| H7 some panel targets are non-discriminative and explain the silent antiepileptics | **REFUTED** | none of 35 targets ranks below AUROC 0.70; the cause is the operating point, with median deployed sensitivity 0.83 and 3 targets under 0.50 |
| H8 engaged targets are independent observations | **REFUTED** | 38 targets fire across approved drugs but span only 15 independent directions; 5 homologous pairs correlate above 0.5 |


## Table 8. Recovery of approved clinical indications, by condition

Ground truth is ChEMBL's drug_indication table restricted to phase 4, mapped to the panel through a keyword list fixed before any prediction was computed. 'Ranking only' removes the reporting threshold, separating the ability to rank conditions from the decision to stay silent. 'Silent' counts drugs for which no condition reached the reporting threshold.


| Indication | Drugs | Recovered in top 3 | Ranking only | Silent |
|---|---|---|---|---|
| Alzheimer's disease | 10 | 0.200 | 0.800 | 5 |
| Parkinson's disease | 39 | 0.154 | 0.538 | 16 |
| Depression / anxiety | 113 | 0.611 | 0.779 | 30 |
| Psychosis / schizophrenia | 73 | 0.644 | 0.740 | 14 |
| Addiction | 17 | 0.529 | 0.588 | 5 |
| ADHD | 31 | 0.226 | 0.258 | 16 |
| Chronic pain | 144 | 0.229 | 0.250 | 83 |
| Sleep / wakefulness | 27 | 0.296 | 0.296 | 11 |
| Epilepsy | 58 | 0.103 | 0.224 | 47 |


## Table 9. Targeted expansion: candidates, audit and outcome

Candidates were selected because the clinical-indication test identified epilepsy and chronic pain as the conditions served worst, and a systematic ChEMBL query found these mechanisms clear both a volume bar of 800 measured activities and a source-diversity bar. Each was then audited on its own fetched data before training, and again on deployed specificity afterwards.


| Candidate | Actives | Scaffolds | Measured inactives | AUROC vs inactives | Sensitivity | Outcome | Reason if rejected |
|---|---|---|---|---|---|---|---|
| a4b2nAChR | 796 | 278 | 146 | 0.937 | 0.899 | deployed |  |
| a3b4nAChR | 398 | 153 | 187 | 0.973 | 0.955 | deployed |  |
| Nav1_6 | 681 | 187 | 45 | 0.634 | 0.565 | deployed |  |
| Nav1_8 | 501 | 163 | 39 | 0.876 | 0.728 | deployed |  |
| Cav3_2 | 633 | 201 | 33 | 0.971 | 0.907 | withdrawn after training | active band compressed near zero, calibrated threshold 0.065, atenolol scores 0.084 |
| GABAA_a5 | 672 | 284 | 4 | n/a | n/a | not trained | only 4 measured inactives, cannot set a threshold honestly |
| CGRP | 761 | 333 | 26 | 0.985 | 0.997 | deployed | only 26 measured inactives, cannot set a threshold honestly |
| DHODH | 1421 | 360 | 145 | 0.964 | 0.967 | deployed |  |
| RIPK1 | 2349 | 826 | 719 | 0.983 | 0.969 | deployed |  |

### Prospective validation: sensitivity under a scaffold hold-out

Cross-validated performance answers how well a model interpolates within the chemistry it has seen.
The question that matters for a triage tool is how it behaves on a chemical series it has never
seen. We therefore withheld 20% of Bemis-Murcko scaffolds for each target, retrained all 39
trainable binder models on the remaining scaffolds, and recalibrated each decision threshold using
only held-out negatives and an independent background sample, so no held-out compound influenced
either the model or the threshold that scores it. The evaluation set comprises **16,874 held-out
compounds across 6,435 withheld scaffolds**.

Pooled recall is **0.790 (95% CI 0.783 to 0.796)**; the median per-target recall is 0.807 and 19 of
36 targets reach 0.80 or better (Figure 8, Table 6). Recall varies substantially and informatively
across the panel: mTOR (0.977), histamine H3 (0.926), PDE4B (0.923) and the dopamine transporter
(0.903) generalise well to unseen scaffolds, whereas five targets fall below 0.50 (SIRT1, mGluR5,
melatonin MT1, KEAP1 and AMPA GluA2), all of them endpoints with small measured-active sets.

Three targets, orexin OX2, LRRK2 and NLRP3, produced thresholds that collapsed to the permitted
floor, meaning the model could not separate its actives from background chemistry at any cut. Their
nominally high recall is an artefact of a threshold near zero rather than evidence of
discrimination, and they are excluded from the pooled estimate and flagged in Table 6 rather than
being allowed to flatter the headline figure.

### Prospective validation: specificity on 1000 non-CNS compounds

A tool that fires indiscriminately can achieve high recall while being useless, so specificity was
measured on a scale sufficient to bound it. We sampled 1000 compounds from the measured library that
carry no recorded activity at any modelled target, and scored them through the deployed pipeline.
**875 of 1000 returned no actionable disease signal, a specificity of 0.875 (95% CI 0.853 to 0.894)**
and a false-positive rate of 12.5% (Figure 9). An independently drawn replicate sample gave 12.3%,
so the estimate is stable. False positives concentrate in the neurodegenerative classes, with median
scores of 0.35 to 0.40, only just above the 0.30 actionable threshold, and 80 of the 125 fired on a
single condition rather than producing a diffuse profile.

Two caveats bound the interpretation. These compounds are presumed inactive because no activity is
recorded, not proven inactive by measurement, so 0.875 is a lower bound on the true specificity.
Second, all sampled compounds fall inside the applicability domain because they are drawn from the
training library, so this test does not probe specificity on distant chemistry.

Taken together the two tests place the deployed operating point at a sensitivity of 0.790 and a
specificity of 0.875, a balanced accuracy of 0.832 (Figure 10).

### Validation on reference compounds

Across a panel of well-characterised CNS drugs the top predicted condition matched the known
indication for donepezil (Alzheimer's disease), selegiline (Parkinson's disease), fluoxetine
(depression), haloperidol (psychosis), morphine (chronic pain), methylphenidate (addiction and
attention deficit) and resveratrol (neuroprotection). Two acetylcholinesterase inhibitors of atypical
chemotype, rivastigmine and galantamine, were under-called to low magnitude rather than mis-assigned,
and were correctly flagged as low-confidence by the applicability-domain layer, which illustrates the
value of reporting uncertainty alongside every prediction. The decoy-aware receptor models recover
expected selectivity, for example haloperidol as a D2 and 5-HT2A binder and fluoxetine as a selective
serotonin-transporter binder, while non-CNS controls return no distinctive engagement.

### Falsification analysis

Every result above was produced by investigators who wanted the tool to work. To counter that, each
central claim was restated so that it could fail and paired with a null model able to produce the
same apparent success by accident. Where predictive power was at issue, scoring used scaffold
hold-out models that never saw the compounds they scored. The analysis is versioned separately from
the validation results so that a refutation can never be mistaken for a confirmation (Table 7).

**The disease layer carries information.** Across held-out compounds the correct condition appears in
the top three predictions for 80.5% of cases, against 18.8% for a null that shuffles which disease
each target maps to (p = 0.005) and 48.2% for always answering with the three commonest conditions.

**The curated edge weights do not.** Recomputing the same predictions with uniform weights gives
0.8036 and with weights randomly permuted across edges 0.8025, against 0.8045 for the curated values,
a spread of 0.002. The ablation was repeated after the graph gained two conditions and three targets,
with new curated weights written for each, and the conclusion did not move. The information lies in which target connects to which condition, not in how
strongly. The weights are therefore described in this manuscript as a mechanistic prior and as graph
structure, not as tuned parameters, and the graph would be no less accurate without them.

**Blood-brain-barrier gating is a filter, not a discriminator.** The gate multiplies every condition
by the same probability and therefore cannot change which condition ranks first. It determines
whether anything is reported at all. Earlier wording implying that it sharpens the disease call
overstated its role and has been corrected throughout.

**Specificity transfers to distant chemistry.** On compounds drawn by random PubChem identifier,
independent of every set used to build the tool, the false-positive rate is 0.051 among those most
distant from training chemistry, against 0.125 measured on library compounds. The suspicion that the
headline specificity was an artefact of range restriction is refuted.

**Read-across is validated only in its intended regime.** It recovers the true target for 97.0% of
held-out compounds against 6.0% for a frequency baseline, with the query and any identical structure
excluded. This does not show that read-across works for a target class the index does not contain,
and the figure is not quoted as if it did.

### Validation against clinical indications

The preceding test asks whether the disease layer recovers the condition a compound's target maps to
under this project's own map, which establishes internal consistency rather than external truth. A
stronger test uses ChEMBL's drug_indication table restricted to phase 4, mapped to the panel through
a keyword list fixed before any prediction was computed and deliberately narrow, so that an unmatched
heading is discarded rather than coerced. Auditing that mapping rather than trusting it caught
Wolff-Parkinson-White syndrome, a cardiac condition, matching the substring "parkinson"; it was
excluded.

On 467 drugs the top-3 accuracy is 0.379. On the 162 whose exact structure appears nowhere in the
training chemistry it is 0.352, statistically indistinguishable from the 0.387 achieved on structures
the models have seen, which establishes that the result is prediction rather than memorisation. Both
exceed a permutation null of 0.145 decisively (p = 0.001). Neither exceeds a frequency null of 0.654,
because chronic pain, depression and psychosis account for most approved central-nervous-system
indications and a constant answer naming those three is right about two-thirds of the time. That
constant answer carries no information about any individual compound, but it is a real bar and the
tool does not clear it on this metric. Removing the reporting threshold and judging the ranking alone
raises accuracy to 0.497, which locates much of the gap in the decision to stay silent rather than in
the ranking itself. Per-indication recovery is reported in full (Table 8) because the aggregate
conceals a wide spread, from 0.644 for psychosis to 0.103 for epilepsy.

### Deployed specificity and the withdrawal of an endpoint

Two gaps in the evidence above motivated a further test. The background-specificity calibration bounds
the false-positive rate on a random sample of the measured training library, which is drug-like
medicinal chemistry and therefore a weak negative set; a model can pass it and still fire on a sugar.
The scaffold hold-out analysis used genuinely random structures but scored the hold-out twins rather
than the models the server delivers. Every deployed binder model was therefore scored at its deployed
threshold against 600 random PubChem structures and against molecules no central-nervous-system
target plausibly binds.

The median false-positive rate across the deployed panel is 0.0017. Four endpoints failed. Nav1.1
assigned binder probabilities of 0.806 to glucose, 0.809 to urea, 0.802 to acetate and 0.816 to
glycine against a threshold of 0.796, at an overall false-positive rate of 0.098 on random chemistry.
The failure is not one of calibration: every trivial control lies within a band of 0.015, and the
threshold that would restore 5% specificity still calls glycine a binder. This compressed-probability
pathology also explains that endpoint's deployed sensitivity of 0.32 and its provenance, 78% of
measurements from a single assay. Nav1.1 was withdrawn from the panel. Cav3.2 was rejected before
deployment for the same pathology inverted, its calibrated threshold falling to 0.065 while atenolol
scored 0.084. SIRT1 and alpha3beta4 were re-thresholded at negligible cost to sensitivity. After
these decisions every deployed endpoint holds below 5% on random chemistry and none fires on a
trivial molecule.

### Targeted expansion, and a limit that more targets cannot remove

Because the clinical-indication test identified epilepsy and chronic pain as the conditions served
worst, ChEMBL was queried systematically for the mechanisms those conditions require. Six cleared
both a volume bar of 800 measured activities and a source-diversity bar; four survived audit and
training and are deployed: alpha4beta2 and alpha3beta4 nicotinic receptors, Nav1.6 and Nav1.8
(Table 9). Nicotine and varenicline, previously invisible to the tool, now register their known
nicotinic engagement, and recovery of approved indications for addiction rose from 0.412 to 0.529.

Epilepsy did not improve, and the reason is instructive. In ChEMBL's own measurements carbamazepine
scores pChEMBL 4.49 at Nav1.7 and lamotrigine 4.77, both labelled inactive; mexiletine reaches 4.93
and lacosamide 6.74, below the binder threshold of 7. The classic antiepileptic drugs are not
high-affinity sodium-channel ligands; they are use-dependent, state-dependent blockers acting at tens
of micromolar. Silence on carbamazepine is therefore the correct response to the measured data, and
no binding-affinity endpoint will recognise that pharmacological class. The epilepsy gap is a
mismatch between what the panel measures and how a class of drugs works, not a missing target, and it
is reported as such rather than as a coverage deficiency that further expansion would close.

A related caution emerged from calibrating the new endpoints. Nav1.6 ranks its held-out actives
against random chemistry at an AUROC of 0.997 with a false-positive rate of 0.000, which appears to
license relaxing its threshold of 0.9605 and would raise sensitivity from 0.591 to 0.989. Tested
against hard negatives the proposal fails: compounds assayed against Nav1.6 and found inactive carry
a median probability of 0.532 and a 90th percentile of 0.951, overlapping a binder distribution whose
10th percentile is 0.907, and the relaxed threshold would call 73% of tested-inactive compounds
binders. A low false-positive rate on unrelated chemistry does not license a lower threshold, and
discrimination against random chemistry can be near-perfect while discrimination against the chemistry
that matters is poor. The dual constraint, measured inactives and background together, performs work
that neither term performs alone.

### Panel redundancy: how many independent mechanisms does an engaged compound show?

A panel that reports the number of engaged targets invites the reading that three engaged targets are
three independent pieces of evidence. For homologous proteins that reading is wrong, and the
magnitude of the error was measured rather than assumed. Across 400 approved drugs, dopamine D2 and
D3 fire together with a phi correlation of 0.813, and D3 fires for 78% of compounds engaging D2; the
mu and kappa opioid receptors correlate at 0.791, 5-HT2A and 5-HT7 at 0.699, and the dopamine and
noradrenaline transporters at 0.642, with the noradrenaline transporter firing for 84% of compounds
engaging the dopamine transporter. Of 52 targets, 36 fire at least once on that set, but a
correlation analysis of the firing pattern resolves only 14 independent directions. A raw count
therefore overstates the evidence by roughly a factor of two and a half.

Co-firing is not itself an error. A genuinely promiscuous ligand should engage both members of a
homologous pair, and that is a true fact about the compound. Presenting it as corroboration is the
error. The server therefore groups engaged targets into six homology families, reports the number of
independent mechanisms alongside the number of targets, and states the measured correlation wherever
two members of a family both fire. Haloperidol, for example, is reported as five engaged targets but
approximately three independent mechanisms, with dopamine D2 and D3 and the two serotonin receptors
each labelled with their correlation. The correction propagates to the batch table and to both
machine-readable exports, so a stored result cannot be read as better corroborated than it is. The
nicotinic endpoints added in the targeted expansion are the least correlated family measured, at
0.35, so that expansion contributed more independent information than the panel average.

The same analysis corrected an error in its own first formulation, which is recorded because it
bears on how such rates should be reported. Defining engagement as any positive target signal gave an
apparent rate of 60% of random compounds engaging at least one target. That figure is an artefact of
the definition rather than a property of the server: the six base-rate-enrichment endpoints return a
positive signal whenever the predicted probability merely exceeds the training base rate, which
occurs for about 12% of random compounds at each endpoint, whereas the 43 binder endpoints require a
calibrated threshold to be crossed and fire for about 0.6%. Pooling the two counts an ordinary
above-average probability as an event. The quantity that reaches the user is whether any condition
crosses the reporting threshold, which occurs for 11.5% of 600 random structures, against 55.8% of
approved drugs. That rate sits close to the 12.5% measured on library compounds and above the 5.1%
measured on chemistry most distant from training, as expected for a sample drawn without regard to
similarity. Correlation among the endpoints slightly reduces rather than inflates the aggregate: at
least one binder endpoint fires for 22.2% of random structures, against 23.9% for a hypothetical
panel with the same per-endpoint rates firing independently.

### Two conditions opened, and a gating rule corrected

The coverage audit distinguished conditions absent for want of data from conditions absent for want
of an endpoint, and found two of the latter. The calcitonin-gene-related-peptide receptor (1,578
measured activities across 26 sources) and dihydroorotate dehydrogenase (2,558 across 55) were added,
bringing migraine and multiple sclerosis into the panel and taking it to 52 targets across sixteen
conditions. Amyotrophic lateral sclerosis was already scored but through seven mechanisms shared with
other conditions and none specific to it; RIPK1-mediated necroptosis (5,291 activities across 55
sources) was added to close that, and it also contributes to multiple sclerosis, where its inhibitors
are likewise in trial. Two further ALS candidates were rejected on data rather than judgement: SARM1
has 85 measured activities from 3 sources and KCNQ2 has 103 from 12.

Prospectively, on compounds sharing no scaffold with anything their models saw, the CGRP endpoint
reaches an AUROC of 1.000 at a sensitivity of 0.993 and RIPK1 an AUROC of 0.995 at 0.950, both with a
false-positive rate of 0.000 on random chemistry. DHODH ranks at 0.989 but fires for only 0.330 of
its held-out actives at its deployed threshold, and is reported as a low-sensitivity endpoint.

Adding these targets exposed an error in the gating rule that had been invisible while every modelled
condition was central. The blood-brain-barrier term encodes an assumption that the target lies behind
the barrier. For these two conditions it does not. The CGRP receptor acts at the trigeminal ganglion
and dural vasculature, which are outside the barrier, and that is precisely why the gepants and the
anti-CGRP antibodies are effective without central penetration; dihydroorotate dehydrogenase
inhibition acts on proliferating peripheral lymphocytes. Applying the gate suppressed correct calls:
rimegepant engages the receptor at 0.99 and was silenced by a predicted barrier probability of 0.33.
Migraine and multiple sclerosis are therefore exempt from gating, after which all three marketed
gepants report migraine and brequinar reports multiple sclerosis, while peripheral negative controls
remain silent. RIPK1 is deliberately not exempt, its relevance being central.

Two incidental observations support the underlying models rather than the interface. Leflunomide, an
inactive prodrug, scores 0.06 at dihydroorotate dehydrogenase while its active metabolite
teriflunomide scores 0.91, so the model separates them. And teriflunomide's own measured affinity is
pChEMBL 6.45, below the value of 7 that defines a binder here, so declining to call it is
definitionally correct rather than a failure of sensitivity.

### A limitation quantified rather than asserted: state-dependent block

The inability to represent use-dependent sodium-channel blockers was stated above as a boundary. It
was then tested, because the scientifically correct quantity for that class is not the potency, which
depends on the voltage protocol and is therefore not a function of structure, but the shift in
potency between a depolarised and a hyperpolarised holding potential, which is a property of the
molecule and is what use-dependence means.

Across five sodium-channel subtypes and 13,919 activities carrying a pChEMBL value, 2,086 report a
holding potential that can be parsed from the assay description, spread across fifteen distinct
values, and 89 compounds are measured at two or more. Against the 800 measured activities every
deployed endpoint had to clear, that is not a trainable endpoint.

The same 89 compounds show why the question was worth asking and why pooling protocols is not an
option. The median absolute shift is 0.82 log units, the 90th percentile 1.97 and the maximum 3.19,
so a compound can appear sixfold to a thousandfold more or less potent depending only on the protocol
used to measure it. That is the same magnitude as the total error of a good regression model, which
is the quantitative form of the objection: a model trained on pooled protocols would have an error
budget consumed entirely by the experimenter's choice of holding potential. The effect is real, it
matters clinically, and public data cannot support learning it.

## The BrainSafe AI server

**Input.** The user types a compound name or pastes a SMILES string. Names are resolved through
PubChem; structures are standardised with RDKit (salt stripping, largest-fragment selection).

**Output.** Results are laid out as a single scrollable report (Figure 2): a summary card with the
rendered structure, the free-brain-exposure verdict and headline metrics; the mechanistic map
(Figure 3); a brain-relevance panel giving the exposure-scaled score for each condition with its driving
mechanism; a target-engagement table reporting each calibrated probability, its base rate and a
base-rate-aware call; a receptor-binding table; an ADME table; a physicochemical profile; and an
applicability and confidence card giving global and per-endpoint domain flags with the nearest
measured analogue. An About tab documents the methods, the model-selection and validation results,
and an explicit coverage panel stating which mechanisms and diseases the tool can and cannot yet
assess, so a null result is read as an honest unknown rather than as inactivity. A target-engagement
profile plots every modelled mechanism against its own reporting threshold; because each endpoint
carries a different base rate and a different cut, raw probabilities are not comparable between
endpoints, and what is plotted is the distance above each threshold, which is the quantity the
disease layer consumes. Sub-threshold targets are shown as a percentage of their own cut, so a
compound that engages nothing still yields an interpretable figure rather than an empty panel.

Where the server reports nothing, it states the two reasons a null result can arise: the compound may
act through a mechanism outside the panel, or through a modelled target but below its reporting
threshold. The second is quantified for the user, since thresholds are set for precision, holding the
false-positive rate near 0.2% on random chemistry, and the measured cost is a median sensitivity of
0.77 falling to 0.26 at the strictest endpoints. Silence is not evidence of inactivity and the
interface says so with the number attached.

**Export.** Every result is downloadable in four forms generated from a single tidy table, so that
they cannot disagree with the screen or with each other: a data table with units and training context
for every endpoint; a self-contained report that inlines its figures, structure image and styling and
therefore opens offline and prints to a portable document; a structured record carrying thresholds,
screening mode, provenance and the caveats that must travel with the numbers; and the mechanistic map
as a vector figure for direct use in a manuscript.

**Batch screening.** Up to 300 compounds may be pasted or uploaded as delimited or plain text and are
returned as one row each, giving exposure, engaged mechanisms, reported conditions and
applicability-domain status, with the whole table downloadable. This is the mode intended for
prioritising a library before laboratory work is committed.

**Deployment.** The server is distributed as a container image running as an unprivileged process
with a health endpoint for orchestration. A pre-flight script verifies that every declared dependency
resolves at its pinned version, that all model artefacts load, that the knowledge graph is internally
consistent, that chemically unrelated compounds yield distinct and directionally correct profiles,
and that every export format is well formed and self-contained; it exits non-zero on any failure and
is intended to gate a release.

## Discussion

BrainSafe AI combines BBB penetration, a broad measured-data target panel, free-brain exposure and
safety into a single calibrated, applicability-aware and mechanistically interpretable readout for
arbitrary structures. Its central design choice is honesty about uncertainty: enrichment over base rate
rather than raw probability, per-endpoint applicability domains, an explicit coverage statement, and a
variance decomposition that quantifies how far each endpoint depends on chemotype.

The limitations are inherited from the training data and are stated plainly. The most important concerns
prospective use. Scaffold-grouped cross-validation gives a mean classifier AUROC of 0.919, but on
compounds published after the training cutoff the aggregate falls to 0.752. The analysis above shows
that this aggregate conceals two regimes: inside the applicability domain the classifiers retain 0.827
AUROC and the potency models a rank correlation of 0.56, while outside it performance is close to
chance. The practical consequence is that BrainSafe AI should be read as a reliable triage instrument
for chemistry related to measured data, and as a hypothesis generator only for novel chemotypes; the
server states which case applies for every query. Quantitative potency values remain weak priors rather
than affinity predictions, which is why the deployed receptor models are binder classifiers and the
interface shows a predicted pKi only for compounds already classified as binders. The antioxidant
regression that drives the neuroprotection axis has an in-domain prospective rank correlation of 0.59
but essentially no out-of-domain signal, so that axis is a qualitative flag outside the domain.

Target coverage, though broad, is finite, and the boundaries have been measured rather than assumed.
A systematic query of ChEMBL for the mechanisms the panel omits found that the reasons differ in kind
and that only one of them is a matter of quantity. Some targets have too few measured ligands for any
model: kainate receptors (244, 139 and 34 activities for GluK1, GluK2 and GluK3), the vesicular
monoamine transporter VMAT2 (149), SOD1 (29), TDP-43 (9), and C9orf72, which has no ChEMBL target
record at all. Some have volume without diversity, which is a more insidious failure: tau carries
95,345 potency values, more than any deployed endpoint, but a 1,000-activity sample draws on a single
document and 86% of it is one thioflavin-S displacement campaign, so a scaffold-split model would
learn that screening library rather than the protein; huntingtin is worse at 98% from one assay.
Some are obscured by target annotation rather than by data: ChEMBL assigns activity to the protein
rather than to the binding site, so the phencyclidine channel site used by ketamine and the GluN2B
allosteric site used by ifenprodil are pooled despite having unrelated structure-activity
relationships, which is why the latter is modelled and the former is not. Finally, some conditions are
not targets at all. Migraine and multiple sclerosis are reachable through defined mechanisms not yet
included (the calcitonin-gene-related-peptide receptor, 1,578 activities, and dihydroorotate
dehydrogenase, 2,558), whereas stroke and cerebral ischaemia offer no comparable small-molecule
target set, consistent with four decades of failed neuroprotection trials.

Expanding the panel is not costless, which bears on how these gaps should be closed. Each added
endpoint is a further opportunity to fire spuriously, so the probability that a compound engaging
nothing nevertheless receives a reported finding rises with panel size; it stands at 11.5% on random
chemistry. Added endpoints also overlap: 36 of 52 targets fire at least once across approved drugs
but span only 14 independent directions. Coverage should therefore be extended where a measured
weakness demands it, as was done here for addiction, rather than pursued for completeness, and
preference given to mechanisms that are not homologous to those already present.

An operational consequence follows for the coverage statement itself. The list of what the server can
and cannot assess is read as an audit, and an audit that is wrong about its own models is worse than
none. A hand-maintained version of that list drifted: it advertised an endpoint that had been
withdrawn, denied two families that had been added, and quoted sensitivities of 0.54 and 0.58 for
endpoints whose deployed values were 0.774 and 0.660. The modelled list and every sensitivity figure
are now generated from the deployed model registry when the page loads, and a pre-flight check
asserts that no deployed endpoint is missing from the list, that no withdrawn endpoint is advertised,
that no mechanism declared absent is in fact deployed, and that every quoted sensitivity equals the
deployed value.

Three target endpoints carry thin negative classes, and an attempt to enlarge them failed in a way
worth recording. The CGRP receptor, Nav1.6 and AMPA GluA2 have 26, 45 and 32 measured inactives
respectively, and PubChem BioAssay holds 198, 196 and 516 assays for the corresponding genes, which
appeared to offer a remedy. It does not: a scan of 120 assays per target returned no measured
inactive compounds at all, and querying the largest six assays of each directly returns active
compound identifiers against none inactive. These targets have never been through a large public
primary screen; their assays are medicinal-chemistry dose-response series in which every compound
tested was already of interest. The route remains valid for targets that have been screened at scale,
which is why curated inactive sets exist for acetylcholinesterase and GSK-3-beta, but for these three
the thin negative class cannot be remedied from public data and their sensitivity figures stand as
reported.

Two ADME endpoints, hepatocyte clearance and plasma protein binding, are weak under the scaffold
split and are reported as such. Decoy-based validation gives an optimistic bound because decoys are
presumed rather than measured inactives, which is why the near-miss figure and the background
false-positive rate are reported together, and why deployed specificity is additionally measured on
random chemistry and on molecules nothing should bind. The target-to-disease weights in the knowledge
graph are expert-curated rather than learned; the falsification analysis shows they carry no
measurable predictive content beyond the graph topology, and they are retained as a mechanistic prior
and versioned in the repository so that they can be inspected and revised. Predictions concern
molecular target engagement and physicochemical properties, not clinical efficacy, and do not
distinguish agonism from antagonism. The tool is intended for research prioritisation and hypothesis
generation and is not for medical, diagnostic or treatment decisions.

One limitation deserves separate statement because it is not a coverage gap and cannot be closed by
adding targets. The panel defines engagement as high-affinity binding, and an entire pharmacological
class works otherwise. Use-dependent, state-dependent channel blockers act at concentrations one to
two orders of magnitude weaker than the binder threshold, and ChEMBL's own measurements label the
classic antiepileptic drugs inactive at the sodium channels through which they are understood to act.
The server is therefore silent on carbamazepine, phenytoin and lamotrigine, and that silence is
correct with respect to the measured data while being unhelpful with respect to the clinical
question. Extending the tool to such mechanisms would require a different class of endpoint, modelling
use-dependent block or a phenotypic outcome, not a further binding target.

## Data availability

All code, the curated knowledge graph, per-fold validation artefacts and the scripts that regenerate
every table and figure in this manuscript are available at https://github.com/krishna-g-999/brainsafe-ai
under the MIT licence. The trained estimators (0.78 GB) and the raw API responses retrieved from
ChEMBL, BindingDB and PubChem are too large for version control and are deposited separately; both
archives carry a manifest, committed with the code, recording the SHA-256 of the archive and of every
file inside it, so a download is verified rather than trusted. The models archive supersedes
doi:10.5281/zenodo.21858576, which holds the pre-audit models and describes different bytes from the
ones this manuscript reports; the manifest records that superseded identifier explicitly rather than
pointing at a record that would serve the wrong models.

Reproducing the panel from the repository requires only the endpoint tables and
`src/brainsafe/models/train_rf.py`; every random seed is fixed at 42, and the adversarial suite
includes a check that retraining an endpoint reproduces its reported score, which it does to three
decimal places (Figure 6D).

**[TO BE SUPPLIED BEFORE SUBMISSION]** the public server URL, the deposit DOI for version 1.1, the
author list, the corresponding author's address, and the funding statement. These are the only
placeholders in this manuscript; every other value is computed from an artefact in the repository.

## Funding

[TO BE SUPPLIED BEFORE SUBMISSION]

## Figure legends

![Figure 1](figures/Figure1_architecture.png)

**Figure 1.** How a prediction is assembled. (A) A query structure is reduced to its largest organic
fragment, sanitised, and represented as one fixed 1,036-column vector. That vector is scored by four
model families: nine exposure and ADME endpoints, twelve target potency and activity endpoints, the
49-endpoint binder panel, and two auxiliary regressions. Every target score is then admitted only in
proportion to the predicted probability that the compound reaches the brain, and the surviving scores
are ranked by base-rate enrichment rather than by raw probability, so an endpoint that fires often
across the library cannot dominate. Every reported score carries a calibrated probability, a
conformal interval and an applicability-domain distance. (B) The counts in (A) are deployed
estimators. Each was preceded by twenty fits that never serve a prediction and exist only to measure
how the twenty-first behaves on compounds withheld from it: ten random folds and ten scaffold folds.
Across the panel that is 1,420 cross-validation fits over 71 cross-validated endpoints, behind 72
deployed estimators and eight isotonic calibrators.

![Figure 2](figures/Figure2_feature_vector.png)

**Figure 2.** The model input, computed for donepezil by the same featuriser the models use. (A) The
structure, after standardisation. (B) All 1,024 fingerprint bits drawn as a 32 x 32 grid, 47 of them
set. Two properties of this representation bound what the models can do and are stated rather than
left implicit: folding means a set bit reports that some substructure environment hashing to that
index is present, not which one, and excluding chirality means two enantiomers produce byte-identical
rows. The second is why rows identical in feature space are collapsed before any split rather than
left to fall on both sides of one. (C) The twelve descriptors with the values this molecule has.
They are unscaled, because a random forest splits on thresholds and is unchanged by any monotone
rescaling, so no scaler is fitted and none can leak across a split.

![Figure 3](figures/Figure3_cv_design.png)

**Figure 3.** Two cross-validation schemes and the distance between them. (A) The same compounds
partitioned two ways. Under a random split the held-out fold is scattered through every scaffold
class, so a test compound usually has a close analogue in training and the score reports
interpolation. Under a scaffold split the held-out fold is a whole Bemis-Murcko class absent from
training, so the score reports generalisation to chemistry the model has no near neighbour for.
(B) Both scores for every classifier endpoint; the median cost of withholding a scaffold class is
0.027 AUROC. (C) All ten folds behind each mean, because a mean over ten tight folds and a mean over
ten dispersed folds are not the same claim.

![Figure 4](figures/Figure4_pools_and_thresholds.png)

**Figure 4.** Why the decision thresholds are measured on a pool they were not set on. (A) The
background library is partitioned into three pools by a stable hash of the canonical structure, so a
compound's pool is a property of the molecule and never depends on run order. Decoys are drawn from
the first, thresholds are set on the second, and the false-positive rate is measured on the third; no
compound appears in more than one. (B) A threshold chosen as a quantile of a sample and then scored
on that same sample returns the quantile it was given, and cannot exceed it. Measured on the disjoint
evaluation pool, three endpoints exceed the 0.05 target, which under the previous procedure was
arithmetically impossible. That the number can now disagree with its target is the evidence that it
is a measurement. (C) Each deployed endpoint carries two operating points, a sensitive triage
threshold and a stricter screening threshold, with sensitivity measured on actives withheld by
scaffold.

![Figure 5](figures/Figure5_negative_class.png)

**Figure 5.** Recovering the measured negative class. (A) A censored bound settles a label whenever
the whole interval it defines falls on one side of the activity cut. "IC50 > 10 uM" places the true
potency strictly below pChEMBL 5.0 and is a measured non-binder; "IC50 > 100 nM" spans both classes
and is discarded as undecidable rather than guessed at. (B) Class balance for the 57 endpoints
extended, before and after: 21,994 measured non-binders were added and the endpoints above 90 per
cent active fell from 35 to 13. (C) The effect on cross-validated performance, per endpoint, against
the panel as it stood immediately before the merge. Classification becomes slightly harder (median
-0.0040) and regression improves (median +0.0202). Both are the expected direction: real non-binders
are harder negatives than the decoys they join, and a censored bound is a real low-potency anchor
where a regression previously had nothing. BBB and antioxidant_DPPH are flat at exactly zero, not
missing, because neither draws from a ChEMBL target.

![Figure 6](figures/Figure6_validation.png)

**Figure 6.** Four validations that a cross-validated score cannot substitute for. (A) Expected
calibration error before and after isotonic regression fitted on out-of-fold predictions, so the
calibrator never sees a compound in its own fit; mean ECE falls from 0.0795 to 0.0161. (B) Recall on
whole scaffold classes withheld before training, with 95 per cent Wilson intervals (34) and marker area
proportional to the number of withheld actives, so an interval that is wide because the evidence is
thin looks thin. (C) Specificity on chemistry the server should stay quiet about, and external
discrimination on approved drugs absent from the training source. (D) The adversarial suite, in which
each check is written so that it can fail. Five of six pass. The sixth is shown at the same size as
the rest: the applicability-domain flag does not separate non-drug-like chemistry from unseen drugs,
and it is reported rather than retuned until it passes.

![Figure 7](figures/Figure7_binder_panel.png)

**Figure 7.** The binder panel, all 49 endpoints. (A) Each endpoint placed by what it discriminates,
AUROC against compounds measured and found inactive at the same target, and what it recovers,
sensitivity on actives withheld by scaffold. Marker area is the number of measured actives, and
performance tracks it. (B) Every endpoint named, so a reader can look up a target rather than accept
a panel average. The two endpoints withdrawn after specificity testing and those below the
reliability gate are marked, not omitted. An endpoint is withdrawn when its probability band is too
compressed for any threshold to separate real ligands from trivial metabolites; Nav1.1 scored glucose,
urea, acetate and glycine between 0.80 and 0.82 against a threshold of 0.796.

**[TO BE SUPPLIED BEFORE SUBMISSION]** a screenshot of the server interface for a worked example, and
the mechanistic map for a single compound, both of which require the deployed instance.

## References

Every entry was resolved by a live query against CrossRef or Europe PMC and accepted only on a title match, or, where the identity is known and the registered title is a short form, by resolving the DOI and confirming the title and first author. The requested title, the matched title and the score are recorded in `manuscript/references_verified.json`, so the list can be re-checked mechanically. None is written from memory.

1. Zdrazil B, Felix E, Hunter F et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Research. 2024. doi:10.1093/nar/gkad1004.
2. Gilson M, Liu T, Baitaluk M et al. BindingDB in 2015: A public database for medicinal chemistry, computational chemistry and systems pharmacology. Nucleic Acids Research. 2016. doi:10.1093/nar/gkv1072.
3. Meng F, Xi Y, Huang J et al. A curated diverse molecular database of blood-brain barrier permeability with chemical descriptors. Scientific Data. 2021. doi:10.1038/s41597-021-01069-5.
4. Loryan I, Reichel A, Feng B et al. Unbound Brain-to-Plasma Partition Coefficient, Kp,uu,brain—a Game Changing Parameter for CNS Drug Discovery and Development. Pharmaceutical Research. 2022. doi:10.1007/s11095-022-03246-6.
5. Rogers D, Hahn M. Extended-Connectivity Fingerprints. Journal of Chemical Information and Modeling. 2010. doi:10.1021/ci100050t.
6. Breiman L. Random Forests. Machine Learning. 2001. doi:10.1023/a:1010933404324.
7. Bemis G, Murcko M. The Properties of Known Drugs. 1. Molecular Frameworks. Journal of Medicinal Chemistry. 1996. doi:10.1021/jm9602928.
8. Niculescu-Mizil A, Caruana R. Predicting good probabilities with supervised learning. Proceedings of the 22nd international conference on Machine learning  - ICML '05. 2005. doi:10.1145/1102351.1102430.
9. Mysinger M, Carchia M, Irwin J et al. Directory of Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking. Journal of Medicinal Chemistry. 2012. doi:10.1021/jm300687e.
10. Alavijeh M, Chishty M, Qaiser M et al. Drug metabolism and pharmacokinetics, the blood-brain barrier, and central nervous system drug discovery. NeuroRX. 2005. doi:10.1602/neurorx.2.4.554.
11. Wager T, Hou X, Verhoest P et al. Moving beyond Rules: The Development of a Central Nervous System Multiparameter Optimization (CNS MPO) Approach To Enable Alignment of Druglike Properties. ACS Chemical Neuroscience. 2010. doi:10.1021/cn100008c.
12. Saify Z, Sultana N. Role of Acetylcholinesterase Inhibitors and Alzheimer Disease. Drug Design and Discovery in Alzheimer's Disease. 2014. doi:10.1016/b978-0-12-803959-5.50007-6.
13. Decourt B, Macias M, Sabbagh M et al. BACE1 Inhibitors: Attractive Therapeutics for Alzheimer’s Disease. Drug Design and Discovery in Alzheimer's Disease. 2014. doi:10.1016/b978-0-12-803959-5.50010-6.
14. Dezsi L, Vecsei L. Monoamine Oxidase B Inhibitors in Parkinson's Disease. CNS & neurological disorders drug targets. 2017. doi:10.2174/1871527316666170124165222.
15. West A. Achieving neuroprotection with LRRK2 kinase inhibitors in Parkinson disease. Experimental Neurology. 2017. doi:10.1016/j.expneurol.2017.07.019.
16. Riemann D, Spiegelhalder K. Orexin receptor antagonists: a new treatment for insomnia?. The Lancet Neurology. 2014. doi:10.1016/s1474-4422(13)70311-9.
17. Duan Y, Kelley N, He Y. Role of the NLRP3 inflammasome in neurodegenerative diseases and therapeutic implications. Neural Regeneration Research. 2020. doi:10.4103/1673-5374.272576.
18. Yamazaki H, Tanji K, Wakabayashi K, Matsuura S, Itoh K. Role of the Keap1/Nrf2 pathway in neurodegenerative diseases. Pathology international. 2015. doi:10.1111/pin.12261.
19. Kovalenko M, Erdin S, Andrew M et al. Histone deacetylase knockouts modify transcription, CAG instability and nuclear pathology in Huntington disease mice. eLife. 2020. doi:10.7554/elife.55911.
20. RG M, JD M, M L et al. Riluzole for amyotrophic lateral sclerosis (ALS)/motor neuron disease (MND). Amyotrophic Lateral Sclerosis and Other Motor Neuron Disorders. 2003. doi:10.1080/14660820310002601.
21. Ritter J. Cardiac safety, drug‐induced QT prolongation and torsade de pointes (TdP). British Journal of Clinical Pharmacology. 2012. doi:10.1111/j.1365-2125.2012.04193.x.
22. Huang K, Fu T, Gao W, Zhao Y, Roohani Y, Leskovec J, Coley CW, Xiao C, Sun J, Zitnik M. Artificial intelligence foundation for therapeutic science. Nature chemical biology. 2022. doi:10.1038/s41589-022-01131-2.
23. Wu Z, Ramsundar B, Feinberg E et al. MoleculeNet: a benchmark for molecular machine learning. Chemical Science. 2018. doi:10.1039/c7sc02664a.
24. Pedregosa F, Varoquaux G, Gramfort A et al. Scikit-learn: Machine Learning in Python. Journal of Machine Learning Research. 2011;12:2825-2830. https://www.jmlr.org/papers/v12/pedregosa11a.html
25. Chen T, Guestrin C. XGBoost: A Scalable Tree Boosting System. Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. 2016. doi:10.1145/2939672.2939785.
26. Lin H, Lin C, Weng R. A note on Platt’s probabilistic outputs for support vector machines. Machine Learning. 2007. doi:10.1007/s10994-007-5018-6.
27. Norinder U, Carlsson L, Boyer S et al. Introducing Conformal Prediction in Predictive Modeling. A Transparent and Flexible Alternative to Applicability Domain Determination. Journal of Chemical Information and Modeling. 2014. doi:10.1021/ci5001168.
28. Jaworska J, Nikolova-Jeliazkova N, Aldenberg T. QSAR Applicability Domain Estimation by Projection of the Training Set in Descriptor Space: A Review. Alternatives to Laboratory Animals. 2005. doi:10.1177/026119290503300508.
29. Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research. 2000. doi:10.1093/nar/28.1.27.
30. Milacic M, Beavers D, Conley P, Gong C, Gillespie M, Griss J, Haw R, Jassal B, Matthews L, May B, Petryszak R, Ragueneau E, Rothfels K, Sevilla C, Shamovsky V, Stephan R, Tiwari K, Varusai T, Weiser J, Wright A, Wu G, Stein L, Hermjakob H, D'Eustachio P. The Reactome Pathway Knowledgebase 2024. Nucleic acids research. 2024. doi:10.1093/nar/gkad1025.
31. Harding SD, Armstrong JF, Faccenda E, Southan C, Alexander SPH, Davenport AP, Spedding M, Davies JA. The IUPHAR/BPS Guide to PHARMACOLOGY in 2024. Nucleic acids research. 2024. doi:10.1093/nar/gkad944.
32. Streamlit: an open-source app framework. https://streamlit.io
33. RDKit: Open-source cheminformatics. https://www.rdkit.org
34. Wilson E. Probable Inference, the Law of Succession, and Statistical Inference. Journal of the American Statistical Association. 1927. doi:10.1080/01621459.1927.10502953.
