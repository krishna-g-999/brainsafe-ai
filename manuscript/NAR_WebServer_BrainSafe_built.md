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
integrates 63 models trained on measured public bioactivity data (ChEMBL, BindingDB and the B3DB
blood-brain-barrier database): 50 molecular-target endpoints spanning blood-brain-barrier
penetration, the principal neurodegenerative, psychiatric, neuroinflammatory and sleep-related target
classes, and two cardiac safety liabilities, together with a nine-endpoint ADME and exposure layer
that includes a directly modelled unbound brain-to-plasma partition coefficient. Every endpoint is
validated under 10-fold cross-validation in two regimes, a random split and a scaffold-grouped split
that holds out entire chemical series. Predictions are probability-calibrated, carry an
endpoint-specific applicability-domain flag with the nearest measured analogue, and are combined into
per-disease relevance scores, filtered by blood-brain-barrier exposure and traced through a curated
target-to-pathway-to-disease knowledge graph spanning fourteen brain conditions, including
Alzheimer's disease, Parkinson's disease, amyotrophic lateral sclerosis, Huntington's disease and
epilepsy. The measured-label classifier panel reaches a mean scaffold-split AUROC of 0.92. The binder
panel is validated not against the decoys used to train it but against compounds experimentally
tested on the same target and found inactive, giving a mean AUROC of 0.96 across 38 targets; each
target's decision threshold is calibrated on held-out measured inactives to hold the false-positive
rate near 10%, yielding a mean sensitivity of 0.90. A variance decomposition shows that the wider scaffold
error bars are dominated by genuine chemotype heterogeneity rather than statistical noise, which
motivates the applicability-domain layer. The server returns an auditable, mechanistically
interpretable brain-relevance profile rather than a single opaque score. BrainSafe AI is freely
available at [URL].

---

## Introduction

Central-nervous-system drug discovery has a distinctive failure profile: a candidate can be potent
against its intended target yet never reach the brain, or reach it but carry an unacceptable safety
liability, or engage unintended targets that reshape its clinical profile. Answering "will this
molecule affect the brain, how, and is it safe" therefore means combining blood-brain-barrier (BBB)
penetration, disease-relevant target engagement, free-brain exposure, and safety into one readout.
Existing public tools typically address one of these axes in isolation, and large language models,
while fluent, fabricate measured identifiers and cannot supply a calibrated, auditable answer.

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
*Target engagement*: the panel covers the cholinergic axis (AChE, BChE, alpha-7 nicotinic receptor)
and the amyloid and tau axes (BACE1, GSK-3beta) central to Alzheimer's disease; monoamine oxidase B
and LRRK2 for Parkinson's disease, the latter being the most common genetic cause; the serotonergic,
dopaminergic, noradrenergic, opioid, cannabinoid, histaminergic, adenosine and sigma-1 systems that
underlie depression, anxiety, psychosis, addiction, attention deficit, chronic pain, and sleep
regulation. *Safety*: hERG blockade is a leading cause of late-stage cardiovascular attrition and is
modelled as an explicit liability. *Developability*: solubility, lipophilicity, plasma protein binding
and hepatocyte clearance determine whether an achievable dose sustains exposure, and a measured
antioxidant endpoint captures the oxidative-stress axis common to neurodegeneration.

In total the server comprises 59 endpoints (50 molecular-target endpoints and 9 ADME endpoints),
realised as 63 trained models because four receptor targets are represented both as potency
regressions and as binder classifiers.

### Training data

Protein-target activity is pooled from ChEMBL (pChEMBL values) and BindingDB at the compound level;
BBB penetration uses the B3DB database augmented with FDA-curated approved drugs; the antioxidant
endpoint uses measured DPPH pIC50 values; and the nine ADME endpoints use measured sets from
Therapeutics Data Commons, MoleculeNet, B3DB and ChEMBL. The core target panel is trained on 64,474
measured records across 61,317 unique compounds. No value is imputed and no source overrides a
measurement. Per-endpoint compound counts, scaffold counts and class balance are given in Tables 1 to 3.

### Methods compared and model selection

Molecules are represented by a 1024-bit ECFP-4 fingerprint plus twelve interpretable physicochemical
descriptors (molecular weight, cLogP, topological polar surface area, hydrogen-bond donors and
acceptors, rotatable bonds, aromatic rings, fraction sp3, ring count, heavy-atom count, formal charge,
and QED drug-likeness). We benchmarked five model families under identical 10-fold cross-validation on
both split regimes (Figure 4a). Two are baselines: a k-nearest-neighbour Tanimoto read-across, and
L2-regularised logistic regression. Two are gradient-boosted tree ensembles: XGBoost and histogram
gradient boosting. The fifth is a random forest. Separately, we trained a graph isomorphism network
(GIN) directly on molecular graphs to test whether a learned representation would beat fixed
descriptors (Figure 4b).

On the scaffold split the random forest gave the highest mean classifier AUROC (0.914), ahead of
XGBoost (0.905), histogram gradient boosting (0.901), the k-nearest-neighbour read-across (0.867) and
logistic regression (0.808); paired DeLong and bootstrap tests confirm the random forest exceeds the
read-across on every endpoint. The graph neural network did not outperform the random forest on any
tested endpoint (Figure 4b), consistent with the fingerprint-plus-descriptor representation being
sufficient at this data scale. We therefore deploy a probability-calibrated random forest for every
endpoint: it gives the best scaffold-split accuracy, calibrates stably, is interpretable through
feature importance, and needs no GPU, which keeps the web server lightweight.

### Cross-validation design

Every endpoint is evaluated by 10-fold cross-validation under two regimes (Figure 5a, 5b). In the
**random** regime folds are drawn at random, so a test compound will usually have close analogues in
the training folds; this measures interpolation within known chemistry. In the **scaffold-grouped**
regime folds are formed by GroupKFold over generic Bemis-Murcko scaffolds, so every compound sharing a
scaffold is assigned to the same fold and each fold holds out entire chemical series unseen during
training; this measures extrapolation to new chemotypes. Fold sizes are equal to within one compound in
both regimes, so the two are directly comparable. All reported means and standard deviations are across
the ten folds, and per-fold values for every endpoint are provided in the repository
(`results/tables/manuscript_T2_per_fold.csv`).

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
property-matched decoys sampled from a 75,000-compound background library with a maximum ECFP-4
Tanimoto below 0.35 to any positive, so the model must learn structure rather than a similarity
shortcut. Each classifier is a compact random forest with prefit sigmoid calibration. Validated
against held-out near-miss decoys (Tanimoto 0.35 to 0.55) rather than the easy training decoys, these
models reach AUROC 0.86 to 0.99 with a background false-positive rate of 0.6% to 9.8% at a 0.5
threshold (Table 2, Figure 6b). We report the near-miss figure as the honest metric because the
easy-decoy AUROC (0.99 throughout) is inflated by the dissimilarity of the negatives.

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

Every classifier is isotonically calibrated (mean expected calibration error 0.072 to 0.012). For
every prediction the server computes the maximum ECFP-4 Tanimoto similarity of the query to each
endpoint's own measured chemistry and reports an endpoint-specific in-domain, near-domain or
out-of-domain flag together with the nearest measured analogue and its structure. The ADME layer
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
hsa05012), the Reactome KEAP1-NFE2L2 oxidative-stress response (R-HSA-9755511) and IUPHAR Guide to
Pharmacology associations. A per-disease relevance score is the strongest engaged target for that
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

The server is a single-page application built with Streamlit and RDKit. The user submits either a
compound name, resolved to a structure through PubChem, or a SMILES string. Models are loaded once and
cached; a complete profile across all 40 models, including both applicability-domain calculations,
returns in approximately 3.3 seconds on a single CPU core, of which the domain calculations against the
75,000-compound reference library account for under 0.05 seconds. The interface is self-contained and
requires no installation or GPU.

## Results

### Per-endpoint performance

Tables 1 to 3 report, for every endpoint, the measured data size, the number of distinct scaffolds,
class balance where applicable, the training and test set size per fold, and the mean and standard
deviation of the appropriate metric under both cross-validation regimes. Figure 6 summarises the same
results graphically.

The eight target classifiers reach a mean AUROC of 0.960 (random) and 0.919 (scaffold); BACE1 is the
most robust to chemotype change (0.967 to 0.956) and MAO-A the least (0.947 to 0.868). An external
test of the BBB model on 306 FDA-curated approved drugs absent from training gives AUROC 0.774. The
receptor potency regressions reach R2 0.60 to 0.68 (random) and 0.39 to 0.58 (scaffold). The 18
decoy-aware binder classifiers reach 0.86 to 0.99 against near-miss decoys. Within the ADME layer,
performance ranges from strong (solubility R2 0.76, P-gp inhibition AUROC 0.94 under the scaffold
split) to weak and explicitly disclosed (hepatocyte clearance R2 0.19, K_p,uu R2 0.35).

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
| logbb | Total brain/plasma ratio (logBB) | regression | 1,058 | R2 | 0.577 ± 0.081 | 0.455 ± 0.145 | Bulk brain distribution |
| caco2_permeability | Caco-2 permeability | regression | 897 | R2 | 0.734 ± 0.051 | 0.593 ± 0.126 | Passive membrane permeability |
| pgp_substrate | P-glycoprotein substrate | classification | 1,371 | AUROC | 0.858 ± 0.032 | 0.808 ± 0.054 | Active efflux out of the brain |
| pgp_inhibition | P-glycoprotein inhibition | classification | 1,212 | AUROC | 0.955 ± 0.018 | 0.937 ± 0.024 | Efflux-mediated drug interactions |
| solubility | Aqueous solubility (logS) | regression | 9,573 | R2 | 0.804 ± 0.017 | 0.763 ± 0.066 | Formulation and absorption |
| lipophilicity | Lipophilicity (logD) | regression | 4,200 | R2 | 0.639 ± 0.028 | 0.564 ± 0.054 | Permeability/promiscuity balance |
| plasma_protein_binding | Plasma protein binding | regression | 1,797 | R2 | 0.434 ± 0.090 | 0.374 ± 0.104 | Determines free fraction |
| clearance_hepatocyte | Hepatocyte clearance | regression | 1,020 | R2 | 0.230 ± 0.104 | 0.193 ± 0.048 | Metabolic stability, exposure duration |


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
| AChE | Acetylcholinesterase | 2020 | 3,392 | 839 | AUROC | 0.785 |
| BChE | Butyrylcholinesterase | 2021 | 2,015 | 527 | AUROC | 0.737 |
| BACE1 | Beta-secretase 1 | 2017 | 6,497 | 1,604 | AUROC | 0.908 |
| GSK3B | Glycogen synthase kinase-3 beta | 2021 | 3,270 | 417 | AUROC | 0.657 |
| MAO_A | Monoamine oxidase A | 2020 | 1,667 | 446 | AUROC | 0.611 |
| MAO_B | Monoamine oxidase B | 2020 | 2,697 | 755 | AUROC | 0.781 |
| hERG | hERG potassium channel | 2019 | 4,570 | 1,248 | AUROC | 0.785 |
| D2 | Dopamine D2 receptor | 2019 | 6,081 | 1,396 | R2 | 0.042 |
| A2A | Adenosine A2A receptor | 2021 | 4,479 | 1,057 | R2 | 0.338 |
| HT2A | Serotonin 5-HT2A receptor | 2020 | 4,276 | 895 | R2 | 0.182 |
| SERT | Serotonin transporter | 2015 | 3,380 | 999 | R2 | 0.100 |
| antioxidant_DPPH | Radical-scavenging capacity | 2016 | 2,340 | 522 | R2 | 0.009 |


Classifier endpoints: mean AUROC 0.752 (range 0.611 to 0.908). Regression endpoints: mean R2 0.134 (range 0.009 to 0.338).


## Table 6. Prospective validation under a scaffold hold-out

Twenty per cent of Bemis-Murcko scaffolds were withheld per target and every model retrained on the remainder, so no held-out compound shares a scaffold with anything its model saw. Thresholds were recalibrated on held-out negatives and an independent background sample. Targets marked excluded produced a threshold at the permitted floor, meaning no separation from background chemistry, and do not contribute to the pooled estimate.


| Target | Train actives | Held-out actives | Held-out scaffolds | Threshold | Recall | 95% CI | Note |
|---|---|---|---|---|---|---|---|
| OX2 | 2,890 | 962 | 233 | 0.050 | 0.992 | [0.984, 0.996] | excluded |
| LRRK2 | 906 | 267 | 90 | 0.050 | 0.985 | [0.962, 0.994] | excluded |
| mTOR | 2,421 | 563 | 216 | 0.377 | 0.977 | [0.961, 0.986] |  |
| H3 | 2,615 | 597 | 308 | 0.758 | 0.926 | [0.903, 0.945] |  |
| PDE4B | 919 | 209 | 99 | 0.148 | 0.923 | [0.879, 0.952] |  |
| DAT | 1,130 | 207 | 93 | 0.315 | 0.903 | [0.855, 0.937] |  |
| OX1 | 2,601 | 711 | 178 | 0.955 | 0.895 | [0.870, 0.915] |  |
| OPRM1 | 3,004 | 714 | 303 | 0.838 | 0.895 | [0.870, 0.915] |  |
| NLRP3 | 186 | 36 | 27 | 0.053 | 0.889 | [0.747, 0.956] | excluded |
| PDE10A | 3,212 | 872 | 326 | 0.980 | 0.877 | [0.854, 0.897] |  |
| SERT | 2,293 | 735 | 190 | 0.573 | 0.861 | [0.834, 0.884] |  |
| HT6 | 2,023 | 718 | 167 | 0.961 | 0.858 | [0.830, 0.882] |  |
| CSF1R | 1,543 | 241 | 152 | 0.726 | 0.846 | [0.796, 0.887] |  |
| D3 | 3,337 | 790 | 351 | 0.980 | 0.846 | [0.819, 0.869] |  |
| Nav1_1 | 40 | 25 | 9 | 0.325 | 0.840 | [0.653, 0.936] |  |
| OPRK1 | 2,497 | 562 | 235 | 0.948 | 0.838 | [0.805, 0.866] |  |
| HDAC1 | 2,414 | 641 | 292 | 0.900 | 0.836 | [0.806, 0.863] |  |
| CB1 | 2,201 | 479 | 200 | 0.809 | 0.827 | [0.790, 0.858] |  |
| NET | 1,307 | 265 | 90 | 0.646 | 0.819 | [0.768, 0.861] |  |
| Nav1_7 | 2,250 | 494 | 179 | 0.941 | 0.818 | [0.781, 0.849] |  |
| HT1A | 2,902 | 801 | 327 | 0.980 | 0.814 | [0.786, 0.839] |  |
| P2X7 | 2,432 | 506 | 139 | 0.979 | 0.800 | [0.763, 0.833] |  |
| A2A | 3,428 | 924 | 364 | 0.994 | 0.777 | [0.749, 0.803] |  |
| A1 | 1,514 | 429 | 177 | 0.972 | 0.769 | [0.727, 0.807] |  |
| a7nAChR | 251 | 86 | 35 | 0.648 | 0.767 | [0.668, 0.844] |  |
| GABA_A | 184 | 17 | 16 | 0.089 | 0.765 | [0.527, 0.904] |  |
| HT7 | 1,166 | 312 | 133 | 0.954 | 0.721 | [0.669, 0.768] |  |
| HDAC6 | 2,960 | 755 | 335 | 0.990 | 0.699 | [0.666, 0.731] |  |
| Nav1_5 | 198 | 44 | 27 | 0.197 | 0.682 | [0.534, 0.800] |  |
| Sigma1 | 1,571 | 357 | 186 | 0.961 | 0.675 | [0.625, 0.722] |  |
| GluN2B | 652 | 252 | 35 | 0.999 | 0.619 | [0.558, 0.677] |  |
| HT2A | 3,161 | 790 | 317 | 0.987 | 0.618 | [0.583, 0.651] |  |
| COX2 | 918 | 264 | 67 | 0.964 | 0.587 | [0.527, 0.645] |  |
| D2 | 2,880 | 796 | 341 | 0.989 | 0.580 | [0.546, 0.614] |  |
| SIRT1 | 134 | 28 | 20 | 0.529 | 0.464 | [0.295, 0.642] |  |
| mGluR5 | 975 | 225 | 106 | 0.995 | 0.444 | [0.381, 0.510] |  |
| MT1 | 545 | 152 | 48 | 0.996 | 0.283 | [0.217, 0.359] |  |
| KEAP1 | 97 | 25 | 15 | 0.978 | 0.240 | [0.115, 0.434] |  |
| GluA2 | 74 | 23 | 9 | 0.816 | 0.217 | [0.097, 0.419] |  |


Pooled recall 12,325/15,609 = 0.790; median per-target 0.807; 19 of 36 targets at or above 0.80.

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
assess, so a null result is read as an honest unknown rather than as inactivity.

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

Target coverage, though broad, is finite: ionotropic glutamate and GABA-A receptors,
protein-aggregation and neuroinflammation phenotypes, and epigenetic targets are not yet modelled, and
conditions outside the current mechanistic panel (for example amyotrophic lateral sclerosis,
Huntington's disease and stroke) are not scored. Two ADME endpoints, hepatocyte clearance and plasma
protein binding, are weak under the scaffold split and are reported as such. Decoy-based validation
gives an optimistic bound because decoys are presumed rather than measured inactives, which is why the
near-miss figure and the background false-positive rate are reported together. The target-to-disease
weights in the knowledge graph are expert-curated rather than learned, and are versioned in the
repository so that they can be inspected and revised. Predictions concern molecular target engagement
and physicochemical properties, not clinical efficacy, and do not distinguish agonism from antagonism.
The tool is intended for research prioritisation and hypothesis generation and is not for medical,
diagnostic or treatment decisions.

## Data availability

All code, trained models, the curated knowledge graph, per-fold validation artifacts and the scripts
that regenerate every table and figure in this manuscript are available at
https://github.com/krishna-g-999/brainsafe-ai under [license]. The server is freely accessible at
[URL] with no login requirement.

## Funding

[To be completed.]

## Figure legends

![Figure 1](figures/Figure1_endpoint_rationale.png)

**Figure 1.** Endpoint selection rationale. Left, the four sequential questions a CNS candidate must
satisfy and the endpoints modelled in each layer, with measured compound counts. Right, the eleven
conditions the panel informs, the targets that drive each, and the mechanisms not yet modelled.

**Figure 2.** The BrainSafe AI report for a query compound (screenshot to be inserted): summary card,
mechanistic map, brain-relevance panel, target and receptor tables, ADME panel, and the applicability
and confidence card.

![Figure 3](figures/Figure2_mechanism.png)

**Figure 3.** The mechanistic map, shown for haloperidol. A four-tier diagram (compound, target,
pathway, disease) in which only engaged targets are drawn, connector weight encodes engagement
strength, node order is arranged to minimise crossings, and pathways are anchored to KEGG and Reactome
identifiers.

![Figure 4](figures/Figure3_model_selection.png)

**Figure 4.** Model selection. (a) Mean scaffold-split AUROC across the classifier panel for the five
model families compared: logistic regression, k-nearest-neighbour read-across, histogram gradient
boosting, XGBoost, and the deployed random forest (gold). (b) Held-out comparison of a graph
isomorphism network against the random forest per endpoint.

![Figure 5](figures/Figure5_cv_design_and_errorbars.png)

**Figure 5.** Cross-validation design and the origin of the error bars. (a) Random 10-fold: every fold
contains all chemical series, so test compounds have close training analogues. (b) Scaffold-grouped
10-fold: each fold holds out whole chemical series, so test compounds are unseen chemotypes.
(c) Decomposition of the between-fold standard deviation into sampling noise and genuine chemotype
heterogeneity for each endpoint (left bar of each pair, random split; right bar, scaffold split).
Percentages give the heterogeneity share of the scaffold variance.

![Figure 6](figures/Figure6_all_endpoints.png)

![Figure 7](figures/Figure7_temporal_by_domain.png)

**Figure 7.** Prospective performance is governed by the applicability domain. Models are trained on
compounds published before a cutoff year and tested on those published after it, with the future
compounds stratified by maximum Tanimoto similarity to the training set. (a) Classifier AUROC per
endpoint. (b) Rank correlation for the potency and antioxidant models, the decision-relevant metric for
triage. (c) Summary across endpoints: predictive power is retained inside the domain and lost outside
it. The server reports which regime each query falls into.

**Figure 6.** Complete per-endpoint performance. (a) Target classifiers under both cross-validation
regimes (bars, mean; whiskers, standard deviation across the ten folds). (b) Decoy-aware binder
classifiers evaluated against near-miss decoys. (c) ADME regression endpoints under the scaffold split;
the dotted line marks R2 = 0.3, below which an endpoint is reported as weak.

![Figure 8](figures/Figure8_scaffold_holdout.png)

**Figure 8.** Prospective sensitivity under a scaffold hold-out. Twenty per cent of Bemis-Murcko
scaffolds were withheld per target and all 39 trainable binder models retrained on the remainder, so
no held-out compound shares a scaffold with anything its model saw. (a) Recall per target with 95%
Wilson intervals, coloured green at 0.80 or above, amber between 0.50 and 0.80, red below 0.50; the
dashed line is the pooled estimate. (b) Distribution across targets. Three targets whose thresholds
collapsed to the floor are excluded.

![Figure 9](figures/Figure9_specificity.png)

**Figure 9.** Specificity on 1000 compounds with no recorded activity at any modelled target. (a)
Distribution of the highest disease score, with the 0.30 actionable threshold marked. (b) Proportion
silent against proportion firing, with the 95% Wilson interval. (c) Conditions to which the false
positives are assigned.

![Figure 10](figures/Figure10_performance.png)

**Figure 10.** (a) The deployed operating point, sensitivity from the scaffold hold-out against
specificity from the non-CNS set, with 95% intervals on both axes. (b) Distribution of AUROC against
held-out measured inactives across the binder panel.

## References

[To be completed. Anchor citations: ChEMBL; BindingDB; B3DB; RDKit; Therapeutics Data Commons;
KEGG; Reactome; IUPHAR/BPS Guide to Pharmacology; scikit-learn; Streamlit; Bemis and Murcko scaffolds;
isotonic calibration; DeLong test.]

Each entry was resolved by exact-title query against CrossRef or Europe PMC and accepted only above a normalised title-similarity of 0.82. The requested title, the matched title and the similarity score are recorded in `references_verified.json`, so every entry can be re-checked mechanically. None is written from memory.

1. Wilson E. Probable Inference, the Law of Succession, and Statistical Inference. Journal of the American Statistical Association. 1927. doi:10.1080/01621459.1927.10502953
2. DeLong E, DeLong D, Clarke-Pearson D. Comparing the Areas under Two or More Correlated Receiver Operating Characteristic Curves: A Nonparametric Approach. Biometrics. 1988. doi:10.2307/2531595
3. Bemis G, Murcko M. The Properties of Known Drugs. 1. Molecular Frameworks. Journal of Medicinal Chemistry. 1996. doi:10.1021/jm9602928
4. Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research. 2000. doi:10.1093/nar/28.1.27
5. Niculescu-Mizil A, Caruana R. Predicting good probabilities with supervised learning. Proceedings of the 22nd international conference on Machine learning  - ICML '05. 2005. doi:10.1145/1102351.1102430
6. Jaworska J, Nikolova-Jeliazkova N, Aldenberg T. QSAR Applicability Domain Estimation by Projection of the Training Set in Descriptor Space: A Review. Alternatives to Laboratory Animals. 2005. doi:10.1177/026119290503300508
7. Wager T, Hou X, Verhoest P et al. Moving beyond Rules: The Development of a Central Nervous System Multiparameter Optimization (CNS MPO) Approach To Enable Alignment of Druglike Properties. ACS Chemical Neuroscience. 2010. doi:10.1021/cn100008c
8. Mysinger M, Carchia M, Irwin J et al. Directory of Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking. Journal of Medicinal Chemistry. 2012. doi:10.1021/jm300687e
9. Saify Z, Sultana N. Role of Acetylcholinesterase Inhibitors and Alzheimer Disease. Drug Design and Discovery in Alzheimer's Disease. 2014. doi:10.1016/b978-0-12-803959-5.50007-6
10. Decourt B, Macias M, Sabbagh M et al. BACE1 Inhibitors: Attractive Therapeutics for Alzheimer’s Disease. Drug Design and Discovery in Alzheimer's Disease. 2014. doi:10.1016/b978-0-12-803959-5.50010-6
11. Riemann D, Spiegelhalder K. Orexin receptor antagonists: a new treatment for insomnia?. The Lancet Neurology. 2014. doi:10.1016/s1474-4422(13)70311-9
12. Yamazaki H, Tanji K, Wakabayashi K, Matsuura S, Itoh K. Role of the Keap1/Nrf2 pathway in neurodegenerative diseases. Pathology international. 2015. doi:10.1111/pin.12261
13. Gilson M, Liu T, Baitaluk M et al. BindingDB in 2015: A public database for medicinal chemistry, computational chemistry and systems pharmacology. Nucleic Acids Research. 2016. doi:10.1093/nar/gkv1072
14. Dezsi L, Vecsei L. Monoamine Oxidase B Inhibitors in Parkinson's Disease. CNS & neurological disorders drug targets. 2017. doi:10.2174/1871527316666170124165222
15. Wu Z, Ramsundar B, Feinberg E et al. MoleculeNet: a benchmark for molecular machine learning. Chemical Science. 2018. doi:10.1039/c7sc02664a
16. Genuer R, Poggi J. Random Forests. Use R!. 2020. doi:10.1007/978-3-030-56485-8_3
17. Zdrazil B, Felix E, Hunter F et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Research. 2024. doi:10.1093/nar/gkad1004
18. Milacic M, Beavers D, Conley P, Gong C, Gillespie M, Griss J, Haw R, Jassal B, Matthews L, May B, Petryszak R, Ragueneau E, Rothfels K, Sevilla C, Shamovsky V, Stephan R, Tiwari K, Varusai T, Weiser J, Wright A, Wu G, Stein L, Hermjakob H, D'Eustachio P. The Reactome Pathway Knowledgebase 2024. Nucleic acids research. 2024. doi:10.1093/nar/gkad1025
19. Harding SD, Armstrong JF, Faccenda E, Southan C, Alexander SPH, Davenport AP, Spedding M, Davies JA. The IUPHAR/BPS Guide to PHARMACOLOGY in 2024. Nucleic acids research. 2024. doi:10.1093/nar/gkad944
20. extended connectivity fingerprints. The IUPAC Compendium of Chemical Terminology. 2025. doi:10.1351/goldbook.11443

## Software

- RDKit: Open-source cheminformatics. https://www.rdkit.org
- Streamlit: an open-source app framework. https://streamlit.io

Requested but not resolved above the similarity threshold, and therefore not cited: b3db, sklearn, tdc, kpuu, xgboost, gin_gnn, platt, conformal, herg_pred, bbb_ml, cns_attrition, lrrk2_pd, nlrp3_neuro, hdac_hd, riluzole_als.
