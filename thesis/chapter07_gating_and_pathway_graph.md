# Chapter 7. Exposure gating and the pathway graph

> Every quantitative statement in this chapter was computed during this session from the artefact
> named beside it. The falsification suite was re-run on 29 and 30 August 2026 and several figures
> moved; the values here are those of the current artefacts, and section 7.10 records which changed.

---

## 7.1 What this layer is for, and what it is not

Everything to this point produces per-target probabilities. This chapter is about the rule that turns
seventy independent numbers into one answer, and it is the part of the system that has been most
misdescribed, including by its own authors.

The disease layer takes the engaged targets, weights them through a curated graph, gates them by
predicted exposure, and returns a ranked list of conditions. It is a **route from a mechanism to the
conditions that mechanism touches**. It is not an indication prediction, and Chapter 10 records the
structural reason it cannot be: 27 of the 51 targets in the graph drive more than one of its 16
conditions, and what selects among them, dose, regimen and patient population, is not present in a
structure.

Three components need separate justification, and the falsification suite refuted the stated purpose
of two of them.

## 7.2 The engagement signal: three maps onto one scale

The endpoints are not commensurable. Each has its own base rate, its own threshold, and one has a
continuous readout. A probability of 0.6 at hERG and a probability of 0.6 at BACE1 are not the same
statement, and averaging or ranking them directly would be meaningless.

Three maps carry them onto a common [0, 1] scale.

**Where the negatives are measured**, with base rate $b_t$, the signal is the signed enrichment over
prevalence, clipped at zero:

$$E_t(p) = \begin{cases}\dfrac{p - b_t}{1 - b_t}, & p \ge b_t \\[1.5ex] \dfrac{p - b_t}{b_t}, & p < b_t\end{cases} \in [-1, 1], \qquad s_t(x) = \max\{0,\ E_t(\hat g_t(\hat q_t(x)))\}.$$

**For a decoy-trained binder** the threshold takes the role the base rate takes above:

$$s_t(x) = \max\left\{0,\ \frac{\hat g_t(\hat q_t(x)) - \tau_t}{1 - \tau_t}\right\},$$

so a sub-threshold call, which by Chapter 5's construction lies within the measured background
false-positive rate, contributes exactly zero.

**For the antioxidant assay** the signal is the percentile of the predicted value against the measured
training distribution.

### Why enrichment rather than probability

The base rates make the case concretely. Across the eight core classifiers
(`submission_package/07_MODELS/endpoint_base_rates.json`) they run from **0.236 at hERG to 0.862 at
BACE1**, a span of nearly four-fold:

| Endpoint | n | Base rate |
|---|---:|---:|
| BACE1 | 8,962 | 0.862 |
| GSK-3β | 5,640 | 0.750 |
| BBB | 7,807 | 0.635 |
| AChE | 5,318 | 0.596 |
| BChE | 3,386 | 0.543 |
| MAO-B | 4,806 | 0.505 |
| MAO-A | 3,789 | 0.238 |
| hERG | 10,276 | 0.236 |

A calibrated probability of 0.60 is therefore **strong evidence of activity at hERG**, where only 24
per cent of the measured set is active, and **evidence of inactivity at BACE1**, where 86 per cent is.
Ranking the raw probabilities would put those two compounds in the same place. The enrichment map
puts them on opposite sides of zero, which is the pharmacologically correct answer.

The negative branch matters as much as the positive one. A probability below an endpoint's base rate
is evidence of inactivity, not weak evidence of activity, and mapping it to a small positive number
would misstate it. Clipping at zero means such a target contributes nothing to any condition rather
than contributing a little.

The interface always shows the untransformed calibrated probability beside the signal, because only
$s_t$ is comparable across endpoints and only the probability is interpretable on its own.

## 7.3 The pathway graph

The graph is curated, versioned and anchored to KEGG [31], Reactome [32] and IUPHAR [33]. Read directly
from `app.py`, it holds:

| | |
|---|---:|
| Targets | **51** |
| Conditions | **16** |
| Edges | **83** |
| Conditions whose mechanism is peripheral | 2, migraine and multiple sclerosis |
| Targets driving more than one condition | **27 of 51** |
| Targets per condition | median 5.5, from 1 to 11 |
| Conditions per target | median 2, maximum 3 |

Each edge carries a pathway name, a KEGG identifier, a condition and a weight, for example
`AChE → ('Cholinergic synapse', 'hsa04725', "Alzheimer's disease", 1.0)`. The best-connected condition
is depression and anxiety with 11 targets; the least connected has one.

The weights take **8 distinct values from 0.40 to 1.00**, with a mean of 0.780 and a median of 0.80,
distributed as 1.00 (16 edges), 0.90 (13), 0.85 (2), 0.80 (18), 0.70 (17), 0.60 (9), 0.50 (4) and
0.40 (4). Section 7.6 reports what they contribute.

## 7.4 Aggregation is a maximum, not a sum

A condition is scored by its single strongest engaged mechanism:

$$S_d(x) = \max_{(t,w)\in G(d)} w\cdot s_t(x).$$

Engaging three of a condition's targets therefore scores exactly as engaging its strongest. That looks
like discarding information and is a deliberate choice with a measured justification.

The justification is H8. Counted across approved drugs
(`inversion/results/H8_panel_independence.csv`), **37 targets fire at least once but span only 16
independent directions**. Targets fire in correlated families, and five homologous pairs correlate
above φ = 0.5:

| Family | Pair | φ |
|---|---|---:|
| Opioid | μ-opioid, κ-opioid | **0.798** |
| Monoamine transporters | SERT, NET | 0.732 |
| Dopamine | D2, D3 | 0.720 |
| Serotonin GPCR | 5-HT2A, 5-HT7 | 0.671 |
| Serotonin GPCR | 5-HT1A, 5-HT7 | 0.583 |

A sum over engaged targets would count one observation several times. A promiscuous ligand engaging
both opioid receptors has provided one piece of evidence about the opioid system, not two, and at
φ = 0.798 the second is nearly determined by the first. The maximum is immune to this by construction,
and it has the further property that one target explains the score and can be named to the user.

Co-firing is not itself an error: a promiscuous ligand *should* engage both homologues. Presenting it
as corroboration is the error, and the interface now groups engaged targets by homology family and
quotes the measured correlation wherever two members fire.

The firing pattern is otherwise sparse. A random compound fires a mean of **0.242** binder endpoints
and an approved drug **1.49**, so for most queries the maximum is taken over a set with one non-zero
member, and the choice between maximum and sum rarely arises at all.

## 7.5 The exposure gate, and why H3 is refuted by construction

The gate is the system's central design claim:

$$\tilde S_d(x) = \gamma_d(x)\, S_d(x), \qquad \gamma_d(x) = \begin{cases}1, & d \in P\\ \hat q_{\mathrm{BBB}}(x), & \text{otherwise,}\end{cases}$$

with $P$ the two peripheral conditions. Potency at a target the compound cannot reach contributes
essentially nothing. This is the one place the models are combined multiplicatively and it encodes a
pharmacological fact rather than a statistical convenience.

**H3 asked whether gating discriminates between diseases. It does not, and the proof is algebraic
rather than empirical** (`inversion/results/H3_gating.csv`). Because $\gamma_d$ takes the same value
for all 14 non-peripheral conditions, multiplying every one of their scores by it is a common positive
scaling. A common positive scaling is rank-invariant. The gate therefore cannot change which condition
is ranked first, second or third; it can only move all of them together, and hence decide whether
anything clears the reporting threshold.

The correct description is that **the gate is an exposure filter over a mechanism ranking**. It
determines *whether* to report, never *which* condition to report.

Three things follow, and it is worth separating them because a refutation is easy to over-read.

**The gate is not thereby useless.** Deciding whether to speak is the system's most-used behaviour:
Chapter 4 records that 925 of 1,000 non-CNS compounds receive no disease call, and the gate is a
principal reason. A filter that suppresses 92.5 per cent of irrelevant chemistry is doing substantial
work.

**But it cannot be credited with the disease ranking.** Any claim that exposure gating sharpens the
disease call is arithmetically false, and the manuscript wording was changed because of this result.

**The two peripheral conditions are the exception that shows the rule.** Migraine and multiple
sclerosis are exempted from the gate because their mechanisms act outside the barrier. For those two,
$\gamma_d = 1$ while the other 14 are scaled by $\hat q_{\mathrm{BBB}}$, so the gate *does* change the
ranking between a peripheral and a central condition. H3's refutation is exact for the 14 and does not
extend to the boundary between the two groups. No document in the project states this qualification,
and it should.

## 7.6 The curated weights add nothing measurable

H2 replaced the curated edge weights with uniform weights and with randomly permuted ones, leaving the
topology untouched, and re-scored the disease layer
(`inversion/results/H2_weight_ablation.csv`):

| Weights | Top-3 accuracy |
|---|---:|
| Curated | 0.7901 |
| Uniform, all 1.0 | 0.7897 |
| Randomly permuted | 0.7874 |

The spread from curated to permuted is 0.0027. On the 7,008 compounds the evaluation uses, that is
the difference between about 5,537 and about 5,518 correct answers: **replacing hand-assigned weights
with uniform ones changes the top-3 outcome for roughly 3 compounds in 7,008, and permuting them
changes it for roughly 19.**

The mechanism is visible in section 7.4. Under a maximum, a weight affects the result only by scaling
the winning term or by changing which term wins. The weights span 0.40 to 1.00 with three quarters of
them at 0.70 or above, while the engagement signals span the full range and are exactly zero for
every sub-threshold target. When a condition has one engaged target, which the sparse firing pattern
makes the common case, the weight is a pure scale factor on a single term. The information lies in
**which target connects to which condition**, not in how strongly.

The consequence for the thesis is a change in how the graph is described. The weights are **structure,
not tuned parameters**, no claim is made for them, and the graph would be simpler and no less accurate
without them. That is a refutation the project acted on rather than argued with.

## 7.7 What the disease layer does establish

Three hypotheses tested the layer itself, and they answer different questions. Taken together they
support a narrow claim rather than a broad one.

**H1: the layer carries information about the target-to-disease map.** Scored with hold-out models
only, and restricted to compounds that are training actives of no other panel target, top-3 accuracy
is **0.7901** with a 95 per cent interval of 0.7804 to 0.7995 on 7,008 compounds. Shuffling the
target-to-disease map gives **0.1625** (interval 0.0595 to 0.2919 over 200 permutations, p = 0.005).
Always answering with the three commonest conditions gives **0.5508**. The layer beats both nulls.

The caveat is in the phrase "this project's own map". H1 asks whether the layer recovers the condition
a compound's target maps to under the graph the project itself wrote. That is a real test of internal
consistency and it is not a test against the world.

**H6: against real clinical indications, the layer is far above chance and below a constant answer.**
Ground truth is ChEMBL's phase-4 drug indications, mapped to the panel through a keyword list fixed
before any prediction was made. On the **162 drugs whose exact structure appears nowhere in training**,
top-3 accuracy is **0.3519** (interval 0.2825 to 0.428) against a permutation null of **0.1449**
(p = 0.0005) and a frequency null of **0.6543**.

Both readings are true and both belong in the record. The layer beats the permutation null decisively,
so its output depends on the compound and is not memorisation: drugs whose structure appears nowhere
in training score much as drugs that do, 0.3519 against 0.3934. It does not beat the frequency null,
because pain, depression and psychosis account for most approved CNS indications and a constant answer
naming those three is right about 65 per cent of the time.

Removing the reporting threshold and judging the ranking alone raises the never-seen figure from
0.3519 to **0.4506**, which locates a substantial part of the gap in the decision to stay silent
rather than in the ranking.

**H9: on the two metrics a constant predictor cannot pass, the layer beats chance.** Top-3 accuracy
against a frequency null is the wrong comparison, because that null answers with the same three
conditions for every compound and cannot rank one molecule against another, which is the only thing a
triage tool does. Two metrics immune to that were computed
(`inversion/results/H9_disease_discrimination*.csv`):

- **Mean per-indication AUROC 0.6158**, where a constant predictor scores 0.500 by construction
  whatever its top-k accuracy. The layer beats chance on **7 of 9** indications.
- **Macro-averaged top-3 recall 0.3583 against 0.3333**, averaging per indication so that naming only
  the common conditions cannot carry it.

The spread matters more than the mean:

| Indication | Drugs | AUROC |
|---|---:|---:|
| Depression / anxiety | 31 | **0.7981** |
| Psychosis / schizophrenia | 18 | 0.7647 |
| Alzheimer's disease | 9 | 0.6565 |
| Chronic pain | 64 | 0.6456 |
| Parkinson's disease | 15 | 0.6220 |
| Addiction | 7 | 0.5659 |
| ADHD | 12 | 0.5006 |
| Sleep / wakefulness | 10 | 0.4993 |
| Epilepsy | 15 | **0.4898** |

The layer responds to the compound decisively for depression and psychosis, marginally for ADHD at
0.5006, and not at all for sleep and epilepsy, both at or just below chance. Publishing the mean
without this table would conceal that two of the nine conditions are not being predicted at all.

Epilepsy is the informative failure, and Chapter 5 explains it: most approved antiepileptics are
small, simple and low-affinity, precisely the chemistry that sits below a strict cut, and the
endpoints that would carry them are among the six firing for fewer than half their own actives.

## 7.8 Silence, and the reporting threshold

Below a gated disease score of **0.30** nothing is reported. Silence means no modelled mechanism
cleared its threshold. It is not a claim of inactivity, and Chapter 4 gives the number that makes it
interpretable: the recall the panel achieves at that compound's distance from training chemistry.

Two measured rates bound what the reporting rule does. A user-visible finding is returned for
**7.17 per cent** of random chemistry and **51.5 per cent** of approved drugs
(`inversion/results/H8_panel_independence.csv`). The first is the false-lead rate a further expansion
of the panel would have to beat; the second says the layer stays silent on about half of a set of
compounds that are, by construction, pharmacologically active somewhere.

## 7.9 What this chapter establishes, honestly summarised

The gating and disease layer is the part of the system with the weakest claims, and stating them
precisely is more useful than defending them.

- **Exposure gating works as a filter and cannot work as a discriminator.** Refuted by construction
  for the 14 non-peripheral conditions; the exemption of migraine and multiple sclerosis is a
  qualification the project has not stated.
- **The curated edge weights carry no measurable information.** Uniform weights change roughly 3
  answers in 7,008. The graph's content is its topology.
- **The layer carries real information about mechanism**, at 0.7901 top-3 accuracy against a
  permutation null of 0.1625 under the project's own map.
- **It ranks conditions better than chance for most of them and not for all**, at a mean
  per-indication AUROC of 0.6158, from 0.7981 for depression down to 0.4898 for epilepsy.
- **It does not predict indication**, does not beat a constant answer on top-3 accuracy, and is
  correctly described as a route from a mechanism to the conditions it touches.

The maximum-not-sum rule is the one component of this chapter whose stated justification survived
testing, and it survived because H8 supplied the evidence for it: 37 firing targets spanning 16
independent directions, with five homologous pairs above φ = 0.5.

## 7.10 A note on what changed under re-running

The falsification suite was re-run on 29 and 30 August 2026, after Chapters 1 to 6 of this thesis were
drafted, and four figures moved. They are recorded here because a thesis that quotes a suite designed
to be re-runnable should say what re-running it did.

| Figure | Before | Now |
|---|---:|---:|
| Specificity on non-CNS chemistry | 0.949 | **0.925** |
| H9 mean per-indication AUROC | 0.6035 | **0.6158** |
| H9 macro-averaged top-3 recall | 0.3847 | **0.3583** |
| H8 targets firing across approved drugs | 36 | **37** |

The direction of the specificity change looks alarming: the false-positive rate on non-CNS chemistry
rose from 0.051 to 0.075, and H8's user-visible finding rate on random chemistry rose from 0.0467 to
0.0717 over the same period. **It has been traced, and the system did not become less specific. The
old figure was inflated by a false-negative bug.**

`models_rf/endpoint_context.json` supplies the base rate each core classifier's engagement signal is
computed against. It had been built before the neutralisation retrain, and every base rate except
the barrier model's was too high. All eight have since been corrected:

| Endpoint | Stale base rate | Correct base rate | Change |
|---|---:|---:|---:|
| GSK-3β | 0.9312 | 0.7504 | **-0.1809** |
| hERG | 0.4133 | 0.2363 | -0.1770 |
| MAO-B | 0.6750 | 0.5046 | -0.1705 |
| MAO-A | 0.4071 | 0.2383 | -0.1688 |
| BChE | 0.7043 | 0.5434 | -0.1609 |
| AChE | 0.7240 | 0.5959 | -0.1281 |
| BACE1 | 0.9104 | 0.8621 | -0.0483 |
| BBB | 0.6348 | 0.6348 | 0.0000 |

Because engagement is the enrichment $(p - b)/(1 - b)$, **an inflated base rate mechanically
suppresses real signal**. The specificity run of 29 August is the first to use correct base rates, so
0.949 was measured by a system that was under-reporting and 0.925 is the honest figure.

The arithmetic reproduces exactly. Comparing the two prediction files, the barrier probability, the
hERG probability and the applicability-domain distance are **identical on all 1,000 compounds**, and
only the disease score moved, on 373 of them, 312 upward with a median rise of 0.06. Of the 24
compounds that newly fire, **20 are explained by the MAO-B correction alone**, MAO-B being the
weight-1.0 driver of Parkinson's disease in the graph, which is why 20 of the 24 new false positives
are assigned to that one condition.

Taking stilbene as the worked case: its MAO-B calibrated probability is 0.7363 and its barrier
probability 0.8906. Under the stale base rate the enrichment is $(0.7363 - 0.6750)/0.3250 = 0.1885$,
and gated it gives 0.1679. Under the correct base rate it is $(0.7363 - 0.5046)/0.4954 = 0.4676$, and
gated 0.4165. The artefact recorded 0.168 before and 0.417 after. One compound, one probability, two
base rates.

Two things follow. The system's specificity did not degrade, so nothing needs fixing; the number
needed correcting, and it has been. And the residual 7.5 per cent should be read with the caveat the
artefact already carries: these compounds are presumed inactive because nothing is recorded about
them, and several of the new false positives, stilbene and the isothiocyanates among them, belong to
chemotypes with genuine reported monoamine oxidase activity. The specificity of 0.925 remains a lower
bound, and the correction has made it a more honest one.

---

## Outstanding items for this chapter

1. **The H3 refutation has an unstated exception.** Gating is rank-invariant across the 14
   non-peripheral conditions but not between those and the two peripheral ones, where $\gamma_d = 1$.
   Every document states the refutation without the qualification.
2. ~~The specificity change from 0.949 to 0.925 is untraced.~~ **Traced.** It is not a degradation:
   stale base rates in `endpoint_context.json` were suppressing real engagement signal, and commit
   `fd7aaa0` corrected them. Section 7.10 gives the mechanism and reproduces it on a worked example.
3. The graph's 83 edges carry weights that measurably do nothing. Removing them would simplify the
   artefact and lose no accuracy, and the project has chosen to retain and describe them instead;
   either is defensible, but the choice should be recorded as a choice.
4. H1 is a test against the project's own target-to-disease map. It should never be quoted without
   that qualification, and H6 is the test against the world.
