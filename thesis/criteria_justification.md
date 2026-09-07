# Justification of the criteria, from collection to testing

> **Provenance.** Every figure here was read from a file in this repository during the session in
> which this was written. The quantitative claims are pinned by `thesis/verify_chapter_numbers.py`
> under the key `criteria`.

## The standard this document holds itself to

A methods section that says "we used a random forest with 300 trees" has justified nothing. A
justification has to say why 300 and not 100, what the alternative would have cost, and whether the
answer was measured or assumed. The honest difficulty is that in any real pipeline the criteria fall
into four very different classes, and a document that presents them all in the same voice is
misleading whatever its individual sentences say.

So every criterion below is labelled:

| | |
|---|---|
| **[MEASURED]** | Tested in this project against the alternative. The number is given. |
| **[CONVENTION]** | Standard practice in the field, adopted deliberately, not tested here. |
| **[FORCED]** | The data left no choice. |
| **[DECLARED]** | A design decision with no evidence behind it. Stated as a choice, defensible, but not a finding. |

There are more **[DECLARED]** entries than a reader might expect, and section 9 collects them
deliberately rather than leaving them scattered. The most important structural fact about this
pipeline is stated there: **no hyper-parameter search of any kind was run.** That has a cost and a
benefit, and both are argued.

---

# Part 1. Data collection

## 1.1 Two independent measured sources, pooled at the compound level — **[MEASURED]**

ChEMBL and BindingDB are pooled for the protein targets. Each source contributes a per-compound
median potency, and the two are combined by InChIKey of the desalted parent. Nothing is imputed and
neither source overrides the other.

*Why two.* One source is a single point of failure for curation error, and the review asked for more
data. ChEMBL was verified near-complete for these targets, live totals within about 3 per cent of
what was held, so genuine growth required an independent source rather than a deeper query.

*Why this is not simply "more data is better".* BindingDB's affinity export returns actives, so adding
it could have inflated performance by enriching the positive class. That was audited rather than
assumed: against the ChEMBL-only baseline the scaffold-split headline metric moved by a mean of
**−0.0002**, range −0.015 to +0.013. The addition neither inflates nor degrades. A2A gained most,
+0.013 R² on 1,238 extra compounds.

*What it cost.* hERG remained ChEMBL-only because BindingDB rate-limited the harvest. That asymmetry
is real and is declared rather than smoothed over.

## 1.2 Potency measurements only: IC50, Ki, Kd, EC50 — **[CONVENTION]**

These four are the measurement types that can be placed on one comparable affinity scale. Percentage
inhibition at a single concentration, thermal shift and similar readouts cannot, because they depend
on assay conditions that are not recorded uniformly. Including them would mean pooling quantities
that are not the same quantity.

## 1.3 The pChEMBL scale — **[CONVENTION]**

All potencies are handled as $-\log_{10}$ of the molar value. Two reasons, and the second is the one
usually left unsaid. It makes IC50, Ki, Kd and EC50 numerically comparable; and affinity is
log-distributed, so a linear scale would let a handful of nanomolar compounds dominate any mean or
any regression loss.

## 1.4 Median across replicates — **[CONVENTION]**

A compound measured many times against one target contributes its per-source median. The median
rather than the mean because public bioactivity data contains transcription errors and unit mistakes,
which are outliers, and the median is insensitive to them.

## 1.5 Censored bounds recovered, with an asymmetric rule — **[DECLARED]**, and the reasoning is the justification

A compound tested and found inactive is often deposited as an inequality, "IC50 > 10 µM", rather than
a number. The conventional query for a numeric potency discards exactly those rows, which throws away
the measured negatives and leaves a model learning actives against decoys.

Two rules, and the asymmetry between them is deliberate:

- **A bound settles a label only when the whole interval falls one side of the activity cut.** One
  that spans the cut is undecidable and is discarded. This is the conservative choice; the
  alternative, taking the bound at face value, would manufacture labels the experiment did not
  support.
- **Bounds aggregate by minimum where exact values aggregate by median**
  (`rebuild_endpoints.py:136` against `:110`). A bound is an upper limit on potency, so across
  several bounds the weakest is the only one every measurement agrees with.

*The evidence that this matters*: **14,420 rows sit at pChEMBL exactly 5.0 with label 0**. 5.0 is
10 µM, the standard inactive bound. Those are measured negatives that a numeric-only query would have
lost.

## 1.6 Measured labels only, never curated annotation — **[MEASURED]**, and this one is a scar

Every supervised endpoint is trained on measured bioactivity. Not on qualitative annotation, not on
disease association, not on anything derived from the same knowledge that built the pathway graph.

*Why the rule exists.* An earlier prototype trained on curated annotation scores. Feature ablation
showed it was reading the answer back out of the disease-association features rather than learning
from structure: with structure alone its performance **collapsed to near zero**. The rule is the
scar tissue from that, and it is the reason the disease layer in this system is a separate,
non-learned graph rather than a trained model. If the graph were trained on the same annotations that
define the conditions, every disease result would be circular.

## 1.7 Bulk high-throughput inactives tested and **rejected** — **[MEASURED]**

This is the strongest single justification in the document, because the decision went against the
obvious move.

Adding 4,276 measured PubChem inactives to GSK-3β, the worst-skewed endpoint at 93 per cent active,
corrected the predicted DrugBank base rate from 71.6 per cent to 16.4 per cent, which is a genuine
improvement. It also raised scaffold AUROC from **0.9369 to 0.9891**.

That second number is why the change was reverted. The added negatives have a **median Tanimoto of
0.288 to the actives**: they are chemically unlike inhibitors, so the model is separating "looks like
a kinase inhibitor" from "does not", which is an easy-negative artefact and not discrimination. The
honest estimate stays 0.9369. The base-rate skew is handled instead by calibration and the
applicability flag, and disclosed as a data limitation.

*The principle generalised*: a change that improves a metric is not thereby an improvement. The
correct future approach is stated in the log as similarity-matched hard negatives only, which is
precisely what the binder panel's hybrid scheme (section 5.4) later implemented.

---

# Part 2. Labelling criteria

## 2.1 Active at pChEMBL ≥ 6.0 — **[CONVENTION]**

pChEMBL 6.0 is 1 µM. It is the conventional boundary at which a compound is treated as active in
medicinal chemistry, and it is the point below which a hit is not usually worth progressing.

## 2.2 Inactive below pChEMBL 5.0 — **[CONVENTION]**

5.0 is 10 µM, the concentration at which a compound is conventionally called inactive and the ceiling
of most screening assays, which is why so many censored bounds sit exactly there.

## 2.3 The 5–6 grey zone discarded — **[DECLARED]**

`label_from` returns −1 for anything in [5.0, 6.0) and those rows are dropped.

*The argument.* A compound at 300 nM is neither confidently active nor confidently inactive. Forcing
it into either class teaches the model a boundary the data does not support, and because
assay-to-assay variability on public data is itself of the order of half a log unit, the grey zone is
roughly the width of the measurement error.

*The cost, stated honestly.* Discarding a band of the data reduces the training set and removes the
examples nearest the decision boundary, which are the informative ones for a margin-based learner.
This is a trade and not a free choice. The alternative — a soft or ordinal label — was not tried.

## 2.4 The binder panel's stricter cut at pChEMBL ≥ 7.0 — **[DECLARED]**

`train_binders_hybrid.py:57` selects actives at 7.0, 100 nM, not 6.0.

*Why stricter.* The core classifiers learn actives against **measured** inactives. The binder panel
learns actives against negatives that are largely **presumed**, being decoys drawn from background
chemistry. When the negative class is presumed rather than measured, the positive class has to carry
more of the burden of being right, so only high-confidence actives are used.

This is the correct reading of the two cuts, and it is a question to expect. They are not
inconsistent; they belong to two different learning problems.

## 2.5 Classification for eight endpoints, regression for five — **[FORCED]**

D2, A2A, 5-HT2A and SERT are modelled as potency regression rather than binary classification,
because their ChEMBL sets are **96 to 98 per cent active** — only binders get reported for a
well-studied receptor. A binary task on such a set is ill-posed: Matthews correlation ran 0.21 to 0.44
and failed the deployment quality gate. Regression on pChEMBL is the appropriate task for that data
shape. This is not a preference; the data admits nothing else.

---

# Part 3. Processing and representation

## 3.1 Desalting and neutralisation before featurisation — **[CONVENTION]**

Only the largest fragment is kept and charges are neutralised, so a hydrochloride salt and its free
base give the same vector. They are the same molecule for the purposes of a binding prediction and
the deposited form is an artefact of how the compound was crystallised.

## 3.2 ECFP-4 folded to 1,024 bits — **[CONVENTION]**, and untested here

Circular fingerprints of diameter 4 are the default substructure representation for ligand-based
activity prediction, and 1,024 or 2,048 bits are the usual folding widths.

*What was not done.* **No sweep over radius or bit width was run.** A 2,048-bit vector would halve
the collision rate and cost nothing but memory. The applicability-domain measure does use 2,048 bits,
so the project already holds both, and the fact that the model uses the narrower one is a choice
inherited rather than justified by measurement.

## 3.3 The twelve descriptors — **[MEASURED]**, and the result is deflationary

A block ablation over 13 endpoints compared fingerprint alone, descriptors alone, and both together:

| Comparison | Mean | Range | Wins |
|---|---:|---|---:|
| Combined − fingerprint only | **+0.0014** | −0.0104 to +0.0149 | 9 of 13 |
| Combined − descriptors only | **+0.1049** | +0.0088 to +0.2401 | 13 of 13 |

Read it plainly. **The fingerprint carries essentially all of the signal, and the twelve descriptors
add about 0.001.** They are retained because they cost nothing, because they make the vector
interpretable in a way folded bits are not, and because the exposure layer needs them anyway. They
are not retained because they were shown to help, and the manuscript should not imply otherwise.

The converse comparison is the one that justifies the fingerprint: descriptors alone lose a tenth of
an AUROC on average and up to 0.24 on the antioxidant endpoint.

## 3.4 Deduplication on the feature vector, not the identifier — **[FORCED]**

The featuriser is stereo-blind, so stereoisomers, salts and protonation variants give byte-identical
vectors. Neither the SMILES string nor the InChIKey collapses them, and both are unique for every BBB
row. Deduplicating on an identifier would therefore leave copies of one compound on both sides of a
fold, and the model would be scored on rows it had memorised. There is no choice here: the level at
which rows are indistinguishable is the vector, so that is the level at which they must be collapsed.

## 3.5 Label conflicts dropped rather than resolved — **[DECLARED]**

When a duplicate group disagrees on the label, the whole group is discarded rather than majority-voted
or averaged. The same input carrying both labels cannot be learned from, and picking one is an
arbitrary decision dressed as data. Regression takes the group median, which is defensible because a
potency is a quantity and a median of a quantity is meaningful, where a median of a contradiction is
not.

---

# Part 4. Data separation

## 4.1 Scaffold-grouped cross-validation as the primary estimate — **[MEASURED]**

Random 10-fold gives a mean AUROC of **0.9575** across the eight core classifiers; scaffold-grouped
gives **0.9252**. The gap of about 0.032 is the analogue leakage a random split permits.

*Why the harder number is the headline.* The published medicinal-chemistry record is series: a paper
reports twenty analogues of one scaffold. A random split therefore places close analogues of nearly
every test compound in the training set, and answers a question no user ever asks. A user's compound
is novel by construction, which is what a scaffold split simulates.

## 4.2 Ten folds — **[MEASURED]**

The fold count was compared and it barely matters: mean classification AUROC 0.912 at scaffold-5
against 0.958 at random-5 and 0.964 at random-10. **The split type changes the result by about 0.05;
the fold count changes it by about 0.006.** That comparison is the justification for spending the
argument on the split and not on the fold count.

## 4.3 Three disjoint background pools, assigned by hash — **[FORCED]** in principle, **[DECLARED]** in proportion

The pipeline asks the background library for three things: decoys to train against, a sample to set
each threshold as a quantile, and a sample to measure the false-positive rate at that threshold.

*Why they must be disjoint.* If the threshold is the 95th percentile of a sample, the false-positive
rate on that same sample is 5 per cent whatever the model does. It is the quantile restated, not a
measurement. Drawing decoys from either of the others is worse: the model was trained to score those
compounds near zero and is then congratulated for it. This part is forced by logic.

*Why a hash and not a shuffle.* `blake2b(salt + canonical SMILES) mod 100` gives two properties a
shuffle does not: the assignment does not depend on the order the library happens to be in, and
adding compounds later leaves every existing assignment untouched, so a threshold set today stays
comparable with a rate measured next year.

*The 60/20/20 split itself is declared*, not derived. Decoys need the most because they run at three
per active; the other two need only enough for a stable quantile and a stable rate.

## 4.4 The binder panel's own hold-outs — **[DECLARED]**

A fifth of active scaffold groups is withheld from training so sensitivity is measured on
structurally distinct compounds, and measured inactives are split in half, one half training as hard
negatives and the other setting the threshold. Both halves are written to
`models_rf/holdout/<T>_binder_holdout.json` **so a later step cannot reach past the hold-out by
re-reading the endpoint table**.

That last clause is not hypothetical. It is the defect Chapter 5 documents: a downstream script did
re-read the whole table and set the threshold on 40 of 49 deployed endpoints from data the model had
seen. The file-based hold-out exists because the convention alone failed.

## 4.5 The temporal split — **[MEASURED]**, and it needed a control

Models were refitted before a date cutoff and scored on compounds measured afterwards. **The
size-matched random control is the part that makes this an experiment rather than an anecdote**: a
time split trains on less data as well as none of the future, so without a control at the same
training-set size a drop cannot be attributed to either cause.

The control earned its place. In aggregate the time split looks like decay, sensitivity 0.4886
against 0.8715. Read band by band it is mostly composition: **83.3 per cent of the random split's
test compounds are close analogues of its own training set against 16.3 per cent for the time split**.

---

# Part 5. Training

## 5.1 A random forest as the primary estimator — **[MEASURED]**, and the result is honest rather than flattering

Five families were compared on the same features and the same folds: random forest, XGBoost,
histogram gradient boosting, L2-regularised logistic regression, and a five-nearest-neighbour
Tanimoto read-across.

*Against the baselines the forest is decisively better.* Over the 8 classification endpoints under
scaffold split it beats the read-across on 8 of 8, median +0.0381, p = 0.00781, and logistic
regression on 8 of 8, median +0.0818.

*Against the other ensembles it is not distinguishable.* Over all 13 endpoints it beats histogram
gradient boosting on 8, median +0.0026, **p = 0.89258**, and XGBoost on 8, median +0.0033,
**p = 0.73535**. On the five regression endpoints it loses to XGBoost on 5 of 5, which the test
cannot certify: at five pairs the signed-rank floor is 0.0625, so no result there can reach 0.05.

*So why the forest.* Not accuracy. It is kept for its calibration behaviour, its out-of-bag
structure, and because an untuned forest is a more honest default than a tuned booster. **A candidate
should say this rather than claim the forest won.**

## 5.2 Hyper-parameters — **[DECLARED]**, every one of them

`n_estimators=300, min_samples_leaf=2` (4 on the binder panel), `class_weight="balanced"`,
`random_state=42`. XGBoost at 400 estimators, depth 6, learning rate 0.05; histogram gradient boosting
at 400 iterations, learning rate 0.06.

**None of these was searched.** Section 9 argues why that is defensible and what it forfeits.

The reasoning behind each is real even though it is not empirical here. 300 trees because the variance
term that shrinks with tree count is negligible well before 300, while the correlation term does not
shrink at all, so more trees buy nothing. Leaf size 4 rather than 2 on the binder panel because those
tables are smaller and noisier and a leaf holding one compound is memorisation. A fixed seed so a
rerun reproduces, which is what allowed a corrected artefact to be compared against its predecessor
and the difference attributed to the correction rather than to noise.

## 5.3 Balanced class weights — **[DECLARED]**, and applied inconsistently

Several endpoints are heavily imbalanced and without reweighting a tree gains more impurity reduction
from the majority class.

**An inconsistency that was here and has been fixed.** Histogram gradient boosting was the only
ensemble competing unweighted, while the forest and logistic regression use
`class_weight="balanced"` and XGBoost uses `scale_pos_weight`. On endpoints whose actives outnumber
inactives four to one that is not a neutral difference, and a conclusion of "the forest is not
distinguishable from the boosters" could not rest on one booster having been handicapped.

It now carries `class_weight="balanced"` and the comparison was refitted. The correction is real and
its effect is not: only the histogram-boosting rows moved, by between -0.0013 and +0.0023, and every
verdict is unchanged. The forest is still not distinguishable from either booster over all 13
endpoints. That is a better position than before, because the conclusion now survives a comparison
that is genuinely like-for-like rather than one that flattered it.

## 5.4 Hybrid negatives: measured inactives **and** decoys — **[MEASURED]**

*The failure this fixes.* Decoy-only training saturates. For melatonin MT1 the measured inactives
scored a median binder probability of **0.972** against an active median of **0.997**: the decision
boundary sat where tiny probability differences swing the result, and sensitivity at a controlled
false-positive rate collapsed to 0.52 at MT1 and 0.33 at KEAP1.

The negatives were structurally dissimilar by construction, so the model never learned what a hard
negative looks like. The fix is to give it the ones it never saw: measured inactives split in half,
one half training alongside the decoys.

This is section 1.7's principle applied. Easy negatives inflate a metric and teach nothing; hard
negatives are what a decision boundary needs.

## 5.5 Three decoys per active, at Tanimoto below 0.35 — **[CONVENTION]**

Decoys are drawn only from the decoy pool, three per active, each with maximum Tanimoto below 0.35 to
every active, and the target's own measured inactives are excluded from decoy eligibility so a
compound cannot be both a trained negative and an unseen one.

*The ratio and the ceiling are conventions and were not swept.* Both directions have a cost worth
stating: a lower ratio leaves the negative class thin, a higher one drives the base rate down and the
model toward always answering no; a higher similarity ceiling risks admitting genuine unmeasured
actives as decoys, a lower one makes the negatives too easy, which is exactly the failure of
section 1.7.

---

# Part 6. Calibration

## 6.1 Calibrating at all — **[MEASURED]**

Random-forest votes cluster toward the middle and a raw 0.72 is not a 72 per cent chance. Across the
eight core classifiers, isotonic calibration takes expected calibration error from a mean of
**0.0801 to 0.0147**.

## 6.2 Isotonic on the core, Platt on the panel — **[DECLARED]**

Isotonic regression is non-parametric, so it can repair a badly shaped calibration curve, and it is
fitted on **out-of-fold** predictions with an inner 5-fold because a non-parametric calibrator fitted
on training scores would overfit. Platt scaling has two parameters, cannot repair a non-sigmoidal
distortion, and is far more stable on the small calibration sets the binder endpoints have.

*The consequence, which belongs in the same breath.* The binder endpoints that carry a measured
calibration error average **0.0762**, about five times the core figure, and **nine of the 47 carry no
measured calibration error at all**.

---

# Part 7. Thresholds and the reliability gate

## 7.1 Two constraints, the stricter binding — **[MEASURED]**

Each binder threshold is the stricter of: at most **0.10** of the target's held-out measured inactives
called a binder, and at most **0.05** of the threshold pool.

*Why the second constraint exists*, and the example is the justification: a prospective test found
Nav1.1 assigning glucose a binder probability of 0.806 while scoring 0.979 against its own assay's
negatives. **A model judged only on its own chemistry can be confidently wrong about everything
else.** Nav1.1 was withdrawn on that evidence.

## 7.2 The rate targets themselves — **[DECLARED]**

0.10, 0.05 and the screening set's 0.01 are design choices about how much false positivity a user
should meet, not measurements. What is measured is the price: panel mean sensitivity **0.7638**,
median 0.835, from 0.303 to 0.993.

*The trade is stated in the direction that does not flatter.* Higher sensitivity is available at any
time by relaxing to a background-only threshold. For three endpoints the binding constraint is
separation from compounds an experiment has already reported as non-binders, so relaxing would mean
calling those compounds binders. Higher sensitivity, less truthful server.

## 7.3 The reliability gate at 0.50 sensitivity and 0.75 AUROC — **[DECLARED]**

Six of 47 deployed endpoints fail it and carry a flag on every result. The thresholds are choices.
What justifies the arrangement is not the numbers but that the gate now lives in **one** place,
`panel.py`, imported by five writers, the figure and the app. It previously lived in six places and
disagreed with itself in three, which meant the About page named ten endpoints while six carried the
marker.

---

# Part 8. Testing

## 8.1 Reporting both splits — **[DECLARED]**

Both random and scaffold figures are reported for every endpoint, always. Reporting only the scaffold
number would look rigorous and would hide the size of the analogue effect, which is itself a finding.

## 8.2 The permutation null — **[MEASURED]**

With labels permuted, the same pipeline on the same folds returns **0.4959** random and **0.5026**
scaffold, all sixteen values within **0.0200** of chance.

*Why it is not a formality.* It tests two things at once: that the cross-validated figures are not
obtainable without the labels, and that scaffold grouping alone confers no advantage. The second is
the specific route a reviewer suspects, because whole scaffold classes could in principle carry enough
class-frequency information for a label-free model to exploit the grouping. They do not.

## 8.3 External validation only for the barrier model — **[FORCED]**

There is one genuine external set: FDA-curated approved drugs absent from B3DB, AUROC **0.7934** on
the 241 also novel in feature space. For the target panel no such set exists, and the reason is
structural rather than lazy: for most of those proteins the public record *is* the training set.
Claiming an external validation for the panel would require inventing one.

*Why the 241 and not the 306.* Excluding overlap by InChIKey is not enough, because the InChIKey
separates stereoisomers and salts and the featuriser does not. The 65 compounds that are
feature-identical to a training compound score **0.7102**; they are memorisation, and the difference
between 0.7645 and 0.7934 is the size of it.

## 8.4 The falsification suite — **[DECLARED]**, and it is the methodological core

Ten hypotheses, each stated so that it could fail and each paired with a null model capable of
producing the same apparent success by accident. **Four were refuted and two weakened.**

*Why this is a criterion and not a result.* Every other number in the project was produced by people
who wanted the tool to work. The suite is the only part designed by someone trying to break it, it is
read-only by construction, and acting on a finding is a separate commit so a wording change can never
be mistaken for evidence. A suite in which nothing is ever refuted is a suite of tests that cannot
fail.

---

# Part 9. The criteria with no evidence behind them

Collected here rather than left scattered, because a reader is entitled to see them together.

| Criterion | Value | Status |
|---|---|---|
| Fingerprint radius and width | ECFP-4, 1,024 bits | Convention; no sweep |
| Trees per forest | 300 | Declared; no sweep |
| Minimum leaf size | 2 core, 4 binder | Declared; no sweep |
| Boosting settings | 400 estimators, lr 0.05 / 0.06, depth 6 | Declared; no sweep |
| Decoy ratio | 3 per active | Convention; no sweep |
| Decoy similarity ceiling | Tanimoto < 0.35 | Convention; no sweep |
| Active hold-out share | 0.20 of scaffold groups | Declared |
| Background pool shares | 60 / 20 / 20 | Declared |
| False-positive targets | 0.10, 0.05, 0.01 | Declared |
| Reliability gate | sensitivity 0.50, AUROC 0.75 | Declared |
| Reporting threshold | gated score 0.30 | Declared |
| Conformal error rate | ε = 0.10 | Declared |
| Applicability bands | 0.30 and 0.50 | Convention |
| Grey zone | pChEMBL 5 to 6 discarded | Declared |

## The one that matters most: nothing was tuned

**No grid search, no randomised search, no sweep over any hyper-parameter exists anywhere in this
project.** Every setting above is a default or a declared choice.

*What that forfeits.* No claim can be made that these settings are optimal, or even good. A tuned
gradient booster would very likely beat this forest on several endpoints, and section 5.1 already
shows XGBoost beating it on all five regression endpoints untuned.

*What that buys, and it is not small.* There is **no tuning-set leakage**, because there was no tuning
set. The commonest way a cross-validated figure is inflated in this field is selecting
hyper-parameters on the same folds the final number is reported from, and that failure mode is
structurally absent here. Combined with the permutation null and the scaffold split, the
cross-validated figures are about as free of optimistic bias as this design allows.

*The honest summary sentence.* This project chose to spend its effort on the correctness of the
evaluation rather than on the performance of the estimator, and the results should be read that way:
as a defensible measurement of an untuned model, not as a demonstration of the best achievable
accuracy.

---

## Outstanding items for this document

1. **A bit-width comparison is nearly free and would close section 3.2.** The 2,048-bit fingerprint is
   already computed for the applicability domain, so scoring the panel at both widths on the existing
   folds needs no new data.
2. **The decoy ratio and similarity ceiling have never been varied.** A sweep of either would either
   justify the current values or improve them, and section 5.5 currently defends both by argument
   alone.
3. ~~Histogram gradient boosting runs unweighted in the family comparison~~ **Closed.** It now uses
   `class_weight="balanced"` and the comparison was refitted; only its own rows moved, by at most
   0.0023, and no verdict changed. Section 5.3 records both the correction and its size.
