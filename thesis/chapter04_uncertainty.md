# Chapter 4. Calibration, conformal prediction and the applicability domain

> Every quantitative statement in this chapter was computed during this session from the artefact
> named beside it. Two of those artefacts, `applicability_bbb_validation.csv` and
> `applicability_coverage.csv`, were found stale during the writing of this chapter and were
> regenerated before any figure here was taken from them; section 4.12 records what changed.
> Citation numbers refer to `manuscript/references.md`.

---

## 4.1 Three different questions, and why one answer will not do

A prediction of 0.84 is not yet useful. Three separate questions stand between it and a decision, and
this chapter exists because they are genuinely different questions with genuinely different answers.

**Is 0.84 a probability?** A random forest returns a vote share, not a probability. Averaging the
votes of 300 trees produces a number that behaves like a probability, clusters toward the middle, and
is not one. Calibration is the map that makes a stated 0.84 correspond to being right 84 per cent of
the time. It is a statement about a population of predictions, not about this compound.

**How confident is the model about *this* compound?** Calibration cannot answer this. A well
calibrated model can be well calibrated and uniformly uncertain. Conformal prediction converts the
question into a coverage guarantee that is measured rather than assumed.

**Has the model seen anything like this compound?** Neither of the above answers this either. Both
are computed on held-out chemistry drawn from the same distribution as the training set. A compound
from outside that distribution can receive a confident, well-calibrated, tightly-covered prediction
that is worthless. The applicability domain is the third instrument, and section 4.11 shows it
answers a different question from the one it was built for.

The three are reported together for every value the server returns, and a reader who takes only one
of them has not taken the honest summary.

## 4.2 Calibration: the construction

Core classifiers are calibrated by **isotonic regression** fitted on **out-of-fold predictions**. The
second half of that sentence carries the weight. Fitting a calibrator on the same predictions the
model made about its own training data would produce a calibrator that corrects an optimism the
deployed model does not have, and the reported improvement would be an artefact of the procedure. The
out-of-fold construction means no compound contributes to the calibrator that scores it.

Isotonic regression fits a free monotonic step function by pool-adjacent-violators, minimising
squared error subject only to being non-decreasing. Two consequences follow from monotonicity alone
and both matter.

Because the map cannot reorder compounds, **AUROC is invariant under calibration**. Every
discrimination figure in this thesis is therefore unaffected by anything in this section, which is
why calibration can be reported as a separate axis rather than folded into the headline.

And because the map is free rather than parametric, it can correct an arbitrary distortion, at the
cost of needing enough held-out data to fit a step function without overfitting it. That cost is why
the binder panel does something else, and section 4.4 reports what it costs.

The **binder models use Platt scaling** instead, a two-parameter logistic fit by maximum likelihood
on a stratified fifth withheld from the forest's own fit. A logistic curve is a far more constrained
object than a free step function. It is used because the positive class at one binder endpoint is
often too small for isotonic regression to be anything but overfitting.

## 4.3 What calibration achieves on the core classifiers

Measured over the eight measured-label classifiers, on out-of-fold predictions
(`results/tables/calibration.csv`):

| Endpoint | Brier raw | Brier calibrated | ECE raw | ECE calibrated |
|---|---:|---:|---:|---:|
| BBB | 0.1289 | **0.1358** | 0.0657 | 0.0412 |
| AChE | 0.0810 | 0.0702 | 0.0942 | 0.0110 |
| BChE | 0.0708 | 0.0647 | 0.0735 | 0.0127 |
| BACE1 | 0.0397 | 0.0358 | 0.0497 | **0.0049** |
| GSK-3β | 0.0734 | 0.0669 | 0.0847 | 0.0175 |
| MAO-A | 0.0709 | 0.0609 | 0.0940 | 0.0089 |
| MAO-B | 0.0848 | 0.0770 | 0.0881 | 0.0148 |
| hERG | 0.0717 | 0.0632 | 0.0911 | 0.0069 |
| **Mean** | 0.0801 | | **0.0801** | **0.0147** |

Expected calibration error falls from **0.0801 to 0.0147**, a factor of 5.4, and it falls on **all
eight** endpoints. Per endpoint the calibrated error runs from **0.0049 at BACE1 to 0.0412 at the
barrier model**, which is worst by a factor of eight over the best and remains the panel's weakest
calibration.

**One endpoint gets worse on the other metric, and it is worth stating rather than burying.** The
Brier score improves on seven of the eight and **degrades at the barrier model, from 0.1289 to
0.1358**. That is not a contradiction: the Brier score is a proper scoring rule combining calibration
and discrimination, so a monotone map that improves the calibration component can still cost more in
the sharpness component than it gains. The honest reading is that isotonic calibration buys a real
improvement in the barrier model's reliability at a small cost in its resolution, and that the
barrier endpoint is the one where the calibration machinery works least well by both measures.

## 4.4 What the headline figure does not cover

The figure 0.0801 to 0.0147 appears in the manuscript abstract, in the technical report, in the
evidence map and in Chapter 1 of this thesis. **It describes eight of the seventy deployed
estimators.** No document in the project says so, and it should.

The other 38 scored targets are Platt-scaled, and their calibration is measured separately, on
held-out compounds (`results/tables/integrity_calibration_per_target.csv`):

| | Core classifiers, isotonic | Binder panel, Platt |
|---|---:|---:|
| Endpoints | 8 | 38 |
| Mean ECE | **0.0147** | **0.0762** |
| Median ECE | 0.0118 | 0.0675 |
| Best | 0.0049, BACE1 | 0.030, H3 |
| Worst | 0.0412, BBB | **0.178, TAAR1** |
| Above 0.10 | 0 of 8 | **8 of 38** |

The binder panel's calibrated expected calibration error, at 0.0762, is **five times the core
classifiers' 0.0147 and almost exactly the core classifiers' figure before calibration, 0.0801**. On
this evidence Platt scaling on these endpoints buys very little.

Three things follow, and the chapter states all three rather than the convenient one.

The **honest panel-wide statement** is not that expected calibration error is 0.0147. It is that it
is 0.0147 on the eight core classifiers and 0.0762 across the 38 binder endpoints, and a user reading
a binder probability is reading the second number, not the first.

The **cause is diagnosable and is the one the design anticipated**. Platt scaling was chosen because
the withheld positive class at one binder endpoint is often too small to fit an isotonic step
function. TAAR1, the worst calibrated at 0.178, has 403 rows in its endpoint table and 82 in the
fitted binder set. The constraint was real; the price is now measured.

The **remedy is specific**. The endpoints with enough held-out positives to support isotonic
calibration should have it, and the rest should carry their measured calibration error beside the
probability rather than inheriting a headline figure computed on other endpoints. Neither is done,
and this is the most actionable piece of future work in the chapter.

## 4.5 Conformal prediction: the construction

A calibrated probability describes a population. **Inductive Mondrian conformal prediction** [14]
describes a compound, and it does so with a guarantee that is measured on held-out data rather than
assumed.

At significance ε = 0.10, on a calibration set disjoint from training, the nonconformity of a
labelled compound is α_i = 1 − q̂(y_i | x_i), the calibrated probability the model assigned to the
label that turned out to be correct. For class c with n_c calibration points, the threshold τ_c is
the k-th smallest such α within that class, where k = ⌈(n_c + 1)(1 − ε)⌉. The prediction set for a
new compound is every class whose nonconformity does not exceed its own threshold:

Γ(x) = { c : 1 − q̂(c | x) ≤ τ_c }.

Two features of that construction do real work.

**Mondrian, meaning class-conditional.** Thresholds are computed within each class rather than
pooled. Under exchangeability this delivers coverage of at least 1 − ε *within each class*, which is
what the pooled form does not. The distinction is not academic on this data: Chapter 2 records a
median active fraction of 0.825 and endpoints reaching 0.962, and on such a set marginal coverage of
0.90 can be satisfied while the rare class fails almost completely. The class-conditional form makes
that failure visible instead of averaging it away.

**The prediction set can be empty or contain both classes.** An empty set means the compound conforms
to neither class: it is unlike the actives and unlike the inactives. A set containing both means it
separates neither. Both are informative, both are displayed, and neither is available from a point
probability. This is the specific respect in which the conformal layer says something calibration
cannot.

## 4.6 What coverage the panel actually achieves

Measured on 10,000 held-out compounds across the eight core classifiers
(`results/tables/rf_conformal.csv`):

| Endpoint | Test compounds | Target | Empirical coverage | Mean set size |
|---|---:|---:|---:|---:|
| BBB | 1,561 | 0.90 | 0.904 | 1.013 |
| AChE | 1,064 | 0.90 | **0.921** | 1.019 |
| BChE | 678 | 0.90 | 0.917 | 1.009 |
| BACE1 | 1,793 | 0.90 | 0.906 | 1.007 |
| GSK-3β | 1,128 | 0.90 | 0.908 | 1.030 |
| MAO-A | 758 | 0.90 | **0.889** | 1.012 |
| MAO-B | 962 | 0.90 | 0.899 | 1.025 |
| hERG | 2,056 | 0.90 | 0.904 | 1.079 |

Empirical coverage runs **0.889 to 0.921** against a 0.90 target. **Six of the eight sit at or above
target.** The two that fall short are the monoamine oxidases: MAO-B at 0.899, one thousandth below,
and MAO-A at 0.889, eleven thousandths below on 758 test compounds. Neither shortfall is larger than
the sampling error at those set sizes, and the direction is worth noting rather than explaining away,
because MAO-A and MAO-B are also the two endpoints Chapter 8 records as degrading most under a
scaffold split. Conformal coverage is a guarantee under exchangeability, and these are the endpoints
whose chemistry is least scaffold-diverse.

The uncertainty statement is nonetheless honest in the aggregate: when the server says the truth is in
the set nine times in ten, it is, to within about one part in a hundred.

The **mean set size runs 1.007 to 1.079** on a two-class problem, where 1.0 would mean every
prediction is a confident single label and 2.0 that none is. Between 0.7 and 7.9 per cent of
compounds therefore receive an ambiguous or empty set, hERG being the most equivocal at 1.079. Small
average set sizes are what make the guarantee useful rather than vacuous: a predictor that always
returned both classes would achieve perfect coverage and say nothing.

Two limitations travel with these figures. They cover the **same eight core classifiers** as section
4.3, so the binder panel has no measured coverage statement at all. And the guarantee holds **under
exchangeability**, which is exactly the assumption a compound from a new chemical series violates.
Chapter 8 shows what happens when that assumption fails, and the answer is that recall degrades in a
way the coverage statement does not predict.

## 4.7 The applicability domain: definition and bands

The applicability domain measure is the **maximum Tanimoto similarity of the query to that endpoint's
own measured chemistry**, computed on 2,048-bit Morgan fingerprints, reported with the nearest
measured analogue and its structure so a user can look at it.

The interface reports three bands:

- **In domain**, T ≥ 0.50: the query sits among measured chemistry.
- **Near domain**, 0.30 ≤ T < 0.50.
- **Out of domain**, T < 0.30.

**A naming inconsistency worth fixing.** The evaluation scripts use a single cut, `AD_THRESHOLD =
0.30`, and the artefacts they write label everything above it `in_domain`. So
`applicability_coverage.csv` reporting a barrier "in-domain fraction" of 0.7171 means *in or near
domain* in the interface's language. A reader comparing that figure with the green "In domain" badge
the server shows is comparing two different things.

## 4.8 Choosing the measure

Four candidate measures were compared on their ability to separate chemistry genuinely absent from
the reference library from approved drugs the models had not seen
(`results/tables/applicability_measures.csv`):

| Measure | Median, unseen drugs | Median, non-drug-like | p | Separates | Threshold at 10% drug loss | Distant chemistry caught |
|---|---:|---:|---:|---|---:|---:|
| **Maximum similarity** | 0.5909 | 0.4606 | 0.0072 | yes | 0.3462 | 0.375 |
| Mean of top 5 | 0.5224 | 0.3781 | 0.0075 | yes | 0.3179 | 0.375 |
| 5th nearest neighbour | 0.4444 | 0.3287 | 0.0171 | no | 0.2921 | 0.500 |
| Density within 0.4 | 8.0 | 1.5 | 0.0536 | no | 0.0 | 0.375 |

Two of the four separate at the stated criterion of p < 0.01, and **no alternative beats the deployed
one**. The maximum has the largest separation and the smallest p-value, so the simplest measure is
also the best available here. That is a mildly deflating result and it is the right one to report:
the choice was not tuned, it was tested against three alternatives and won.

The fifth-nearest-neighbour measure catches half of distant chemistry against the maximum's 0.375,
which looks better until one notices it fails the separation criterion. Catching more at a threshold
chosen post hoc is not the same as discriminating.

## 4.9 What the flag was built for, and does badly

The adversarial check asks whether the flag rates chemistry genuinely absent from the reference
library worse than approved drugs it has never seen. It **passes**
(`results/tables/inversion_validation.csv`): median maximum similarity 0.59 for unseen drugs against
0.47 for non-drug-like chemistry, over 25 controls, one-sided Mann-Whitney p = 1.11 × 10⁻³.

Chapter 9 records the history in full, because the check previously failed and the fix was to the
controls rather than to the criterion: 28 of the original controls were measured compounds inside the
flag's own reference library, where calling them in domain is truthful. The passing criterion was not
moved.

**Passing is not the same as being useful, and the gap is large.** At a threshold that rejects a
tenth of genuine drugs, the flag catches **0.375** of genuinely distant chemistry. At the deployed
cut of 0.30 it catches **0.20**. Four out of five non-drug-like structures are waved through. The
flag separates in the aggregate and discriminates poorly on the individual case, and no amount of
restating the p-value changes that.

The technical report's conclusion is right and is worth repeating here: the conformal interval and
the nearest-analogue distance are the stronger statements of confidence, and the interface says so.

## 4.10 What the flag actually predicts

The applicability domain was built to detect non-drug-like chemistry, which it does weakly. It turns
out to be a good predictor of something else entirely, and this is the most useful result in the
chapter.

**Recall is a monotone function of the domain distance.** Across three test sets built by unrelated
rules, withheld by publication date, at random, and by curator
(`results/tables/external_novelty_strata.csv`), recall against maximum Tanimoto to the training
actives runs:

| Band | By date | At random | By curator |
|---|---:|---:|---:|
| below 0.40, different chemotype | 0.1612 | 0.1179 | 0.0494 |
| 0.40 to 0.55, related series | 0.5532 | 0.5232 | 0.4641 |
| 0.55 to 0.70, same series | 0.7403 | 0.7683 | 0.8264 |
| 0.70 and above, close analogue | 0.8616 | 0.9334 | 0.9011 |

The quantity the server already computes for every query, and already displays, predicts the recall
that query will receive. A compound the flag places near the edge of the domain is precisely the
compound whose activity the panel is most likely to miss. The consequence is the one Chapter 1 called
the third absent capability: **the expected recall for the compound in hand is knowable at the moment
of the query**, not only in aggregate, and silence can be reported with a number attached.

So the flag earns its place in a role it was not designed for and was not being credited with. It is
a poor detector of alien chemistry and a good predictor of when the panel will fall silent on a real
active.

## 4.11 The domain flag does not predict discrimination

The natural next assumption is that a compound outside the domain also receives a worse-ranked
prediction. On the one external set where this can be tested, **it does not**
(`results/tables/applicability_bbb_validation.csv`, regenerated for this chapter):

| Subset | n | Mean max Tanimoto | AUROC | Accuracy |
|---|---:|---:|---:|---:|
| In domain, T ≥ 0.30 | 258 | 0.646 | **0.7597** | 0.7326 |
| Out of domain, T < 0.30 | 48 | 0.248 | **0.8059** | 0.7292 |
| All | 306 | 0.584 | 0.7645 | 0.7320 |

The barrier model **ranks the 48 out-of-domain compounds better than the 258 in-domain ones**, by
0.046 AUROC, while accuracy is flat across the two. The result should not be over-read: 48 compounds
is a small set, no interval is computed here, and the two subsets are not matched on anything but
distance. But the direction is the opposite of the assumption, and the assumption is the one a
referee will bring.

The reconciliation with section 4.10 is the point of reporting both. **Recall and ranking are
different quantities.** The domain distance predicts recall, which is a property of where the
decision threshold sits relative to the score distribution of true actives. It does not predict
AUROC, which is threshold-free. A model can rank a distant compound correctly relative to other
distant compounds while placing all of them below a threshold calibrated on near chemistry. That is
precisely the behaviour Chapter 8 documents: on novel chemistry the panel does not become wrong, it
becomes quiet.

The practical statement for a user is therefore narrower than the flag's name suggests. A low domain
distance means *this endpoint is likely to miss a real activity here*. It does not mean *this
endpoint's prediction is unreliable in the ranking sense*, and the two should not be conflated.

## 4.12 How much of drug space is inside the domain

The final quantity is coverage, and it is sobering. Scoring the 11,723 DrugBank approved and
investigational drugs against each endpoint's measured chemistry
(`results/tables/applicability_coverage.csv`, regenerated for this chapter) gives the fraction at or
above the 0.30 cut, and the median distance:

| Endpoint | Fraction at T ≥ 0.30 | Median max Tanimoto |
|---|---:|---:|
| BBB | **0.7171** | 0.393 |
| hERG | 0.4830 | 0.296 |
| AChE | 0.4042 | 0.280 |
| BChE | 0.3423 | 0.267 |
| GSK-3β | 0.3293 | 0.266 |
| MAO-B | 0.3245 | 0.264 |
| MAO-A | 0.3064 | 0.258 |
| BACE1 | **0.2189** | 0.250 |

For every endpoint except the barrier model, **most approved drugs sit below the 0.30 cut**, and at
BACE1 nearly four in five do. The median approved drug is at Tanimoto 0.25 to 0.30 from the nearest
compound ever measured at these targets. Combined with section 4.10, that says the panel will miss a
substantial fraction of genuine activities among approved drugs at the target endpoints, and it is
the quantitative basis for the claim in Chapter 1 that silence is the weakest of the system's
outputs.

The barrier model is the exception because its training set, B3DB augmented with FDA-curated approved
drugs, is drawn from approved-drug chemistry in the first place. Its 0.7171 is a statement about the
overlap of two drug collections, not evidence that the barrier model generalises better.

**A note on provenance.** Both tables in sections 4.11 and 4.12 were built on 13 August against the
pre-retrain barrier model and had never been rebuilt. The stale
`applicability_bbb_validation.csv` reported AUROC 0.7608 on the 306 approved drugs where
`external_bbb_validation.csv`, generated from the same model on the same compounds, reported 0.7645.
Two artefacts described one test set and disagreed, and nothing complained, because neither was
declared in the freshness graph. They were regenerated for this chapter, now agree exactly, and have
been declared. This is the third artefact in the project found stale for want of a declaration rather
than for want of a rebuild.

## 4.13 What a user should read, in order

The three instruments answer three questions and a user who takes one has taken the wrong summary.

**Read the domain distance first.** It tells you whether the rest is worth reading, and it tells you
the recall the panel achieves at that distance. Below 0.40, recall is about 0.16 and a silence means
almost nothing.

**Read the conformal set second.** An empty set or a set containing both classes is a statement that
a point probability cannot make, and it is available for the eight core classifiers only.

**Read the calibrated probability last, and read it against its base rate.** It is well calibrated
for the eight core classifiers, at 0.0147 expected calibration error, and about five times less well
calibrated for the 38 binder endpoints, at 0.0762. Chapter 7 explains why the enrichment over the
base rate, rather than the probability itself, is the quantity the disease layer consumes.

---

## Outstanding items for this chapter

1. **The headline calibration figure should be qualified everywhere it appears.** The manuscript
   abstract, the technical report, the evidence map and Chapter 1 of this thesis all quote 0.0801 to
   0.0147 without saying it covers eight of the seventy deployed estimators. The binder panel's
   figure of 0.0762 belongs beside it.
2. **The binder endpoints have no conformal coverage statement.** The uncertainty stack is complete
   for the eight core classifiers and two-thirds absent for the 38 targets a user is most likely to
   query.
3. **Isotonic calibration should be extended to those binder endpoints with enough held-out
   positives to support it**, and the rest should carry their measured per-endpoint calibration
   error rather than inheriting a headline computed elsewhere.
4. **`in_domain` in the artefacts means T ≥ 0.30, which is the interface's "in or near domain".**
   The naming should be reconciled in one direction or the other.
5. The out-of-domain result in section 4.11 rests on 48 compounds and carries no interval. A Wilson
   or bootstrap interval on both subsets would settle whether the reversal is real or noise, and is
   cheap.
