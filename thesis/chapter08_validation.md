# Chapter 8. Validation, external and prospective, and the composition finding

> Every quantitative statement in this chapter was computed during this session from the artefact
> named beside it. Six of the artefacts used here predate the base-rate correction of 25 August
> described in Chapter 7; each was checked and none depends on the base rates, because all of them
> score model probabilities rather than engagement enrichment.

---

## 8.1 A programme, not a number

Validation here is arranged so that each test answers a question the previous one cannot. A single
figure, however large, would answer only the easiest of them.

| Test | The question it answers | The question it cannot |
|---|---|---|
| Random cross-validation | can the model interpolate within known chemistry? | anything about new chemistry |
| Scaffold-grouped cross-validation | can it generalise to a withheld structural class? | whether the class was withheld realistically |
| Temporal refit | what would it have said about chemistry that did not exist? | whether a drop is the wall or the smaller training set |
| Size-matched random control | which of those two it was | the composition of the test set |
| Novelty stratification | whether the effect is date or distance | absolute performance on distant chemistry |
| Cross-provenance refit | does the signal survive a change of curator? | anything, on more than three endpoints |
| Non-CNS specificity | does it stay quiet when it should? | behaviour on proven inactives |
| Adversarial suite | can any check fail? | what nobody thought to check |

The chapter is organised in that order, and section 8.11 states what the whole programme still does
not establish.

## 8.2 Cross-validation, two regimes

Across the eight measured-label classifiers (`results/tables/rf_cv_summary.csv`), mean AUROC is
**0.958** under a random 10-fold split, ranging 0.899 to 0.976, and **0.925** under a
scaffold-grouped split that withholds entire Bemis-Murcko classes, ranging 0.878 to 0.965.

The distance between the two is the honest statement of how far a model travels, and Chapter 3
records what it costs. The scaffold column is the one this thesis quotes.

## 8.3 The one genuine external set

For the barrier model an external set exists: FDA-curated approved drugs absent from B3DB, scored by
the deployed model (`results/tables/external_bbb_validation.csv`).

| Set | n | Permeable | AUROC | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|
| Not in B3DB by InChIKey | 306 | 203 | 0.7645 | 0.7734 | 0.6505 |
| Also distinguishable from training in feature space | 241 | 171 | **0.7934** | 0.7661 | 0.6571 |
| Of which feature-identical to a training compound | 65 | 32 | 0.7102 | 0.8125 | 0.6364 |

**The second row is the one that supports an external claim.** Excluding overlap by InChIKey is not
enough to call a compound unseen: the InChIKey separates stereoisomers, salts and protonation states,
and the featuriser does not, so a compound can pass that check and still be one the model has
memorised. The third row is the memorisation the first contains, and reporting all three is what lets
a reader see the difference.

An AUROC of 0.793 on genuinely unseen approved drugs is a substantial drop from the 0.878 the
scaffold split reports, and it is the most honest single number in this thesis about the barrier
model.

## 8.4 Why the target panel has no external set

For the target panel no external set of comparable quality exists, and the reason is structural
rather than negligent: **for most of these targets the public measured chemistry is the training
set.** Buying an independent set would mean running new assays.

Two kinds of independence can still be constructed from what exists, and both were, with every model
refitted rather than merely scored. Scoring the deployed models on any of it would have measured only
what they had already been fitted to.

## 8.5 The prospective simulation

Every measured row carries the year of the document it came from. Freezing the data at a cutoff,
fitting on what was known then, and testing on compounds first published afterwards reproduces the
position a user is actually in.

Four things make it a real test rather than a comfortable one.

**The models are refitted**, not scored. **The operating point is frozen too**: thresholds are
quantiles of pre-cutoff held-out measured inactives, so the sensitivity reported is what a user in the
cutoff year would have got. **Decoys are matched to pre-cutoff actives only**, and every post-cutoff
compound is barred from decoy eligibility, so a test compound can never appear as a training negative.
And each endpoint is fitted a second time on a **size-matched random split**, because a time split
trains on less data as well as none of the future, and without a control at the same n a drop cannot
be attributed to either.

**39 of the 47 deployed endpoints** carry enough dated chemistry on both sides of the wall to be
assessed, giving **15,854 post-cutoff test actives**. The other eight are excluded for want of
post-cutoff compounds, never for a poor result.

## 8.6 The aggregate result, which looks like decay

| | Time split | Size-matched random | Difference |
|---|---:|---:|---:|
| AUROC against measured inactives (n = 32) | **0.8231** | 0.9507 | **-0.1276** |
| AUROC against background chemistry (n = 39) | 0.9131 | 0.9752 | -0.0621 |
| Sensitivity at the frozen threshold (n = 39) | **0.4886** | 0.8715 | **-0.3829** |
| False-positive rate on background (n = 39) | 0.0373 | 0.0376 | -0.0003 |

Two rows are reassuring and two are not, and the difference between them is the whole result.

**Specificity transfers essentially unchanged**, 0.0373 against 0.0376. Whatever a temporal wall
costs, it does not make the panel reckless. **Ranking against background chemistry largely survives**,
falling from 0.975 to 0.913.

**Ranking against measured inactives at the same target falls by 0.128**, and **sensitivity at the
frozen threshold falls by 0.383**. Read at face value those two rows say the deployed figures overstate
what a user submitting new chemistry should expect, and the first draft of this analysis said exactly
that.

It is the wrong reading.

## 8.7 The composition finding

A time split flatters nothing; it changes who is being tested. Medicinal chemistry publishes analogue
series, so if the compounds appearing after a cutoff are close relatives of those before it, a high
prospective score demonstrates interpolation within a series. That was the reason for measuring, for
every test compound, its maximum Tanimoto similarity to the training actives of its own split, and
recomputing performance within bands of that distance. **Both classes are binned**, so each band
compares actives against measured inactives at a comparable distance; binning only the actives would
produce bands whose positives are novel and whose negatives are not, and their AUROC would measure
that mismatch rather than the model.

The stratification was built to check whether the time split was too easy. It showed something else.

**The two splits are not testing comparable populations.**

| Band | Time split | Random split |
|---|---:|---:|
| below 0.40, different chemotype | **27.5%** | **1.3%** |
| 0.40 to 0.55, related series | 33.0% | 2.9% |
| 0.55 to 0.70, same series | 23.2% | 12.5% |
| 0.70 and above, close analogue | 16.3% | **83.3%** |

A random split of medicinal-chemistry data is barely a test: **83.3 per cent of its test compounds
are close analogues of something in its own training set**, because the published record is series and
a random draw keeps most of a series on both sides. The time split's test set is genuinely different
chemistry, with a median similarity to training of 0.479 against 0.800, and 4,366 actives below
Tanimoto 0.40 against 212.

**Read band by band, most of the gap disappears.**

| Band | AUROC, time | AUROC, random | Recall, time | Recall, random |
|---|---:|---:|---:|---:|
| below 0.40 | 0.6711 | 0.7516 | 0.1612 | 0.1179 |
| 0.40 to 0.55 | 0.7645 | 0.8217 | 0.5532 | 0.5232 |
| 0.55 to 0.70 | 0.7871 | 0.8379 | 0.7403 | 0.7683 |
| 0.70 and above | **0.8617** | 0.8548 | 0.8616 | 0.9334 |

The aggregate difference of 0.128 in AUROC falls, within a band, to at most **0.081**, and the time
split is the better of the two in **one of the four bands**. The aggregate difference of 0.383 in
recall falls to at most **0.072**, and the time split is better in **two of the four**.

**The panel's accuracy is a function of chemical distance, not of date.** Conditioned on how far a
compound sits from the training chemistry, it makes no difference whether it was withheld by year or
at random. What time changes is the input distribution.

That distinction is not a quibble, because the two readings imply opposite things. Temporal decay
would mean the models need periodic refitting and that their published scores expire. A stable
distance-dependence means the scores do not expire, are conditional on a quantity the server already
measures for every query, and can be quoted for the compound in hand rather than for the average
compound.

**A residual remains and is not explained away.** Coarse bands would produce one if the time split's
compounds sat further out within a band, but they do not: inside the widest band both splits have a
median similarity to training of 0.319. What can honestly be said is that the residual cannot be
attributed with this data, because the random split produces very little novel chemistry to test on,
212 actives spread across 39 endpoints, so no single endpoint carries enough of both splits inside
that band for a paired comparison. **Chemical distance accounts for most of the apparent temporal
effect and not demonstrably all of it.**

## 8.8 Three test sets, one curve

The strongest evidence in this chapter is convergent rather than singular. Three test sets were built
by three unrelated rules: withhold the later quarter by publication date, withhold a quarter at
random, and withhold everything one database curated and the other did not. They share no
construction principle. If recall were a property of the split they would disagree; if it is a
property of chemical distance they should trace one curve.

| Band | By date | At random | By curator |
|---|---:|---:|---:|
| below 0.40 | 0.1612 | 0.1179 | 0.0494 |
| 0.40 to 0.55 | 0.5532 | 0.5232 | 0.4641 |
| 0.55 to 0.70 | 0.7403 | 0.7683 | 0.8264 |
| 0.70 and above | 0.8616 | 0.9334 | 0.9011 |

They trace the same curve. A cross-provenance test set, built by a rule that knows nothing about dates
and nothing about the random seed, reproduces the recall-versus-distance relationship of both others.
That is worth more than any single one of the three, because the ways in which they could individually
mislead do not overlap.

## 8.9 A different curator

Time is not the only kind of independence. A model can be prospectively accurate while depending on
the habits of one curation pipeline: which assays ChEMBL selects, how it derives pChEMBL, which papers
it abstracts. Compounds deposited in BindingDB and absent from ChEMBL are curated by different people
from different papers (`results/tables/external_cross_source.csv`).

This test cannot report AUROC, and the reason belongs before the numbers. **Every BindingDB-only row
at every endpoint is an active**, because BindingDB deposits affinities, so a compound absent from
ChEMBL and present there is one somebody measured and found to bind. There are no independently
curated negatives to rank against, and manufacturing them from the background pool would test the
background pool. What is reported instead is recall at the frozen operating point, paired with the
false-positive rate that same threshold produces on background chemistry, because recall alone is
meaningless: a model answering "active" to everything scores 1.0.

| Endpoint | External actives | Recall | Background FPR | Median novelty |
|---|---:|---:|---:|---:|
| A2A | 758 | 0.6332 | 0.0332 | 0.5145 |
| D2 | 86 | 0.6163 | 0.0690 | 0.6257 |
| 5-HT2A | 459 | 0.3508 | 0.0657 | 0.5195 |

Across **3 endpoints and 1,303 independently curated actives**, of which 1,264 are distinguishable
from training in feature space, mean recall is **0.5334** at a mean background false-positive rate of
**0.0560**, against a mean deployed sensitivity of **0.7160** on the same endpoints.

Read alone that looks like a heavy cost for a change of curator. It is not: the median distance from
training of these test sets is 0.5195, the middle of the range, and the recall observed is what
section 8.7 predicts for compounds at that distance. **A change of curator costs almost nothing once
distance is accounted for.** What it changes is which compounds are available to test on, and the
compounds one database holds and another does not are, unsurprisingly, not the ones both already
agreed about.

## 8.10 Specificity, adversarial checks, leakage

**Specificity.** On 1,000 compounds with no recorded activity at any modelled target, **925 returned
no actionable disease signal**, a specificity of **0.925** with a 95 per cent interval of 0.907 to
0.9397 (`results/tables/noncns_specificity_summary.csv`). Those compounds are presumed inactive
because nothing is recorded about them, so the artefact itself labels the paired false-positive rate
an upper bound. Chapter 7 traces why this figure was 0.949 before 25 August and why the change is a
correction rather than a degradation.

**Adversarial checks.** Six checks, each written so that it could fail, and **6 of 6 pass**
(`results/tables/inversion_validation.csv`). Chapter 4 records the one that previously failed and the
fact that its controls were corrected rather than its criterion moved.

**Leakage.** 52 targets were checked for scaffold overlap between training and test and **none shares
a scaffold** (`results/tables/integrity_audit.csv`). Repeated prediction is bit-identical and the
disease score is reproducible, both recorded as passing.

**Prospective scaffold hold-out.** Whole scaffold classes withheld before training, with recall
measured on them at the deployed threshold, over 41 targets: median **0.8180**, mean 0.7362, from
0.000 to 0.970 (`results/tables/scaffold_holdout_results.csv`). The minimum of zero is one of the six
endpoints Chapter 6 identifies as firing for fewer than half their own actives.

## 8.11 What the programme establishes, and what it does not

**It establishes four things.**

The panel's **specificity is robust** to both kinds of independence tested: the false-positive rate on
background chemistry moves by 0.0003 across a temporal wall. Its **ranking against background
chemistry is largely preserved** on chemistry it has not seen, at a cost of 0.062. Its **recall is a
function of chemical distance rather than of publication date**, traced by three test sets built on
unrelated principles, so published scores do not expire and can be conditioned on a quantity the
server already computes. And **sensitivity is the weakest of its claims** where the chemistry is
genuinely novel.

**Four things it does not establish, and the last is a defect rather than a limitation.**

Recall on genuinely distant chemistry is **poor in absolute terms**, at 0.161 below Tanimoto 0.40, and
no analysis here improves it. The finding is that the poor number is predictable, not that it is
better than it looked.

**None of this is prospective in the strict sense.** No compound here was predicted before it was
measured. The dates are real; the analysis is retrospective.

**The cross-provenance arm rests on three endpoints**, because the number of compounds one database
holds at pChEMBL 7 or above and the other lacks is small.

And **the null models have no artefact.** The technical report states that with labels permuted the
same pipeline returns a mean AUROC of 0.4938 random and 0.4921 scaffold, and that an independent
re-run reproduced all 26 core values to within 4.7 × 10⁻⁵. Both figures are **hard-coded in the report
generator**, and no file in the repository holds them. The label-permutation null is the single
strongest piece of evidence that the cross-validation is not inflated by leakage, and it is the one
claim in this chapter a reader cannot check. It should be regenerated into an artefact and declared,
and until it is, this thesis quotes it as reported rather than as verified.

---

## Outstanding items for this chapter

1. **The label-permutation null and the independent-reproduction figures have no artefact.** They are
   hard-coded in `build_technical_report.py`. This is the same defect class as the censored-bound
   count found in Chapter 3, and it affects a more load-bearing claim.
2. **The technical report's cross-provenance section is stale.** It states a mean deployed sensitivity
   of 0.868 on the three endpoints where the regenerated artefact now gives **0.7160**, the
   difference being the sensitivity correction of Chapter 5 propagating into the artefact but not the
   prose.
3. The scaffold hold-out is now 41 targets with a median recall of 0.8180; the technical report says
   39 targets and 0.815.
4. The 8 endpoints excluded from the prospective arm should be listed in the thesis with their
   reasons, as they are in the artefact, rather than only counted.
