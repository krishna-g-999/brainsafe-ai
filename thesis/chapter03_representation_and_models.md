# Chapter 3. Representation and models

> Every quantitative statement in this chapter was computed during this session from the artefact
> named beside it. Where a figure elsewhere in the project describes an earlier engine, that is said
> rather than glossed. Citation numbers refer to `manuscript/references.md`.

---

## 3.1 One representation for every endpoint

A compound enters this system once and leaves as a fixed vector of **1,036 columns**, verified
directly from `src/brainsafe/features/featurize.py`:

- **1,024 bits** of a folded ECFP-4 fingerprint [20], Morgan radius 2. Bit *j* is set when some atom
  environment of radius at most 2 hashes to *j* modulo 1,024. Folding is not injective, so a set bit
  says that *some* environment hashing there is present, not which one.
- **12 descriptors**, in this order: molecular weight, cLogP, topological polar surface area,
  hydrogen-bond donors, hydrogen-bond acceptors, rotatable bonds, aromatic rings, fraction sp³, ring
  count, heavy atoms, formal charge, and QED.

Nothing is shared between endpoints except this map. Two models that disagree about a compound
disagree about the same 1,036 numbers, which is what makes the panel comparable at all, and it is why
a defect in the featuriser would be a defect in seventy models at once.

The choice of a fingerprint plus a handful of physicochemical descriptors, rather than a learned
representation, is defended in section 3.9 by measurement rather than by preference.

## 3.2 Standardisation, and why neutralisation belongs in the representation

Before featurisation a structure is reduced to its largest organic fragment, sanitised, and
neutralised. The third step is the one that needs arguing, because it looks like housekeeping and is
not.

A drug and its salt must give the same answer. A user who pastes haloperidol hydrochloride and a user
who pastes haloperidol are asking one question. Stripping the counter-ion alone does not achieve
this: it leaves the parent carrying the charge the salt gave it, so the two forms differ in the
`formal_charge` descriptor and, through the hashed environments around the protonated nitrogen, in
the fingerprint as well.

I verified the current pipeline directly rather than taking this on trust. Submitting both forms of
haloperidol gives **byte-identical feature vectors**, and the deployed models return **identical
probabilities: 0.9888 at the barrier and 0.8216 at hERG for both**. The standardiser does what it
claims.

The distinction the rule has to preserve is between a protonation state and a permanent charge. A
protonated amine deposited as a hydrochloride is an artefact of how the compound was crystallised and
is undone. A quaternary ammonium's charge is a property of the molecule and is precisely what stops
it crossing the barrier, so it is kept. Both behaviours were confirmed: neostigmine retains a formal
charge of +1 through standardisation, while a protonated trialkylamine salt is returned neutral.

The rule fires rarely and matters when it does. Across the 228,198 parseable training structures,
neutralisation changes **207 rows, 0.091 per cent**, touching 28 of the 63 endpoints
(`results/tables/uncharging_impact_by_endpoint.csv`). DHODH is the most affected at 65 of 2,079 rows,
3.13 per cent. A rule that alters one row in a thousand is easy to dismiss, and the reason it is not
dismissed here is that its failures are concentrated exactly where a user is most likely to submit a
salt form, which is an approved drug.

**One correction to the project's own documents.** The technical report and the manuscript both
motivate this rule with the figures 0.613 against 0.993 for the barrier and 0.295 against 0.914 for
hERG. Those numbers describe the behaviour of an **earlier engine, before neutralisation was part of
the representation**. The technical report says so; the manuscript does not, and reads as though the
deployed server still behaves that way. It does not. The honest form of the claim is that the defect
was real, was measured at those values on the models then deployed, and is now absent, with the
current pipeline returning identical vectors and identical probabilities for both forms.

## 3.3 What the representation deliberately discards

Chirality is excluded from the hash, so two enantiomers produce identical rows. This is the sharpest
criticism available against a CNS panel and Chapter 10 bounds what it costs, at 0.19 per cent of the
panel where a stereochemical distinction could change a class call.

It has an immediate structural consequence that belongs here rather than in the limitations. Because
enantiomers collide, rows identical in feature space must be collapsed before any split is drawn, or
a model is tested on a compound it was trained on. The scale is visible in the model metadata: the
barrier endpoint enters with 7,805 rows, sheds **3,773 duplicate rows** and **131 conflicting groups**
whose members disagree about the label, and leaves with **3,901** (`models_rf/BBB_meta.json`).
Conflicting groups are dropped rather than resolved by vote, because the featuriser cannot tell their
members apart and a vote would be inventing an answer. Chapter 2 gives the panel-wide figures.

The barrier endpoint loses nearly half its rows to this. That is not a defect in the data; it is the
price of a stereo-blind representation, paid openly at the featurisation boundary rather than
concealed inside an inflated cross-validation score.

## 3.4 Do both blocks earn their place?

The 1,036 columns are two blocks of very different character, and each was tested alone
(`results/tables/feature_block_ablation.csv`).

Across the thirteen endpoints, the combination is best on **9**, the fingerprint alone on **4**, and
the descriptors alone on **none**. But the margin is thin. Over the eight classification endpoints
the combined representation beats the fingerprint alone by a mean of **+0.0026** AUROC. Over the five
regressions it is **worse** by a mean of **-0.0004**.

The fingerprint carries almost all of the signal. The descriptors are close to free and are retained,
but a reader should not be told they are doing heavy work.

Where they do earn their place is more specific than the project's documents currently say. The
useful comparison is not how much the descriptors add to the fingerprint but how close they come to
it alone, because that measures whether physicochemistry is the mechanism. At the barrier endpoint
the descriptors alone reach 0.9498 against 0.9556 for the fingerprint, a gap of **0.006**. At GSK-3β
the same comparison is 0.8779 against 0.9644, a gap of **0.087**. Twelve numbers nearly suffice to
predict whether a compound crosses the blood-brain barrier, and come nowhere near predicting whether
it inhibits a kinase. That is the expected result if barrier penetration is governed by
physicochemistry and target binding by substructure, and it is a cleaner statement of the same
intuition.

**A correction.** The technical report says the descriptors "earn their place mainly on the exposure
endpoints". Under the reading above, that is right. Under the more natural reading, how much the
descriptors add to the fingerprint, it is not: the largest gain is at hERG, at +0.0072, ahead of
MAO-B at +0.0037 and the barrier model at +0.0030. The sentence should say which quantity it means.

Two caveats travel with this table and both are already recorded in the technical report. It uses
five-fold random cross-validation, not the ten-fold random and scaffold-grouped regimes of the
headline figures, so its absolute values are not comparable with Chapter 8. And it predates
deduplication, which inflates the barrier endpoint specifically: measured against the deployed panel,
seven of the eight classifiers agree to within 0.011 AUROC while the barrier model reads 0.060 high,
that endpoint being where the feature-identical duplicates concentrate. The comparison the table
exists to make, between blocks on identical data, is unaffected, because whatever inflates one block
inflates all three.

Per-descriptor permutation importance across thirteen endpoints
(`results/tables/feature_descriptor_importance.csv`) puts topological polar surface area at the top,
informative at 12 of the 13, followed by cLogP and hydrogen-bond donors at 11 each. Formal charge is
informative at only 3, which is consistent with section 3.2: it matters intensely for the few
compounds that carry a permanent charge and not at all for the rest. No descriptor is uninformative
everywhere, which is the weak justification for keeping all twelve.

## 3.5 The estimator

A **random forest** [22] is fitted per endpoint, with hyperparameters read from the model metadata
rather than from the prose: **300 trees**, `min_samples_leaf` **2** for the core classifiers and 4 for
the binder panel, `class_weight="balanced"`, `random_state=42`. Trees are grown on bootstrap
resamples considering ⌊√1036⌋ = 32 features per split.

A forest averages the votes of many decision trees, each grown on a resample of the rows and each
choosing among a random subset of features at every split. The variance a single deep tree would have
is averaged away. Three properties suit this problem, and only the first is about accuracy.

**It handles the representation as it is.** A sparse binary fingerprint sits beside twelve continuous
descriptors on wildly different scales, and a tree ensemble needs no scaling, no imputation and no
encoding to use both.

**It does not extrapolate.** A forest's prediction is an average of training labels, so it cannot
return a value outside the range it has seen. For a system whose central honest claim is that it
behaves predictably at the edge of its applicability domain, refusing to extrapolate is the correct
behaviour rather than a weakness. A boosted model with a linear leaf, or a neural network, will
happily produce confident nonsense off the manifold.

**It is exactly explainable.** TreeSHAP [24] computes exact Shapley values for a tree ensemble in
polynomial time. For most other model families the same explanation is an approximation, and an
explanation a user cannot check is a liability in a tool of this kind.

Class weighting is not decorative. Chapter 2 records a median active fraction of 0.825 across the
deployed endpoints, reaching 0.962 at P2X7. Without weighting, answering "active" at P2X7 scores
0.962 accuracy, and the forest would have every incentive to do so.

## 3.6 The model-family comparison

Five families were compared under identical five-fold cross-validation on both split regimes, over
thirteen endpoints (`results/tables/model_comparison.csv`). Two are baselines a reader is entitled to
demand: a five-nearest-neighbour Tanimoto read-across, which is what a medicinal chemist does by eye,
and L2-regularised logistic regression. Three are ensembles: the random forest, XGBoost [23], and
histogram gradient boosting.

Under the scaffold-grouped split, the family means are:

| Family | 8 classification endpoints, AUROC | 5 regression endpoints, R² |
|---|---:|---:|
| Random forest | **0.9212** | 0.5190 |
| XGBoost | 0.9156 | 0.5413 |
| Histogram gradient boosting | 0.9153 | **0.5424** |
| kNN read-across | 0.8829 | 0.4745 |
| L2 logistic regression | 0.8352 | 0.1397 |

The result splits cleanly by task, and saying so is the only honest way to report it. The
histogram-gradient-boosting figure is stated here with `class_weight="balanced"` set, matching the
forest's own weighting and XGBoost's `scale_pos_weight`; an earlier run of this comparison left it
unweighted, the only family competing without one, which read 0.9149 and is corrected in the current
artefact.

**On classification the forest leads**, best on 7 of the 8 endpoints. It loses AChE to histogram
gradient boosting, 0.9148 against 0.9247.

**On regression the forest leads nothing.** It is best on **0 of the 5**. Histogram gradient boosting
wins three and XGBoost two. The gap in the family mean is 0.023 R².

The read-across is the baseline that matters most, because a model that cannot beat nearest-neighbour
lookup is an expensive lookup table. The forest exceeds it on **all thirteen endpoints**, by margins
running from +0.0095 at D2 to +0.1099 on the antioxidant assay. That it wins everywhere, but by a
margin of only about 0.038 in the classification mean, is itself informative: much of the signal in
this problem genuinely is similarity, and the forest's contribution is the part that is not.

## 3.7 Is the margin real?

Every document in this project has quoted those deltas as point estimates. The decisions log makes
precisely the right criticism of an earlier comparison, that "point deltas alone do not establish
significance", and answers it with DeLong tests, but **those tests were run on the superseded
ensemble and were never re-run for the deployed forests**. The significance table they cite reports
ensemble AUROCs of 0.921 at the barrier and 0.867 at MAO-A, which are not this panel's numbers.

The gap is now closed, at least partly
(`results/tables/model_family_significance.csv`, generated by
`src/brainsafe/evaluation/model_family_significance.py`). Pairing the families by endpoint and
applying a Wilcoxon signed-rank test to the per-endpoint means under the scaffold split gives:

| Comparison | Forest higher on | Median Δ | Mean Δ | p | Verdict |
|---|---:|---:|---:|---:|---|
| vs kNN read-across, all 13 | 13 of 13 | +0.0382 | +0.0407 | 0.00024 | distinguishable |
| vs L2 logistic regression, all 13 | 13 of 13 | +0.1030 | +0.1988 | 0.00024 | distinguishable |
| vs XGBoost, all 13 | 8 of 13 | +0.0033 | **-0.0051** | 0.735 | **not distinguishable** |
| vs histogram gradient boosting, all 13 | 8 of 13 | +0.0026 | **-0.0054** | 0.893 | **not distinguishable** |
| vs XGBoost, 8 classification only | 8 of 8 | +0.0048 | +0.0056 | 0.0078 | distinguishable |
| vs histogram gradient boosting, 8 classification only | 7 of 8 | +0.0049 | +0.0059 | 0.078 | not distinguishable |

Both columns are given because against the boosting methods they disagree in sign, and the
disagreement is the result. The forest sits marginally above them on 8 of the 13 endpoints, which is
what the positive median reports, and below them by larger amounts on the five regressions, which is
what drags the mean to -0.0051. A summary that quoted only the median would say the forest is
narrowly ahead; one that quoted only the mean would say it is narrowly behind. Neither is a real
difference, which is what the p-values say.

The picture that emerges is sharper than any point estimate. Against the two baselines the forest is
decisively better and the result is not close. Against the two boosting methods, pooled over all
thirteen endpoints, **it is statistically indistinguishable**, because the classification wins and
the regression losses cancel. Restricted to classification it beats XGBoost on all eight endpoints
significantly, and its lead over histogram gradient boosting does not reach the 0.05 level.

Two limitations of this test must travel with it. It is **not a DeLong test**: DeLong compares two
ROC curves on identical compounds and accounts for their correlation, which is the stronger
instrument, but it applies per endpoint on classification only and needs fold-level predictions
rather than the summary table. That remains the better analysis and has not been done for this panel.
And the regression subset is **underpowered**: with five endpoints, the smallest attainable Wilcoxon
p-value is 0.0625, so "the forest loses all five" is a clean count and not a significant result by
this instrument.

The correct summary is therefore that the forest was **not selected because it is more accurate than
boosting**, since on this evidence it is not. It was selected on three grounds that are not about the
mean: it does not extrapolate, TreeSHAP is exact for it, and it is the most stable of the three under
hyperparameter choice. A reader is entitled to say that boosting would have served about as well on
classification and rather better on the receptor regressions, and the cost of the uniformity is
stated at 0.023 R² rather than left to be discovered.

## 3.8 A learned representation, tested

A graph isomorphism network was trained on the raw molecular graph and compared against the forest on
an identical scaffold hold-out, over four endpoints (`results/gnn/gnn_vs_rf.csv`).

| Endpoint | Metric | GIN | Random forest |
|---|---|---:|---:|
| BBB | AUROC | 0.8868 | **0.9235** |
| BACE1 | AUROC | 0.9243 | **0.9571** |
| MAO-A | AUROC | 0.7363 | **0.8100** |
| A2A | R² | 0.4673 | **0.5485** |

The forest wins all four, and the margin is not marginal: 0.074 at MAO-A and 0.081 R² at A2A. This is
consistent with the literature on data scale, since graph networks generally need far more data per
task, or pretraining, to overtake fingerprints. It is four endpoints and not thirteen, so it bounds
rather than settles the question, and fine-tuning a pretrained chemical language model remains the
most defensible thing an unconvinced reader could ask for next. It has not been tried.

## 3.9 Would more data have helped?

`results/tables/learning_curve.csv` fits four endpoints at 10, 25, 50, 75 and 100 per cent of their
training data. The technical report reads this as showing that "the curves flatten well before the
full training set", and uses that as the argument against gathering more data and for having spent
the effort on the negative class instead.

**That reading does not survive the table, and this is the clearest overstatement I have found in the
project's documents.** Measuring the gain from half the data to all of it:

| Endpoint | Metric | Half data | Full data | Gain |
|---|---|---:|---:|---:|
| BBB | AUROC | 0.9276 | 0.9325 | +0.0049 |
| BACE1 | AUROC | 0.9293 | 0.9481 | +0.0188 |
| MAO-A | AUROC | 0.8410 | 0.8764 | +0.0354 |
| A2A | R² | 0.5369 | 0.6223 | **+0.0854** |

Only the barrier model has flattened. MAO-A gains 0.035 AUROC on its final doubling, which is larger
than the entire margin separating the forest from the read-across baseline. A2A gains 0.085 R² and
shows no sign of saturating: at a quarter of its data it reaches only 75.0 per cent of its full-data
score, where the three classifiers reach 92 to 96 per cent.

The defensible conclusion is the opposite of the one recorded, and it is more useful. The
**classification** endpoints are approaching saturation, so more measured actives would buy little
there, and the decision to spend the effort on recovering the negative class instead was right for
them. The **receptor regressions are data-limited**, and A2A is the endpoint that gained most from
the BindingDB expansion in Chapter 2, at +0.0129 R² on 1,238 added compounds. Two independent lines
of evidence therefore point the same way: potency regression on these receptors would improve with
more data, and that is a specific, costed piece of future work rather than a general wish.

## 3.10 Why seventy models and not one

The obvious alternative is a single multi-task network with one output per endpoint. It is rejected
for reasons about the data rather than about preference, and three of the four are established in
Chapter 2.

**The label matrix is almost entirely missing.** 228,200 measurements over 63 endpoints and roughly
169,000 compounds fill under two per cent of a compound-by-endpoint matrix. A compound measured at
AChE has almost never been measured at OX2. Multi-task learning shares strength across tasks; at this
sparsity it mostly shares absence, and a missing measurement is not a negative result.

**The negative class means different things per endpoint.** Where censored bounds were recovered a
negative is a measurement; where they were not it is a property-matched decoy, which is an
assumption. One loss over both silently averages the two.

**Base rates are incompatible.** A shared output layer propagates the 0.962 active fraction at P2X7
into another endpoint's decision boundary.

**Failure stays local.** Five endpoints were trained, tested and withdrawn. In an independent panel a
withdrawal removes one output. Under shared weights the same data would have shaped every other
endpoint's representation, and withdrawing it cleanly would not be possible.

The cost is real and is not hidden: a multi-task model can borrow statistical strength for a small
endpoint from a large related one, and the smallest deployed set here has 387 compounds and would
plausibly benefit. That benefit is forgone deliberately, because the price is a shared representation
in which no endpoint's threshold, negative class or withdrawal is independently defensible.

---

## Outstanding items for this chapter

1. ~~The haloperidol figures in the manuscript describe a fixed defect without saying so.~~
   **Resolved for the full manuscript**, which now reads "on the models this server *previously
   deployed* it returned a barrier probability of 0.613 against 0.993", correctly scoping the figures
   to the superseded engine. `NAR_condensed_draft.md` states the same pair more tersely, "haloperidol
   hydrochloride *then* scored" the salt figure "against 0.993 for the free base", which carries the
   same past-tense qualification in fewer words; neither document claims 0.993 is the free base's
   current score, which now stands at 0.9888 after the retrain this chapter records.
2. **The learning-curve conclusion in `docs/TECHNICAL_REPORT.md` section 7.3 is not supported by its
   own table** and should be replaced by the split reading in section 3.9 above.
3. **The descriptor claim in section 7.2 of the technical report is ambiguous** between two
   quantities that give different answers, and should name which it means.
4. **Two evidence pointers in `docs/decisions_log.md` are broken.** The entries of 2026-07-20 cite
   `results/tables/STable14_significance.csv` and `results/tables/STable15_cv_comparison.csv`. Both
   files exist at `supplementary/`, not `results/tables/`. Both also describe the superseded
   ensemble, which the entries should state.
5. A DeLong test on the deployed forests, per classification endpoint against the read-across
   baseline, remains undone and is the stronger version of section 3.7.
