# Viva preparation

> **Provenance.** Every number in this document was read from a file in this repository during the
> session in which it was written, and the quantitative claims are pinned by
> `thesis/verify_chapter_numbers.py` under the chapter key `viva`, so a number that drifts fails a
> check rather than surviving into an answer.

## How to use this

Three things get people through a viva on a computational thesis, in this order.

**Know what every symbol means.** An examiner who suspects the candidate is reciting will stop and
ask what a term is. Part I builds each piece from nothing, with a worked example in the thesis's own
numbers, so that "what is expected calibration error" has an answer that starts from what a
probability is.

**Know the spread, never just the mean.** Every headline in this thesis has a range, and the range is
usually the interesting part. Part II is the table to memorise. A candidate who says "AUROC 0.925" is
weaker than one who says "0.925 on average across the eight core classifiers, from 0.8777 at the
barrier model to 0.9648 at BACE1, under scaffold-grouped cross-validation".

**Know where it is weak, and say so first.** Part IV is the five questions that could go badly. In
every one of them the strongest move is to state the limitation before the examiner does, with the
number attached. This thesis is unusually well placed for that: four of its ten hypotheses were
refuted and it says so on the front page of its own falsification chapter.

---

# Part I. The mathematics

## 1. From a molecule to a vector

A model cannot consume a molecule. It consumes a fixed-length list of numbers, and everything the
model can possibly learn is a consequence of what that list captures and what it throws away.

**The fingerprint.** An extended-connectivity fingerprint of diameter 4, ECFP-4, is built by the
Morgan algorithm:

1. Give every atom an initial integer identifier from its own properties: element, degree, charge,
   whether it is in a ring, number of attached hydrogens.
2. For two rounds, replace each atom's identifier by a hash of its own identifier together with the
   sorted list of its neighbours' identifiers and the bond orders joining them.
3. The set of all identifiers seen at any round is the fingerprint: an unbounded set of integers,
   each standing for one circular substructure up to two bonds from a centre, hence diameter 4.

That set is then **folded** into 1,024 bits by taking each identifier modulo 1,024 and setting that
bit. Folding is what makes the vector fixed-length and is also the first thing to admit under
questioning: two different substructures can land on the same bit. This is a **hash collision**, it
is not detectable from the vector, and it is the price of a fixed length.

**Be ready for:** *what does bit 512 mean?* Nothing on its own. A folded bit is the union of every
substructure whose hash falls there, so individual bits are not interpretable and the thesis never
interprets one. What is interpretable is a *distance* between two whole vectors, which is section 2.

**The descriptors.** Twelve whole-molecule numbers are appended: molecular weight, cLogP, topological
polar surface area, hydrogen-bond donors, hydrogen-bond acceptors, rotatable bonds, aromatic rings,
fraction of sp3 carbon, ring count, heavy atoms, formal charge, and QED. Together: 1,024 + 12 =
**1,036 columns**, identical for every endpoint.

**What is thrown away, and this matters.** The featuriser is computed on the desalted, neutralised
parent, and ECFP is a topological descriptor: it encodes connectivity, not geometry. So the vector is
**stereo-blind**. Two enantiomers give byte-identical vectors and therefore identical predictions.
The thesis states this rather than hiding it, and it is the reason the external barrier set has two
numbers instead of one, because excluding overlap by InChIKey is not enough when the InChIKey
separates stereoisomers and the featuriser does not.

## 2. Tanimoto similarity

For two binary vectors $A$ and $B$,

$$T(A,B) \;=\; \frac{|A \cap B|}{|A \cup B|} \;=\; \frac{c}{a + b - c}$$

where $a$ is the number of bits set in $A$, $b$ the number set in $B$, and $c$ the number set in
both. It runs from 0, no shared bits, to 1, identical bit patterns.

**Worked example.** $A$ has 40 bits on, $B$ has 50, and 30 are shared.
$T = 30 / (40 + 50 - 30) = 30/60 = 0.50$.

**Why this and not Euclidean distance.** These vectors are sparse and binary. Euclidean distance
counts the zeros, and two molecules agree on thousands of absent substructures whatever they are, so
almost every pair would look close. Tanimoto counts only the bits that are on in at least one of
them, which is the comparison a chemist actually makes.

**Where it is used here.** Twice, and on a *different* fingerprint from the model's: applicability
distance uses 2,048-bit Morgan fingerprints, unfolded to 2,048 rather than 1,024, so the domain
measure is not limited by the same collisions the model suffers.

Rough calibration for the viva: **T ≥ 0.7** is a close analogue, one or two atoms different;
**0.4–0.55** is a related series; **below 0.4** is a different chemotype, and is where this system's
recall collapses (Part II).

## 3. The random forest

**One decision tree** splits the data repeatedly on single features, choosing at each node the split
that most reduces impurity. It has almost no bias and enormous variance: change a few training rows
and the tree changes completely.

**Bagging** averages that variance away. Draw $B$ bootstrap samples of the training set, fit a tree
to each, and average their predictions. If the trees were independent with variance $\sigma^2$, the
average would have variance $\sigma^2/B$. They are not independent, so the actual result is

$$\operatorname{Var}\!\left(\bar{f}\right) \;=\; \rho\sigma^2 \;+\; \frac{1-\rho}{B}\,\sigma^2$$

where $\rho$ is the pairwise correlation between trees. **The whole design of a random forest follows
from that formula.** The second term vanishes as you add trees; the first does not. So the way to
improve a forest is to reduce $\rho$, and that is what the *random subspace* step does: at each
split, consider only a random subset of features. It makes each tree individually worse and the
ensemble better.

**The settings in this thesis, and why.** 300 trees; minimum 2 samples per leaf on the core
classifiers and 4 on the binder panel; `class_weight="balanced"`; `random_state=42`.

- *300 trees*: the $\,(1-\rho)\sigma^2/B\,$ term is negligible well before 300; more is wasted time.
- *Minimum leaf size 4 on the binder panel*: those endpoints are smaller and noisier, so leaves are
  held larger to stop single compounds forming their own leaf.
- *Balanced class weights*: several endpoints are heavily imbalanced, and without reweighting the
  tree gains more impurity reduction from the majority class.
- *A fixed seed*: so a rerun reproduces, which is what let the thesis compare a corrected artefact
  against its predecessor and attribute the difference to the correction rather than to noise.

**Why a forest gives a probability at all.** The fraction of trees voting for the positive class is
a number in $[0,1]$, but it is not automatically a probability in the sense of "60 per cent of the
compounds I score 0.6 are active". Making it one is calibration, section 10.

**Be ready for:** *why a forest rather than a neural network or gradient boosting?* The honest answer
is in the thesis and is not flattering: the forest is **not** measurably more accurate. Over the
8 classification endpoints under scaffold split it beats XGBoost on 8 of 8 with a median margin of
0.0048 and Wilcoxon p = 0.00781, but over all 13 endpoints including regression it wins on only 8
and the difference is not distinguishable, median delta 0.0033, p = 0.73535. Against histogram
gradient boosting on all 13 it is likewise indistinguishable, p = 0.89258. It is kept for its
calibration behaviour and its out-of-bag structure, not because it wins.

## 4. What a probability means, and why it is not comparable

Each endpoint's model is fitted on that endpoint's own table, and those tables have very different
**base rates** — the proportion of the table that is active. Across the eight core classifiers the
base rate runs from **0.2363 at hERG to 0.8621 at BACE1, a spread of 0.6258**.

The consequence is the single most important idea in the scoring layer: a probability of 0.60 at
hERG, where under a quarter of the table is active, is strong evidence of activity; the same 0.60 at
BACE1, where six sevenths of the table is active, is evidence *against*. A single threshold applied
across the panel would be asking a different question at every endpoint.

## 5. The enrichment map

The fix is to measure how far a probability sits from its own endpoint's base rate, signed, and
scaled so that the answer always lands in $[-1, 1]$:

$$E_t(p) \;=\; \begin{cases} \dfrac{p - b}{1 - b}, & p \ge b \\[2ex] \dfrac{p - b}{b}, & p < b \end{cases}$$

**Read the algebra.** Above the base rate, the largest possible excess is $1-b$, so dividing by $1-b$
maps $p=1$ to $+1$. Below it, the largest possible shortfall is $b$, so dividing by $b$ maps $p=0$ to
$-1$. Both branches give 0 at $p=b$, so the function is **continuous**, and it is **piecewise linear
with a kink at $p=b$** — the two arms have different slopes, $1/(1-b)$ above and $1/b$ below,
which are equal only when $b = 0.5$.

**Worked example, from the thesis.** Donepezil at MAO-A: $p = 0.0237$, $b = 0.2383$. Since
$p < b$, $E = (0.0237 - 0.2383)/0.2383 = -0.2146/0.2383 = -0.90$. At AChE: $p = 1.0000$, $b =
0.5959$, so $E = (1 - 0.5959)/(1 - 0.5959) = +1.00$.

**The clip.** The engagement signal is $\max(0, E_t)$: a probability below base rate is evidence of
*inactivity*, and evidence of inactivity should contribute nothing to a disease score rather than
contributing negatively. That is a design decision, and an examiner may reasonably push on it.

**Be ready for:** *is this not just a linear rescaling?* It is two linear rescalings joined at the
base rate, and the point is the joint. A single linear map could not put the neutral point at $b$ and
the extremes at $\pm 1$ simultaneously for every endpoint.

## 6. Aggregation: the maximum, not the sum

A condition is reached from several targets through a curated graph. Each engaged target contributes
$w \cdot \max(0, E_t)$, where $w$ is the graph edge weight. Those contributions are combined by a
**maximum**, not a sum.

**The argument.** A sum is correct if the engaged targets are independent observations, because three
independent pieces of evidence are stronger than one. They are not. Measured across 400 approved
drugs, **37 targets fire at least once but span only 16 independent directions**, and five homologous
pairs co-fire at $\varphi$ above 0.5: μ and κ opioid at 0.798, SERT and NET at 0.732, D2 and D3 at
0.720. A sum would count the same pharmacology two or three times.

**How "16 independent directions" is obtained.** Build the correlation matrix of the binary firing
indicators over those drugs, take its eigenvalues, and count how many exceed 1. This is Kaiser's
rule. A panel of $N$ perfectly correlated targets has one eigenvalue of $N$ and the rest zero, giving
1; a panel of $N$ independent targets has all eigenvalues 1, giving $N$. It is a blunt instrument and
the thesis says so, but it answers the question asked: how many independent things is the panel
actually measuring?

The maximum **understates** the evidence rather than overstating it, which is the safer direction for
a tool whose purpose is to avoid a wasted experiment.

## 7. The gate, and the proof that it cannot rank

The final score for condition $d$ on compound $x$ is

$$\tilde S_d(x) \;=\; \gamma_d(x)\cdot\max_{(t,w)\in G(d)} w\cdot\max\!\bigl(0,\,E_t(\hat q_t(x))\bigr),
\qquad
\gamma_d(x) = \begin{cases}1, & d \in P\\ \hat q_{\mathrm{BBB}}(x), & \text{otherwise.}\end{cases}$$

$P$ is the set of peripheral conditions, migraine and multiple sclerosis, whose mechanisms act
outside the barrier.

**The proof you must be able to give.** For the 14 non-peripheral conditions, $\gamma_d(x)$ takes the
*same* value $\hat q_{\mathrm{BBB}}(x)$ for every $d$. Multiplying every element of a list by one
positive constant $c$ preserves order: if $S_1 > S_2$ then $cS_1 > cS_2$ for any $c>0$. Therefore the
gate **cannot change which condition ranks first**. It changes only the absolute values, and hence
only whether anything clears the reporting threshold of 0.30.

This is hypothesis H3 and it is **refuted by construction** — refuted by algebra, not by an
experiment. The correct description is that the gate is an **exposure filter over a mechanism
ranking**: it decides *whether* to speak, never *which* condition to name.

**The exception, which you should raise yourself.** The proof holds for the 14. It does not hold
across the boundary to the two exempt conditions, where $\gamma_d = 1$ while the others are scaled.
On teriflunomide, a DHODH inhibitor used in multiple sclerosis, the barrier probability is 0.559, so
the exemption is worth a factor of 1.79 on that condition's score. The gate therefore *does* reorder
a peripheral condition against a central one. No document in the project stated this before the
thesis did.

## 8. AUROC

**Definition.** Sweep the decision threshold from 1 to 0. At each threshold compute the true-positive
rate $\mathrm{TPR} = \mathrm{TP}/(\mathrm{TP}+\mathrm{FN})$ and the false-positive rate
$\mathrm{FPR} = \mathrm{FP}/(\mathrm{FP}+\mathrm{TN})$. Plot TPR against FPR. AUROC is the area under
that curve.

**The interpretation to give in a viva**, because it is the one that shows understanding:

$$\mathrm{AUROC} \;=\; \Pr\bigl(\hat q(x^{+}) > \hat q(x^{-})\bigr)$$

for a randomly chosen active $x^{+}$ and a randomly chosen inactive $x^{-}$, with ties counted as
half. It is the probability that the model ranks a random active above a random inactive. 0.5 is
chance; 1.0 is perfect ranking.

**Equivalently** it is the normalised Mann-Whitney $U$ statistic, $\mathrm{AUROC} = U/(n_+ n_-)$,
which is why a rank test and an AUROC are two views of the same quantity.

**What AUROC does not tell you, and this is the pivot of Chapter 8.** It is *threshold-free*. It
measures ranking only. A model can rank beautifully and still fire for almost none of the actives at
its deployed threshold. That is exactly what happens here: on the prospective split the AUROC against
measured inactives falls from 0.9507 to 0.8231, a drop of 0.1276, while sensitivity at the frozen
threshold falls from 0.8715 to 0.4886, a drop of 0.3829 — **three times larger**. Ranking survives
what a fixed operating point does not.

**Be ready for:** *why not AUPRC?* Under heavy class imbalance, precision-recall is the more
informative curve and its absence is a fair criticism. The defence is that the operating-point
metrics reported alongside AUROC — sensitivity, and the false-positive rate on background — carry the
same information in the form a user actually meets, and that the thresholds are set on a disjoint
pool so those rates are measured rather than fitted.

## 9. Sensitivity, specificity, and the operating point

$$\text{sensitivity (recall)} = \frac{\mathrm{TP}}{\mathrm{TP}+\mathrm{FN}}, \qquad
\text{specificity} = \frac{\mathrm{TN}}{\mathrm{TN}+\mathrm{FP}} = 1 - \mathrm{FPR}$$

Sensitivity answers "of the compounds that really are active, what fraction do I catch?". Specificity
answers "of the compounds that really are not, what fraction do I correctly leave alone?".

**The trade is unavoidable.** Lowering a threshold raises sensitivity and lowers specificity. This
project chose specificity: thresholds hold the background false-positive rate at or under 0.05 on
every deployed endpoint, maximum 0.0500, and the price is a panel mean sensitivity of **0.7638**,
median 0.835, from **0.303 at GABA-A to 0.993 at CGRP**.

**The consequence a user meets**: silence is not evidence of inactivity. Six of the 47 deployed
endpoints fire for fewer than half their own held-out actives, and the server reports the expected
recall at a compound's own distance from training chemistry beside every negative result.

## 10. Calibration

A model is **calibrated** if, among all compounds it scores 0.7, about 70 per cent really are active.
Ranking well and being calibrated are different properties: a monotone transform of the scores leaves
AUROC untouched and can destroy calibration entirely.

**Expected calibration error.** Partition the predictions into $M$ bins $B_1,\dots,B_M$ by predicted
probability and take the weighted average gap between confidence and accuracy:

$$\mathrm{ECE} \;=\; \sum_{m=1}^{M} \frac{|B_m|}{n}\,\bigl|\operatorname{acc}(B_m) - \operatorname{conf}(B_m)\bigr|$$

**State the weakness without being asked:** ECE depends on the binning. Change the number of bins or
whether they are equal-width or equal-frequency and the number moves. It is a summary, not a test.

**Isotonic regression** is the calibrator used on the eight core classifiers. It fits the
best-fitting *non-decreasing* step function from raw score to observed frequency, minimising squared
error subject to monotonicity, and it is solved exactly by the **pool-adjacent-violators algorithm**:
walk the sequence, and wherever a value is lower than the one before it, merge the two into a block
carrying their mean; repeat until the sequence is non-decreasing. It is non-parametric, so it can fix
a badly shaped curve, and for that reason it can overfit on small samples, which is why it is fitted
on **out-of-fold** predictions rather than on the training scores.

**Platt scaling**, used on the binder panel, is the parametric alternative: fit
$\Pr(y=1\mid s) = 1/\bigl(1 + \exp(As+B)\bigr)$ by maximum likelihood. Two parameters, so it is far
more stable on small samples and cannot fix a non-sigmoidal distortion.

**The numbers.** Over the eight core classifiers ECE falls from a mean of **0.0801 to 0.0147**, over
a range of **0.0049 at BACE1 to 0.0412 at BBB**. Two things to volunteer here. First, the barrier
model is the **worst calibrated of the eight**, at nearly three times the mean, and its probability is
not merely reported: it multiplies every gated condition, so its error propagates into all of them.
Second, the headline covers eight estimators. The 38 binder endpoints that carry a measured
calibration error average **0.0762**, about five times the core figure, and nine of the 47 carry no
measured calibration error at all.

## 11. Conformal prediction

Calibration says a probability is honest on average. Conformal prediction gives something stronger: a
**prediction set** with a finite-sample coverage guarantee.

**The construction, inductive (split) conformal.** Hold out a calibration set the model never trained
on. Define a nonconformity score $\alpha$ measuring how badly an example fits a proposed label — for
a classifier, $\alpha = 1 - \hat q(\text{label})$ serves. Compute $\alpha$ for every calibration
point. For a new compound and each candidate label, compute its $\alpha$ and include the label in the
prediction set if its $\alpha$ is below the $\lceil (n+1)(1-\varepsilon) \rceil$-th smallest
calibration score.

**The guarantee.** $\Pr(y \in \Gamma^{\varepsilon}(x)) \ge 1 - \varepsilon$, and it requires only
**exchangeability** of calibration and test data, not any distributional assumption and not that the
model is any good. A useless model still achieves coverage; it does so by returning both labels.

**Mondrian, or class-conditional, conformal** — which is what this thesis uses — applies the quantile
*within each class* separately. Marginal conformal can hit 90 per cent overall while covering the
minority class far worse; Mondrian prevents that, at the cost of needing enough calibration data per
class.

**Why set size is reported beside coverage.** Coverage alone is gameable: always returning every
label gives perfect coverage and no information. Average set size near 1 means the interval commits
to a single label. Here, at a target of 0.90, empirical coverage runs **0.889 to 0.921** with set
sizes from 1.007 to 1.079.

**Volunteer this:** two endpoints undershoot, MAO-A at 0.889 and MAO-B at 0.899. Both are marginal,
and the guarantee is asymptotic in the calibration-set size, but the honest phrasing is "coverage
holds at six of eight and is marginally short at two", not "coverage holds".

**And the bigger gap:** conformal covers the 8 core classifiers only. None of the 47 deployed binder
endpoints has a conformal statement of any kind.

## 12. The applicability domain

The domain measure is the **maximum Tanimoto similarity** from the query to any compound in a
reference library, on 2,048-bit Morgan fingerprints. Bands: **T ≥ 0.50** in domain, **0.30 to 0.50**
near domain, **below 0.30** out of domain. In the evaluation code the flag `in_domain` means
T ≥ 0.30, which is the interface's "in or near", and the thesis lists that naming mismatch as an
outstanding item.

**Be ready for the uncomfortable finding.** The flag does **not** predict discrimination. On an
external set the out-of-domain AUROC of 0.8059 was *higher* than the in-domain 0.7597, on 48
compounds and with no interval. What the distance does predict, strongly, is **recall** — section
Part II. So the honest description is that the applicability band is a statement about how likely the
panel is to miss a real activity, not about whether a returned answer can be trusted.

## 13. Cross-validation, and why the split is the experiment

**$k$-fold cross-validation** partitions the data into $k$ parts, trains on $k-1$ and tests on the
held-out one, $k$ times. Here $k = 10$.

**Random $k$-fold is the wrong experiment for chemistry.** The published record is *series*: a paper
reports twenty analogues of one scaffold. Split at random and close analogues of every test compound
sit in the training set, so the model is asked to interpolate within a series it has already seen.

**Bemis-Murcko scaffold grouping** is the fix. Strip a molecule to its ring systems plus the linkers
joining them, discarding all side chains: that is its scaffold. Group compounds by scaffold and keep
every group whole within a fold. A test compound is then guaranteed to have no training compound
sharing its scaffold.

**The numbers.** Random 10-fold gives a mean AUROC of **0.9575** across the eight core classifiers,
from 0.8990 to 0.9764. Scaffold-grouped gives **0.9252**, from **0.8777 at BBB to 0.9648 at BACE1**.
The gap of about 0.032 is the honest cost of asking the harder question, and the scaffold figure is
the one to quote.

**And neither is obtainable without the labels.** With labels permuted, the same pipeline on the same
folds returns **0.4959** random and **0.5026** scaffold, every one of the sixteen values within
**0.0200** of chance. That is the evidence that the cross-validated numbers are not an artefact of
leakage, and it also shows that scaffold grouping alone confers nothing.

## 14. The statistical tests used, and why each one

**Wilson score interval**, for a proportion. Preferred to the normal approximation
$\hat p \pm z\sqrt{\hat p(1-\hat p)/n}$, which misbehaves near 0 and 1 and can produce bounds outside
$[0,1]$. Wilson inverts the score test instead:

$$\frac{\hat p + \dfrac{z^2}{2n} \;\pm\; z\sqrt{\dfrac{\hat p(1-\hat p)}{n} + \dfrac{z^2}{4n^2}}}{1 + \dfrac{z^2}{n}}$$

*Worked, and it reproduces the thesis exactly.* Specificity is 925 of 1,000, so $\hat p = 0.925$,
$n = 1000$, $z = 1.96$. Denominator $1 + 3.8416/1000 = 1.00384$. Centre
$(0.925 + 0.00192)/1.00384 = 0.92337$. Margin
$1.96\sqrt{0.925\times0.075/1000 + 3.8416/4\,000\,000}/1.00384 = 0.01637$. Interval
$0.9070$ to $0.9397$, which is what `noncns_specificity_summary.csv` records.

**Paired bootstrap**, for the difference between two models on the same compounds. Resample
*compounds* with replacement, recompute both models' AUROC on each draw, and take the distribution of
the difference. Resampling compounds rather than predictions preserves the pairing, which is what
makes the comparison powerful. Used in H10 with 2,000 resamples.

**Wilcoxon signed-rank test**, for comparing two model families across endpoints. Non-parametric,
paired by endpoint. Used instead of a $t$-test because 8 or 13 endpoints is too few to trust
normality, and instead of DeLong because DeLong compares two AUROCs on *one* dataset whereas the
question here is whether one family beats another *across* endpoints. The thesis is explicit that
this is not a DeLong test.

**Permutation test**, for a null with no closed form. Shuffle the labels, recompute the statistic,
repeat; the $p$-value is the fraction of shuffles at least as extreme as the observed value. Used for
the leakage null and for H1's target-to-disease map, where the smallest attainable $p$ with 200
permutations is $1/201 = 0.005$.

**Fisher's exact test**, for a $2\times2$ table with small counts, where the chi-squared
approximation fails. Relevant to H4, where the distant stratum has **one** false positive in 61
compounds.

---

# Part II. The numbers to know cold

## The headline set

| Quantity | Value | Spread | Source |
|---|---|---|---|
| Core AUROC, random 10-fold | 0.9575 | 0.8990 (BBB) to 0.9764 (BACE1) | `rf_cv_summary.csv` |
| Core AUROC, scaffold-grouped | **0.9252** | **0.8777** (BBB) to **0.9648** (BACE1) | `rf_cv_summary.csv` |
| Permuted-label null | 0.4959 random, 0.5026 scaffold | all 16 within 0.0200 of chance | `permutation_null.csv` |
| Panel sensitivity, held-out | **0.7638** | median 0.835, 0.303 (GABA-A) to 0.993 (CGRP) | `binder_modes.json` |
| Background FPR, deployed | mean 0.0234 | maximum 0.0500 | `binder_modes.json` |
| Core ECE after isotonic | 0.0147 | 0.0049 (BACE1) to 0.0412 (BBB) | `calibration.csv` |
| Binder ECE, the 38 measured | 0.0762 | median 0.0675, 0.0300 to 0.1780 | `integrity_calibration_per_target.csv` |
| Conformal coverage, target 0.90 | 0.889 to 0.921 | set size 1.007 to 1.079 | `rf_conformal.csv` |
| Specificity, non-CNS library | **0.925** | 0.907 to 0.9397, an upper bound | `noncns_specificity_summary.csv` |
| External barrier AUROC | **0.7934** on 241 unseen drugs | 0.7645 on 306, 0.7102 on the 65 memorised | `external_bbb_validation.csv` |
| Panel size | 47 deployed of 52 | 6 fail the reliability gate | `binder_modes.json` |
| Graph | 16 conditions, 51 targets | | `GRAPH_FINGERPRINT.json` |

## Recall against chemical distance — the most important table in the thesis

| Nearest-training Tanimoto | Time split | Random split | Cross-provenance |
|---|---:|---:|---:|
| Below 0.40, different chemotype | **0.1612** | 0.1179 | **0.0494** |
| 0.40 to 0.55, related series | 0.5532 | 0.5232 | 0.4641 |
| 0.55 to 0.70, same series | 0.7403 | 0.7683 | 0.8264 |
| 0.70 and above, close analogue | 0.8616 | 0.9334 | 0.9011 |

Three test sets built by unrelated rules — held out by date, at random, by curator — trace one curve.
**Recall is a function of chemical distance, not of how the set was held out.** On a genuinely novel
chemotype the tool finds roughly one true activity in six, and across a database boundary roughly one
in twenty.

## The falsification suite

| | Verdict |
|---|---|
| H1 the disease score is informative | SUPPORTED |
| H2 the curated edge weights add value | **REFUTED** |
| H3 BBB gating discriminates between diseases | **REFUTED**, by construction |
| H4 specificity transfers to novel chemistry | SUPPORTED |
| H5 read-across beats a frequency baseline | SUPPORTED |
| H6 the disease scores match clinical indications | WEAKENED |
| H7 some panel targets cannot rank at all | **REFUTED** |
| H8 engaged targets are independent observations | **REFUTED** |
| H9 the disease layer discriminates between compounds | SUPPORTED |
| H10 the barrier model earns its place over a rule | WEAKENED |

**Four refuted, two weakened.** What the refutations cost: the graph is now described as structure
rather than tuned parameters; the gate is described as a filter rather than a discriminator; the
interface groups homologous targets and quotes their correlation; and four endpoints were withdrawn
for firing on glucose and urea.

## H10, the newest and the one most likely to be probed

| | AUROC on 241 unseen approved drugs | Margin against the deployed forest |
|---|---:|---|
| Deployed forest, 1,036 features | **0.7934** | — |
| Random forest, 12 descriptors only | 0.7701 | −0.0233, interval −0.0553 to **+0.0106**, p = 0.083 |
| Logistic regression, 12 descriptors | 0.7420 | −0.0514, interval −0.0940 to −0.0101 |
| TPSA alone, no fitting | 0.7440 | −0.0494, interval −0.0938 to −0.0047 |
| CNS rule: TPSA ≤ 90 and MW ≤ 400 | 0.6369 | −0.1565, interval −0.2108 to −0.1004 |

On the scaffold hold-out the fingerprint's advantage over twelve descriptors is 0.0338 with an
interval that excludes zero. **On unseen approved drugs the interval crosses zero**, so the verdict is
WEAKENED. The fingerprint's contribution to barrier prediction is demonstrated on the training
distribution and not demonstrated on chemistry the model has never seen.

---

# Part III. Questions, and how to answer them

## On the data

**What is your training data and where does it come from?** Public bioactivity records, principally
ChEMBL, assembled per endpoint, with the negative class partly recovered from censored bounds. Across
the endpoint tables there are **170,619 distinct SMILES**, of which 170,617 parse to a desalted
parent.

**What is a censored bound and why does it matter?** A compound tested and found inactive is often
deposited not as a number but as an inequality, "IC50 > 10 μM". A conventional query for a numeric
potency discards exactly those rows, which throws away the measured negatives and leaves the model
learning from actives against decoys. The rule used here: a bound settles a label when the whole
interval falls one side of the activity cut, and is discarded when it spans both.

**How do you deduplicate?** By InChIKey of the desalted parent. Be ready to say why that is not
sufficient for the external test: the InChIKey separates stereoisomers and salt forms, and the
featuriser does not, so a compound can pass an InChIKey exclusion and still be one the model has
memorised. That is why the external barrier result is reported twice, 0.7645 on all 306 and 0.7934 on
the 241 that are also distinguishable in feature space.

## On the models

**Why 300 trees?** See section 3: the variance term that shrinks with $B$ is negligible well before
300, and the correlation term does not shrink at all.

**Did you tune hyper-parameters, and on what?** Be honest and precise about what was and was not
tuned, and note that the thresholds — the parameters that most affect what a user sees — are set on
a pool disjoint from the one that measures them.

**Why is the barrier model your weakest core classifier?** Scaffold AUROC 0.8777 against a panel mean
of 0.9252, worst calibration of the eight at 0.0412, and H10 shows its advantage over twelve
descriptors is not established on unseen drugs. Three independent measurements agreeing that the
barrier model is the weakest link, and it sits at the gate. Say this before the examiner assembles it.

## On validation

**What is your external validation?** One genuine external set: FDA-curated approved drugs absent
from the barrier model's training database, AUROC 0.7934 on the 241 that are also novel in feature
space. For the target panel no external set of this kind exists, and the reason is structural: for
most of those proteins the public record *is* the training set.

**Your prospective simulation shows decay. Is the model getting worse over time?** No, and this is the
best answer in the thesis. In aggregate the time split looks like decay: sensitivity 0.4886 against
0.8715 for a size-matched random control. But **83.3 per cent of the random split's test compounds
are close analogues of its own training set, against 16.3 per cent for the time split**. The two
splits are not testing comparable populations. Read band by band, the AUROC gap of 0.128 falls to at
most 0.081 and the recall gap of 0.383 falls to at most 0.072.

**Then is it entirely a composition effect?** No, and do not claim it is. Chemical distance accounts
for most of the effect and not demonstrably all of it, because the random split yields only 212 novel
actives across 39 endpoints, too few for a paired within-band comparison.

## On the disease layer

**What does the disease score actually establish?** Three hypotheses, three different questions.
Against the project's own target-to-disease map, top-3 accuracy 0.7901 on 7,008 held-out compounds
against a permutation null of 0.1625. Against real clinical indications, 0.3519 on 162 drugs never
seen in training, which beats a permutation null of 0.1449 decisively and does *not* beat a frequency
null of 0.6543. Against metrics a constant answer cannot pass, mean per-indication AUROC 0.6158,
beating chance on 7 of 9 indications, from 0.7981 for depression down to 0.4898 for epilepsy.

**So does it predict indication?** No. It ranks mechanisms and the conditions those mechanisms touch.
Two of the nine conditions are not being predicted at all.

## On limitations

**What is the single biggest limitation?** Accuracy is a function of chemical distance. Recall runs
from 0.1612 to 0.8616 across four novelty bands, a factor of five, and on a genuinely novel chemotype
the tool finds about one true activity in six. No analysis in the thesis improves that; what the
thesis establishes is that the figure is *predictable*, which is why the server reports it beside a
negative result.

**What would you do next?** In order of what it would change: extend the uncertainty stack to the
target panel, which needs no new data and is the largest single improvement available; report
intervals on every per-endpoint figure; and run a genuinely prospective test, since every result in
the thesis is retrospective — the dates are real, but no compound was predicted before it was
measured.

---

# Part IV. The five that could sink it

**1. "Your model is stereo-blind. How can you predict CNS activity without stereochemistry?"**

Concede it immediately and completely: two enantiomers give byte-identical vectors and identical
predictions, and for targets where activity is enantiospecific this is a real limitation, not a
simplification. Then give the consequence the thesis actually measured: it is why the external
barrier set is reported twice, and the 65 compounds that are feature-identical to a training compound
score 0.7102 against 0.7934 for the genuinely novel ones. The system is a triage tool over
mechanisms, not a stereochemical potency predictor.

**2. "Your headline AUROC of 0.9252 is a mean over eight models. Is that not cherry-picking?"**

The mean is reported with its range every time, 0.8777 to 0.9648, and the weakest member is named:
the barrier model. Go further and volunteer that the eight core classifiers are *not* the endpoints a
user most often queries — the 47 binder endpoints are — and those carry Platt scaling only, at a
calibration error five times the core figure, with no conformal statement at all. The thesis states
that the uncertainty machinery is most complete for the endpoints a user is least likely to query.

**3. "You claim 0.925 specificity. On what?"**

On 1,000 non-CNS compounds **presumed** inactive because nothing is recorded about them, not proven
inactive. It is an upper bound and the artefact labels it as one. They are also drawn from within the
applicability reference, so the figure does not bound behaviour on genuinely distant chemistry. H4
tested exactly that and is the suite's weakest supported verdict: on 61 distant compounds the
false-positive rate is 0.0164, but that is **one** false positive, the 95 per cent interval runs
0.0029 to 0.0872 and contains the 0.0750 comparator. Specificity does not degrade with distance; that
it improves is not established.

**4. "Four of your ten hypotheses were refuted. Does your system work?"**

This is the question to be pleased to get. A refutation is the most valuable output a falsification
suite can produce, and each one cost the project something concrete: a claim about the knowledge
graph, a claim about the gate, an interface change, and four withdrawn endpoints that were firing on
glucose and urea. What survived testing is stated narrowly and holds: the component models rank
compounds well and not by leakage; recall is a measured function of chemical distance; the disease
layer carries real information about mechanism and does not predict indication. A suite in which
nothing is ever refuted is a suite of tests that cannot fail.

**5. "How do I know any of these numbers are real?"**

Every figure in the thesis is read from an artefact at build time, `tools/check_freshness.py` fails
if any declared artefact is older than its inputs, and `thesis/verify_chapter_numbers.py` recomputes
each pinned claim from its artefact and fails if the chapter states something the artefact does not
support. Then volunteer the counter-example, because it is in the thesis: the project's own recurring
defect was numbers that lived in a generator rather than an artefact, and Chapter 10 lists six
instances, of which the label-permutation null and the natural-product coverage figures were given
artefacts during the writing.

---

# Part V. If you do not know

Three sentences that are always better than bluffing.

*"I do not know that off the top of my head, but it is in `<artefact>` and I can tell you the method:
…"* Then give the method. Examiners are testing understanding, not recall.

*"That is a limitation I have not measured. What I would do to measure it is …"* Naming the
experiment is most of the credit.

*"I think you may be right, and here is what would change if you are."* A candidate who can trace the
consequence of being wrong is demonstrating exactly the skill the thesis is about. The project did
this to itself repeatedly: the sensitivity figure that turned out to be measured on training
compounds, the co-firing correlations the interface had drifted from, the barrier model whose
advantage does not survive its own interval.

---

## Outstanding items for this document

1. The mock-viva drill has not been run. Sitting the questions in Part III aloud, against a clock, is
   worth more than reading them.
2. Part I assumes the examiners ask about the methods actually used. It does not prepare answers on
   methods deliberately *not* used, such as graph neural networks beyond the four-endpoint comparison,
   transformer embeddings, or docking, and at least one examiner usually asks why not.
