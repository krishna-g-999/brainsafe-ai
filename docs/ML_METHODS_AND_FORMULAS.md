# ML methods and formulas, with live worked examples

Every formula below is the one actually implemented in the deployed pipeline, not a textbook
restatement of it, and every numeric example is computed live by `tools/build_ml_formulas.py` from
the deployed models on donepezil (`COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2`), the same worked example the manuscript and
graphical abstract use, so the figures here cannot disagree with the ones there. Generated 2026-09-17.

## 1. Feature representation

Each compound is reduced to its largest organic fragment, neutralised and sanitised
(`src/brainsafe/features/featurize.py:parent_mol`), then represented as a fixed **1,036-column
vector**:

**1,024-bit folded ECFP-4 (Morgan) fingerprint, radius 2.** For each atom, the algorithm enumerates
its circular substructure environment out to radius 2 bonds, hashes each environment to an integer,
and sets bit `hash(environment) mod 1024` to 1. Folding is lossy by design: a set bit means *some*
environment hashed there, not which one, which is why a single bit cannot be read as "this
substructure is present" on its own.

**Twelve physicochemical descriptors**, each RDKit's standard implementation: molecular weight,
Crippen logP, topological polar surface area, hydrogen-bond donor and acceptor counts, rotatable
bond count, aromatic ring count, fraction of sp3 carbons, ring count, heavy-atom count, formal
charge, and QED (a composite 0-1 drug-likeness score).

**Worked example, donepezil:** 47 of the 1,024 fingerprint bits are set. Descriptor values:
molecular weight 379.50, cLogP 4.36, TPSA 38.77
Å², QED 0.747. The complete 1,036-value vector for this molecule, and for every
compound behind every model, is `results/tables/master_feature_vectors.csv`.

## 2. Random forest prediction

A random forest is an average of many decision trees, each fitted on a bootstrap resample of the
training rows and, at every split, choosing the best split among a random subset of features
(`max_features = sqrt(1036)` ≈ 32 by scikit-learn's classification default). A forest's raw vote for
class 1 is the mean of every one of its trees' votes:

$$\hat{p}_{\text{forest}}(x) = \frac{1}{n_{trees}}\sum_{i=1}^{n_{trees}} \mathbb{1}[\text{tree}_i(x) = 1]$$

Core target-potency and exposure classifiers use 300 trees, `min_samples_leaf=2`; binder classifiers
use 200 trees, `min_samples_leaf=6` (a coarser leaf, since binder training sets are smaller and more
imbalanced). Both use `class_weight="balanced"` and a fixed `random_state=42`.

For the eight endpoints with a `*_calibrated.joblib` (section 3), this raw forest vote is not the
final answer: it is computed separately inside each of five calibration folds, each fold's vote is
isotonically calibrated on its own held-out data, and the five calibrated values are then averaged.
Section 3 below carries this same worked example through that second step and confirms the two
together reproduce the number the server actually reports.

**Worked example, donepezil against fold 0 of the deployed BBB model's five internal forests:**
individual tree votes in this one fold range from 0.000 to
1.000 across its 300 trees; this fold's raw forest vote is
their mean, 0.9170. Tree 0 of this fold, truncated to depth 3, is Figure S: it
first splits on molecular weight, then QED, before it needs any fingerprint bit; the complete tree is
exported as text in `results/tables/decision_tree_example_full.txt`.

## 3. Calibration

A raw forest vote is not a calibrated probability: a set of compounds a forest votes 0.7 for will
not, in general, turn out to be active 70% of the time. Two methods are used, chosen by how much
calibration data each endpoint's held-out set can support.

**Isotonic regression** (the eight core target/exposure classifiers). Two related but distinct
uses of it exist and are kept separate rather than conflated:

- *Measuring calibration quality.* `src/brainsafe/models/calibrate.py` recalibrates the already-saved
  random-split out-of-fold predictions with a single isotonic regressor fitted via 5-fold
  `cross_val_predict`, so no compound calibrates the function that scores it, and reports the honest
  before/after error. This is what turns the panel's mean expected calibration error from 0.0801
  (raw) to 0.0147 (calibrated).
- *The deployed model.* Separately, `CalibratedClassifierCV(forest, method="isotonic", cv=5)` is
  fitted fresh on all data and saved as `<endpoint>_calibrated.joblib`. This wraps **five** internal
  random forests, one per calibration fold, each with its own isotonic calibrator fitted on that
  fold's held-out rows. At prediction time every fold's forest votes, every vote is passed through
  that fold's own calibrator, and the deployed probability is the mean of the five calibrated
  values:

$$\hat{p}(x) = \frac{1}{5}\sum_{k=1}^{5} g_k\big(\hat{p}_{\text{forest},k}(x)\big)$$

  where $g_k$ is fold $k$'s isotonic calibrator.

**Worked example, donepezil, all five BBB folds:** calibrated per-fold probabilities
0.9829, 0.9834, 0.9897, 0.9974, 1.0000; mean 0.9907. `app.py`'s own reported
probability for this query is 0.9907 — the two agree to the digit, because they are the same
computation checked twice, once inside the deployed object and once by this script independently
loading its five internal forests and calibrators.

**Sigmoid (Platt) scaling** (the binder panel): a raw score $s$ is mapped through a fitted logistic
function $\hat{p} = 1 / (1 + \exp(As + B))$, with $A, B$ fitted by
`CalibratedClassifierCV(method="sigmoid")`. Used for binders specifically because the held-out
calibration set for a single receptor is often too small to fit a non-parametric step function
without overfitting it.

## 4. Conformal prediction

Alongside a calibrated probability, the eight core classifiers report a **Mondrian conformal
prediction set**: not a single number but a statement of which classes are plausible at a stated
confidence level (`src/brainsafe/evaluation/rf_conformal_temporal.py`).

At significance $\varepsilon = 0.10$ (target coverage 90%), on a calibration set disjoint from
training, the *nonconformity* of a calibration compound truly labelled $c$ is
$\alpha_i = 1 - \hat{p}(y_i = c \mid x_i)$: how far the model's own probability for the true
class falls short of certainty. For each class $c$ separately (Mondrian: one threshold per class,
not one threshold overall), the threshold is the

$$k\text{-th smallest calibration nonconformity, } k = \lceil (n_c + 1)(1 - \varepsilon) \rceil$$

among the $n_c$ calibration compounds of that class. At prediction time, class $c$ is included in a
compound's prediction set if its nonconformity $1 - \hat{p}(c \mid x)$ falls at or below that
class's threshold; the reported set can therefore hold zero, one, or both classes.

**BBB model, measured performance**: empirical coverage 0.896 against
a 0.9 target, mean prediction-set size
1.215 on a two-class problem where 1.0 is a single
confident label (`results/tables/rf_conformal.csv`).

## 5. Base-rate enrichment

Training sets are frequently active-heavy (the active fraction reaches over 90% for some deployed
endpoints), so a raw calibrated
probability is not by itself evidence of engagement: an endpoint whose training set is 96% active
will output a high probability for almost anything. Every target score reported to a user is
therefore an **enrichment over that endpoint's own training base rate**
(`app.py:enrichment`), a signed quantity in $[-1, 1]$:

$$e(x) = \begin{cases} \dfrac{\hat{p}(x) - b}{1 - b} & \hat{p}(x) \geq b \\[6pt]
\dfrac{\hat{p}(x) - b}{b} & \hat{p}(x) < b \end{cases}$$

where $b$ is the endpoint's training base rate. $e = 0$ at exactly the base rate, $e \to +1$ as
$\hat{p} \to 1$, $e \to -1$ as $\hat{p} \to 0$; a target is reported as engaged at $e \geq 0.2$.

**Worked example, donepezil at hERG:** calibrated probability 0.7342, training base rate
0.2363. Since 0.7342 $\geq$ 0.2363: $e = (0.7342 - 0.2363) /
(1 - 0.2363) = $0.6520, matching the deployed value to the digit.

**At AChE:** probability 1.0000, base rate 0.5959, enrichment 1.0000.

## 6. Exposure gating and the disease layer

A target's engagement signal is never reported on its own; every disease score is the strongest
engaged target for that condition, gated multiplicatively by predicted barrier penetration
(manuscript Methods, "Exposure gating and the disease layer"):

$$\tilde{S}_d(x) = \gamma(x) \cdot \max_{(t,w) \in G(d)} w \cdot s_t(x)$$

where $s_t(x)$ is the engagement signal (the positive part of enrichment) at target $t$, $G(d)$ is
the set of (target, edge-weight) pairs the knowledge graph connects to disease $d$, and $\gamma(x)$
is predicted barrier penetration ($=1$ for the two peripherally-acting conditions in the graph).
Because $\gamma$ multiplies every condition identically, it can silence a call but cannot re-rank
which condition is reported: this was posed as a falsifiable hypothesis (H3) and confirmed, refuted
by construction, rather than assumed.

**Worked example, donepezil:** barrier penetration $\gamma(x) = $0.9907, driving target
AChE at engagement signal
1.0000, so
$\tilde{S}_{\text{Alzheimer's disease}}(x) = 0.9907 \times
1.0000 = 0.9907$.

## 7. Applicability domain

The applicability-domain distance is the **maximum Tanimoto similarity** between the query's
ECFP-4 fingerprint and every fingerprint in the 158,890-compound background reference library:

$$T(A, B) = \frac{|A \cap B|}{|A \cup B|}, \qquad d(x) = \max_{r \in \text{library}} T(x, r)$$

computed as popcount of the bitwise AND over popcount of the bitwise OR of the two folded
fingerprints. This is reported beside every prediction with the nearest analogue's structure, so a
user can tell interpolation from extrapolation directly rather than trust a single number.

**Worked example, donepezil:** maximum similarity 1.0000 against
158,890 reference compounds (donepezil is itself a training compound, hence the exact match),
tier "In domain", expected recall at this distance 0.8616
(from 2,586 compounds previously measured at a comparable distance).

## 8. Evaluation metrics

**AUROC** (area under the receiver-operating-characteristic curve): the probability that a randomly
chosen active scores higher than a randomly chosen inactive,
$P(\hat{p}(x^+) > \hat{p}(x^-))$, estimated by scikit-learn's `roc_auc_score` over every
active/inactive pair in the held-out fold.

**Sensitivity** (recall): $TP / (TP + FN)$ at the endpoint's deployed decision threshold, measured
on actives withheld by scaffold and never seen during that fit.

**Wilson score interval** (used for recall confidence intervals throughout, in preference to the
normal approximation because it stays inside $[0, 1]$ even at small $n$ or extreme proportions):
for $k$ successes in $n$ trials, at $z = 1.96$ (95%),

$$p = k/n, \quad d = 1 + z^2/n, \quad c = \frac{p + z^2/(2n)}{d}, \quad
h = \frac{z\sqrt{p(1-p)/n + z^2/(4n^2)}}{d}, \quad \text{CI} = [c - h,\ c + h]$$

(`src/brainsafe/evaluation/noncns_specificity.py:wilson`).

---
*Generated by `tools/build_ml_formulas.py` from the deployed models and
`results/tables/rf_conformal.csv`, `models_rf/binder_modes.json`. Regenerate after any change to the
panel with the same command.*
