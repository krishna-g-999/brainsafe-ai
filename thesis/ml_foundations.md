# The mathematical and machine-learning foundations

> Written for a reader who needs to understand this system from the ground up: what a sample is, what
> the features are, where the labels came from, what each algorithm actually computes, and what was
> optimised. Every number was read from a file in this repository on 7 September 2026 at commit
> `f339b5a`. Nothing here is quoted from a figure or from memory.
>
> Read it alongside `script_compendium.md`, which says which file produces which artefact.

---

## 1. The sample space

### 1.1 What one sample is

The underlying population is small drug-like organic molecules. A single **observation** in this
project is not a molecule on its own: it is a **(molecule, endpoint) pair with a measured
experimental outcome**.

Formally, let $\Omega$ be the set of small molecules. For an endpoint $t$ (a protein target, the
barrier, or a physicochemical property) there is an unknown function

$$f_t : \Omega \rightarrow \mathcal{Y}_t$$

that we are trying to approximate, where $\mathcal{Y}_t = \{0, 1\}$ for a classification endpoint and
$\mathcal{Y}_t \subseteq \mathbb{R}$ for a regression endpoint. We never observe $f_t$. We observe a
finite sample of noisy measurements of it, drawn from the published assay literature.

This matters because it defines what the models can and cannot claim. **The sample is not drawn
uniformly from $\Omega$.** It is drawn from what medicinal chemists chose to synthesise and publish,
which is heavily clustered around a small number of chemical series per target. Every design decision
downstream, above all the scaffold split, exists to deal with that non-uniformity.

### 1.2 The size of the sample space

**170,619 distinct SMILES** appear across the endpoint tables, of which **170,617** parse as a
desalted parent (`library_sp3_coverage.csv`). That is the chemistry the whole system has ever seen.

Its composition constrains what can be claimed: median fraction sp3 is **0.3448**, and only
**9.23 per cent** of it is both sp3-rich (fraction sp3 at least 0.55) and has at most one aromatic
ring. That is the region where terpenoid and steroidal natural-product chemistry sits, so the library
is dominated by flat, aromatic, synthetic-medicinal-chemistry space.

### 1.3 Per-endpoint sample sizes, after deduplication

These are the numbers the models are actually fitted and scored on, after byte-identical feature
vectors have been collapsed (section 3.4).

**Core classification endpoints:**

| Endpoint | n | positives | negatives | positive rate | AUROC (random) | AUROC (scaffold) |
|---|---:|---:|---:|---:|---:|---:|
| hERG | 9,933 | 2,370 | 7,563 | 23.9 % | 0.9565 | 0.9270 |
| BACE1 | 8,207 | 7,096 | 1,111 | 86.5 % | 0.9764 | 0.9648 |
| GSK-3β | 5,439 | 4,056 | 1,383 | 74.6 % | 0.9649 | 0.9425 |
| AChE | 5,125 | 3,014 | 2,111 | 58.8 % | 0.9659 | 0.9212 |
| MAO-B | 4,534 | 2,299 | 2,235 | 50.7 % | 0.9630 | 0.9170 |
| BBB | 3,901 | 2,473 | 1,428 | 63.4 % | 0.8990 | 0.8777 |
| MAO-A | 3,585 | 857 | 2,728 | 23.9 % | 0.9619 | 0.9059 |
| BChE | 3,278 | 1,760 | 1,518 | 53.7 % | 0.9724 | 0.9451 |

**Core regression endpoints** (metric is $R^2$):

| Endpoint | n | $R^2$ (random) | $R^2$ (scaffold) |
|---|---:|---:|---:|
| D2 | 7,905 | 0.6403 | 0.5311 |
| A2A | 6,743 | 0.7231 | 0.6205 |
| 5-HT2A | 6,075 | 0.6996 | 0.5560 |
| SERT | 4,479 | 0.6897 | 0.4612 |
| antioxidant (DPPH) | 2,782 | 0.6589 | 0.4153 |

Two things to notice, because a viva will ask about both. **Class balance varies enormously**, from
23.9 per cent positive (MAO-A, hERG) to 86.5 per cent (BACE1); this is why class weighting is not
optional. And **every endpoint drops under the scaffold split**, by between 0.0116 (BACE1) and 0.0560
(MAO-A) for classification, and much more for regression, up to 0.2285 for SERT. That gap is the
honest measure of how much of the random-split performance was analogue recall.

### 1.4 The binder panel and the ADME layer

The **binder panel** holds **52 registered endpoints, 47 of them deployed**. Deployed sensitivity
against held-out actives has mean **0.7638**, median **0.835**, and ranges from **0.303** to
**0.993**. AUROC against measured inactives has median **0.9470**, range 0.719 to 0.985. **Six of the
47 fail the reliability gate** and are marked as such rather than removed, so a user is told the call
is unreliable instead of silently receiving it.

The **ADME layer** carries eight further models. Their sample sizes are much smaller and their
performance correspondingly weaker:

| Endpoint | n | $R^2$ (random) | $R^2$ (scaffold) |
|---|---:|---:|---:|
| solubility | 9,573 | 0.8008 | 0.7263 |
| lipophilicity | 4,200 | 0.6389 | 0.5659 |
| plasma protein binding | 1,797 | 0.4336 | 0.3645 |
| logBB | 1,058 | 0.5836 | 0.4131 |
| hepatocyte clearance | 1,020 | 0.2302 | 0.2062 |
| Caco-2 permeability | 897 | 0.7359 | 0.5818 |
| **Kp,uu** | **566** | **0.4056** | **0.3523** |

Plus two classifiers: P-gp inhibition (n = 1,212, scaffold AUROC 0.9346) and P-gp substrate
(n = 1,371, scaffold AUROC 0.8070).

**Kp,uu is the weakest link and deserves to be flagged in any defence.** It is trained on 566
compounds and reaches $R^2 = 0.3523$ under a scaffold split, yet it produces the headline
"Predicted Kp,uu" verdict the interface shows most prominently. That verdict should be read as a
coarse three-way banding, not a quantitative prediction.

---

## 2. What a "feature" is here

A machine-learning model cannot read a molecule. It reads a vector of numbers. **Featurisation** is
the map

$$\phi : \Omega \rightarrow \mathbb{R}^{1036}$$

from a molecular structure to a fixed-length numeric vector. Every core and binder model in this
project uses the identical $\phi$, defined in `src/brainsafe/features/featurize.py`. That is a
deliberate design choice: it means a comparison between two endpoints, or between two algorithms, is
never confounded by a difference in representation.

### 2.1 The two blocks

$$\phi(m) = \big[\underbrace{b_0, b_1, \ldots, b_{1023}}_{\text{ECFP-4, binary}},\ \underbrace{d_1, \ldots, d_{12}}_{\text{descriptors, real}}\big] \in \mathbb{R}^{1036}$$

**Block 1: a 1,024-bit ECFP-4 fingerprint** (Morgan fingerprint, radius 2). Each $b_i \in \{0, 1\}$.

**Block 2: twelve physicochemical descriptors**, in this fixed order:

| # | Name | Meaning |
|---:|---|---|
| 1 | `mw` | molecular weight |
| 2 | `clogp` | calculated lipophilicity (Crippen) |
| 3 | `tpsa` | topological polar surface area |
| 4 | `hbd` | hydrogen-bond donors |
| 5 | `hba` | hydrogen-bond acceptors |
| 6 | `rotatable_bonds` | conformational flexibility |
| 7 | `aromatic_rings` | count of aromatic rings |
| 8 | `fraction_csp3` | fraction of carbons that are sp3 (three-dimensionality) |
| 9 | `ring_count` | total rings |
| 10 | `heavy_atoms` | non-hydrogen atom count |
| 11 | `formal_charge` | permanent charge |
| 12 | `qed` | quantitative estimate of drug-likeness |

These twelve are the classical determinants of passive membrane permeability, which is why they are
present at all: the barrier endpoint is the one place where physicochemistry, not substructure, is
the mechanism.

### 2.2 How the fingerprint is computed, step by step

The Morgan / ECFP algorithm is iterative and purely combinatorial. There is no learning in it.

**Step 0.** Each atom $a$ receives an initial integer identifier $\mathrm{id}^{(0)}_a$, hashed from
its own properties: element, degree, charge, attached hydrogens, ring membership.

**Step $k$.** Each atom's identifier is updated from itself and its immediate neighbours:

$$\mathrm{id}^{(k)}_a = H\Big(\mathrm{id}^{(k-1)}_a,\ \big\{(\beta_{ab},\ \mathrm{id}^{(k-1)}_b) : b \in \mathcal{N}(a)\big\}\Big)$$

where $H$ is a hash function, $\mathcal{N}(a)$ the bonded neighbours, and $\beta_{ab}$ the bond type.
The neighbour set is sorted before hashing so the result does not depend on atom ordering.

**Radius 2** means this is done twice, so each final identifier encodes the atom together with
everything within two bonds of it: a circular substructure. "ECFP-**4**" refers to the diameter
(2 × radius), which is why radius 2 and ECFP-4 are the same thing.

**Folding.** The set of identifiers produced across a whole library is unbounded. To get a
fixed-length vector, each identifier $h$ is reduced modulo the bit width:

$$b_i = 1 \iff \exists\, h \in \mathcal{I}(m) \text{ such that } h \bmod 1024 = i$$

### 2.3 The consequence of folding: collisions are guaranteed

This is the single most important caveat about the representation, and it is measured rather than
asserted (`fingerprint_collisions.csv`).

On a sample of 20,000 structures from this project's own tables, **52,882 distinct atomic
environments map onto 1,024 bits**. By the pigeonhole principle collisions are unavoidable once the
number of environments exceeds the bit width, and in practice:

- **all 1,024 bits carry more than one environment**;
- the median bit carries **52** distinct environments, the maximum **73**, the minimum **33**.

Therefore **bit $i$ does not mean "substructure $i$ is present"**. It means "at least one of roughly
52 different substructures is present". The practical consequence is confined to interpretation: a
bit that a fitted tree ranks as important names a *set* of environments, so a SHAP attribution to a
single bit cannot be read as a chemical substructure. It does not invalidate prediction, because the
model only needs the mapping to be consistent, and it is.

A second consequence: **the fingerprint is stereo-blind**. `MurckoScaffoldSmiles` is called with
`includeChirality=False` and the fingerprint ignores stereochemistry, so enantiomers produce identical
vectors. Combined with desalting and neutralisation, this means salt forms, protonation states and
stereoisomers of one compound all collapse to one vector. Section 3.4 explains why that forces a
deduplication step.

### 2.4 Standardisation before featurisation

`parent_mol` desalts, neutralises and sanitises before any fingerprint is computed. This is not
cosmetic. Before neutralisation was adopted, haloperidol hydrochloride written the way PubChem serves
it scored BBB **0.613** against **0.993** for the free base, and hERG **0.295** against **0.914**. A
user pasting the salt form lost a cardiac liability flag on a compound that has one.

Permanent charges are deliberately preserved. Choline and neostigmine keep their quaternary nitrogen
and their +1, because that charge is real chemistry and is exactly the property that stops such
compounds crossing the barrier; erasing it would teach the model that they should. Measured over the
library: of 170,617 structures, **198 (0.12 per cent) change representation** and **1,155 charged
structures are correctly left untouched**.

---

## 3. The labels: where $y$ came from

### 3.1 The measurement scale

Activities are pooled from ChEMBL and BindingDB and expressed on the **pChEMBL** scale:

$$\text{pChEMBL} = -\log_{10}\big(\text{molar potency}\big)$$

so an IC50 of 1 µM ($10^{-6}$ M) is pChEMBL 6, and 100 nM is pChEMBL 7. It is a logarithmic scale, so
each whole unit is a tenfold change in potency. Working in logarithms is what makes the errors roughly
additive and the regression targets roughly symmetric.

### 3.2 The classification cuts

From `rebuild_endpoints.py:56`:

$$y = \begin{cases} 1 & \text{pChEMBL} \geq 6.0 \quad (\leq 1\ \mu\text{M, active}) \\ 0 & \text{pChEMBL} < 5.0 \quad (> 10\ \mu\text{M, inactive}) \\ \text{discarded} & 5.0 \leq \text{pChEMBL} < 6.0 \end{cases}$$

**The grey zone between 5 and 6 is thrown away, not assigned.** This is a deliberate and costly
choice. A compound at pChEMBL 5.5 is neither usefully active nor confidently inactive, and forcing it
to one side would teach the model a boundary that the assay data does not support. The cost is
discarded data; the benefit is that the two classes are separated by a full log unit of potency, which
is comfortably larger than typical inter-laboratory assay variability.

The **binder panel uses a stricter cut**: `ACTIVE_P = 7.0` (100 nM), because a binder call is a
stronger claim than an activity call.

### 3.3 Replicates and censored values

Where a compound has been measured more than once against the same target, the replicates are
resolved by **median**, not mean. The median is robust to the single badly transcribed outlier that
is common in aggregated assay literature.

**Censored measurements** ("> 10 µM", meaning the assay saw no effect up to its highest tested
concentration) are handled separately. They are a genuine and informative observation of inactivity,
so discarding them would throw away most of the negative class. They are labelled directly as
inactive rather than passed through `label_from`, because a bound is not a point estimate and treating
it as one would be wrong. Where bounds must be aggregated they combine by `min`; exact values combine
by `median`.

### 3.4 Deduplication, and why it must happen on the feature vector

Because $\phi$ is stereo-blind and operates on the desalted, neutralised parent, several distinct
database rows can produce **byte-identical** feature vectors. To the model they are the same input.

Left in place, any splitter puts copies of one compound on both sides of a fold, and the model is then
scored on rows it has memorised. Deduplication must happen at the level of the **feature vector**,
because that is the level at which the rows are indistinguishable: neither the SMILES string nor the
InChIKey collapses them, and both are unique for every BBB row.

The scale of this is not marginal. **BBB goes from 7,807 raw rows to 3,901**: essentially half the
table is duplicate input.

Where duplicated rows **disagree on the label**, the whole group is **dropped**, not voted on. The
same input carrying both labels cannot be learned from, and picking one would be an arbitrary decision
dressed as data. For regression, the group takes the median.

---

## 4. Feature selection: there is none, and that is the answer

The question "from many features, how did you get to these?" has a specific and slightly unusual
answer here, and it should be given directly.

**No feature selection was performed anywhere in this project.** There is no `SelectKBest`, no
`SelectFromModel`, no recursive feature elimination, no variance threshold, no $\chi^2$ or mutual
information filter, and no PCA. I verified this by searching the entire live tree. **All 1,036 columns
go into every model, always.**

### 4.1 Why not

Feature selection performed on the full dataset, before cross-validation, is one of the classic ways
a reported score is inflated: the selector has seen the test fold's labels, so the "held-out" score is
not held out. Doing it correctly requires nesting the selection inside every fold, which multiplies
the compute and complicates every downstream artefact. Given that tree ensembles perform their own
implicit feature selection at each split, the gain would have been small and the leakage risk real.

**Not selecting is therefore a decision, and its justification is that it is the conservative one.**

### 4.2 What was done instead: an ablation, reported honestly

Rather than selecting features, the two blocks were **ablated** and the result reported whichever way
it came out (`feature_block_ablation.csv`, 13 endpoints):

| Comparison | mean change |
|---|---:|
| combined minus **descriptors only** | **+0.1049** |
| combined minus **fingerprint only** | **+0.0014** |

Read these carefully, because the second is against the project's own interest.

The **fingerprint earns its place decisively**: removing it and keeping only the twelve descriptors
costs 0.1049 on average. That is a large effect and it justifies the substructure representation.

The **twelve descriptors add almost nothing**: 0.0014 on average, and they help on only 9 of 13
endpoints. They are retained for two honest reasons, neither of which is predictive performance. They
cost nothing to compute. And the exposure layer needs them independently, since `mw` and `clogp` are
used directly for property-matched decoy selection and the CNS MPO score.

**The correct statement is "the descriptors were retained, not validated."** Saying that first is much
stronger than having it extracted under questioning.

---

## 5. The linear algebra

### 5.1 The design matrix

Stack the featurised training compounds row-wise:

$$\mathbf{X} \in \{0,1\}^{n \times 1024} \times \mathbb{R}^{n \times 12} \subset \mathbb{R}^{n \times 1036}, \qquad \mathbf{y} \in \{0,1\}^n \ \text{or}\ \mathbb{R}^n$$

For BBB, $\mathbf{X}$ is $3901 \times 1036$. Every algorithm in section 6 except the GNN consumes
exactly this object.

### 5.2 The matrix is wide, sparse and badly scaled

Three properties of $\mathbf{X}$ drive most of the modelling decisions.

**It is sparse.** A drug-like molecule sets only a small fraction of the 1,024 bits, so most entries
are zero. Tree ensembles handle this natively; distance-based methods need a similarity measure
appropriate to sparse binary data, which is why the read-across baseline uses Jaccard rather than
Euclidean distance.

**It is badly scaled across blocks.** Columns 1 to 1024 are in $\{0,1\}$. Column 1025 (`mw`) runs to
several hundred. Column 1036 (`qed`) is in $[0,1]$. For a tree this is irrelevant, because trees split
on thresholds and are invariant to any monotone rescaling of a column. For a distance- or
gradient-based method it is fatal, since `mw` would dominate every distance. **This is why logistic
regression is wrapped in a `StandardScaler` pipeline** and the forests are not: it is a correctness
requirement for one family and a no-op for the other.

**It is wide relative to its height.** With $p = 1036$ and $n$ from 566 to 9,933, some endpoints have
fewer samples than features. In that regime an unregularised linear model can fit the training data
exactly and generalise not at all; this is why the logistic regression baseline carries an $L_2$
penalty and why the tree ensembles, which regularise by averaging and by `min_samples_leaf`, are
better suited to the data as it stands.

### 5.3 Similarity as an inner product

Chemical similarity throughout this project is **Tanimoto**, which for binary vectors is the Jaccard
index and has a clean algebraic form. For fingerprints $\mathbf{a}, \mathbf{b} \in \{0,1\}^p$:

$$T(\mathbf{a}, \mathbf{b}) = \frac{|\mathbf{a} \cap \mathbf{b}|}{|\mathbf{a} \cup \mathbf{b}|} = \frac{\mathbf{a}^\top \mathbf{b}}{\|\mathbf{a}\|^2 + \|\mathbf{b}\|^2 - \mathbf{a}^\top \mathbf{b}}$$

since for binary vectors $\mathbf{a}^\top\mathbf{a} = \|\mathbf{a}\|^2$ is just the number of bits set.
$T \in [0,1]$, and $1 - T$ is the Jaccard distance.

Tanimoto is used in four distinct roles, and it is worth keeping them apart:

1. the **applicability domain**, as $\max_j T(\phi(x), \phi(x_j))$ over training compounds;
2. the **decoy similarity ceiling** during binder training, `TAN_MAX = 0.35`;
3. the **kNN read-across baseline**, as a distance metric;
4. the **leakage checks**, where $T = 1.0$ identifies a compound the model cannot distinguish from a
   training row.

Note that roles 1 and 4 use a **2,048-bit** fingerprint while the models use 1,024 bits. Doubling the
width reduces collisions and so gives a sharper similarity estimate; it does not need to match the
model's own width because it is measuring chemistry, not feeding an estimator.

---

## 6. The algorithms, and what each one optimises

### 6.1 Random forest: the deployed model

A random forest is an average of $B$ decision trees, each fitted to a bootstrap resample of the
training set, each choosing its splits from a random subset of features.

**What a single tree optimises.** At each node the tree searches for the feature $j$ and threshold
$\tau$ that most reduce impurity. For classification the criterion is **Gini impurity**:

$$G(S) = 1 - \sum_{c} p_c^2, \qquad p_c = \frac{\text{weighted count of class } c \text{ in } S}{\text{total weight in } S}$$

and the chosen split maximises the decrease

$$\Delta G = G(S) - \frac{|S_L|}{|S|}G(S_L) - \frac{|S_R|}{|S|}G(S_R)$$

For regression the criterion is variance reduction, equivalently mean squared error.

**Why an ensemble.** A single deep tree has low bias and very high variance. Averaging $B$ trees
trained on bootstrap resamples leaves the bias roughly unchanged and reduces the variance. If the
trees were independent with variance $\sigma^2$, the average would have variance $\sigma^2/B$; they
are correlated, so the true reduction is smaller, and this is exactly why each split considers only a
random subset of features, which decorrelates them.

**How a probability is produced.** The forest's output is the mean over trees of each tree's class
fraction at the leaf the compound lands in:

$$\hat{p}(x) = \frac{1}{B}\sum_{b=1}^{B} p_b(x)$$

This is a *vote share*, not a calibrated probability, which is why section 7 exists.

**The deployed settings** (`RF_COMMON` in `train_rf.py`):

| Parameter | Value | Why |
|---|---|---|
| `n_estimators` | 300 | enough for the variance reduction to have plateaued |
| `min_samples_leaf` | 2 (core), 4 (binder) | prevents single-compound leaves, the main overfitting route |
| `class_weight` | `"balanced"` | see below |
| `random_state` | 42 | reproducibility |
| `max_features` | sklearn default, $\sqrt{p} \approx 32$ | decorrelates the trees |

**Class weighting.** With `class_weight="balanced"`, each class $c$ receives weight

$$w_c = \frac{n}{k \cdot n_c}$$

where $n$ is the total sample count, $k$ the number of classes, $n_c$ the count in class $c$. For
MAO-A (857 positives, 2,728 negatives, $n = 3585$, $k = 2$) this gives $w_1 = 3585/(2 \times 857) =
2.09$ and $w_0 = 3585/(2 \times 2728) = 0.657$. Each positive therefore counts about 3.2 times as much
as each negative in the impurity calculation. Without this, a model on hERG (23.9 per cent positive)
could reach 76 per cent accuracy by predicting "inactive" for everything.

### 6.2 Gradient boosting: XGBoost and HistGradientBoosting

Boosting builds an **additive** model, fitting each new tree to the errors of the current ensemble:

$$F_m(x) = F_{m-1}(x) + \eta \, h_m(x)$$

where $\eta$ is the learning rate and $h_m$ the new tree. This is the opposite philosophy to a forest:
a forest averages many independent low-bias, high-variance trees in parallel, while boosting adds many
shallow high-bias trees in sequence, each correcting its predecessors.

**XGBoost's objective.** It minimises a regularised loss using a second-order Taylor expansion. With
$g_i$ and $h_i$ the first and second derivatives of the loss at the current prediction, the objective
for a tree with leaves $j$ is

$$\mathcal{L} \approx \sum_j \left[ G_j w_j + \tfrac{1}{2}(H_j + \lambda) w_j^2 \right] + \gamma T, \qquad G_j = \sum_{i \in j} g_i,\ H_j = \sum_{i \in j} h_i$$

which is quadratic in the leaf weight $w_j$ and so is solved in closed form:

$$w_j^* = -\frac{G_j}{H_j + \lambda}$$

This closed-form leaf value, plus the $\lambda$ and $\gamma$ regularisation terms, is what
distinguishes XGBoost from classical gradient boosting.

**Deployed comparison settings:**

| Family | Settings |
|---|---|
| XGBoost | `n_estimators=400`, `max_depth=6`, `learning_rate=0.05`, `subsample=0.8`, `colsample_bytree=0.8`, `tree_method="hist"`, `eval_metric="logloss"`, `scale_pos_weight` set per endpoint |
| HistGradientBoosting | `max_iter=400`, `learning_rate=0.06`, `class_weight="balanced"` |

Both are class-balanced: XGBoost through `scale_pos_weight` (the ratio of negatives to positives) and
HGB through `class_weight`. **HGB's `class_weight` was added during the audit**; it had been the only
ensemble competing unweighted, which made the comparison unfair in its disfavour.

### 6.3 Logistic regression

The linear baseline. It models the log-odds as a linear function of the features:

$$\log\frac{p}{1-p} = \mathbf{w}^\top \mathbf{x} + b \qquad \Longleftrightarrow \qquad p = \sigma(\mathbf{w}^\top\mathbf{x} + b) = \frac{1}{1 + e^{-(\mathbf{w}^\top\mathbf{x}+b)}}$$

and fits $\mathbf{w}$ by minimising the $L_2$-penalised log-loss:

$$\mathcal{L}(\mathbf{w}) = -\sum_{i=1}^{n} w_{y_i}\Big[ y_i \log p_i + (1-y_i)\log(1-p_i) \Big] + \frac{1}{2C}\|\mathbf{w}\|_2^2$$

Settings: `C=1.0`, `max_iter=2000`, `class_weight="balanced"`, wrapped in `StandardScaler`.

It is a genuine baseline rather than a straw man, and it is worth knowing that it loses badly: on the
scaffold split it trails the forest on all thirteen endpoints, and on the regression endpoints it
produces **negative $R^2$** (SERT $-0.2638$, antioxidant $-0.7424$), meaning it does worse than
predicting the mean. Structure-activity relationships are not linear in fingerprint space.

### 6.4 kNN read-across

The chemist's own method, formalised. For a query $x$, find its five nearest training compounds by
Jaccard distance on the fingerprint block and average their labels. `n_neighbors=5`,
`metric="jaccard"`. Since Jaccard distance is $1 - T$, this is literally "look at the five most
similar measured compounds", which is what a medicinal chemist does by eye.

It is included because it is the honest null for the whole enterprise: **if a fitted model cannot beat
read-across, it has not learned anything a similarity search would not have given.** The forest beats
it on all thirteen endpoints.

### 6.5 The graph neural network (GIN)

The GNN is the one component that does **not** use $\phi$. Instead of a fixed hashed fingerprint it
learns its own representation from the raw molecular graph.

**Input.** Each molecule becomes a graph: atoms are nodes, bonds are undirected edges. Each node
carries a **37-dimensional** feature vector, built from one-hot blocks:

| Block | Dimensions |
|---|---:|
| element (14 common elements + "other") | 15 |
| degree (0 to 5 + other) | 7 |
| hybridisation (SP, SP2, SP3, SP3D, SP3D2 + other) | 6 |
| attached hydrogens (0 to 4 + other) | 6 |
| formal charge, aromaticity, ring membership | 3 |
| **total** | **37** |

**Nothing is hashed.** That is the point of the method: no information is lost to folding.

**The layer update.** A Graph Isomorphism Network layer updates each atom's vector from the **sum** of
its neighbours':

$$h_v^{(k)} = \mathrm{MLP}^{(k)}\Big( (1 + \epsilon^{(k)}) \cdot h_v^{(k-1)} + \sum_{u \in \mathcal{N}(v)} h_u^{(k-1)} \Big)$$

Sum aggregation, rather than mean or max, is what gives GIN its expressive power: it is provably as
discriminative as the Weisfeiler-Lehman graph isomorphism test, which mean and max aggregation are
not. The $\epsilon$ term is learned and lets the layer weight an atom's own state against its
neighbourhood.

After $K = 3$ layers, each atom's vector encodes its 3-bond neighbourhood. The molecule vector is then
the **mean over atoms**, and a two-layer readout MLP maps it to the endpoint:

$$h_G = \frac{1}{|V|}\sum_{v \in V} h_v^{(K)}, \qquad \hat{y} = \mathrm{MLP}_{\text{readout}}(h_G)$$

**Architecture and optimisation:**

| Parameter | Value |
|---|---|
| hidden dimension | 64 |
| layers | 3 |
| dropout (readout) | 0.2 |
| batch normalisation | after every layer |
| optimiser | Adam, learning rate $10^{-3}$, weight decay $10^{-5}$ |
| batch size | 128 |
| max epochs | 120, early stopping with patience 18 |
| loss (classification) | `BCEWithLogitsLoss`, `pos_weight` = negatives/positives |
| loss (regression) | MSE on standardised targets |

Adam maintains per-parameter running estimates of the first and second moments of the gradient and
scales each step accordingly, which is what makes it robust to the very different gradient scales
across a deep network's layers. Weight decay is $L_2$ regularisation. Early stopping on a validation
split is itself a form of regularisation: training halts when validation AUROC has not improved for 18
epochs, and the best state is restored.

**The result. The GIN loses to the random forest on all four endpoints tested**
(`results/gnn/gnn_vs_rf.csv`), on the identical scaffold hold-out:

| Endpoint | metric | GIN | Random forest | winner |
|---|---|---:|---:|---|
| BBB | AUROC | 0.8868 | **0.9235** | RF |
| BACE1 | AUROC | 0.9243 | **0.9571** | RF |
| MAO-A | AUROC | 0.7363 | **0.8100** | RF |
| A2A | $R^2$ | 0.4673 | **0.5485** | RF |

This is the expected result at this data scale and should be presented as such rather than apologised
for. Graph networks have many more parameters than a forest has effective degrees of freedom, and
they typically need tens of thousands of examples per endpoint before their learned representation
beats a well-made fixed fingerprint. With 3,000 to 10,000 compounds, the fingerprint wins.

---

## 7. Was the protocol the same for every algorithm? No, and the differences matter

This is the question most likely to be asked, and the honest answer has three tiers.

| | endpoints | splits | folds | dedup | class-weighted | purpose |
|---|---:|---|---:|---|---|---|
| **Random forest** | 13 | random **and** scaffold | 10 | yes | yes | **deployed** |
| XGBoost | 13 | random and scaffold | 10 | yes | yes (`scale_pos_weight`) | comparison |
| HistGradientBoosting | 13 | random and scaffold | 10 | yes | yes (added in audit) | comparison |
| Logistic regression | 13 | random and scaffold | 10 | yes | yes + scaled | baseline |
| kNN read-across | 13 | random and scaffold | 10 | yes | n/a | null |
| **GIN (graph network)** | **4** | **scaffold only** | **1 hold-out** | **no** | yes (`pos_weight`) | **exploratory** |

**The four comparison families are exactly like-for-like with the forest**: same 13 endpoints, same
features, same folds, same deduplicated rows. Any difference between them is attributable to the
estimator alone. That is what makes the comparison meaningful.

**The GIN is not.** It runs on 4 endpoints, with a single 70/10/20 scaffold hold-out rather than
10-fold cross-validation, and it does not go through `_dedup_features`. It is therefore an
**exploratory demonstration, not a fair benchmark**, and it should be described that way. Two
consequences follow: its numbers carry no fold-level error bars, and its test set may contain
compounds whose graph is a stereoisomer of a training compound. Since the GIN loses anyway, neither
qualification changes the conclusion, but the conclusion should be stated as "a graph network did not
beat the fingerprint in the one exploratory comparison we ran", not as a benchmarked result.

### 7.1 The two cross-validation schemes

Every core model is validated twice.

**Random 10-fold** (`StratifiedKFold` for classification, `KFold` for regression, both shuffled with
seed 42). Compounds are assigned to folds at random, with stratification preserving the class ratio in
each fold. This is the conventional estimate and it is optimistic, because congeneric analogues land
on both sides.

**Scaffold-grouped 10-fold** (`GroupKFold`). Compounds are grouped by **Bemis-Murcko scaffold**, the
ring systems and the linkers between them with side chains stripped, and whole groups are assigned to
folds so that no scaffold appears in both training and test. This measures generalisation to
structurally new chemistry, which is the regime a user is actually in.

Two subtleties in the implementation are worth knowing:

- The scaffold is computed on the **same desalted parent** the featuriser uses, so a salt and its free
  base cannot land in different folds while being identical to the model.
- **Acyclic compounds share one group** rather than each receiving its own. Giving each its own group,
  as an earlier version did, quietly turned the scaffold split back into a random split for that part
  of the set.

The gap between the two schemes is the honest measure of analogue recall, and it is reported
everywhere rather than averaged away.

---

## 8. What was optimised, and what was not

### 8.1 There was no hyperparameter search

**No grid search, no randomised search, no Bayesian optimisation, and no manual sweep exists anywhere
in this repository.** I verified this by searching the entire live tree for `GridSearchCV`,
`RandomizedSearchCV`, `optuna`, `hyperopt`, `BayesSearchCV` and `param_grid`. There are no matches.

Every hyperparameter in section 6 was chosen by convention and fixed: 300 trees and
`min_samples_leaf=2` are standard defaults-plus-convention, and the one deliberate departure,
`min_samples_leaf=4` on the binder panel, is documented with its reason.

**Both halves of this must be said together.**

*What it costs.* No claim of optimality is available. There is direct evidence that performance was
left on the table: an untuned XGBoost scores higher than the deployed forest on **all five** regression
endpoints, and HistGradientBoosting beats it on several classification endpoints too.

*What it buys.* There was no tuning set, therefore there can be no tuning-set leakage. Selecting
hyperparameters on the same folds from which the final number is reported is the single commonest way
a figure in this field is quietly inflated, and this project is structurally immune to it. Combined
with the scaffold split, the deduplication and the permutation null, the cross-validated numbers here
carry about as little optimistic bias as this design admits.

The forest is deployed rather than XGBoost **for its calibration behaviour and its out-of-bag
structure, not for its accuracy**, and the significance testing supports that framing: over all 13
endpoints the forest is not distinguishable from either booster (Wilcoxon $p = 0.89258$ against HGB,
$p = 0.73535$ against XGBoost).

### 8.2 Optimisation at three levels

It helps to separate the three distinct things the word "optimisation" refers to here.

**Level 1, inside a single tree.** A greedy search over (feature, threshold) pairs maximising impurity
decrease. Exact, deterministic given the data and the random feature subset, no iteration.

**Level 2, fitting the ensemble.** The forest has no global objective; it is an average of
independently fitted trees. Boosting does have one, minimised by additive stagewise descent. The GIN
has one, minimised by Adam over 120 epochs with early stopping.

**Level 3, choosing the hyperparameters.** This is the level at which nothing was optimised.

Most published "we optimised our model" claims refer to level 3. Here, only levels 1 and 2 occurred.

---

## 9. The probability layer

### 9.1 What the raw number is, and is not

A forest's output $\hat p(x)$ is the mean vote share across 300 trees. It is a score in $[0,1]$, and
it is **systematically miscalibrated**: forest probabilities cluster toward the middle, because
averaging bounded quantities pulls extremes inward. A raw 0.72 is not a 72 per cent chance.

Measured over the eight core classifiers, the mean **expected calibration error** before calibration is
**0.0801**.

### 9.2 Expected calibration error

Partition $[0,1]$ into $M = 10$ equal-width bins. In each bin compare the mean predicted probability
with the observed frequency of positives, and average by bin occupancy:

$$\mathrm{ECE} = \sum_{m=1}^{M} \frac{|B_m|}{n} \Big| \underbrace{\bar{p}(B_m)}_{\text{confidence}} - \underbrace{\bar{y}(B_m)}_{\text{accuracy}} \Big|$$

A perfectly calibrated model has $\mathrm{ECE} = 0$: among compounds it scores 0.7, exactly 70 per
cent are active.

### 9.3 Isotonic calibration

The fix is **isotonic regression**: fit a non-decreasing step function $g$ mapping raw scores to
calibrated probabilities, minimising squared error subject to monotonicity:

$$\min_{g \text{ non-decreasing}} \sum_i \big(y_i - g(\hat p_i)\big)^2$$

solved exactly by the **pool-adjacent-violators algorithm** (PAVA), which repeatedly merges adjacent
blocks that violate monotonicity and replaces them with their common mean, in $O(n)$ time.

Two properties matter. Because $g$ is monotone, **calibration cannot change the ranking of compounds**,
so AUROC is unchanged; it changes only what the number means. And because it is non-parametric it can
correct an arbitrary distortion, unlike Platt scaling which assumes a sigmoid.

It is fitted on **out-of-fold** predictions under an inner 5-fold, so nothing is calibrated and
evaluated on the same rows. Mean ECE falls from **0.0801 to 0.0147**.

*Caveat, measured during the audit:* this figure is computed under the random split. Under the
scaffold split it is **0.0169**, so the choice of split costs about 0.0018 and the calibration claim
survives essentially unchanged.

The binder panel uses **Platt scaling** (a fitted sigmoid) rather than isotonic, because its
calibration sets are smaller and isotonic overfits on small samples.

### 9.4 Base rates and enrichment

A calibrated probability still cannot be read in isolation, because endpoints have very different
prevalences: 23.9 per cent for hERG, 86.5 per cent for BACE1. A probability of 0.5 is strong evidence
of activity in a rare endpoint and evidence of *inactivity* in a common one.

The **enrichment map** rescales the probability against the endpoint's base rate $b$:

$$E_b(p) = \begin{cases} \dfrac{p - b}{1 - b} & p \geq b \\[2ex] \dfrac{p - b}{b} & p < b \end{cases}$$

It is piecewise linear with a kink at $p = b$, continuous there (both branches give 0), and maps
$E_b(1) = +1$, $E_b(0) = -1$. It is **strictly increasing in $p$**, which is the key property: it
cannot reorder compounds within an endpoint, so it changes interpretation without changing ranking.

### 9.5 Conformal prediction

Calibration gives a probability. Conformal prediction gives a **set** with a coverage guarantee.

**Inductive Mondrian (class-conditional) conformal at $\varepsilon = 0.10$.** Split the data into
training, calibration and test. Define the **nonconformity** of a calibration compound as

$$\alpha_i = 1 - \hat{p}(y_i \mid x_i)$$

that is, one minus the predicted probability of its *true* class. For each class $c$ separately, take
the threshold at rank

$$k = \lceil (n_c + 1)(1 - \varepsilon) \rceil$$

of the sorted $\alpha$ values for that class. The prediction set for a test compound is every class
whose nonconformity falls below its class threshold:

$$\Gamma(x) = \{ c : 1 - \hat{p}(c \mid x) \leq \alpha^{(k)}_c \}$$

The guarantee, which holds under **exchangeability** of calibration and test data, is
$\mathbb{P}(y \in \Gamma(x)) \geq 1 - \varepsilon$. "Mondrian" means the threshold is taken per class,
which prevents the guarantee being met on the majority class alone.

**Reading the set size.** With two classes, a set can have size 0, 1 or 2, and

$$\mathbb{E}[|\Gamma|] = 1 + \mathbb{P}(\text{ambiguous}) - \mathbb{P}(\text{empty})$$

so the excess over 1.0 is ambiguity **minus** emptiness. A size-2 set says "either class is
plausible"; a size-0 set says "neither is", which is a legitimate signal that the compound is
anomalous. Measured coverage runs **0.876 to 0.933** against the 0.90 target and mean set size from
**0.956 to 1.215**. The smaller figure is the one that proves the formula matters: a mean below 1.0
is impossible unless sets are empty, and BACE1 returns nothing at all for **4.4 per cent** of
compounds. BBB is ambiguous on **21.5 per cent** under the random split, rising to **27.0 per cent**
under a scaffold split.

### 9.6 Interval estimates

Proportions are reported with the **Wilson score interval** rather than the normal approximation,
because the latter misbehaves near 0 and 1 and can produce bounds outside $[0,1]$:

$$\frac{\hat{p} + \frac{z^2}{2n} \pm z\sqrt{\frac{\hat{p}(1-\hat{p})}{n} + \frac{z^2}{4n^2}}}{1 + \frac{z^2}{n}}$$

Differences between models on the same compounds use a **paired bootstrap**; differences across
endpoints use the **Wilcoxon signed-rank test**, with the minimum attainable $p$ reported alongside so
that an underpowered comparison is not mistaken for a null. At $n = 5$ that floor is **0.0625**, above
0.05, so the five regression endpoints cannot reach significance however unanimous the result.

---

## 10. What type of classification this is

Three points that are commonly misunderstood about the architecture.

**These are independent binary classifiers, not a multi-class or multi-label model.** Each endpoint has
its own model, its own training table, its own threshold and its own calibrator. There is no shared
representation, no softmax over targets, and no joint loss. A compound predicted active at AChE and
BChE is the output of two entirely separate forests.

**The consequence is that the panel does not enforce consistency.** Nothing prevents the panel
returning a chemically implausible combination, because no component sees more than one endpoint. This
is a real limitation and it is the price of the modularity that lets an endpoint be withdrawn without
retraining anything else.

**The disease layer is not learned at all.** It is a fixed, curated knowledge graph mapping targets to
conditions with weights, verified against KEGG and Reactome. Disease score is
$\max_t (w_t \cdot s_t) \times \text{BBB}$, the strongest engaged mechanism scaled by predicted
exposure. This is deliberate: an earlier prototype that learned the disease layer from curated
annotation was shown by ablation to be reading the answer back out of its own features, with
structure-only performance collapsing to near zero. Making the layer non-learned is what stops the
disease results being circular.

---

## 11. The one-page summary

| Question | Answer |
|---|---|
| What is a sample? | A (molecule, endpoint) pair with a measured assay outcome |
| How many? | 170,619 distinct structures; 566 to 9,933 per endpoint |
| What are the features? | 1,024-bit ECFP-4 + 12 descriptors = **1,036**, identical for every model |
| Are the features interpretable? | The 12 descriptors yes; the 1,024 bits **no**, median 52 environments share each bit |
| How were labels made? | pChEMBL ≥ 6.0 active, < 5.0 inactive, 5.0–6.0 **discarded**; binder panel uses 7.0 |
| How were features selected? | **They were not.** All 1,036 always used; justified by ablation, not selection |
| Do the descriptors help? | +0.0014 on average. Retained, not validated |
| Does the fingerprint help? | +0.1049. Decisively yes |
| Which algorithm is deployed? | Random forest, 300 trees, `min_samples_leaf` 2, `class_weight` balanced, seed 42 |
| Same protocol for all? | Yes for XGBoost, HGB, logistic regression, kNN. **No for the GIN** (4 endpoints, 1 hold-out) |
| Did the GNN win? | No, it lost on all four |
| What was tuned? | **Nothing.** No search of any kind exists in this repository |
| What kind of classification? | Independent binary classifiers per endpoint, plus 5 regressors and 8 ADME models |
| How are probabilities made honest? | Isotonic calibration, ECE 0.0801 → 0.0147, plus conformal sets at 90 % |
