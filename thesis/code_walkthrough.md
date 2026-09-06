# Code walkthrough: how an endpoint is built, split, trained and thresholded

> **Provenance.** Every setting quoted here was read from the source file named beside it during the
> session in which this was written, and the counts were recomputed from the endpoint tables. The
> quantitative claims are pinned by `thesis/verify_chapter_numbers.py` under the key `code`.

## What this document is for

Supervisors asked for the scripts explained: the criteria that define an endpoint, how the random
forest, the boosting models and the graph network were used, and what was done for data separation.
This walks the pipeline in the order it runs, names the file and the line for each decision, and
gives the reason the code has that value rather than another.

Read it with the source open. Every heading is a real file.

---

# 0. Before anything else: a discrepancy, now corrected

This section is kept in full although the defect is fixed, because the reasoning is the answer to a
question a reader of the earlier report will ask, and because a circulated document cannot be
unsent. If your copy of the technical report predates 6 September 2026, its section 3.1.1 shows the
version described below.

The diagram now appears as two, one per stage, and states 6.0 for the table builder and 7.0 for the
binder panel with the decoy rule beside the stage that owns it.
`tests/test_panel_app_consistency.py::TestPipelineDiagramMatchesThePipeline` fails if either cut,
the decoy rule or the reliability gate drifts from the code again, and it checks the 64,419 figure
against the tables rather than trusting the sentence.

**What it used to say.** `docs/TECHNICAL_REPORT.md:373` drew the labelling rule as

> pChEMBL >= 7 : active, pChEMBL <= 5 : inactive, in between : discarded

**The code that builds the tables does not do that.** `src/brainsafe/data/rebuild_endpoints.py:56`:

```python
def label_from(pvalue: float) -> int:
    return 1 if pvalue >= 6.0 else (0 if pvalue < 5.0 else -1)
```

Active at **6.0**, not 7. The tables settle it: **64,419 rows carry label 1 with a pChEMBL between
6.0 and 7.0**, and under the diagram's rule not one of them could be active.

What has happened is a **conflation of two stages that legitimately use different cuts**:

| Stage | File | Active cut | Why |
|---|---|---|---|
| Build the endpoint table | `rebuild_endpoints.py:56` | pChEMBL ≥ 6.0 | The conventional 1 µM boundary, with 5–6 dropped as a grey zone |
| Select actives for the binder panel | `train_binders_hybrid.py:57` | pChEMBL ≥ 7.0 | A stricter bar, because this model's negatives are presumed decoys rather than measured inactives |

The diagram also places "Tanimoto below 0.35" in the table-building stage, which is `TAN_MAX` from
the binder script. So the figure is describing the binder-panel flow while captioned as table
construction. It is a documentation error: **the code and the data agree with each other**, and only
the diagram disagrees with both.

**One apparent contradiction to be able to answer.** `label_from` discards everything in [5, 6), yet
**14,420 rows sit at pChEMBL exactly 5.0 and carry label 0**. They did not come through
`label_from`. pChEMBL 5.0 is 10 µM, the standard censored inactive bound, and those rows enter
through the measured-inactive path where the label is set directly. That is censored-bound recovery,
visible in the data.

---

# 1. The criteria: what makes an endpoint

**`src/brainsafe/data/rebuild_endpoints.py`** — writes `data/endpoints/<target>.csv` with columns
`smiles, label, pchembl, year, source`. There are **63** such tables on disk.

**Sources are pooled at the compound level, not concatenated.** Each of ChEMBL and BindingDB
contributes a per-compound median potency on the shared $-\log_{10}$ molar scale; the two are pooled
by InChIKey of the desalted parent. Nothing is imputed and no source overrides a measurement.

**Only potency types are kept**: IC50, Ki, Kd, EC50.

**The label rule** is `label_from` above. The grey zone from 5 to 6 is **dropped rather than
assigned**, because a compound at 300 nM is neither confidently active nor confidently inactive and
forcing it into a class teaches the model a boundary the data does not support.

**Censored bounds.** A compound tested and found inactive is often deposited as an inequality rather
than a number. The conventional query for a numeric potency discards exactly those rows, which throws
away the measured negatives. The recovery rule, `rebuild_endpoints.py:130-137`:

```python
recs.append({... "pchembl": float(r["pchembl_bound"]) ...})
# The bound is an upper limit on potency, so the weakest bound is the safest summary.
return df.groupby("inchikey").agg(smiles=("smiles", "first"),
                                  pchembl=("pchembl", "min"),
                                  year=("year", "min")).reset_index()
```

Two things to be able to explain. **`min`, not `median`** — the aggregation for exact values is a
median (`:110` and `:172`) but for bounds it is a minimum, because a bound is an upper limit and the weakest
one is the claim you can actually defend. And a bound settles a label **only when the whole interval
falls one side of the cut**; one that spans the cut is undecidable and is discarded.

**The script refuses to run on incomplete inputs.** Every cached source is checked before anything is
touched, every table is rebuilt in memory and checked for plausibility before any file is written,
and a rebuild that would empty or halve a table exits non-zero unless `--allow-shrink` is passed. On
a fresh clone the caches are absent and the correct outcome is a clear error, not sixty-three empty
tables.

---

# 2. The representation: `src/brainsafe/features/featurize.py`

Every endpoint sees the same **1,036 columns**: 1,024 folded ECFP-4 bits plus 12 descriptors
(`mw, clogp, tpsa, hbd, hba, rotatable_bonds, aromatic_rings, fraction_csp3, ring_count,
heavy_atoms, formal_charge, qed`), computed on the desalted, neutralised parent.

The two facts that matter downstream: folding to 1,024 bits means **hash collisions are possible and
undetectable**, and ECFP is topological, so the vector is **stereo-blind**. Enantiomers produce
byte-identical vectors. Section 3.1 is a direct consequence.

---

# 3. Data separation: four mechanisms, four different questions

This is the part supervisors asked about most specifically, and it is four separate things. Confusing
them is the commonest way to misdescribe the pipeline.

## 3.1 Deduplication, on the feature vector — `train_rf.py:_dedup_features`

```python
seen: dict[bytes, list[int]] = {}
for i, vec in enumerate(X):
    seen.setdefault(vec.tobytes(), []).append(i)
```

**Why the feature vector and not the InChIKey.** The featuriser is stereo-blind, so stereoisomers,
salt forms and protonation variants of one compound give identical vectors. Neither the SMILES string
nor the InChIKey collapses them — both are unique for every BBB row. If they are left in place, any
splitter puts copies of one compound on both sides of a fold and the model is scored on rows it has
memorised. Deduplication must happen at the level at which the rows are indistinguishable, which is
the vector.

**Conflicts are dropped, not voted on.** If a duplicate group disagrees on the label, the whole group
is discarded:

```python
labels = {int(y[i]) for i in idxs}
if len(labels) > 1:
    conflicts += 1
    continue
```

The same input carrying both labels cannot be learned from, and choosing one is an arbitrary decision
dressed as data. Regression takes the group median instead.

## 3.2 Scaffold grouping — `train_rf.py:_scaffold_groups`

Bemis-Murcko: strip side chains, keep ring systems and the linkers joining them.

```python
scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
codes.append(mapping.setdefault(scaf or "_ACYCLIC_", len(mapping)))
```

Three decisions in two lines, and be ready for all three.

- **`includeChirality=False`** — the scaffold retains atom and bond types but not stereochemistry,
  consistent with a stereo-blind featuriser. It is not the generic carbon framework, which would need
  `MakeScaffoldGeneric` and would group far more loosely.
- **Computed on the same desalted parent the featuriser uses**, so a salt and its free base cannot
  land in different folds while being identical to the model.
- **Acyclic compounds share one group.** They have no scaffold, and the docstring records what
  happened when they each got their own: it "quietly turned the scaffold split into a random split
  for that part of the set". A fold boundary must never run through them.

## 3.3 The three background pools — `src/brainsafe/models/pools.py`

The binder pipeline asks the background library for three different things, and taking two of them
from one place makes the reported number an arithmetic consequence rather than a measurement.

```python
SHARES = {"decoy": 60, "threshold": 20, "evaluation": 20}

def _band(smiles, salt="brainsafe-background-v1"):
    digest = hashlib.blake2b(f"{salt}:{smiles}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big") % 100
```

- **decoy (60)** — trained on as presumed negatives.
- **threshold (20)** — sets the operating point, as a quantile of the score distribution.
- **evaluation (20)** — reports the false-positive rate at that point.

**The argument to give.** If the threshold is the 95th percentile of a sample, the false-positive rate
on *that same sample* is 5 per cent whatever the model does. It is not a measurement, it is the
quantile restated. Drawing decoys from either of the other two is worse: the model was explicitly
trained to score those compounds near zero and is then congratulated for doing so.

**Why a hash and not a shuffle.** Two properties a shuffle does not give: the assignment does not
depend on the order the library happens to be in, and **adding compounds later leaves every existing
assignment untouched**, so a threshold set today stays comparable with a rate measured next year.

## 3.4 The binder panel's own holdouts — `train_binders_hybrid.py`

```python
ACTIVE_P, DECOY_RATIO, TAN_MAX = 7.0, 3, 0.35
ACTIVE_HOLDOUT = 0.20   # share of active scaffold groups withheld from training
```

- **A fifth of the active scaffold groups never enters training**, so sensitivity is measured on
  structurally distinct compounds rather than on recall of the training set.
- **Measured inactives are split in half**: one half trains as hard negatives, the other sets the
  threshold and validates.
- Both halves are written to `models_rf/holdout/<T>_binder_holdout.json`, **so a later threshold step
  cannot reach past the holdout by re-reading the endpoint table.** That is not hypothetical: it is
  the defect Chapter 5 documents, where a downstream script re-read the whole table and set the
  threshold on 40 of 49 endpoints from data the model had seen.

---

# 4. The random forest — `src/brainsafe/models/train_rf.py`

The primary model. Eight endpoints are classification, five are regression, because those receptors'
measured sets are almost entirely active and a binary split is ill-posed.

```python
SEED = 42
N_SPLITS = 10
RF_COMMON = dict(n_estimators=300, min_samples_leaf=2, n_jobs=-1, random_state=SEED)
```

Classifiers add `class_weight="balanced"`; regressors do not.

**The cross-validation, `_cv`:**

```python
schemes = {
    "random": (StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED)
               if task == "classification" else KFold(...)),
    "scaffold": GroupKFold(N_SPLITS),
}
```

Both schemes are run for every endpoint and both are reported. `GroupKFold` takes the scaffold codes
as `groups`, which is what keeps a scaffold whole within a fold. `GroupKFold` is deterministic and
takes no seed, which is why no `random_state` appears on it.

**Out-of-fold predictions are saved**, not just the metrics:

```python
oof[split_name] = pd.DataFrame({"smiles": smiles, "y_true": y, "fold": fold_id,
                                "prediction": preds, "scaffold_group": groups})
```

This is what later work stands on: calibration is fitted on these rather than on training scores, and
H10 reuses the stored fold assignment so its comparison sits on **identical** folds rather than
similar ones.

**Outputs**: `rf_cv_folds.csv` (one row per endpoint × split × fold), `rf_cv_summary.csv` (mean and
sd), `models_rf/<endpoint>.joblib` (refit on all data), `<endpoint>_meta.json`.

---

# 5. The binder panel — `src/brainsafe/models/train_binders_hybrid.py`

The 47 deployed endpoints. Different problem from the core eight: measured inactives are scarce, so
the negatives are a **hybrid** of measured inactives and property-matched decoys.

```python
RF = dict(n_estimators=300, min_samples_leaf=4, n_jobs=-1,
          random_state=42, class_weight="balanced")
```

**Leaf size 4 rather than 2**, because these tables are smaller and noisier and a leaf of one
compound is memorisation.

**Decoy selection**, `:159-169`:

```python
need = max(DECOY_RATIO * len(act) - len(ina_train), len(act))
...
if max(DataStructs.BulkTanimotoSimilarity(fp, afps)) < TAN_MAX:
    dec.append(str(bg_ok[i]))
```

Three decoys per active (`DECOY_RATIO = 3`), each at maximum Tanimoto **below 0.35** to every active,
drawn **only from the decoy pool**, with the target's own measured inactives excluded from decoy
eligibility so a compound cannot be both a trained negative and an unseen one.

**Why the hybrid exists at all**, from the docstring: decoy-only training saturates. For melatonin MT1
the measured inactives scored a median binder probability of 0.972 against an active median of 0.997,
so the decision boundary sat where tiny probability differences swing the result, and sensitivity at a
controlled false-positive rate collapsed. The fix is to give the model the hard negatives it never
saw.

**Calibration is Platt, not isotonic**, `:184`:

```python
forest = RandomForestClassifier(**RF).fit(Xt, yt)
cal = CalibratedClassifierCV(FrozenEstimator(forest), method="sigmoid").fit(Xc, yc)
```

`method="sigmoid"` is Platt scaling: two parameters, stable on the small calibration sets these
endpoints have. `FrozenEstimator` means the forest is **not** refit inside the calibrator, so the
calibration is fitted on held-out scores from the fitted model rather than on a re-cross-validated
one.

**A related script, `train_measured_label_holdout.py`**, handles SIRT1 and Nav1.5, which have too few
binders at pChEMBL ≥ 7 for the hybrid scheme and are trained on the measured activity label instead
(active ≥ 6, inactive < 5). Its docstring is worth reading aloud in a viva: it records that an earlier
version "claimed more than the code delivered", because every positive still entered training while
sensitivity was reported on all of them.

---

# 6. Calibration of the core eight — `src/brainsafe/models/calibrate.py`

Two separate things happen, and they are easy to conflate.

**Measurement** — how much calibration helps — is done on the **already-saved out-of-fold
predictions**, with an inner 5-fold so nothing is evaluated on data it was fit on:

```python
p_cal = cross_val_predict(IsotonicRegression(out_of_bounds="clip"), p, y, cv=5)
```

**Deployment** — the model actually served — wraps the forest:

```python
cal = CalibratedClassifierCV(base, method="isotonic", cv=5)
```

`out_of_bounds="clip"` matters: isotonic regression is only defined on the range it was fitted over,
and a test score outside that range is clipped to the nearest end rather than extrapolated.

**Output**: `results/tables/calibration.csv`, Brier and ECE before and after, per endpoint.

---

# 7. Do other model families do better? — `src/brainsafe/evaluation/model_comparison.py`

Five families, the same features, the same folds, 5-fold random and scaffold-grouped.

```python
"RandomForest":          RandomForestClassifier(n_estimators=300, min_samples_leaf=2,
                                                class_weight="balanced", random_state=SEED)
"XGBoost":               XGBClassifier(n_estimators=400, max_depth=6, learning_rate=0.05,
                                       subsample=0.8, colsample_bytree=0.8, tree_method="hist",
                                       eval_metric="logloss", scale_pos_weight=pos_weight)
"HistGradientBoosting":  HistGradientBoostingClassifier(max_iter=400, learning_rate=0.06)
"kNN read-across":       KNeighborsClassifier(n_neighbors=5, metric="jaccard")
"LogisticRegression":    make_pipeline(StandardScaler(),
                                       LogisticRegression(max_iter=2000, C=1.0,
                                                          class_weight="balanced"))
```

For regression the same five, with `Ridge(alpha=1.0)` standing in for logistic regression.

**Why five and not three.** From the docstring: "a tree ensemble beating another tree ensemble is a
weak result on its own". Two of the five are baselines a reader is entitled to demand — a
nearest-neighbour read-across, which is what a medicinal chemist does by eye and which any model must
beat to justify itself, and L2-regularised logistic regression, the simplest thing that could work on
this representation.

**Three details worth knowing.**

- **`metric="jaccard"` is 1 − Tanimoto.** The kNN is therefore literally a read-across over the five
  nearest measured analogues, not a generic distance baseline.
- **The read-across sees only the fingerprint** (`FINGERPRINT_ONLY = {"kNN read-across"}`), because
  Jaccard is undefined on continuous descriptors. Every other estimator sees all 1,036 columns.
- **`scale_pos_weight` for XGBoost, `class_weight="balanced"` for the others.** Same intent,
  different API; HistGradientBoosting has neither and is left unweighted, which is a real asymmetry
  in the comparison and should be conceded if asked.

**Two defects this script was written to fix**, and the docstring says so: an earlier version ran only
the three ensembles while the manuscript reported all five, so two quoted numbers had no artefact
behind them; and it split the raw endpoint table where `train_rf.py` splits the deduplicated one,
which compares estimators on a task slightly easier than the one they are deployed on.

---

# 8. The graph neural network — `src/brainsafe/gnn/`

Three files: `graph_features.py` turns a molecule into a graph, `gin_model.py` is the network,
`train_gnn.py` runs the comparison.

**The graph.** Each atom becomes a node with a **37-dimensional** feature vector: one-hot atomic
number over 14 common elements plus "other", one-hot degree 0–5, one-hot hybridisation over five
types, one-hot hydrogen count 0–4, then formal charge, aromaticity and ring membership as scalars.
Each bond becomes **two directed edges**, so messages flow both ways.

**The network** is a Graph Isomorphism Network:

```python
class GINLayer(nn.Module):
    def forward(self, x, edge_index):
        agg = torch.zeros_like(x)
        agg.index_add_(0, dst, x[src])          # sum neighbour messages into each node
        return self.mlp((1 + self.eps) * x + agg)
```

This is the GIN update $h_v' = \mathrm{MLP}\bigl((1+\varepsilon)h_v + \sum_{u \in N(v)} h_u\bigr)$.
**Sum aggregation, not mean or max**, is the whole point of GIN: sum is injective over multisets, so
the layer can distinguish neighbourhoods that mean-pooling would collapse. `eps` is learned.

Hidden width 64, 3 layers, batch norm after each, mean pooling as readout, dropout 0.2.

**The comparison, `train_gnn.py`.** Four endpoints — BBB, BACE1, MAO_A, A2A — on a **single scaffold
hold-out**: 20 per cent of scaffold groups as test, a further slice as validation for early stopping,
`MAX_EPOCHS = 120`, `PATIENCE = 18`, batch 128. The GIN and a random forest are trained on the
identical training compounds and scored on the identical test compounds, so the only thing that
differs is the model.

**What to concede before being asked.** This is four endpoints, not thirteen, and one split, not ten
folds. It **bounds** the question rather than settling it, and the honest next step — fine-tuning a
pretrained chemical language model — has not been tried. The docstring says CPU is fine for this
demonstration and that the full run belongs on a GPU cluster.

---

# 9. Thresholds: four scripts, one file, one correct order

All four write `models_rf/binder_modes.json`, and each depends on the one before it.
`tools/check_freshness.py:123` records the hazard: running any of them alone silently reverts the
later ones.

```
final_thresholds.py → screening_thresholds.py → apply_specificity_decisions.py
                                              → calibrate_background_specificity.py
```

**`final_thresholds.py`** sets each threshold from the **stricter of two independent requirements**:

```python
TARGET_FPR = 0.10        # at most this share of the target's held-out measured inactives
BACKGROUND_FPR = 0.05    # at most this share of the threshold pool
```

Requirement A is selectivity against tested non-binders, read from the holdout JSON — **only the
withheld half counts**. Requirement B is specificity against unrelated chemistry, on the threshold
pool, which is disjoint from the decoys the model trained on and from the evaluation pool the rate is
reported on. The reason B exists is concrete: a prospective test found Nav1.1 assigning glucose a
binder probability of 0.806 while scoring 0.979 against its own assay's negatives. A model judged only
on its own chemistry can be confidently wrong about everything else.

**`screening_thresholds.py`** computes a second, stricter set at `SCREENING_FPR = 0.01`, for
low-prevalence screening where a false positive is expensive. Both sets are stored and neither is a
correction of the other.

**`apply_specificity_decisions.py`** acts on `results/deployed_specificity_audit.csv`, which scores
every binder at its calibrated threshold against 600 random PubChem structures and against molecules
no CNS target plausibly binds — glucose, urea, acetate, ethanol, glycine, lactate, atenolol. Nav1.1
was withdrawn on this evidence. Every decision is **re-derived whenever the models change, never
inherited**, because a withdrawal is a claim about a particular fit.

**`calibrate_background_specificity.py`** runs last and therefore decides. This is the script Chapter 5
is about: it used to select actives as `df.loc[p >= 7, "smiles"]` over the whole endpoint table,
scoring the model on compounds it had been fitted on while the registry still carried
`sensitivity_basis: "held_out_actives_by_scaffold"`. It now scores `models_rf/holdout/` and declares
its basis when it cannot.

**The reliability gate** lives in one place, `src/brainsafe/panel.py:60`:

```python
MIN_SENSITIVITY = 0.50
MIN_AUROC = 0.75
```

Five writers, the figure script and `app.py` all import it. It previously lived in six places and
disagreed with itself in three.

---

# 10. Questions to expect on the code

**"Why deduplicate on the feature vector rather than the InChIKey?"** Because the featuriser is
stereo-blind, so the InChIKey does not collapse rows the model cannot tell apart. Both are unique for
every BBB row.

**"Your scaffold split — what happens to acyclic compounds?"** They share one group. Giving each its
own group turns the scaffold split into a random split for that part of the set, and the code says so
because it happened.

**"Why is the threshold pool separate from the evaluation pool?"** Otherwise the reported
false-positive rate is the quantile restated. It is arithmetic, not a measurement.

**"Why isotonic on the core and Platt on the panel?"** Isotonic is non-parametric and can fix a badly
shaped calibration curve, but it needs data; Platt has two parameters and is stable on the small
calibration sets the binder endpoints have.

**"Why does the binder panel use pChEMBL ≥ 7 when the tables are built at ≥ 6?"** Because its
negatives are presumed decoys rather than measured inactives, so it selects only high-confidence
actives. The technical report's diagram used to state the ≥ 7 cut at the wrong stage; section 0.

**"Is `class_weight='balanced'` applied consistently across the comparison?"** No, and concede it.
Random forest and logistic regression use `class_weight="balanced"`, XGBoost uses `scale_pos_weight`,
and HistGradientBoosting has neither and runs unweighted.

**"Your GNN comparison is four endpoints and one split. Is that enough?"** No, and the docstring says
so. It bounds the question rather than settling it.

---

## Outstanding items for this document

1. ~~`docs/TECHNICAL_REPORT.md:373` states the wrong active cut~~ **Closed.** The diagram is split
   into its two stages, each labelled with the file that implements it and the cut that file applies,
   and the reliability gate it stated as 0.60 now reads the deployed 0.50. Five tests in
   `TestPipelineDiagramMatchesThePipeline` pin it. Anyone holding a copy of the report from before
   6 September 2026 has the conflated version; section 0 records what it said.
2. The comparison in section 7 leaves HistGradientBoosting unweighted while three of the other four
   families are class-balanced. Either weight it or state the asymmetry in the manuscript's
   model-comparison paragraph.
3. `src/brainsafe/models/` contains `_train_nav17.py` and `_train_ox.py`, whose leading
   underscore marks them as one-off scripts, alongside `train_batch2.py`, `train_new_binders.py` and
   `train_receptor_binders.py`. None is described here. Before submission they should either be
   documented as part of the panel's build history or removed.
