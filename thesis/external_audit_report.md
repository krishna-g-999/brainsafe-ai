# External audit of the BrainSafe AI codebase

> Written as an external auditor, on 7 September 2026, against commit `d25245b`. Every number in this
> report was read from a file in this repository during the audit or produced by a script written
> during it; none is quoted from the manuscript, from a figure, or from memory. Where I reproduced a
> published result independently I say so and give both numbers. Where I formed a suspicion that the
> evidence then refuted, I have kept the suspicion in the report rather than deleting it, because a
> list containing only confirmed faults tells the reader nothing about how hard anyone looked.

## 1. Scope and method

The live codebase is **197 Python files** once the environment, the packaged copies under
`submission_package/` and `reviewer_package/`, and the `_ARCHIVE_2026-08-10_Brainsafe/` legacy tree
are excluded. The audit covered all of them by inspection, and concentrated on the scientifically
load-bearing subset: data curation (33 files), featurisation (2), model fitting (18), evaluation (33),
the falsification suite (10), the reproduction ladder (7), the independent evidence checks (5), and
the served application.

The method was not a read-through. For each claim that a script makes about its own correctness I
tried to reproduce the number independently, from the raw inputs, with code written for the audit
rather than imported from the project. Seven such reproductions were run. Reporting, figure and deck
builders were inspected but not audited for scientific content, because they render results rather
than produce them.

Three classes of defect were searched for specifically, because they are the ones that survive review:
a guard applied in one place and not in the parallel place; a cached derived value that its inputs
have since invalidated; and a statistic reported outside the regime in which it means what it is read
as. All three were found, one of each.

## 2. The six previously reported defects

All six were re-checked from the artefacts rather than from the commit messages. **Five are correctly
and completely fixed. One has not been addressed.**

### 2.1 Deduplication in the conformal analysis: fixed, and independently reproduced

`rf_conformal_temporal.py` now calls `_dedup_features` before splitting, and `rf_conformal.csv`
records `n_total` alongside `n_test`. BBB reads `n_total` 3,901 against the 7,807 raw rows, which is
the figure Chapter 3 states.

I re-ran the corrected protocol from the raw endpoint tables with independently written code. The
published random-split numbers reproduce **exactly on all eight endpoints**, to the third decimal, for
both coverage and mean set size:

| Endpoint | published coverage | mine | published set size | mine |
|---|---:|---:|---:|---:|
| BBB | 0.896 | 0.896 | 1.215 | 1.215 |
| AChE | 0.933 | 0.933 | 1.045 | 1.045 |
| BChE | 0.905 | 0.905 | 1.017 | 1.017 |
| BACE1 | 0.887 | 0.887 | 0.956 | 0.956 |
| GSK-3B | 0.876 | 0.876 | 0.991 | 0.991 |
| MAO-A | 0.888 | 0.888 | 1.006 | 1.006 |
| MAO-B | 0.914 | 0.914 | 1.052 | 1.052 |
| hERG | 0.903 | 0.903 | 1.094 | 1.094 |

The new columns are internally consistent: `frac_ambiguous + frac_empty + frac_singleton` sums to 1.0
on every row, and `avg_set_size` equals `frac_singleton + 2 x frac_ambiguous` on every row (BBB:
0.7849 + 2 x 0.2151 = 1.2151 against a stated 1.215).

The decision to put the scaffold arm in a separate file rather than as extra rows is correct and the
reason given is verifiable: six consumers key this table by endpoint alone, and a second row per
endpoint would have been silently overwritten or mixed into a min/max.

### 2.2 The set-size inference: fixed

Chapter 4 now states the identity `1 + P(ambiguous) - P(empty)` explicitly, reports both fractions
from the artefact, and says what the earlier edition claimed and why it was wrong. The test-compound
column sums to 8,803, which is the figure the section states.

### 2.3 The fingerprint collision claim: fixed

The `featurize.py` docstring now says the opposite of what it said, which is to say it now says the
truth, and points at an artefact. I measured the collision rate independently on a different 20,000
structure sample: **52,835 distinct atomic environments, all 1,024 bits ambiguous, median 52 per bit,
maximum 74**, against the artefact's 52,882 / 1,024 / 52 / 73. The small differences are the sample;
the conclusion is identical and is not sample-dependent, since 52,000 environments cannot occupy 1,024
bits without collision under any hash.

### 2.4 The underpowered significance test: fixed

`min_attainable_p` is now computed and reported per comparison. I verified the floor independently:
at n = 5 the smallest attainable two-sided p is **0.0625**, at n = 8 it is 0.00781, at n = 13 it is
0.00024. All eight regression rows now carry the verdict "underpowered: no result at this n can reach
alpha" instead of "not distinguishable", which is the honest reading in both directions: it stops the
project claiming a null against XGBoost, which wins 5 of 5, and equally stops it claiming one against
logistic regression, which it beats 5 of 5.

### 2.5 The handicapped booster: fixed

`HistGradientBoostingClassifier` now carries `class_weight="balanced"`, matching the forest's
`class_weight` and XGBoost's `scale_pos_weight`. The comparison is now like-for-like across all four
weightable families. Logistic regression was already correctly wrapped in a `StandardScaler` pipeline,
so that baseline was never unfair.

### 2.6 The applicability domain: NOT addressed

`applicability_domain.py` is unchanged since commit `79d41de`. Its docstring still states the
expectation "the model should be clearly better where it is in domain", and
`applicability_bbb_validation.csv` still shows the opposite ordering, in-domain AUROC 0.7597 on 258
compounds against out-of-domain 0.8059 on 48.

This is the least urgent of the six and the only one that is arguably not a defect at all, but it
should not be left as it stands, because a docstring that states a failed expectation invites a
reader to conclude the component does not work. My measurement was that the difference is not
significant (bootstrap 95% CI on in minus out **[-0.1965, +0.1247]**, 5,000 resamples, so nothing is
established either way at n = 48), while the correlation between maximum Tanimoto and the absolute
probability error is **-0.1422, p = 0.0128**. The domain measure carries information about how far the
probability is from the truth and none about which side of 0.5 it falls on. The fix is to state that
expectation instead, and to test it.

## 3. New findings

### Finding A. The external validation set is not as novel as the flag says (HIGH)

**This is the one finding in this audit that changes a headline number, and it should be fixed before
submission.**

The external barrier validation is reported three ways: on all 306 FDA-curated drugs absent from B3DB
by InChIKey, on the **241** of those "also distinguishable from training in feature space", and on the
**65** that are feature-identical to a training compound. The 241-compound figure is the one that
supports the claim of performance on genuinely novel chemistry, and it is the number Chapter 8, the
viva pack and the manuscript lean on.

The 241/65 split is read from `novel_to_model` and `feature_identical_to_training`, two columns cached
in `data/external/processed/external_bbb_test.csv` by `integrate_external.py`. **Those flags were
computed against a representation the model no longer uses.**

`featurize.parent_mol` desalts *and neutralises* before fingerprinting, a change adopted after the
flags were written and documented in its own docstring (haloperidol hydrochloride scored BBB 0.613
against 0.993 for the free base). The flags, and the independent check in
`audit/evidence/leak_external.py`, both fingerprint the raw SMILES with `Chem.MolFromSmiles`, without
that step. Under the representation the deployed model actually sees, the overlap is larger:

| criterion | overlap with training, of 306 |
|---|---:|
| full InChIKey, the exclusion criterion used | 0 |
| Tanimoto = 1.0 on the **raw** SMILES, the cached flag | 65 |
| InChIKey skeleton, ignoring stereo, salt and protonation | 78 |
| Tanimoto = 1.0 under the **deployed** featuriser | **79** |

Fourteen compounds are flagged `novel_to_model = True` and are nonetheless feature-identical to a
training row. None goes the other way. Recomputing the published metrics with the deployed
representation:

| set | n | AUROC | accuracy | balanced acc. |
|---|---:|---:|---:|---:|
| as published, stored flag | 241 | **0.7934** | 0.7344 | 0.7116 |
| corrected, novel under the current featuriser | 227 | **0.7666** | 0.7225 | 0.6883 |
| the 14 misclassified as novel | 14 | 0.8333 | **0.9286** | 0.9167 |
| all 306 | 306 | 0.7645 | 0.7320 | 0.7119 |

The fourteen carry the signature of memorisation: 92.9 per cent accuracy against 72.3 per cent on the
truly novel remainder. **The strict-novelty AUROC falls from 0.7934 to 0.7666**, which is a drop of
0.0268 and which places it essentially level with the all-306 figure of 0.7645. The correct statement
after this fix is that novelty filtering does not improve the external result, not that it improves it
by about 0.03.

Three things make this worth more than its 0.0268.

The project's own machinery came close to catching it and did not. `leak_external.py` prints the
skeleton overlap of **78**, within one of the right answer, directly above the Tanimoto figure of 65
that was adopted. The stricter number was computed and not used.

The staleness is invisible to the freshness checker. `external_bbb_test.csv` appears in
`check_freshness.py` only at line 308, as an *input* to `models_rf/ad_reference.pkl`. It is never
declared as a derived artefact with inputs of its own, so nothing requires it to be newer than the
featuriser that determines its contents. This is a genuine gap in an otherwise strong guard: the
freshness system protects results but not the cached data that defines which rows a result is computed
over.

It has spread. `validate_inversion.py` reads `novel_to_model` at two places and
`applicability_measures.py` at one, so the falsification suite and the domain analysis inherit the
same 241-compound set.

The fix is contained: re-run `integrate_external.py` so the flags are computed with `parent_mol`,
declare `external_bbb_test.csv` in `check_freshness.py` as derived from the featuriser and the BBB
table, correct `leak_external.py` to use the deployed representation, and update the affected claims.
The verifier will flag them, since they are pinned.

### Finding B. The binder inactives are split at random beside actives split by scaffold (MEDIUM)

`train_binders_hybrid.py` splits the actives by Bemis-Murcko scaffold and gives the reason in its own
comment: "The split is by scaffold rather than at random so a held-out active is a structurally
different compound, not a close analogue of a training one." Three lines earlier it splits the measured
inactives with `rng.permutation` and takes the first half.

That argument applies to the inactives with equal force and more consequence, because the held-out
inactives set the operating threshold as the 0.90 quantile of their scores, and **42 of the 52
endpoints are bound by that requirement rather than by the background pool**
(`final_thresholds.csv`). If a held-out inactive is a memorised analogue of a training inactive it
scores low, drags the quantile down, and leaves the threshold too permissive.

I tested this on the deployed models without retraining, comparing held-out inactives that are close
analogues of a training inactive (maximum Tanimoto at least 0.7) against those that are structurally
distant (below 0.4), across the 40 endpoints with enough of each:

- the median endpoint has **42.8 per cent** of its held-out inactives sitting at Tanimoto 0.7 or above
  to a training inactive, so the structural overlap is substantial and real;
- in **18 of 40** endpoints the close analogues score significantly lower than the distant ones
  (Mann-Whitney, p < 0.05), so the memorisation is measurable, not hypothetical.

**But the expected consequence does not appear.** The median shift in the 0.90 quantile if only
structurally distant inactives were used is **-0.0315**, and it moves upward on only **16 of 40**
endpoints. If anything the current thresholds are marginally conservative rather than permissive. The
reason is that the threshold is an upper-tail quantile while memorisation acts on the bulk of the
distribution, and the per-endpoint shifts are dominated by noise: the distant subsets are small, and
individual shifts range from -0.694 to +0.695, which is not a signal.

So the honest verdict is narrower than the concern that prompted it. **The design is inconsistent and
the project's own stated reasoning condemns it, but no bias in the delivered thresholds is
demonstrated.** It should be corrected for consistency, and the correction should be expected to change
little. It is not a reason to distrust the current panel.

### Finding C. Two scaffold groups can straddle the conformal split boundaries (LOW)

The new scaffold arm orders compounds by shuffled group rank and cuts at fixed indices
`idx[:a], idx[a:b], idx[b:]`. Groups are contiguous in that ordering, but a group spanning index `a`
or `b` is divided across two partitions. At most two groups are affected per endpoint. The exposure is
bounded and small, but it is avoidable by assigning whole groups by cumulative size, and this project
holds itself to that standard elsewhere.

### Finding D. The deployed binder models are calibrated on a random, non-scaffold split (LOW)

`train_binders_hybrid.py` line 182 uses `train_test_split(..., test_size=0.2, stratify=y,
random_state=42)`, fits the forest on the 80 per cent and the Platt sigmoid on the 20 per cent. Two
consequences worth stating rather than fixing: the deployed binder estimator is fitted on 80 per cent
of its available data, and its sigmoid is calibrated on compounds that are analogues of its training
set, so binder probabilities should be expected to be less well calibrated on novel chemistry than the
core panel's isotonic figures suggest. Chapter 10's coverage matrix already records that no binder
endpoint carries a conformal set or an isotonic calibration, so this is consistent with what the thesis
claims; it is the mechanism behind that gap rather than a new one.

### Finding E. The domain analysis scores an undeployed model (INFORMATIONAL)

`applicability_domain.py` loads `BBB.joblib`, not the `BBB_calibrated.joblib` that `app.py` serves.
AUROC is rank-invariant so it is unaffected, but the accuracy column is computed at a 0.5 cut on
uncalibrated probabilities and is therefore not the deployed model's accuracy.

## 4. Suspicions the evidence refuted

Both of these are natural lines of attack and both should be in the candidate's pocket with the number
attached, because the answer is stronger than a denial.

**That the deterministic partition inflates the scaffold results.** `GroupKFold` is unshuffled and
assigns the largest group first, and there is only ever one partition, so no uncertainty from that
choice is reported anywhere. The concern is concrete: the largest scaffold group falls entirely in
fold 1 and constitutes 60 per cent of BBB's fold-1 test set, which is that endpoint's worst fold at
0.8176 against a 0.8777 mean. Fold-level spread is genuinely wider under the scaffold split (mean
per-endpoint standard deviation 0.0237 against 0.0093 random, and BBB ranges 0.8176 to 0.9178).

I re-ran BBB, MAO-A and AChE under `StratifiedGroupKFold` with ten different shuffles. The pooled
estimate is **stable to the third decimal**: BBB 0.8765 with a standard deviation of 0.0008 and a total
spread of 0.0029 across all ten partitions, against the deployed 0.8779. MAO-A 0.9063 with a standard
deviation of 0.0026, AChE 0.9214 with 0.0013. The fold-to-fold variation is real and averages out.
**The headline scaffold numbers do not depend on the partition.**

**That calibration was measured on the easier split.** Calibration quality is reported from the random
out-of-fold predictions while scaffold out-of-fold predictions exist in the same directory. I computed
both. Mean expected calibration error after isotonic regression is **0.0151 under the random split and
0.0169 under the scaffold split**; raw, it is 0.0801 against 0.0833. The choice costs about two parts
in a thousand and can now be defended with a number rather than conceded.

## 5. State of the verification machinery

Run during this audit, on the current tree:

- **81 tests pass**, with 57 subtests, in 97 seconds.
- **Freshness clean**: every declared artefact is newer than its inputs; 252 model entries verified
  with 0 checksum mismatches and 0 missing.
- **Panel consistent**: 52 endpoints, 47 deployed, registry and fitted models agree.
- **229 quantitative claims** checked across the ten chapters and the three aside documents, 0 failures.

The gap Finding A exposes is not in this machinery's execution but in its coverage: it verifies that
artefacts are current with respect to declared inputs, and `external_bbb_test.csv` was never declared.

## 6. Verdict

The codebase is in unusually good scientific order. The evidence for that is not that I found little,
but the character of what I found: after a deliberate search for leakage, circularity and misapplied
statistics across 197 files, the substantive defects reduce to one stale cached flag, one design
asymmetry with no demonstrated consequence, and three minor items. Every serious guard this project
needs is present somewhere in it, including several that most work of this kind lacks entirely:
feature-level deduplication before splitting, a background pool partitioned three ways so no compound
can serve two roles, thresholds set from data the model never saw, an explicit permutation null, a
falsification suite that has refuted four of its own hypotheses, and an independent reproduction
ladder.

The recurring failure mode is not carelessness but incomplete propagation: a guard is invented, argued
for correctly, and then applied to one of the two places that need it. Finding A, Finding B and the
original conformal defect are all the same shape. That is worth saying to a committee plainly, because
it is a much better answer than treating each as an isolated slip.

**One item, Finding A, must be fixed before submission**, because it changes a number the manuscript
relies on and because the correction is small, contained and in the direction that costs the project
0.0268 of AUROC. Everything else can be scheduled.

---

## Outstanding items

- Finding A: recompute the novelty flags with the deployed featuriser and update the affected claims.
- Finding B: split the binder measured inactives by scaffold and re-derive the 42 affected thresholds.
- Finding 2.6: restate the applicability-domain expectation as calibration rather than discrimination,
  and test it.
- Findings C, D and E: correct at convenience.
- `check_freshness.py`: declare `external_bbb_test.csv` as a derived artefact.
