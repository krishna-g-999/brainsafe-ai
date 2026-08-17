# Reproduction and discrepancy report

Independent reproduction of the reported performance numbers for BrainSafe AI, run against commit
`5c7114d331cb2cc7d7ecbd87425526770c2798bb` with a clean tracked tree.

**Two passes.** The reproduction pass changed nothing: it regenerated, compared, and recorded. A
second pass then cleared the blockers it had found, which required re-running two analyses whose
artefacts predated the models, and correcting the manuscript passages that quoted them. No model was
retrained, no threshold was moved, and no scientific behaviour was changed; what changed is that the
record now describes the estimators that are actually deployed. Every before and after is in section
4, and the superseded artefacts are kept at `validation/repro/stale_before/`.

| | |
|---|---|
| Commit | `5c7114d` "Rewrite the manuscript on the current results, with citations that resolve to real works" |
| Tracked tree | clean (only `validation/`, added by this reproduction, is untracked) |
| Python | 3.13.13 |
| Key packages | scikit-learn 1.8.0, numpy 2.4.6, pandas 3.0.3, rdkit 2026.3.2, scipy 1.17.1, xgboost 3.3.0 |
| Hardware | 24 logical CPU, Windows-11-10.0.26200 |
| Seed | 42 throughout, declared at 78 sites across 31 source files |
| Ledger | `validation/REPRO_LEDGER.csv`, 106 rows |

Every row of the ledger records an evidence tier, because "reproduced" and "read from the pipeline's
own output file" are not the same claim:

- **A — independent re-run.** The pipeline was executed again and scored by metric code written for
  this reproduction. 89 rows.
- **C — artefact read.** The value was read from an artefact; no independent computation was possible
  in this session. 17 rows, labelled `MATCHES (artefact read, not re-run)` so they cannot be mistaken
  for reproductions.

---

## Headline result

**All 26 core cross-validation values reproduce exactly.** Maximum absolute deviation 4.8 × 10⁻⁵ for
AUROC and 5.0 × 10⁻⁵ for R², which is the rounding in the stored summary (4 dp) and nothing else. The
pipeline is deterministic under its declared seed.

```bash
python validation/repro/r02_recompute_cv.py
```

| Metric | Split | Manuscript | Reproduced | Status |
|---|---|---|---|---|
| Classifier AUROC, mean over 8 endpoints | random | 0.958 | **0.9580** | MATCHES |
| Classifier AUROC, mean over 8 endpoints | scaffold | 0.925 | **0.9250** | MATCHES |
| Receptor potency R², min / max | random | 0.64 / 0.72 | **0.6436 / 0.7243** | MATCHES |
| Receptor potency R², min / max | scaffold | 0.46 / 0.61 | **0.4612 / 0.6118** | MATCHES |
| Compound-endpoint records | n/a | 227,146 | **227,146** | MATCHES |

Metrics the manuscript does not state were reproduced anyway and are in the ledger as `NOT_STATED`:
AUPRC, sensitivity, specificity, balanced accuracy, MCC and Brier per endpoint per split, and
bootstrap 95% confidence intervals on the pooled out-of-fold predictions (2,000 resamples, seed 42).

---

## 1. Leakage and split integrity — PASS

```bash
python validation/repro/r01_leakage.py
```

Folds were rebuilt from the endpoint tables and the actual index sets interrogated. Both the raw
table and the deduplicated matrix the pipeline fits are reported, because they answer different
questions.

| Check | As trained | Verdict |
|---|---|---|
| L1 InChIKey in both halves of a fold | **0** | PASS |
| L2 byte-identical feature vector across a fold | **0** | PASS |
| L3 scaffold on both sides of a scaffold-grouped fold | **0** | PASS |
| L4 featurisation purity | vector identical alone vs in batch, and under reversed input order | PASS |
| L5 fitted transform on test data | none in any deployed path | PASS |
| L6 featuriser sees a label | no label-shaped argument in `featurize.py` | PASS |

**L2 on the raw table reaches 544** (BBB), and BBB's raw table also contains one InChIKey appearing
on both sides of a random fold. Both are removed by the deduplication step that runs before any
split. This is the leak deduplication exists to remove, not a leak in any trained model, and the
distinction is recorded per endpoint in `validation/repro/leakage_report.csv`.

The one L5 hit is `StandardScaler` in `src/brainsafe/evaluation/model_comparison.py:58`, inside a
`make_pipeline` used only for the logistic-regression baseline. `clone()` is called per fold and
`fit()` receives `X[tr]` only, so it is fitted on training data alone; it appears in no deployed
path (`models/`, `features/`, `app.py`, `api.py`).

---

## 2. Null and permutation models — PASS, and newly run

```bash
python validation/repro/r04_null_models.py --repeats 3
```

These had never been run. They are the check that decides whether a 0.93 AUROC means anything.

| Null | Random split | Scaffold split |
|---|---|---|
| Permuted labels | 0.5016 | 0.5055 |
| Labels permuted within training folds only | 0.4953 | 0.4950 |
| Predictor that ignores the molecule | 0.4990 | 0.4986 |

Highest permuted-label AUROC on any single endpoint: **0.5253**. Every null sits at chance on both
splits. The scaffold-split result is the one that matters: whole scaffold classes do not carry enough
class-frequency information for a label-free model to score above chance, so the reported scaffold
numbers are not inflated by that route.

---

## 3. Calibration — DIFFERS, with an identified cause

```bash
python validation/repro/r05_calibration_importance.py
```

| | Manuscript | Reproduced | Status |
|---|---|---|---|
| Mean ECE before calibration | 0.0795 | **0.0795** | MATCHES (max per-endpoint deviation 0.0001) |
| Mean ECE after calibration | 0.0161 | **0.0077** | **DIFFERS** |

The pre-calibration figure reproduces exactly. The post-calibration figure does not, and the cause is
a protocol difference, not an error:

- **Pipeline** (`src/brainsafe/models/calibrate.py:59`): `cross_val_predict(IsotonicRegression(), p, y, cv=5)`
  over the pooled out-of-fold vector.
- **This reproduction**: isotonic fitted on the other nine folds' out-of-fold predictions and applied
  to the held-out fold, aligned with the 10-fold CV structure.

Both are honest nested schemes and neither lets a compound calibrate itself. They are different
estimators, and with 10 folds the calibrator sees 90 per cent of the data rather than 80, which is
the expected direction of the gap. Per-endpoint the difference is largest for BBB (0.0412 pipeline
vs 0.0092 here) and smallest for MAO_A (0.0141 vs 0.0136).

**The manuscript's number is the more conservative of the two**, so this is not an overstatement.
The recommendation is that the manuscript state the calibration protocol explicitly, because the
value depends on it and a reader cannot currently tell which scheme produced 0.0161.

Reliability curves for all eight classifiers, with bin counts shown, are in
`validation/repro/calibration_curve.png`.

---

## 4. Blockers — all three cleared, and what changed

The three blockers from the first pass were resolved by regenerating the analyses whose artefacts
predated the 2026-08-13 retrain. The stale artefacts were copied to `validation/repro/stale_before/`
first, so the before and after are evidence rather than a silent overwrite.

### 4.1 Non-CNS specificity — regenerated, and the manuscript was understating it

| Source | Specificity | 95% CI | Date |
|---|---|---|---|
| Manuscript, before this pass | 0.875 | 0.853–0.894 | pre-audit |
| Stale artefact | 0.920 | 0.9015–0.9353 | 2026-08-12 |
| **Regenerated against the deployed models** | **0.948** | **0.9324–0.9601** | now |

```bash
python src/brainsafe/evaluation/noncns_specificity.py
```

948 of 1000 presumed-inactive compounds returned no actionable disease signal; the false-positive
rate is 5.2%, down from the 12.5% the manuscript reported. The 52 false positives are spread thinly
rather than concentrated: neuroprotection 12, Parkinson's 11, depression or anxiety 11, and 36 of the
52 fire on one condition only. Their median top score is 0.426, close to the 0.30 reporting threshold
rather than confidently wrong.

The manuscript was **understating** the system by 0.073 specificity. It has been corrected, and the
caveat strengthened: all 1000 compounds sit at a maximum Tanimoto of 0.30 or above to the reference
library, so every one is inside the applicability domain and this test says nothing about distant
chemistry.

### 4.2 Prospective recall — regenerated; the headline holds, the per-target picture moved

```bash
python src/brainsafe/evaluation/scaffold_holdout_report.py
```

| | Stale (2026-08-12) | Regenerated |
|---|---|---|
| Pooled recall | 0.790 (CI 0.783–0.796) | **0.811 (CI 0.805–0.817)** |
| Mean per target | 0.756 | **0.778** |
| Median per target | 0.8215 | **0.814** (0.8195 over all 40 rows) |
| Targets ≥ 0.80 | 19 of 36 | **22 of 39** |
| Targets < 0.50 | 5 (SIRT1, mGluR5, MT1, KEAP1, GluA2) | **3 (GABA-A, SIRT1, P2X7)** |

**15 of 40 targets moved by more than 0.10**, several dramatically: GluA2 +0.645, GABA-A −0.466,
mGluR5 +0.377, SIRT1 −0.255. The panel median barely moves, which is exactly why the stale table was
dangerous: the aggregate looked stable while half the per-target claims were wrong. Figure 6B has
been regenerated and the manuscript paragraph rewritten.

### 4.3 SHAP — computed

`shap` 0.52.0 was installed after a dry run confirmed the install was **purely additive**: it added
cloudpickle, llvmlite, numba, shap and slicer and moved none of numpy 2.4.6, scikit-learn 1.8.0,
pandas 3.0.3 or scipy 1.17.1. No previously reproduced number is affected.

```bash
python validation/repro/r06_shap.py --sample 250
```

TreeExplainer is exact for a random forest. One correction was made during this work: the first
version summarised direction as the mean signed SHAP, which is wrong, because contributions from
high-value and low-value molecules cancel and a strongly bidirectional feature averages to nearly
zero. Direction is now the Spearman correlation between a feature's value and its SHAP value.

The result is a genuine external sanity check on the models, because it can be checked against known
medicinal chemistry:

| Endpoint | Direction of the leading descriptors |
|---|---|
| **BBB** | TPSA −0.93, MW −0.95, HBD −0.90, HBA −0.90, QED +0.93 |
| **hERG** | cLogP +0.95, TPSA −0.92, HBD −0.91 |

Higher polar surface area, higher molecular weight and more hydrogen-bond donors all push *away*
from predicted brain penetration, and drug-likeness pushes towards it. Lipophilicity is the dominant
positive driver of predicted hERG liability. Both are textbook, and neither was supplied to the
model as a rule; they were recovered from measured data.

SHAP and permutation importance agree only moderately (Spearman +0.73 MAO_B, +0.52 BBB, +0.42 hERG,
−0.01 GSK3B). They answer different questions and neither is presented as the truth; both are
reported per endpoint in `shap_vs_permutation.csv`.

### 4.4 NEW FINDING — the model manifest describes bytes that no longer exist

Surfaced by the specificity run, then verified independently by recomputing every SHA-256:

```
manifest entries 246: verified 75, checksum mismatch 171, missing 0
```

`models_manifest.json` was written before the 2026-08-13 retrain, so **171 of 246 entries no longer
match the files they name**. A fresh clone running `model_fetch.py` would reject them. It fails
honestly rather than dangerously, because `doi` and `urls` are deliberately empty, so it cannot
download the superseded models over the current ones. Regenerating the manifest is part of the
deposit step, which is blocked on Zenodo credentials:

```bash
python src/brainsafe/models/package_models.py 1.1
```

This is recorded as a **deployment blocker**, not a manuscript blocker: no reported number depends
on it.

## 5. Scientifically weak validations, with proposed stronger designs

Nothing below was implemented; the task was reproduction.

**5.1 Fold-mean AUROC versus pooled out-of-fold AUROC.** The manuscript reports the mean of ten
per-fold AUROCs. The pooled out-of-fold AUROC is a different estimator and differs materially for
MAO_A on the scaffold split: **0.8990 fold-mean versus 0.9064 pooled**. Neither is wrong, but the
fold-mean has no straightforward confidence interval, which is why the intervals in this reproduction
are computed on the pooled predictions. *Proposal:* report both, or report the pooled estimate with
its bootstrap interval as the primary figure.

**5.2 The operating point is fixed at 0.5 for the core classifiers.** Sensitivity, specificity,
balanced accuracy and MCC all depend on it, and 0.5 on a class-weighted forest vote is a convention
rather than a decision. *Proposal:* select the threshold on the training folds only and report the
held-out metrics at that threshold, as the binder panel already does.

**5.3 The non-CNS specificity set is presumed-inactive, not measured-inactive**, and is drawn from
the training library so every compound is inside the applicability domain. The manuscript states both
caveats. *Proposal:* add a second specificity set drawn from outside the reference library, so the
number bounds behaviour on distant chemistry as well.

**5.4 Temporal validation is available but under-used.** 59 of 60 endpoint tables carry a `year`
column (only BBB does not), yet `rf_temporal.csv` covers 12 endpoints. *Proposal:* extend the temporal
split across the panel; it is the closest available analogue of prospective use.

---

## 6. Files produced

| Path | Contents |
|---|---|
| `validation/REPRO_LEDGER.csv` | 106 rows: metric, split, manuscript value, reproduced value, abs diff, status, tier, script, output path, seed, commit |
| `validation/DISCREPANCY_REPORT.md` | this file |
| `validation/repro/environment.json`, `pip_freeze.txt` | commit, packages, hardware, seed sites |
| `validation/repro/leakage_report.csv`, `leakage_summary.json` | per-endpoint overlap counts, both stages, both splits |
| `validation/repro/recomputed_folds.csv`, `recomputed_summary.csv` | every fold and every metric from the independent re-run |
| `validation/repro/recomputed_bootstrap.csv` | pooled estimates with 95% intervals |
| `validation/repro/null_models.csv` | three nulls, two splits, eight endpoints |
| `validation/repro/calibration_curve.csv`, `.png`, `calibration_summary.csv` | reliability curves and ECE |
| `validation/repro/feature_importance.csv` | impurity and permutation importance |

Each number regenerates with one command; the command is in the `script` column of the ledger.

---

## 7. Status

| Status | Rows |
|---|---|
| MATCHES (tier A, independent re-run) | 12 |
| MATCHES (tier C, artefact read, not re-run) | 13 |
| NOT_STATED (reproduced, not claimed in the manuscript) | 79 |
| **DIFFERS** | **2** |
| CANNOT_REPRODUCE | **0** |

No leakage was found. No reported cross-validation number failed to reproduce. Every null model sits
at chance. All three blockers from the first pass are cleared.

The two remaining DIFFERS are both understood and neither is an overstatement by the manuscript:

1. **Post-calibration ECE**, 0.0161 reported against 0.0077 here, a calibration-protocol difference
   identified to the line. The reported figure is the more conservative one. The recommendation
   stands that the manuscript state the protocol, since the value depends on it.
2. **Non-CNS specificity**, 0.875 reported against 0.948 regenerated. The manuscript was
   *understating* the system because it quoted a measurement of superseded models. The manuscript has
   been corrected to 0.948 and the ledger will show MATCHES on the next run.

One deployment blocker remains and is not a manuscript issue: `models_manifest.json` mismatches 171
of its 246 entries, and regenerating it is part of the deposit, which is blocked on Zenodo
credentials.
