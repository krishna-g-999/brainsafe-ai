# Audit fix log

Running record of Critical findings from [`AUDIT_REPORT.md`](AUDIT_REPORT.md) as they are fixed.
Each entry records what changed, the evidence that it works, and the before/after effect on any
scientific output. Fixes are validated in an isolated copy of the tree, never against the live data.

---

## BS-C-12 — the documented rebuild command silently destroyed the training data

**Status: FIXED** · `src/brainsafe/data/rebuild_endpoints.py`

### The defect

On a fresh clone the API caches are absent (they are gitignored). `chembl_compound_level` returned an
**empty DataFrame** instead of raising, and `main()` wrote every result unconditionally, so the
documented command overwrote all 11 core endpoint CSVs with header-only files and **exited 0**.

### What changed

| Change | Rationale |
|---|---|
| `missing_inputs()` preflight in `main()` | Every required cache is checked **before** anything is touched. Exits 1 with the list of missing files and an explicit statement that nothing was written. |
| `chembl_compound_level` / `bindingdb_compound_level` now raise `FileNotFoundError` | A missing source is an error, not an empty result. |
| `check_plausible()` with `MIN_ROWS = 100` | A rebuild producing fewer than 100 rows is a broken cache, not a result. The smallest legitimate table is BChE at 2,621 rows. |
| `check_plausible()` with `SHRINK_LIMIT = 0.5` plus `--allow-shrink` | Refuses to replace a table with one less than half its current size unless the operator says the reduction is intended. |
| All tables built and checked in memory before any write | A failure part way through can no longer leave `data/endpoints/` half rebuilt and half stale. |
| Backup is now unconditional and timestamped | The old `if not backup.exists()` guard silently skipped on every run after the first, and the fixed name `endpoints_chembl_only_2026-07-21` described contents it no longer held. |

### Validation

All tests run in isolated copies under the scratchpad. The live tree was never used as a target.

**1. Fresh-clone scenario (the defect).**
```
$ python <isolated>/src/brainsafe/data/rebuild_endpoints.py     # no caches present
rebuild_endpoints.py rewrites the training tables and cannot run without the cached
source responses. Missing:
  data/_chembl_cache/  (the whole cache directory is absent)
  data/_bindingdb_cache/  (the whole cache directory is absent)
...Nothing has been written and data/endpoints/ is untouched.

REAL EXIT CODE = 1
tables emptied: 0
```
Before: 11 tables emptied, exit 0. After: 0 tables emptied, exit 1.

**2. Partial-cache scenario (`AChE_y.json` truncated from 8,189 to 800 records).**
```
[AChE] the rebuild produced 1758 rows against 4387 already on disk, a drop of more than
50 per cent. Refusing to overwrite. Re-run with --allow-shrink if the reduction is intended.
Nothing was written.
EXIT=1
AChE unchanged on disk: True
backup dir created: False
```
The refusal happens before the backup step, so a rejected run leaves no trace at all.

**3. Happy path is behaviour-preserving.** The unmodified script from `HEAD` and the fixed script
were run on the same caches in two separate isolated trees:

```
hERG   live=211005d25ad445f2  ORIGINAL=51c3560e745aecf2  MINE=51c3560e745aecf2   original==mine: True
AChE   live=c03dedbe4fc8a4ca  ORIGINAL=c03dedbe4fc8a4ca  MINE=c03dedbe4fc8a4ca   original==mine: True
D2     live=08b4489fcd4e9d50  ORIGINAL=08b4489fcd4e9d50  MINE=08b4489fcd4e9d50   original==mine: True
```
**Byte-identical output to the pre-fix script on every endpoint.** Row counts unchanged: AChE 4,387;
BChE 2,621; BACE1 8,501; GSK3B 4,958; MAO_B 3,665; MAO_A 2,228; D2 7,734; A2A 6,785; HT2A 5,989;
SERT 4,572; hERG 5,875.

### Scientific behaviour: unchanged

No label, potency or row changed. The guards are refusals only.

---

## Incidental fix — `.round()` silently skipped for hERG

**Status: FIXED** (same file, separate commit) · **New finding, not in the original report**

Found while validating BS-C-12. `data/_bindingdb_cache/hERG_labelled.csv` is header-only (31 bytes),
because BindingDB contributes nothing for hERG. Reading it yields a `pchembl` column of dtype
**object**; `pd.concat` propagates that to the pooled frame; `groupby().median()` keeps it; and
`.round({"pchembl": 3})` on an object column is a **silent no-op**.

Confirmed directly under pandas 3.0.3:
```
empty bd rows: 0   bd.pchembl dtype: object
after concat, pchembl dtype: object
after groupby, pchembl dtype: object
value after .round(3): 4.6899999999999995        <- round did nothing
if cast to float first, .round(3) gives: 4.69
```

hERG was therefore the one endpoint of eleven written with unrounded potencies:

| File | Rows | `pchembl` values with more than 3 decimal places |
|---|---|---|
| live `AChE`, `D2`, `BACE1`, `MAO_A` | — | **0** |
| live `hERG` | 5,875 | **66** |
| `hERG` rebuilt by the pre-fix script | 5,875 | **76** |
| `hERG` rebuilt after the fix | 5,875 | **0** |

**Fix:** coerce `pchembl` and `year` with `pd.to_numeric` after the concat, so every target is pooled
and rounded on the same numeric dtype.

**Scientific impact: none.** Comparing the live hERG table with the rebuilt one:
```
rows live/rebuilt: 5875 5875     same SMILES set: True     labels identical: True
max |pchembl difference|: 1.776e-15
rows where the label would change: 0
```
This is text formatting, not chemistry.

### Open provenance discrepancy (recorded, not resolved)

The committed `data/endpoints/hERG.csv` carries **66** unrounded values; the committed code run on the
committed caches produces **76**. **Neither the pre-fix nor the post-fix script reproduces the shipped
hERG table byte-for-byte**, so that file was written by a version of the pipeline, a pandas release, or
a cache state that is no longer present. The other ten tables reproduce exactly.

This is a concrete instance of **BS-C-13** (the retrieval scripts were deleted and the caches are not
committed) and should be resolved when BS-C-13 is addressed: after regenerating the caches from a
restored fetcher, `hERG.csv` should be rebuilt so that the shipped file and the code agree.

---

## BS-C-13 — the core tables could not be regenerated

**Status: FIXED, apart from the deposit itself** · `src/brainsafe/data/fetch_{endpoints,antioxidant,clinical}.py`,
`src/brainsafe/data/package_caches.py`, `.gitignore`, `data/raw/measured_endpoints_SOURCE.md`

### The defect

`measured_endpoints_SOURCE.md` cited `BS_fetch_endpoints.py`, `BS_fetch_antioxidant.py` and
`BS_fetch_clinical.py` as the retrievers for the ChEMBL activities, B3DB, the DPPH set and the
clinical reference. All three were deleted in commit `fea5029`, and the caches they produce were
gitignored, so no reader could get from a public database to a training table.

### What changed

**1. The three scripts are restored** from `fea5029^` and adapted to the current layout: paths derive
from the repository root rather than `os.chdir` to the script's own directory, which would now
resolve inside `src/`. Queries, kept standard types, labelling and aggregation are unchanged.

`fetch_endpoints.py` deliberately stops at the cache for the eleven targets. The archived version
wrote `data/endpoints/<target>.csv` from ChEMBL alone; those tables are now produced by
`rebuild_endpoints.py`, which pools ChEMBL with BindingDB, so writing them here would replace the
pooled tables with ChEMBL-only ones. BBB is still written here, B3DB being its only source. Each
script refuses to overwrite a populated output with an empty or halved one, the page and assay caps
warn instead of truncating silently, and failed request batches are counted rather than passed over.

**2. The caches are committed.** 77 files, 31 MB, previously ignored. They are the only record of
what the three sources returned on the recorded retrieval dates; ChEMBL serves the current release,
so refetching after the next one returns different rows and nothing would reveal the drift. `.git`
grew from 121 MB to 119 MB (the JSON compresses).

**3. They are also packaged for deposit.** `package_caches.py` mirrors `package_models.py`: one
archive plus `caches_manifest.json` recording the archive SHA-256 and the size and checksum of every
file. Because the caches are also tracked, `--verify` doubles as an integrity check on the working
tree.

### Validation

```
$ python src/brainsafe/data/package_caches.py
wrote dist/brainsafe_source_caches_v1.0.tar.gz (4.5 MB, 77 files)
wrote caches_manifest.json: archive sha256 408cb275e331b99d...

$ python src/brainsafe/data/package_caches.py --verify
OK: every cached response matches its published checksum.          EXIT=0

# after flipping a single bit in AChE_y.json, in an isolated copy:
  changed: 1
    data/_chembl_cache/AChE_y.json
FAILED: the caches on disk are not the caches that were published. EXIT=1
```
Offline path checked: `fetch_target_activities("AChE", "CHEMBL220")` reads 8,189 records from cache
with no network call, and the record schema (`smiles`, `pchembl`, `year`) is what
`rebuild_endpoints.py` consumes.

### Scientific behaviour: unchanged
No table, model or number was regenerated.

### Remaining
Depositing `dist/brainsafe_source_caches_v1.0.tar.gz` and recording its DOI and URLs in the `doi`
and `urls` fields of `caches_manifest.json`. That needs credentials this session does not have.

---

## BS-C-14 — BindingDB censored measurements were treated as exact

**Status: FIXED in code; shipped tables not yet regenerated** · `src/brainsafe/data/fetch_bindingdb.py`

### The defect

`parse_affinity` applied `lstrip("><~=")` before converting, so `">10000"` became precisely 10 uM and
`"<1"` precisely 1 nM. 1,141 of 24,014 cached records (4.75 per cent) are censored.

### What changed

The relation is returned with the value. Only exact measurements contribute to the pooled potency.
Censored records go to a companion `_censored.csv` with their bound and relation, marked decisive
where the bound falls outside the grey zone. `~` is dropped.

### Quantified impact, and it is smaller than the defect suggests

| Measure | Result |
|---|---|
| Censored records previously read as exact | **1,141 of 24,014 (4.75%)** |
| Compounds whose **label** changes | **0** |
| Compounds present both ways | 21,791 |
| Of those, pooled potency shifted at all | **6 (0.03%)**, two by more than 0.5 log unit |
| Compounds leaving the BindingDB pool (censored-only) | 1,131 |
| Censored records that are decisive | 849, **all actives**, zero inactives |

Endpoint tables, rebuilt in an isolated tree:

| Endpoint | Before | After | Δ | % |
|---|---|---|---|---|
| GSK3B | 4,958 | 4,577 | **−381** | −7.68 |
| A2A | 6,785 | 6,421 | **−364** | −5.36 |
| HT2A | 5,989 | 5,836 | −153 | −2.55 |
| SERT | 4,572 | 4,507 | −65 | −1.42 |
| D2 | 7,734 | 7,671 | −63 | −0.81 |
| MAO_B | 3,665 | 3,616 | −49 | −1.34 |
| BACE1 | 8,501 | 8,489 | −12 | −0.14 |
| AChE | 4,387 | 4,380 | −7 | −0.16 |
| BChE | 2,621 | 2,615 | −6 | −0.23 |
| MAO_A | 2,228 | 2,224 | −4 | −0.18 |
| hERG, BBB | unchanged | | 0 | 0.00 |
| **Total** | **65,122** | **64,018** | **−1,104** | **−1.70** |

### Correction to the audit's own framing

`AUDIT_REPORT.md` presented BS-C-14 alongside BS-C-15 as losing the measured negative class. For
BindingDB that is **not** what the data show: the censored records are overwhelmingly `<` bounds,
compounds too potent to measure rather than too weak. Honouring the relation recovers no inactives
here, and **the shipped labels were never wrong because of this**. The defect is real and worth
fixing, but its practical severity is far below Critical. BS-C-15, which concerns the ChEMBL
`pchembl_value` filter, is a separate and still-open matter.

### Awaiting approval
The shipped `data/endpoints/*.csv` are untouched. Regenerating them costs 1,104 rows and would
require retraining every affected model.

---

## BS-C-06 — feature-identical duplicates split across cross-validation folds

**Status: FIXED in code; shipped models not retrained** · `src/brainsafe/models/train_rf.py`,
`src/brainsafe/features/featurize.py`

### The defect

Nothing deduplicated the endpoint tables before splitting. The featuriser folds a molecule to a
stereo-blind 1024-bit fingerprint over its desalted parent, so stereoisomers, salt forms and
protonation variants of one compound become byte-identical vectors, and `StratifiedKFold` has no way
to keep them together. Invisible at the SMILES level: all 7,807 BBB strings are distinct, with zero
string-level duplicates and zero string-level conflicts.

### What changed

`_dedup_features` collapses rows by feature vector before the split **and before the deployed
refit**. Classification groups whose members disagree on the label are dropped (288 across the
panel); regression takes the group median. Two related repairs in the same path: the scaffold is now
computed on the same desalted parent the featuriser uses, and acyclic compounds share one group
instead of each receiving their own (which had turned the scaffold split into a random split for
that part of the set). `featurize._mol_from_smiles` became public as `parent_mol` so the two agree by
construction.

### Quantified impact: all 13 endpoints, both splits

| | Before | After |
|---|---|---|
| Random 10-fold AUROC (8 classifiers) | 0.947–0.969, mean **0.9604** | 0.899–0.979, mean **0.9534** |
| Scaffold 10-fold AUROC | 0.868–0.956, mean **0.9186** | 0.874–0.967, mean **0.9168** |
| Random 10-fold R² (5 regressors) | mean 0.6381 | mean 0.6377 |
| Scaffold 10-fold R² | mean 0.4743 | mean **0.4842** |

Per endpoint, random split:

| Endpoint | n before | n after | Before | After | Δ |
|---|---|---|---|---|---|
| **BBB** | 7,805 | **3,901** | 0.9605 | **0.8990** | **−0.0615** |
| BACE1 | 8,501 | 7,832 | 0.9667 | 0.9787 | **+0.0120** |
| hERG | 5,875 | 5,706 | 0.9541 | 0.9516 | −0.0025 |
| BChE | 2,621 | 2,533 | 0.9684 | 0.9661 | −0.0023 |
| AChE | 4,387 | 4,216 | 0.9626 | 0.9619 | −0.0007 |
| MAO_A, MAO_B, GSK3B | | | | | ≤0.0005 |
| SERT (R²) | 4,572 | 4,142 | 0.6016 | 0.6216 | **+0.0200** |

Scaffold split: BBB 0.9197 → **0.8777** (−0.0420); every other classifier moves by ≤0.0114 and most
move **up**; all four receptor regressions improve.

### The honest reading

**The duplication problem is BBB, and only BBB.** Excluding it, the panel after deduplication is
0.946–0.979 mean **0.9612** random and 0.874–0.967 mean **0.9224** scaffold, which is at or above
what was published. BACE1 improves because 41 contradictory groups had been pulling against it.

So the README's "Random 10-fold 0.95–0.97 (mean 0.96)" and "Scaffold 0.87–0.96 (mean 0.92)" remain
defensible **for the target panel**. What cannot stand is the BBB figure, and BBB is the gate that
multiplies into every per-disease score.

### Awaiting approval
No model retrained, no results table rewritten.

---

## BS-C-16 — expansion fetchers deduplicated on the raw SMILES string

**Status: FIXED in code; shipped tables not regenerated** · `build_compound_library.py` and seven fetchers

### What changed

`add_parent_key` standardises to the desalted parent and attaches its InChIKey; `fetch_batch2-5`,
`fetch_new_targets`, `fetch_readacross_targets` and `fetch_pka` now group on that instead of the raw
SMILES. It deliberately does **not** merge stereoisomers: enantiomers are different compounds whose
potencies can differ by orders of magnitude. Verified on a constructed case: caffeine and caffeine
HCl merge; two alanine enantiomers stay separate.

### Correction to the audit's own claim

`AUDIT_REPORT.md` called this **the root cause of BS-C-06**. That is wrong. Applying parent-InChIKey
keying to the shipped tables removes **411 rows of 198,499 (0.21%)** across 31 endpoints, not the
13,846 duplicates BS-C-06 concerns. The 4,012 salted SMILES mostly carry a parent appearing nowhere
else, so collapsing them creates few duplicates.

The two findings are independent. BS-C-06 exists because the featuriser is stereo-blind while the
tables correctly keep stereoisomers distinct, and no InChIKey policy would change that. Severity here
should be read as **Major**, not Critical. The report has been corrected.

---

## BS-C-01 to BS-C-05 — the threshold and specificity family

**Status: FIXED in code; shipped models and thresholds not regenerated** ·
`models/pools.py` (new), `train_binders_hybrid.py`, `final_thresholds.py`,
`train_measured_label_holdout.py`, `noncns_specificity_fast.py`, `build_reviewer_package.py`,
`_train_nav17.py`, `_train_ox.py`, `_train_small.py`, `.gitignore`

### The defect, as one thing rather than five

The pipeline asked the background library for three different jobs and took all three from the same
place. Decoys were **trained on** as negatives; the threshold was a **quantile** of a sample from the
same library; the false-positive rate was then **reported on that same sample**, which makes it the
quantile restated. Separately, `final_thresholds.py` re-read the endpoint table and set the threshold
from *every* measured inactive, including the half training had withheld, on **40 of 49** deployed
endpoints. And every reported sensitivity was measured on positives the model had been fitted to,
beneath five docstrings asserting the opposite.

### What changed

`models/pools.py` splits the applicability-domain library once into **disjoint** `decoy` (95,515),
`threshold` (31,694) and `evaluation` (31,681) pools, by a stable hash of the structure so the
partition survives library growth and does not depend on file order. Verified: pairwise overlap 0.

| Fix | Where |
|---|---|
| Decoys drawn only from the decoy pool, and the target's own measured inactives excluded from decoy eligibility | `train_binders_hybrid.py` |
| A fifth of **active scaffold groups** withheld and never trained on; sensitivity reported on those | `train_binders_hybrid.py`, `train_measured_label_holdout.py` |
| Both halves of every split persisted to `models_rf/holdout/<T>_binder_holdout.json` | both training scripts |
| Threshold read from the persisted holdout, never from the endpoint table; endpoints with no record are named in a warning rather than silently thresholded on seen data | `final_thresholds.py` |
| False-positive rate measured on the **evaluation** pool, disjoint from decoys and from the threshold pool; the in-sample value retained under an honest name | `final_thresholds.py`, `train_binders_hybrid.py` |
| Non-CNS negatives drawn from the 8,418 approved drugs **absent** from the AD reference (overlap verified 0) instead of from the reference itself | `noncns_specificity_fast.py` |
| Full exclusion set canonicalised instead of an order-dependent 40,000 slice | `noncns_specificity_fast.py` |
| Three verbatim clones reduced to wrappers, and un-ignored so they exist in a clone at all | `_train_nav17/_train_ox/_train_small.py`, `.gitignore` |

### Quantified on three endpoints, retrained in an isolated tree

| | COX2 | OPRM1 | LRRK2 |
|---|---|---|---|
| Inactives setting the threshold | 693 → **347** | 366 → **183** | 18 → **9, below minimum** |
| Actives train / holdout | all → 967 / **215** | all → 3003 / **715** | all → 988 / 185 |
| **Sensitivity** | 0.752 → **0.391** | 0.972 → **0.906** | 0.974 → **cannot report** |
| AUROC vs measured inactives | 0.923 → **0.767** | 0.981 → **0.951** | — |
| Threshold | 0.909 → 0.920 | 0.849 → 0.921 | 0.872 → **0.40 fallback** |
| FPR on the threshold set (in-sample) | 0.101 → 0.101 | 0.104 → 0.104 | — |
| **FPR on the disjoint pool** (n=31,681) | not measured → **0.0087** | not measured → **0.041** | — |

### The honest reading: sensitivity was inflated, specificity was not

The **sensitivity** claims do not survive. COX2 loses 0.36, and LRRK2 turns out not to support a
data-driven threshold at all: only 9 of its 18 measured inactives are held out, below the minimum, so
it falls back to 0.40 and can report no sensitivity, where it previously reported **0.974** from a
threshold set on all 18.

The **specificity** claims essentially do. Measured on the disjoint evaluation pool, COX2 is 0.0087
and OPRM1 0.041, against 0.0087 and 0.0413 as shipped. Sampling 3,000 compounds from a library of
158,890 overlapped the decoy set too little for the reported rate to have been materially inflated.
BS-C-03 and BS-C-04 were real defects of method with, on this evidence, close to no numerical
consequence; **BS-C-01 and BS-C-02 had a large one.**

Note the drop combines two effects that this protocol cannot separate: measuring on unseen actives
rather than trained ones, and training on 80 per cent of active scaffolds instead of all of them.
Both are properties of an honest protocol, but the second is a real reduction in training data.

### Awaiting approval, and still open
Three endpoints of 49 were retrained, to measure rather than to replace. Nothing shipped was
regenerated. `train_receptor_binders.py`, `train_new_binders.py` and `train_batch2.py` still report
sensitivity on training positives and still draw decoys from the whole library; they need the same
treatment before the panel is consistent.

---

## BS-C-02 / BS-C-04 — the remaining three training scripts

**Status: FIXED in code; shipped models not retrained** · `train_receptor_binders.py`,
`train_new_binders.py`, `train_batch2.py`, `models/pools.py`

### Why this mattered beyond the three files

These twelve endpoints were left on the old basis after the rest of the panel moved off it: every
active trained on and then scored, decoys from the whole background library, and the false-positive
rate sampled from that same library. A panel mean that mixes two bases is worse than either alone,
because no single sentence describes it.

### What changed

All three now draw decoys from the `decoy` pool, measure the background rate on the disjoint
`evaluation` pool, withhold a fifth of the active scaffold groups, report sensitivity on those, and
persist the holdout so `final_thresholds.py` honours it.

`scaffold_holdout` and `write_holdout` moved into `pools.py` instead of being written three more
times. `scaffold_holdout` withholds whole groups and returns an empty mask on degenerate scaffold
structure, rather than silently withholding everything or nothing.

### A further defect, fixed here and not previously recorded

All three computed the hard-decoy AUROC from `Xc[yc == 1]`, under a comment describing it as
held-out actives:

```python
pa = cal.predict_proba(Xc[yc == 1])[:, 1]     # "held-out actives"
```

`Xc` is the calibration split and `cal` was fitted on it, so those positives were not held out in any
sense. Every `auroc_hard_decoys` in the panel is optimistic for this reason, independently of
BS-C-02. It now reads the withheld actives.

### Validation

Exercised end to end against a shrunken background pool, so the logic is tested in seconds rather
than after featurising 95,000 decoys:

```
n_active                 513          n_active_holdout    86
auroc_hard_decoys        1.0          sensitivity_at_0.5  1.0
sensitivity_basis        held_out_actives_by_scaffold
background_fpr_basis     held_out_evaluation_pool
```
The holdout is written and read back, the metric fields are emitted, and the evaluation pool is the
one sampled. **Those metric values come from a truncated pool and are not real numbers** — the test
proves the plumbing, not the performance.

The test earned its place: it caught a bug introduced in this very change.
`background_pools(with_fingerprints=True)` returns `(smiles, fingerprints)`, and all three scripts
initially indexed the tuple rather than unpacking it, which would have failed on the first real run.

### Still open
The panel is now internally consistent in code, but **no model has been retrained**, so every shipped
`sensitivity_at_threshold` and `auroc_hard_decoys` still rests on the old basis. On the three
endpoints measured under BS-C-01 to BS-C-05, sensitivity fell by between 0.07 and 0.36. Table 2 and
the abstract's panel-mean sensitivity will move once the full panel is regenerated.
