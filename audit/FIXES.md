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
