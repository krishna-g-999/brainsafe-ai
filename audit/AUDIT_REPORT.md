# BrainSafe AI — Ranked Audit Report

**Audit date:** 2026-08-10
**Repository:** `D:\BRAINSAFE_AI`
**Commit audited:** `5d338d6ef8958c909ce8a82e1ed7cf898a83f5c4` (branch `main`)
**Working tree:** 2 modified, uncommitted files (`manuscript/NAR_WebServer_BrainSafe_built.md`, `..._draft.md`)
**Purpose:** pre-submission diagnostic for the Nucleic Acids Research Web Server Issue
**Scope:** diagnostic only. **No pre-existing repository file was modified.** Everything under `audit/` is new.

A one-page summary of the Critical and Major findings is in
[`audit/EXECUTIVE_SUMMARY.md`](EXECUTIVE_SUMMARY.md).

---

## 0. Ground truth: environment, install, smoke run

Every number below was regenerated. Nothing is quoted from memory, and nothing is taken from the
repository's own summaries without independent recomputation.

### 0.1 Environment

| Item | Value |
|---|---|
| Platform | `Windows-11-10.0.26200-SP0`; Intel64 Family 6 Model 198 Stepping 2; 24 logical CPUs |
| Interpreter for all runs | `brainsafe_env/Scripts/python.exe`, Python 3.13.13 |
| RDKit / scikit-learn / NumPy / pandas / SciPy | 2026.03.2 / 1.8.0 / 2.4.6 / 3.0.3 / 1.17.1 |
| Full package list | [`evidence/venv_freeze.txt`](evidence/venv_freeze.txt) (104 packages) |
| Seed | 42 (`src/brainsafe/models/train_rf.py:46`) |

`brainsafe_env` matches `requirements.txt` exactly. The machine's **base** interpreter does not
(scikit-learn 1.9.0, pandas 2.3.3, SciPy 1.18.0, RDKit 2026.03.4). Loading a deployed model under the
base interpreter raises `InconsistentVersionWarning`; under `brainsafe_env` it loads cleanly. All
results in this report were produced under the pinned environment.

### 0.2 Clean install — **PASSES**

```
python -m venv cleanenv
cleanenv/Scripts/python.exe -m pip install -r requirements.txt      → EXITCODE=0
→ rdkit 2026.03.2 · sklearn 1.8.0 · numpy 2.4.6 · pandas 3.0.3 · scipy 1.17.1 · xgboost 3.4.0
```

### 0.3 Smoke run — **PASSES**

```
brainsafe_env/Scripts/python.exe src/brainsafe/evaluation/app_health.py
→ EXIT=0 · "app.py is healthy: all checks passed" · 20 checks
```
Captured output: [`evidence/health_full.txt`](evidence/health_full.txt).

**There is no Critical finding arising from install or smoke run.** The pipeline installs, the
application starts, all models load, and the probe compounds produce directionally correct
pharmacology. That is a healthy starting position and should be stated before anything else.

### 0.4 Deployed-model integrity — **PASSES**

All 195 files declared in `models_manifest.json` are present with matching byte sizes (0 missing,
0 mismatched).

### 0.5 Reproducibility of the headline metrics — **EXACT**

| Endpoint / split | Published | Reproduced | Δ | Verdict |
|---|---|---|---|---|
| BBB random 10-fold | 0.9605 | 0.9605 | −0.0000 | **MATCHES** |
| AChE random 10-fold | 0.9626 | 0.9626 | −0.0000 | **MATCHES** |
| MAO_A random 10-fold | 0.9467 | 0.9467 | −0.0000 | **MATCHES** |
| hERG random 10-fold | 0.9541 | 0.9541 | +0.0000 | **MATCHES** |
| BBB scaffold 10-fold | 0.9197 | 0.9197 | 0.0000 | **MATCHES** |

Evidence: [`evidence/repro_cv.py`](evidence/repro_cv.py), [`repro_cv.csv`](evidence/repro_cv.csv),
[`repro_scaffold.py`](evidence/repro_scaffold.py).

**This is the strongest thing in the repository.** The numbers are real, the seed is honoured, the
pipeline is deterministic. **Nothing in this report alleges fabrication.** Every figure I tested
regenerated exactly. The findings are about what the numbers *measure*.

### 0.6 Also verified as MATCHING

Independently recomputed and confirmed against `README.md`:

- Mean ECE 0.0721 → 0.0127 (claimed 0.072 → 0.012; the second rounds to 0.013, see BS-N-06)
- Temporal classifier AUROC range 0.611–0.908 (claimed 0.61–0.91)
- External BBB AUROC 0.7741 on n=306 (claimed 0.774)
- Conformal empirical coverage 0.890–0.918 at target 0.90 (claimed 0.89–0.92)
- 914 live files / 118 scripts against `REPOSITORY_INVENTORY.csv` and `SCRIPT_INDEX.csv`

---

## Severity definitions

| Severity | Meaning |
|---|---|
| **Critical** | Undermines a central scientific claim, or a reported number does not measure what it is said to measure. Must be fixed before submission. |
| **Major** | Materially misleads, breaks a documented workflow, or is a compliance gap that will draw a revision request. |
| **Minor** | Real defect, limited blast radius. |
| **Trivial** | Cosmetic, hygiene, drift. |

**Confidence tags.** `[verified]` = I ran the command or read the code myself and reproduce the
finding. `[reported]` = surfaced by a component audit, code quoted, but not independently re-executed
by me. Every Critical finding is `[verified]`.

---

# CRITICAL FINDINGS

The Criticals fall into three families: **thresholds and specificity are measured on data they were
derived from** (C-01 to C-05), **compounds leak between train and test** (C-06 to C-08), and **the
checks that should have caught this cannot fail** (C-09, C-10). C-11 is a submission blocker.

---

## BS-C-01 · Deployed decision thresholds are set using the measured inactives the model was trained on `[verified]`

**Component:** threshold calibration

`train_binders_hybrid.py` carefully splits measured inactives in half and sets the threshold on the
**held-out** half. `final_thresholds.py` then re-reads the endpoint CSV, takes **all** inactives, and
overwrites that threshold — silently undoing the holdout.

**Evidence.** `src/brainsafe/models/final_thresholds.py:69-73`
```python
ina = df.loc[df["label"] == 0, "smiles"].astype(str).tolist()
if len(ina) >= 15:
    Xi, _ = featurize(ina)
    pi = mdl.predict_proba(Xi)[:, 1]
    thr_in = float(np.quantile(pi, 1.0 - TARGET_FPR))
```
No holdout, no exclusion, no reference to `ina_hold`.

Proof from the shipped artefact — `models_rf/binder_modes.json`, endpoint `A1`:
```
n_measured_inactive_train    = 109      <- the half used for training
n_measured_inactive_holdout  = 110      <- the honest half
n_measured_inactive_for_threshold = 219 <- 109 + 110: BOTH halves were used
threshold       = 0.9735
threshold_basis = measured_inactives_and_background
```
`data/endpoints/A1.csv` contains exactly **219** rows with `label == 0`. The shipped threshold was
derived from every one of them, including the 109 the model was fitted on.

Panel-wide `threshold_basis` distribution (recomputed from the shipped JSON):
```
measured_inactives_and_background : 40      <- written only by final_thresholds.py
held_out_measured_inactives       :  5
random_chemistry_specificity      :  3
global_fallback                   :  1      <- CGRP ships at the untested fallback 0.40
```

**40 of 49 deployed endpoints carry a threshold contaminated this way.**

**Severity:** Critical. The threshold is what converts a probability into the "engaged / not engaged"
call the interface reports and every disease score consumes.

**Fix.** Persist `ina_hold` per endpoint at training time and have `final_thresholds.py` read it
instead of `df.label == 0`. Re-derive every shipped threshold and every downstream number.

---

## BS-C-02 · Every reported `sensitivity_at_threshold` is measured on training positives, while five docstrings state the opposite `[reported, code verified]`

**Component:** binder training and threshold calibration

The docstring of `src/brainsafe/models/train_binders_hybrid.py:15-16` says:

> *"Metrics are reported on the held-out measured inactives only, so no compound used for training is
> used to judge the model."*

The code trains on `act` and then evaluates on the same `act`:
```python
# train_binders_hybrid.py:115   -- act goes INTO the training set
smiles = act + dec + ina_train
# :139                          -- aX is featurised from that same act
pa = cal.predict_proba(aX)[:, 1]
# :147
"sensitivity_at_threshold": round(float((pa >= thr).mean()), 3)
```

Only the **negative** class is held out. The positive class is 100% training data, so
`sensitivity_at_threshold` is fully in-sample and `auroc_vs_measured_inactives` is half in-sample
(training positives against held-out negatives).

The same pattern appears in `_train_nav17.py:128,136`, `_train_ox.py:128,136`,
`train_measured_label_holdout.py:49,65,77`, `_train_small.py:49,65,77`, and is re-measured the same
way by `final_thresholds.py:66-82`, `calibrate_background_specificity.py:79-88`,
`screening_thresholds.py:70-76` and `calibrate_binder_thresholds.py:57-77`.

`train_measured_label_holdout.py:5-8` is the sharpest case: its docstring exists specifically to claim
it *fixes* leakage, and it prints `"(previously reported ..., leaked)"` at line 82, while carrying the
identical defect.

**Severity:** Critical. `sensitivity_at_threshold` propagates into `binder_modes.json`,
`final_thresholds.csv`, `screening_thresholds.csv`, the application's audit page, and the manuscript's
panel-mean sensitivity.

**Fix.** Hold out a stratified fraction of positives at training time and evaluate on it. Every
sensitivity figure in the paper needs recomputation. **Correct the five docstrings first** — a
reviewer who reads them and then reads the code will not trust anything else in the submission.

---

## BS-C-03 · Reported false-positive rates are in-sample and bounded by construction `[verified]`

**Component:** threshold calibration

The threshold is chosen as a quantile of a background sample, and the false-positive rate is then
reported on **that same sample**.

**Evidence.** `src/brainsafe/models/final_thresholds.py:58` then `:83`
```python
thr_bg = float(np.quantile(pbg, 1.0 - BACKGROUND_FPR))   # BACKGROUND_FPR = 0.05
...
bgfpr  = float((pbg >= thr).mean())                       # same pbg
```
Identically in `calibrate_binder_thresholds.py:68,76` and `screening_thresholds.py:57,79`.

Recomputed from the shipped `binder_modes.json`:
```
background_fpr_at_threshold : n=49  min=0.0    median=0.0273  max=0.0500
screening_background_fpr    : n=41  min=0.001  median=0.0100  max=0.0330
fpr_at_threshold            : n=45  min=0.015  median=0.1050  max=0.2500
```

`background_fpr_at_threshold` **never exceeds 0.0500**, which is exactly `BACKGROUND_FPR`, and
`screening_background_fpr` sits at exactly `SCREENING_FPR = 0.01` at the median. These values cannot
exceed their targets: the threshold is the corresponding quantile of the very sample being scored.
(They fall *below* the target when the measured-inactive threshold `thr_in` exceeds `thr_bg`, which is
why the median is 0.0273 rather than a flat 0.05 — so this is an in-sample bound rather than a strict
identity. The distinction does not rescue the number: it carries no information about out-of-sample
behaviour.)

**Severity:** Critical. These are the specificity figures that justify the panel.

**Fix.** Draw two disjoint background samples: one to set the quantile, one to report the rate.
`screening_thresholds.py:39` already uses a different seed (`default_rng(11)` vs `default_rng(7)`) —
make that split explicit and reuse it.

---

## BS-C-04 · The "independent" background sample is drawn from the decoy pool `[reported, code verified]`

**Component:** threshold calibration

`final_thresholds.py:8-9` states:

> *"This sample is independent of every target's training set, so it is the component that bounds false
> discovery."*

Both the decoys and the background sample come from the same file, with no mutual exclusion:
```python
# train_binders_hybrid.py:69-70, 111  -- decoys
with (M / "ad_reference.pkl").open("rb") as fh:
    bg_smiles, bg_fps = pickle.load(fh)
    ...
    dec.append(str(bg_ok[i]))

# final_thresholds.py:44-47           -- "independent" background
with (M / "ad_reference.pkl").open("rb") as fh:
    bg_smiles, _ = pickle.load(fh)
idx = rng.choice(len(bg_smiles), size=min(N_BACKGROUND, len(bg_smiles)), replace=False)
```

With up to 12,941 decoys per endpoint drawn from a ~74k library and `N_BACKGROUND = 3000`, hundreds of
background compounds are expected to be compounds the model was explicitly trained to score as zero.
Same defect in `calibrate_background_specificity.py:49-51` and `screening_thresholds.py:43-45`.

**Severity:** Critical, and it compounds BS-C-03.

**Fix.** Partition `ad_reference.pkl` once into a decoy pool and a reserved specificity pool, and
never let a compound cross.

---

## BS-C-05 · The headline specificity of 0.875 is measured on the applicability-domain reference set itself `[verified]`

**Component:** specificity evaluation

`noncns_specificity_fast.py:79-81` draws its 1,000 "non-CNS" negatives from `ad_reference.pkl` — the
very library `app.assess_domain()` measures similarity *against*. Every sampled compound is therefore
its own nearest neighbour.

**Evidence.** Recomputed from the shipped `results/tables/noncns_specificity_predictions.csv` (n=1000):
```
ad_max_tanimoto:  min 1.000   median 1.000   fraction == 1.0  =  1.000
```
**All 1,000 compounds have maximum Tanimoto exactly 1.000 to the reference set.**

The consequence is visible in the output. `results/tables/noncns_specificity_summary.csv` contains:
```
False-positive rate (any actionable disease call), 125, 1000, 0.125
"False-positive rate, in domain (T>=0.5)",         125, 1000, 0.125
```
The "in domain" stratum is **numerically identical** to the unstratified row, and the "near domain"
and "out of domain" rows are absent because those subsets are empty. Nothing in the script flags this.

This figure reaches the manuscript as *"875 of 1000 returned no actionable disease signal, a
specificity of 0.875 (95% CI 0.853 to 0.894)"*, and propagates into
`results/tables/MASTER_validation_summary.csv` as a balanced accuracy of 0.8328 (by averaging it with a
scaffold-held-out sensitivity — two different populations) and into the PPV/NPV decision analysis in
`integrity_audit.py:120-129`.

**Severity:** Critical. It is a statement about compounds the model has seen, presented as specificity
on novel chemistry.

**Fix.** Draw the negatives from a pool disjoint from `ad_reference.pkl`, or retire 0.875 and promote
the inversion suite's random-PubChem rate (H4: 0.1117 overall, n=600) with its stratification intact.
Do not average a specificity and a sensitivity measured on different populations.

---

## BS-C-06 · Feature-identical duplicates are split across random CV folds; BBB is 48% duplicates `[verified]`

**Component:** training pipeline / feature generation

No deduplication happens anywhere between loading an endpoint table and splitting it. The featuriser
leaves `includeChirality` at its default of `False` (`featurize.py:30`) and strips salts to the largest
fragment (`featurize.py:60-62`), so stereoisomers, salt forms and protonation variants produce
**byte-identical feature vectors**. `train_rf.py:101-107` then splits with
`StratifiedKFold(shuffle=True)`, which cannot keep them together.

`grep -n "drop_duplicates\|dedup" src/brainsafe/models/train_rf.py` → **no matches.**

**Measured across all 60 endpoint tables** ([`evidence/leak_internal.py`](evidence/leak_internal.py),
[`duplicate_audit.csv`](evidence/duplicate_audit.csv)):

- **198,499 rows; 13,846 (6.98%) are feature-identical duplicates of another row**
- **288 duplicate groups carry mutually contradictory labels** (identical input, both 0 and 1)
- **21 of 60 endpoints exceed 5% duplication**

| Endpoint | Rows | Unique feature vectors | Duplicate rows | % | Conflicting groups |
|---|---|---|---|---|---|
| **BBB** | 7,807 | 4,032 | **3,773** | **48.34** | **131** |
| NET | 3,067 | 2,660 | 407 | 13.27 | 3 |
| DAT | 2,692 | 2,343 | 349 | 12.96 | 4 |
| OPRK1 | 5,167 | 4,624 | 543 | 10.51 | 2 |
| SERT | 4,572 | 4,142 | 430 | 9.41 | 3 |
| BACE1 | 8,501 | 7,873 | 628 | 7.39 | 41 |

This is **invisible to any SMILES-string check**: all 7,807 BBB SMILES strings are distinct, with zero
string-level duplicates and zero string-level label conflicts.

**Quantified impact** ([`evidence/repro_cv.csv`](evidence/repro_cv.csv)) — republishing the random
10-fold CV with duplicates collapsed:

| Endpoint | Published | Reproduced | Deduplicated | Rows removed | **Inflation** |
|---|---|---|---|---|---|
| **BBB** | 0.9605 | 0.9605 | **0.8990** | 3,904 | **+0.0614** |
| AChE | 0.9626 | 0.9626 | 0.9619 | 171 | +0.0007 |
| MAO_A | 0.9467 | 0.9467 | 0.9462 | 94 | +0.0005 |
| hERG | 0.9541 | 0.9541 | 0.9516 | 169 | +0.0026 |

**Severity:** Critical, but *narrowly so*. For AChE, MAO_A and hERG the effect is negligible (≤0.003)
and the published numbers stand. The damage is concentrated in **BBB** — the gate that multiplies into
every per-disease score, and the endpoint carrying the external-validation claim. The README's
"Random 10-fold 0.95–0.97" is defensible for the target panel but **not for BBB**, whose honest
random-split value is **0.899**.

**Fix.** Deduplicate on the feature vector before splitting; resolve the 288 contradictory groups
explicitly; then either set `includeChirality=True` (defensible for BBB and receptor binding) or state
plainly that the model is stereo-blind. The current state is the worst of both: stereo-blind features
over stereo-distinct rows.

---

## BS-C-07 · The "306 FDA drugs absent from training" external set is 21% memorised `[verified]`

**Component:** external validation

The overlap exclusion uses the **full 27-character InChIKey**, which distinguishes stereoisomers,
salts and protonation states. The model's features do not (BS-C-06). A compound can therefore be
excluded from the "overlap" set and still be an input the model has memorised.

**Evidence.** `src/brainsafe/data/integrate_external.py:173` — `"in_b3db_training": ik in train_keys`,
with `train_keys` built from full InChIKeys at `:142-145`; consumed at
`src/brainsafe/evaluation/external_validation.py:42`.

Recomputed independently ([`evidence/leak_external.py`](evidence/leak_external.py)):
```
BBB training rows: 7807 · external rows: 1683 · flagged novel: 306
unique full InChIKeys in training: 7802
unique InChIKey skeletons in training: 4018      <- 3,784 collapse under stereo/salt/protonation

(a) full-InChIKey overlap    :  0 / 306      <- the repository's own criterion
(b) skeleton-InChIKey overlap: 78 / 306      (25.5%)
(c) Tanimoto == 1.0          : 65 / 306      (21.2%)   <- the decisive one
```

**65 of the 306 "never seen" compounds have a feature vector identical to a training compound.**
Named examples: Warfarin, Spironolactone, Prednisone, Budesonide, Ketoprofen, Minocycline, Miconazole,
Methylprednisolone, Quinidine, Naloxone, Scopolamine, Zopiclone, Estradiol, Cortisone, Levofloxacin,
Betamethasone, Cytarabine, Praziquantel, Pravastatin, Nicardipine (full list in
[`evidence/leak_external.txt`](evidence/leak_external.txt)).

**The direction of the bias is favourable, and this should be said plainly.** A component audit
recomputed the metric on the 241 genuinely novel compounds and obtained **AUROC 0.8015**, against
0.7741 as reported. Removing the contamination *improves* the headline. `[reported]`

**Severity:** Critical — because the claim as written is indefensible, and because it is framed in the
code itself as *"a genuine external test on approved drugs, the population a reviewer cares about
most"* (`external_validation.py:6-7`).

**Fix.** Recompute on the subset novel in **feature space**, report n=241 and both numbers, state the
criterion, and add the fingerprint-identity check to `integrate_external.py`. This finding costs
nothing to fix and makes the paper stronger.

---

## BS-C-08 · A quarter of the inversion suite's "held-out" compounds are training actives of another panel target `[verified]`

**Component:** inversion and falsification suite

`inversion/summarise.py:37-40` writes, and `inversion/REPORT.md:20` publishes, that H1 was *"scored
with hold-out models only, so no compound was seen in training."*

Scaffolds are withheld **per target** (`scaffold_holdout_panel.py:107`), but `inv_disease_layer.py:82-86`
scores every held-out compound with **every** model in the panel. A compound active at two targets is
held out from one model and is a training active of the other.

**Measured independently** over `models_rf/holdout/heldout_actives.json` (46 targets, 17,607 entries)
against the pChEMBL ≥ 7 actives of every other endpoint:

```
held-out entries that are TRAINING ACTIVES of another panel target: 4562 / 17607 = 25.9%

   DAT        heldout=  207   also trained elsewhere=  142  (68.6%)
   D2         heldout=  796   also trained elsewhere=  502  (63.1%)
   a3b4nAChR  heldout=   42   also trained elsewhere=   26  (61.9%)
   SERT       heldout=  735   also trained elsewhere=  396  (53.9%)
   NET        heldout=  265   also trained elsewhere=  139  (52.5%)
   OPRK1      heldout=  562   also trained elsewhere=  256  (45.6%)
   D3         heldout=  790   also trained elsewhere=  351  (44.4%)
```

*Discrepancy recorded rather than resolved:* a component audit reported **31.4%** using a slightly
different active-set definition (and 100% overlap for OX1/OX2). My figure of **25.9%** uses
`pchembl >= 7` over `data/endpoints/*.csv` for the 46 targets in the hold-out file. Both support the
finding; the exact figure should be pinned down before it is quoted.

The affected pairs are exactly the homologous families the suite's own H8 identifies as correlated:
D2/D3 → Psychosis, SERT/NET/DAT → Depression and ADHD, OPRM1/OPRK1 → Chronic pain. The permutation
null does not control for this, because shuffling the target-to-disease map leaves the memorisation
intact and merely sends it to the wrong disease.

**Severity:** Critical. H1's headline top-3 accuracy (0.8045) is materially memorisation, and this is
the claim a reviewer is most likely to test.

**Fix.** Rerun H1 restricted to compounds appearing in no other target's training actives, or withhold
scaffolds **globally** across the panel. Correct the sentence in `summarise.py:37` and regenerate
`REPORT.md`.

---

## BS-C-09 · The inversion duplicate check tests the wrong artefact and cannot fail `[verified]`

**Component:** inversion and falsification suite

The repository advertises *"Validated by inversion (six adversarial checks, all pass)"*
(`README.md:34-35`). One of the six is a duplicate check, and it passes while BS-C-06 is true.

**Evidence.** `results/tables/inversion_validation.csv`, row 2:
```
check : No duplicate compounds in master library
result: PASS
detail: 61,317 rows, 61,317 unique InChIKeys (0 duplicates)
```
against my measurement of the tables actually trained on: **13,846 duplicate rows across 60 endpoint
tables, 3,773 in BBB alone.**

Both statements are true; they are about different files. `validate_inversion.py:56-60` reads
`data/processed/compound_library.csv` — the deduplicated master catalogue — while the models train
from `data/endpoints/*.csv`. It also keys on full InChIKey, which does not collapse stereoisomers.

Two of the other six checks are similarly unfalsifiable `[reported, code verified]`:
- **Check 1 (`validate_inversion.py:43-53`)** asserts that `sklearn.GroupKFold` produces disjoint
  groups — a contract scikit-learn guarantees. It cannot fail for any input, and it does not test that
  `_scaffold_groups` computes the *right* groups, which is the actual defect (BS-N-02, BS-N-03).
- **Check 3 (`:63-77`)** re-executes the identical code path with the identical seed and compares
  against the artefact that path produced. It is a staleness detector, not a reproducibility test.

Checks 4–6 are underpowered by construction: n=200 (an arbitrary `.head(200)` slice), **n=4**
(two CNS versus two peripheral compounds, one of which — donepezil — is BBB training data), and
**n=1** (PFOA), each against a hard-coded pass bar (`0.1`, `pc > pp`, `0.30`) that sits comfortably on
the correct side of the observed value with no recorded pre-registration.

**Severity:** Critical. A falsification suite that cannot fail is worse than none, because it is cited
as positive evidence in the README and the manuscript.

**Fix.** Point the duplicate check at `data/endpoints/*.csv` and test on the **feature vector**;
verify it fails on today's data before trusting it again. Replace check 1 with a chemistry-level
assertion (no test compound within Tanimoto 0.95 of any training compound). Convert checks 5 and 6 to
panels of ≥50 compounds with intervals.

---

## BS-C-10 · `integrity_audit.py` computes calibration on training compounds and labels it "held-out" `[verified]`

**Component:** integrity audit

The section header at `src/brainsafe/evaluation/integrity_audit.py:83` reads
`"B. CALIBRATION against held-out measured inactives"`, and the module docstring (lines 11-14) repeats
the claim. The code holds nothing out:

```python
# integrity_audit.py:92-103
df  = pd.read_csv(f).dropna(subset=["smiles"])
pv  = pd.to_numeric(df.get("pchembl"), errors="coerce")
act = df.loc[pv >= 7, "smiles"].astype(str).tolist()      # the exact positive class trained on
ina = df.loc[df["label"] == 0, "smiles"].astype(str).tolist()   # ALL inactives, incl. ina_train
...
p = np.r_[mdl.predict_proba(Xa)[:, 1], mdl.predict_proba(Xi)[:, 1]]
e, b = ece(y, p), brier_score_loss(y, p)
```

Compare `scaffold_holdout_panel.py:116,138`, which splits inactives in half and trains on
`ina_train`. Neither the actives nor the trained inactive half is excluded here. The reported ECE and
Brier in `results/tables/integrity_calibration_per_target.csv` and `integrity_audit.csv` are
training-set calibration.

**Related, same file:** `integrity_audit.py:139` — `for ep, smis in list(held.items())[:8]:` — the
leakage check silently examines **only the first 8 of 46 targets**, then reports the result as
`"targets checked for scaffold overlap"`. The truncation is undocumented and `json` preserves
insertion order, so the same 8 are always chosen. 38 targets are never checked.

**Severity:** Critical (the mislabelling); Major (the `[:8]`).

**Fix.** Restrict to the held-out inactive half and score actives with the hold-out twins in
`models_rf/holdout/`. Until then, relabel the section "calibration on training compounds
(optimistic)". Remove the `[:8]` or make it a named constant with the fraction reported.

---

## BS-C-11 · The manuscript contains zero in-text citations and an unfilled reference placeholder `[verified]`

**Component:** manuscript / NAR compliance

**Evidence.** A regex sweep of the manuscript body (lines 1–1000) for `[n]`, `(YYYY)` and `et al.`
returns **0 matches**. Twenty numbered references sit at lines 1009–1028 with nothing pointing to
them, so every one is listed but never cited.

The reference section opens with a placeholder (`manuscript/NAR_WebServer_BrainSafe_built.md:1003-1005`):
> `[To be completed. Anchor citations: ChEMBL; BindingDB; B3DB; RDKit; Therapeutics Data Commons;`
> `KEGG; Reactome; IUPHAR/BPS Guide to Pharmacology; scikit-learn; Streamlit; Bemis and Murcko scaffolds;`
> `isotonic calibration; DeLong test.]`

and closes at line 1035 with a self-admitted list of uncited essentials:
> *"Requested but not resolved above the similarity threshold, and therefore not cited: b3db, sklearn,
> tdc, kpuu, xgboost, gin_gnn, platt, conformal, herg_pred, bbb_ml, cns_attrition, lrrk2_pd,
> nlrp3_neuro, hdac_hd, riluzole_als."*

**B3DB — the BBB training source and the basis of the headline external claim — has no reference entry
anywhere.** `README.md:112` gives only a bare prose attribution.

Missing method citations `[reported]`: ECFP-4/Morgan (Rogers & Hahn 2010 — ref 20 substitutes an
authorless IUPAC Gold Book entry), scikit-learn (Pedregosa), random forest (Breiman 2001 — ref 16
substitutes a software textbook), isotonic calibration (Zadrozny & Elkan — ref 5 substitutes an
adjacent empirical study), conformal prediction (Vovk — and the method is absent from the manuscript
entirely despite being claimed in `README.md:14`), PubChem, XGBoost, GIN, the DPPH assay, AqSolDB,
DrugBank, and the FDA-curated drug set behind the external validation. Correctly cited: ChEMBL (17),
BindingDB (13), Bemis–Murcko (3), Wilson (1), DeLong (2), MoleculeNet (15), DUD-E (8), KEGG (4),
Reactome (18), IUPHAR (19), Jaworska (6).

**Severity:** Critical, submission-blocking.

**Fix.** Add in-text citations throughout, complete the reference list, and add primary-source
citations for every method named above.

---

## BS-C-12 · Running the documented rebuild command on a fresh clone silently destroys the training data `[verified by inspection]`

> **STATUS: FIXED.** See [`FIXES.md`](FIXES.md). Validated in an isolated tree: the fresh-clone case now
> exits 1 with 0 tables emptied, and the happy path reproduces the pre-fix output byte-identically.

**Component:** data assembly

`src/brainsafe/data/rebuild_endpoints.py:17` instructs a reviewer to
`Run: python src/brainsafe/data/rebuild_endpoints.py`. On a fresh clone the API caches are absent
(they are gitignored — BS-C-13), and the loader returns an **empty frame instead of raising**:

```python
# rebuild_endpoints.py:46-48
cache = CHEMBL_CACHE / f"{name}_y.json"
if not cache.exists():
    return pd.DataFrame(columns=["inchikey", "smiles", "pchembl", "year"])
```

`main()` then writes the result **unconditionally**:
```python
# rebuild_endpoints.py:115-116
df, p = rebuild_target(name)
df.to_csv(ENDPOINTS / f"{name}.csv", index=False)
```

**All 11 core endpoint CSVs are overwritten with header-only files, with no error and a zero exit
code.** The backup at `:107-109` does not help: it is guarded by `if not backup.exists()`, and
`archive/legacy/` does not exist in this working tree, so on a second run it would preserve nothing.

*I did not execute this command, for obvious reasons.* The finding rests on reading the code path; a
component audit reproduced it in an isolated copy and recorded "rows written on a fresh clone: 0".

**Severity:** Critical. The first thing a reproducibility-minded reviewer does is run the documented
command, and it deletes the data.

**Fix.** Raise on a missing cache; refuse to write an empty table; make the backup unconditional and
timestamped.

---

## BS-C-13 · The core endpoint tables cannot be regenerated: the fetch scripts are deleted and the caches are gitignored `[verified]`

> **STATUS: FIXED** apart from the Zenodo deposit itself. Scripts restored, caches committed and
> packaged with checksums. See [`FIXES.md`](FIXES.md).

**Component:** provenance / reproducibility

`rebuild_endpoints.py:46` reads `data/_chembl_cache/{name}_y.json`. **Nothing in the repository writes
those files.** They were produced by `BS_fetch_endpoints.py`, deleted in commit `fea5029`:

```
$ git log --diff-filter=D --name-only fea5029 | grep fetch
archive/legacy/analysis_scripts/BS_fetch_antioxidant.py
archive/legacy/analysis_scripts/BS_fetch_clinical.py
archive/legacy/analysis_scripts/BS_fetch_endpoints.py
$ ls src/brainsafe/data/ | grep BS_fetch      → (nothing)
```

`data/raw/measured_endpoints_SOURCE.md` still names `BS_fetch_endpoints.py` as the retriever for both
the ChEMBL targets and B3DB, and `BS_fetch_antioxidant.py` for the DPPH set. **Three primary data
sources cite scripts a reviewer cannot obtain from the repository.**

Meanwhile `.gitignore` excludes every raw cache (`data/_bindingdb_cache/`, `data/_chembl_cache/`,
`data/_pubchem_cache/`), so neither the retriever nor its output is distributed. The endpoint CSVs
themselves *are* committed (BS-P-09), which is what saves the situation — but the chain from public
database to training table is unreproducible.

Two further defects visible in the recovered script `[reported]`: `MAX_PAGES, PAGE = 16, 1000` caps
every core target at 16,000 activities with no truncation warning, and it grouped by **raw SMILES**.

**Severity:** Critical for a Web Server submission judged on reproducibility.

**Fix.** Restore all three scripts from `fea5029^` into `src/brainsafe/data/`, and either commit the
caches (~31 MB) or deposit them alongside the Zenodo model archive with a checksum manifest.

---

## BS-C-14 · BindingDB censored measurements are stripped of their relation and used as exact values `[verified]`

> **STATUS: FIXED in code; shipped tables not regenerated.** Quantified in [`FIXES.md`](FIXES.md), and
> the practical impact is far smaller than this entry implies: **no label anywhere flips**, six of
> 21,791 pooled potencies shift, and the censored records are almost all `<` bounds (too potent to
> measure), so **no negative class is recovered here**. Severity should be read as Major, not Critical.

**Component:** data acquisition / label definition

`src/brainsafe/data/fetch_bindingdb.py:85-94`
```python
def parse_affinity(raw: str) -> float | None:
    """Convert a BindingDB affinity string in nM to a -log10(molar) potency value."""
    s = str(raw).strip().lstrip("><~=").strip()
    ...
    return 9.0 - math.log10(v)  # nM -> pX
```

`lstrip("><~=")` discards the censoring operator. `">10000"` becomes exactly 10,000 nM (pX = 5.0);
`"<1"` becomes exactly 1 nM (pX = 9.0, an active). Neither is what the record says.

Measured from the committed cache:
```
leading relation chars across data/_bindingdb_cache/*.json: {'<': 849, '>': 299}
total censored: 1,148
```
A component audit puts this at **1,148 of 24,014 records (4.8%)**, worst for GSK3B (395 of 2,513,
15.7%) and A2A (368 of 3,751, 9.8%). `build_bindingdb_compounds` (`:115`) additionally takes a
**median across a mixture of Ki, Kd, IC50 and EC50** with no distinction.

**The ChEMBL branch is safe, but by accident** `[reported]`. `grep -rn "standard_relation" src/`
returns **zero hits**; the fetchers rely on `pchembl_value__isnull=false`, and ChEMBL assigns a
pChEMBL value only where the relation is `=`. A live API check confirmed zero non-`=` relations among
records carrying a pChEMBL value. This should be made explicit rather than left to convention.

**Severity:** Critical — a well-known QSAR defect that a cheminformatics reviewer will look for.

**Fix.** Retain the relation. Treat `>` as right-censored (inactive if the bound is below the inactive
cut, otherwise drop), `<` as left-censored, exclude `~`. Aggregate within one affinity type, or record
the type as a column.

---

## BS-C-15 · The query design discards the measured negative class, leaving most endpoints near-degenerate `[verified]`

**Component:** data acquisition / label definition

Requiring `pchembl_value__isnull=false` excludes every ChEMBL record whose relation is `>` — and those
records *are* the measured non-binders ("no inhibition up to 10 µM"). Live API counts `[reported]`:

```
AChE   CHEMBL220 : all=20132  with_pchembl=8259   relation_gt=2210   <- 2,210 measured inactives discarded
OX1    CHEMBL5113: all=24348  with_pchembl=10429  relation_gt=662
PDE10A CHEMBL4409: all=14236  with_pchembl=6883   relation_gt=595
```

The consequence, recomputed across all 60 endpoint tables:

```
endpoints audited: 60
endpoints >90% active: 39
endpoints >96% active: 21
endpoints with fewer than 25 inactives: 5
```

| Endpoint | Rows | Active | Inactive | % active |
|---|---|---|---|---|
| OX2 | 6,231 | 6,207 | **24** | 99.61 |
| **GABAA_a5** | 676 | 672 | **4** | **99.41** |
| P2X7 | 4,401 | 4,358 | 43 | 99.02 |
| PDE10A | 5,112 | 5,052 | 60 | 98.83 |
| LRRK2 | 1,478 | 1,460 | 18 | 98.78 |
| HT2A | 5,989 | 5,908 | 81 | 98.65 |

**No target-specific decision threshold can be estimated from 4 negatives.** The project's own
`fetch_pubchem_expansion.py` docstring concedes this for LRRK2, OX2 and GABA-A. The attempted remedy
reached almost nothing: `results/tables/pubchem_merge_provenance.csv` records **one** endpoint (OX1,
33 compounds added), and 18 of the 22 `data/_pubchem_cache/*_inactives.csv` files are 2 bytes.

This interacts directly with BS-C-01 and BS-C-18: thresholds for these endpoints are quantiles of a
handful of points.

**Severity:** Critical.

**Fix.** Issue a second query per target for `standard_relation=>` with `standard_units=nM`, convert
the bound, and admit those compounds as label 0 — the single highest-value data change available.
Report the base rate beside every AUROC, and state plainly which endpoints cannot support a calibrated
decision.

---

## BS-C-16 · Two incompatible deduplication keys; 48 of 59 endpoints deduplicate on the raw SMILES string `[verified]`

**Component:** data assembly — **this is the root cause of BS-C-06**

Only the 11 core targets and the BindingDB pool deduplicate by InChIKey
(`rebuild_endpoints.py:57`, `fetch_bindingdb.py:115`). Every expansion fetcher groups on the **raw
ChEMBL SMILES string**:

```
fetch_batch2.py:92            med = df.groupby("smiles").agg(...)
fetch_batch3.py:69            med = df.groupby("smiles").agg(...)
fetch_batch4.py:120           med = df.groupby("smiles").agg(...)
fetch_batch5.py:69            med = df.groupby("smiles").agg(...)
fetch_new_targets.py:74       med = df.groupby("smiles").agg(...)
fetch_readacross_targets.py:97
fetch_pka.py:98
```

`build_compound_library.standardise` — which strips salts and computes the InChIKey — is imported by
`rebuild_endpoints.py`, `fetch_bindingdb.py`, `fetch_pubchem_inactives.py` and the ADME scripts, and by
**none** of the seven above. No standardisation is applied before these tables are written.

Consequence, measured: **4,012 rows across 43 endpoint files carry multi-fragment (salt or mixture)
SMILES**, plus 3,318 of 44,127 in `data/readacross/`. These are distinct strings that collapse to
identical feature vectors — exactly the duplication quantified in BS-C-06.

`docs/DATA_DICTIONARY.md` asserts *"Deduplication is by full standard InChIKey"* as a property of the
pipeline. That is true only of `compound_library.csv`, not of the endpoint tables the models train on.

**Severity:** Critical.

**Fix.** Standardise to InChIKey (with salt stripping) before writing any endpoint table. Add
neutralisation and tautomer canonicalisation, which are absent repository-wide
(`grep -rn "rdMolStandardize\|Uncharger\|TautomerEnumerator" src/` → nothing).

---

# MAJOR FINDINGS

### Numbers, documentation and provenance

## BS-M-01 · The record and compound counts disagree three ways, and the live server states a superseded figure `[verified]`

"64,474 measured records; 61,317 unique compounds" exists **only as hard-coded prose** — no script
computes it and no results table contains it:
`app.py:8`, `:2267`, `:2666`, `:2687`, `:2825`, `README.md:13`, `docs/BS_BENCHMARK_ANALYSIS.md:12`,
`docs/BS_MODEL_CARD.md:12`, `:215`, `deploy/huggingface/README.md:20`.

| Source | Records | Unique compounds |
|---|---|---|
| Manuscript (`:104`) | **67,984** | **61,226** |
| README (`:13`, `:24`) and the live app | **64,474** | **61,317** |
| `docs/RF_CV_RESULTS.md:16`, `docs/BS_MODEL_CARD.md:26` | 67,982 | 61,317 |
| Exact sum of `rf_cv_summary.csv` / `inversion_validation.csv` | **67,982** | **61,317** |

`docs/BS_MODEL_CARD.md:215` labels its own decomposition authoritative — *"where they differ, these
deployed numbers are correct"* — and it is stale. Recomputed against `rf_cv_summary.csv`:

| Endpoint | Model card *n* | Live *n* | Δ | Card AUROC | Live scaffold AUROC | Δ |
|---|---|---|---|---|---|---|
| AChE | 4,324 | 4,387 | +63 | 0.915 | 0.9212 | +0.0062 |
| BACE1 | 8,067 | 8,501 | **+434** | 0.950 | 0.9556 | +0.0056 |
| GSK3B | 4,044 | 4,958 | **+914** | 0.920 | 0.9369 | **+0.0169** |
| MAO_B | 3,455 | 3,665 | +210 | 0.885 | 0.8895 | +0.0045 |
| hERG | 5,905 | 5,875 | −30 | 0.901 | 0.9208 | **+0.0198** |

**Six of eight AUROC values disagree with the live results table.** The model card carries both 64,474
and 67,982, ninety lines apart.

**Fix.** Compute both counts in a script, emit to `results/tables/`, and have the manuscript build and
`app.py` read them. The repository already does this for `REPOSITORY_MAP.md`; apply the same pattern.

## BS-M-02 · `MASTER_validation_summary.csv` describes a different experiment from the manuscript `[reported]`

| Quantity | Manuscript / `scaffold_holdout_results.csv` | `MASTER_validation_summary.csv` |
|---|---|---|
| Pooled recall | 12,325/15,609 = 0.790 | **11,914/15,069 = 0.7906** |
| Targets pooled | 36 | **40** |
| Threshold-collapsed exclusions | 3: OX2, LRRK2, NLRP3 | **4: CGRP, LRRK2, OX2, RIPK1** |
| Targets ≥ 0.80 recall | 19 of 36 | **22 of 40** |
| Targets < 0.50 recall | five | **six** (adds Nav1_6) |

The file whose name signals it is the authoritative roll-up agrees with neither the manuscript nor the
generating table, and it names **CGRP and RIPK1 — the two flagship new endpoints — as failures**, while
the manuscript reports CGRP at AUROC 1.000 / sensitivity 0.993 and RIPK1 at 0.995 / 0.950. **Reconcile
or withdraw the file before submission.**

## BS-M-03 · Endpoint and model counts contradict each other, inside the manuscript and across the docs `[verified + reported]`

| Claim | Source | Live value |
|---|---|---|
| "Target panel (13 endpoints)" | `README.md:24` | **52** deployed endpoints |
| "all 63 model artefacts load" | `README.md:101` | `app_health.py` prints **69 model objects loaded** |
| "62 endpoints … 69 trained models" | manuscript `:93-96` | — |
| "67 endpoint names / 71 models" | manuscript `:153-155` (the uncommitted edit) | — |
| **"all 40 models"** | manuscript `:298` | — |

Captured (`evidence/health_full.txt`): `69 model objects loaded`; `52 targets, 16 conditions`;
`52 deployed endpoints`. The 63 is a `.joblib` file count in one directory that *includes* two
withdrawn models and *excludes* the entire ADME layer. `deploy/huggingface/README.md:17` says 52,
correctly. The uncommitted manuscript edit introduces a fourth count in a section three paragraphs
from the third.

**Fix.** Generate these counts. Define "endpoint" once. `README.md:9-11,24` currently describes the
retired 13-endpoint panel, so a reviewer opening the repository will not recognise the tool in the
abstract.

## BS-M-04 · Further internal contradictions in the manuscript `[reported]`

- **Weight ablation, two triples for one experiment:** 0.7917/0.7911/0.7899 (`:290`, and `app.py:121-122`)
  versus 0.8045/0.8036/0.8025 (`:41`, `:483`, `:593`). Line 595 explains the rerun; Methods was never updated.
- **Condition count:** sixteen (`:30`, `:730`, correct — `len(DISEASE_ORDER)` = 16) versus fourteen
  (`:283`) versus eleven (`:934`, the Figure 1 legend).
- **Binder sensitivity:** 0.897 (abstract `:34`) versus 0.88 (`:256`) versus Table 2's own mean 0.914;
  and *"no target now falls below the reliability threshold"* (`:260`) against *"8 targets under 0.50"*
  (`:488`).
- **Binder AUROC panel size:** 43 targets (`:33`), 34 (`:245`), 42 (MASTER csv).
- **Hold-out denominator:** 16,874 compounds stated at `:534`, pooled recall computed over 15,609 at
  `:537`, with the 1,265-compound gap never explained.

## BS-M-05 · Calibration degrades BBB, and only the mean is reported `[verified]`

`results/tables/calibration.csv`, recomputed:
```
mean ECE raw 0.0721 → calibrated 0.0127        (README claims 0.072 → 0.012)  MATCHES
BBB:  ECE   0.0403 → 0.0424   (worse)
BBB:  Brier 0.0788 → 0.0947   (worse)
```
BBB is the **only** endpoint that regresses, it started as the best-calibrated, and
`BBB_calibrated.joblib` exists and is loaded preferentially by `app.py:311`. Nothing gates deployment
on improvement. `docs/METHODOLOGY_AUDIT.md:30-31` states this honestly; the manuscript (`:264`) and
README (`:33`) do not.

**The calibration method itself is sound** — see BS-P-03. The likely root cause is BS-C-06: with 48%
duplicate rows, isotonic regression fits a heavily tied distribution.

## BS-M-06 · "Every classifier is isotonically calibrated" is false `[reported, code verified]`

Manuscript `:264`. The **47 binder classifiers use sigmoid (Platt)**, not isotonic:
`train_binders_hybrid.py:128`, `train_receptor_binders.py:120`, `train_new_binders.py:105`,
`train_batch2.py:120`, `train_measured_label_holdout.py:62`, `_train_small.py:62`,
`_train_nav17.py:117`, `_train_ox.py:117` — all
`CalibratedClassifierCV(FrozenEstimator(forest), method="sigmoid")`. Only the eight in
`train_rf.CLASSIFICATION` are isotonic. **Fix:** reword to "the eight property and target classifiers".

## BS-M-07 · The deployed calibrated model is not the model whose ECE is reported `[reported, code verified]`

`calibrate.py:59` measures isotonic-on-out-of-fold-predictions; `calibrate.py:74-77` builds a
*different* estimator (`CalibratedClassifierCV(base, cv=5).fit(X, yv)`) and saves it as
`{ep}_calibrated.joblib`, which `app.py:311-312` loads. The shipped artefact is never evaluated.
Related: `calibrate.py:59` passes `cv=5` to a regressor, giving unshuffled, unstratified `KFold` over
label-ordered OOF files (measured fold-0 test positive fraction 0.84 against 0.322 in fold 3), and ECE
is computed only on the random-split OOF while scaffold OOF files exist unused.

## BS-M-08 · BBB has no temporal validation, and the temporal collapse of the regressions is not surfaced `[verified]`

`data/endpoints/BBB.csv` has columns `smiles,label` only — no `year` — so
`rf_conformal_temporal.py:71-72` silently skips it. Recomputed: `'BBB' in rf_temporal.csv` → **False**.
The README's *"Temporal (future compounds) 0.61–0.91"* row does not cover the BBB gate, and nothing
says so.

The regressions collapse, and the README's headline table shows only classifier AUROC:

| Endpoint | Random-CV R² | Temporal R² | Loss |
|---|---|---|---|
| D2 | 0.6012 | **0.042** | −93% |
| A2A | 0.6824 | 0.338 | −50% |
| HT2A | 0.6362 | 0.182 | −71% |
| SERT | 0.6016 | 0.100 | −83% |
| antioxidant_DPPH | 0.6689 | **0.009** | −99% |

Near-zero temporal R² on four of five regressors is the most important honest limitation in this work,
and it is currently discoverable only by opening a CSV. Declaring it will strengthen the paper.

Related `[reported]`: `fetch_new_targets.py:75` assigns each compound `year=("year","max")`, so a
compound first reported in 2008 and re-measured in 2023 is dated 2023 and lands in the "future" test
set. `min` is the conservative choice; the current setting biases temporal scores upward.

## BS-M-09 · Conformal coverage describes a refit model on a random split, and no conformal set reaches the server `[reported, spot-verified]`

`rf_conformal_temporal.py:45` fits a **fresh** forest on 60% of the data; the deployed
`models_rf/BBB.joblib` is a bare `RandomForestClassifier` fitted on 100%. A conformal guarantee is a
property of a specific predictor and does not transfer. The split is `rng.permutation` (random), not
scaffold — a weaker regime than the manuscript's own headline. Average set size is 1.001–1.108, so at
ε=0.10 "coverage 0.904" is arithmetically the model's accuracy and adds no uncertainty quantification.
`app.py:1131` shows a `"Calibrated + conformal"` badge; no per-compound prediction set is emitted.

**The underlying implementation is correct** — see BS-P-04. The problem is the claim, not the code.

## BS-M-10 · The applicability-domain result contradicts its own stated expectation, silently `[reported]`

`applicability_domain.py:5-7` says the model *"should be clearly better where it is in domain."*
`results/tables/applicability_bbb_validation.csv`:
```
in_domain (T>=0.30),     258,  AUROC 0.7715
out_of_domain (T<0.30),   48,  AUROC 0.8157      <- higher
all,                     306,  AUROC 0.7741
```
Out-of-domain **exceeds** in-domain. The script prints the table and returns — no test, no interval, no
comment. Since the domain flag qualifies every served prediction, this needs a DeLong or bootstrap
comparison and a sentence in the limitations, not silence in a CSV.

## BS-M-11 · The prospective external test is null and is still in the pipeline `[reported]`

`results/tables/external_100_summary.csv`: **top-1 accuracy 1/8 = 0.125, permutation p = 0.368,
Cohen's κ = 0.082**, with per-class recalls at n=1. `scaffold_holdout_panel.py:4-5` acknowledges the
test was underpowered and supersedes it, but `external_100.py`, `external_100_stats.py` and both output
tables remain live and are consumed by `build_review_folder.py`. **Fix:** archive them, or report the
null result explicitly.

### Statistics

## BS-M-12 · Confidence intervals treat clustered observations as independent `[reported]`

Wilson intervals are implemented correctly in seven places and `n` is always reported — genuinely good
practice. But every sample is clustered and none is treated as such. H1 reports
`0.8045 [0.7958, 0.8130]` on n=8,232 rows that are grouped by target (≤250 each), drawn from scaffold
clusters, and include 516 rows from compounds counted more than once. A ±0.9% interval over 46 clusters
is far too narrow. The same objection applies to the pooled scaffold-holdout recall
`0.7906 [0.7841, 0.7971]`. **Fix:** cluster-bootstrap resampling targets or scaffolds; expect the
intervals to widen several-fold.

## BS-M-13 · No multiple-comparison control anywhere `[reported]`

`grep` finds no Bonferroni, Holm or FDR in the repository. `scaffold_holdout_results.csv` reports 39
per-target recalls with 39 intervals; `H6_clinical_per_indication.csv` reports 9; `external_100_stats.py`
emits 12 metrics from n=15; `deployed_specificity_audit.csv` applies a pass/fail rule across 47
endpoints, where roughly two failures are expected by chance at a nominal 5% criterion.

## BS-M-14 · Pass/fail thresholds are hard-coded with no recorded pre-registration `[reported]`

Twelve verdict thresholds are bare literals — `validate_inversion.py:76,85,107` (`0.005`, `0.1`,
`0.30`); `inv_distant_specificity.py:140,145` (`0.125`, `× 1.5`); `summarise.py:33,47,75,91,119`;
`scaffold_holdout_report.py:26,58`; `deployed_specificity_audit.py:41`; `integrity_audit.py:139`.
Several duplicate each other across files and **two pairs actively disagree**:
`MIN_INACTIVES = 25` (`calibrate_binder_thresholds.py:34`) against `len(ina) >= 15`
(`final_thresholds.py:70`), and `MIN_SENS = 0.50` against a `sensitivity >= 0.60` reliability gate.

The sharpest case is H4: its SUPPORTED verdict rests on **n=59, k=3**, with a 95% CI of
[0.0174, 0.1392] that **contains** the comparator 0.125, against a hard-coded `ref_fpr * 1.5`
tolerance — and `ref_fpr = 0.125` is itself the contaminated number from BS-C-05.

**Fix.** One `evaluation/thresholds.py` with a rationale per value and a note on when each was fixed.

### Reproducibility and code health

## BS-M-15 · Six deployed endpoints have no recorded training command `[reported]`

`Nav1_6`, `Nav1_8`, `RIPK1`, `CGRP`, `DHODH`, `Cav3_2` carry `"mode":
"hybrid_decoys_plus_measured_inactives"` in `binder_modes.json` — so `train_binders_hybrid.py` produced
them via `argv` — but appear in **no** `TARGETS` list, and there is no Makefile, shell script or
documented invocation anywhere. **A reviewer cannot regenerate the panel.**

Compounding this, the three underscore-prefixed files are byte-identical clones of their parents apart
from the target list (`_train_nav17.py` and `_train_ox.py` differ from `train_binders_hybrid.py` only
at line 49; `_train_small.py` from `train_measured_label_holdout.py` only at line 31), and the parent
already supports `sys.argv` overrides. They are strictly redundant — but **`_train_nav17.py` is the
sole producer of `Nav1_7_binder.joblib` and `_train_small.py` the sole producer of TAAR1, GluA2 and
Nav1_1**. They are mis-named scratch files doing production work; deleting them would break
reproducibility.

## BS-M-16 · Random-forest hyper-parameters are copy-pasted into four incompatible settings `[reported]`

`train_rf.py:48` is the reference: `n_estimators=300, min_samples_leaf=2`. Across the panel:

| Setting | Files |
|---|---|
| 300 / 2 | `train_rf.py`, and correctly **imported** by `adme/train_adme.py:28` and `calibrate.py:33` |
| 300 / **4** | `train_binders_hybrid.py:47`, `_train_nav17.py:47`, `_train_ox.py:47`, `train_measured_label_holdout.py:33`, `_train_small.py:33` |
| **200 / 6** | `train_batch2.py:43`, `train_receptor_binders.py:112` *and* `:118` (same literal twice), `train_new_binders.py:98` *and* `:103` |
| **500** / 2 | `train_pka.py:42` |

`random_state=42` is written as a literal **14 times**. `train_rf.py:185` records hyper-parameters in
`<ep>_meta.json`, but the binder metas record none, so a reviewer cannot recover which forest produced
which endpoint. **Fix:** one imported constant; record it in every meta JSON.

Related `[reported]`: five binder scripts deploy a forest fit on **80%** of the data
(`train_test_split(test_size=0.2)` then `FrozenEstimator`), whereas `train_rf.py:131-135` and
`train_adme.py:99-100` deploy on 100%. Undocumented and asymmetric.

## BS-M-17 · `binder_modes.json` is accreted in place by six scripts with contradictory semantics and no run stamp `[reported]`

`calibrate_background_specificity.py:86` **ratchets** (`max(thr_old, thr_bg)`);
`final_thresholds.py:79` **recomputes from scratch** — its docstring says this exists precisely so
calibration "cannot ratchet thresholds upward"; `apply_specificity_decisions.py:87` **hard-assigns**.
There is no `produced_by`, no timestamp, and no run log, so the execution order that produced the
shipped panel is unrecoverable. `final_thresholds.py` and `screening_thresholds.py` also lack the
`deployed` guard that `calibrate_background_specificity.py:62-66` has, so re-running either re-includes
the withdrawn `Nav1_1` and `Cav3_2` in the reported medians.

## BS-M-18 · Thresholds set from as few as 15 compounds `[reported]`

`final_thresholds.py:70` admits `len(ina) >= 15`, against the project's own `MIN_INACTIVES = 25`. The
90th percentile of 15 points is the second-largest value. Shipped: MT1 n=15 → 0.9896; GluN2B n=16 →
0.9990; CSF1R n=18 → 0.3784; Nav1_8 n=20 → 0.9350; P2X7 n=22 → 0.9990.

## BS-M-19 · Thresholds hand-tuned so that named compounds do not fire `[reported]`

`apply_specificity_decisions.py:53-63`:
```python
RETHRESHOLD = {
    "a3b4nAChR": (0.450, "ethanol scored 0.437 at the calibrated threshold of 0.392"),
    "SIRT1":     (0.650, "random-chemistry false-positive rate 0.058, above the 5 percent intended"),
    "RIPK1":     (0.500, "fired on acetate, glycine and lactate at the floor threshold of 0.050"),
}
```
`0.450` is `0.437 + ε`. The sensitivities reported afterwards are measured after tuning on those same
compounds. Documented transparently — but it is threshold selection on the evaluation set.
**Fix:** pre-register the trivial-control set, tune on a disjoint half, report on the other.

## BS-M-20 · `noncns_specificity.py` is a stale duplicate of `noncns_specificity_fast.py`, and both write the same files `[reported]`

Both define the same `wilson`, `canon`, `N_SAMPLE = 1000` and `default_rng(31415)`, and both write
`results/tables/noncns_specificity_predictions.csv` and `..._summary.csv` — whichever ran last wins,
and no artefact records which. They sample from **different populations** (`data/adme/*.csv` versus
`ad_reference.pkl`). Determined empirically from the shipped predictions file: 990/1000 sampled SMILES
are in `ad_reference`, so the live 0.875 came from `_fast.py`. **Fix:** archive the duplicate; stamp
the producing module and pool into the summary.

## BS-M-21 · Feature ablation and descriptor importance use random splits while the headline is the scaffold split `[reported]`

`feature_analysis.py:65-66,96` uses `StratifiedKFold`/`KFold` and `train_test_split`. `_scaffold_groups`
is **imported at line 38 and never used** — dead code signalling the intended design. Under a random
split with 1024 ECFP bits, close analogues sit on both sides and the descriptor block's contribution is
systematically understated. The docstring presents both analyses as "statistical, not hand-waving,
justification" without the split caveat.

### Web server, deployment, security, accessibility

## BS-M-22 · Submission blockers: no server URL, no contact, no availability commitment `[verified]`

```
manuscript/NAR_WebServer_BrainSafe_built.md:51   freely available at [URL].
manuscript/NAR_WebServer_BrainSafe_built.md:922  [URL] with no login requirement.
manuscript/NAR_WebServer_BrainSafe_draft.md:51   freely available at [URL].
manuscript/NAR_WebServer_BrainSafe_draft.md:717  [URL] with no login requirement.
CITATION.cff:5                                   family-names: "[To be completed]"
```
Also unfilled: `:3` `[Author list to be finalised]`, `:7` `[email]`, `:921` `under [license]`,
`:926` Funding `[To be completed.]`.

| NAR requirement | Status |
|---|---|
| Freely available, no login | **PRESENT** (`--allow-unauthenticated`; no auth code) |
| Live URL, stated in the abstract | **ABSENT** — four placeholders; no deployment appears to exist |
| Contact / support address | **ABSENT** (README, About tab `app.py:2657-2722`, `api.py:60-87`) |
| Maintenance commitment (≥2 years) | **ABSENT** — repo-wide search returns no hits |
| Interface screenshot | **ABSENT** — `:936` reads *"(screenshot to be inserted)"* |
| Browser compatibility statement | **ABSENT**; `app.py:981-982` uses unguarded CSS `:has()` |
| Manuscript describes the server | **PRESENT** and done well (`:780-824`: input, output, export, batch, deployment) |
| Limitations statement | **PRESENT** and extensive (`:834-915`) |
| Example / tutorial | **PARTIAL** — nine well-chosen one-click examples, no interpretive walkthrough |
| Help pages | **PARTIAL** — strong in-app glossary, nothing static |
| Licence permits academic use | **PRESENT** (MIT) |

## BS-M-23 · Zenodo DOI is absent from the README and the manuscript `[verified + reported]`

`10.5281/zenodo.21858576` is stated **byte-identically** in five places (`REPOSITORY_MAP.md:10`,
`AUDIT_PACKAGE/00_READ_ME_FIRST.md:74`, and three reviewer-package files) — **MATCHES** — but appears
in neither `README.md` nor the manuscript. The manuscript's Data availability (`:919-922`) offers only
the GitHub URL, and **the trained models are not in git**. As written, it points a reader to a
repository that does not contain the models it promises.

## BS-M-24 · The documented `docker run` command produces a dead server `[verified]`

`serve.py:27` — `UI_PORT = os.environ.get("PORT", "7860")`
`Dockerfile:59` — `EXPOSE 8501 8000`; `Dockerfile:64-65` health-checks `localhost:8501`
`Dockerfile` — contains **no** `ENV PORT` (verified: `ENV` appears only at lines 9 and 46)
`README.md:88` — `docker run -p 8501:8501 brainsafe-ai`

Streamlit binds 7860; port 8501 maps to nothing. A reviewer following the README gets a container that
starts, logs cleanly and serves nothing. The healthcheck fails permanently wherever `PORT != 8501`,
including Hugging Face Spaces (`deploy/huggingface/README.md:7` sets `app_port: 7860`), so orchestrators
restart-loop it. `DEPLOY.md:4` states the default is 8501, which is false.

## BS-M-25 · TLS certificate verification is silently disabled on fallback `[verified]`

`app.py:2745-2750`
```python
requests.packages.urllib3.disable_warnings()
...
for verify in (True, False):  # verified first, fall back for SSL-intercepted networks
    try:
        r = requests.get(u, timeout=12, verify=verify)
```
A network attacker returns an arbitrary SMILES and the server renders a complete, confident report
**for a different molecule**, with the warning suppressed one line earlier. For a scientific tool, a
wrong answer that looks right is worse than an outage. **Fix:** delete the fallback; honour
`REQUESTS_CA_BUNDLE`/`SSL_CERT_FILE` if intercepting networks must be supported.

## BS-M-26 · API accepts an unbounded body, with no rate limit and no timeout `[verified]`

`api.py:178-179`
```python
n = int(self.headers.get("Content-Length") or 0)
payload = json.loads(self.rfile.read(n) or b"{}")
```
`n` is attacker-controlled and uncapped. `ThreadingHTTPServer` (`:198`) with
`protocol_version = "HTTP/1.1"` (`:125`) and no `Handler.timeout` — classic slowloris. `MAX_BATCH = 300`
caps compound count, not body size or string length. `DEPLOY.md:151-155` delegates to "an institutional
gateway", but `deploy/gcloud/deploy_cloud_run.md:82` deploys `--allow-unauthenticated` with no Cloud
Armor: the stated mitigation is not the one deployed.

## BS-M-27 · Batch amplifies into PubChem with no aggregate deadline `[verified by inspection]`

`api.py:185-188` and `app.py:2402-2406` resolve up to 300 names serially, each up to two
`requests.get` at `timeout=12`. Worst case: one request holds a thread ~2 hours while issuing 600
outbound calls to NCBI — a self-DoS and an amplification vector against a third party with published
rate limits. Cloud Run's `--timeout 300` kills it first, so the user gets a failure with no partial
results.

## BS-M-28 · Model-fetch failure is swallowed and reported as success `[verified]`

`app.py:296-300`
```python
except Exception:
    return True          # never let the fetch layer break a working local deployment
```
Any exception becomes the assertion "models are present and verified"; `load_models` then dies with an
unhandled `FileNotFoundError` instead of the designed message at `:306-307`. The comment's goal is
already met by `model_fetch.py:89-90`. This defeats the integrity guarantee that BS-P-02 otherwise
provides.

## BS-M-29 · `compress_models.py` reports success even when models fail to load `[reported]`

`compress_models.py:99-103` swallows every exception into a CSV text field and continues; `failed`
collects only prediction-changed models, so `:111-113` prints *"Every recompressed model reproduces its
predictions exactly"* even when models errored. The guarantee also rests on **6 probe molecules**
(`:42-45`) across 60+ models.

## BS-M-30 · Accessibility: contrast, focus indicator, alt text, unlabelled central figure `[reported]`

The palette work is deliberate and well-documented (`app.py:272-281` removes red for
deuteranopia/protanopia and pairs every status with a word). Four concrete failures remain:

| Issue | Evidence | Measurement |
|---|---|---|
| Body text below AA | `MUTE2 = "#94A3B8"` (`:271`), used for **every table column header** (`:972-974`), KPI units (`:970`), SVG axis labels (`:2134`, `:2154`) | **2.56:1**; AA needs 4.5:1 |
| Wordmark | `GOLD = "#F0A500"` (`:1102`, `:2204`) | **2.08:1** |
| Focus ring | `--ring:0 0 0 3px rgba(240,165,0,.34)` (`:873`) with global `outline:none` (`:906-907`) | ~**1.3:1**; SC 1.4.11 needs 3:1, and `box-shadow` does not render in `forced-colors` mode, so High Contrast users get **no focus indicator at all** |
| Missing alt | `:1126`, `:1129` (logos), `:1795` (the nearest-analogue structure, which *is* the evidence for the confidence claim) | no `alt` attribute |
| Unlabelled figure | `build_network_svg` `:1348-1350` — no `role="img"`, no `aria-label`, no `<title>` | `build_profile_svg` `:2125-2127` does it correctly; the pattern exists and is not applied to the main diagram |

The rest of the palette passes (MUTE 5.46:1, AMBER 5.05:1, GREEN 5.02:1, BLUE 5.99:1).

## BS-M-31 · Applicability-domain sweep runs six times per render `[reported]`

`assess_domain()` (`app.py:413-430`) is **not cached** and sweeps ~61,000 fingerprints per call,
invoked at `:1724`, `:1753`, `:2034` (twice), `:2049`, `:2169`. `render_exports` (`:2277-2314`) is the
multiplier: `st.download_button` requires payloads eagerly, so CSV, HTML report, JSON and SVG are all
regenerated on every script run. Streamlit reruns on every widget interaction. Model loading itself
**is** correctly cached (`@st.cache_resource` at `:288, 303, 328, 338, 348, 510, 549, 617, 2431`).

## BS-M-32 · The stated memory rationale for the single-container design is false `[verified by inspection]`

`serve.py:52-57` runs the API in a **thread** of the parent and Streamlit in a **separate process**.
`@st.cache_resource` is per-process, so both load the full ~800 MB model set independently. The claim
appears three times (`serve.py:4-6`, `Dockerfile:56-58`, `DEPLOY.md:3`). `DEPLOY.md:108-110` recommends
"at least 2 GB", which will **OOM**; Cloud Run's `--memory 4Gi` survives by accident.

## BS-M-33 · `coverage_modelled()` cannot list BBB, hERG or A2A `[reported]`

`app.py:1857` filters group members by `m in live`, where `live` derives from `TARGET_KIND` — which
contains **neither BBB, nor hERG, nor A2A**. The groups at `:1859`, `:1871`, `:1873` name them as
members, so all three are silently dropped from the user-facing coverage panel, **including the BBB
gate and hERG, the only safety endpoint**. `app_health.py`'s coverage check shares the same `live` set
and cannot catch it. Separately, `Nav1_5` is in `TARGET_KIND` but absent from `KNOWLEDGE_GRAPH`, so it
is deployed and scorable but contributes to no disease.

## BS-M-34 · Figure numbering is off by one, and the interface figure does not exist `[reported]`

| Manuscript | File embedded | Exists |
|---|---|---|
| Figure 1 (`:930`) | `Figure1_endpoint_rationale.png` | yes |
| **Figure 2 (`:936`)** | **none — "(screenshot to be inserted)"** | **no** |
| Figure 3 (`:940`) | `Figure2_mechanism.png` | yes |
| Figure 4 (`:947`) | `Figure3_model_selection.png` | yes |

Orphans on disk referenced by nothing: `Figure4_validation.png`, `Figure11_decision_analysis.png`. The
Figure 7 legend (`:967-972`) is printed **before** the Figure 6 legend (`:974-977`).

## BS-M-35 · Four documents describe a repository layout that no longer exists `[reported]`

`docs/BS_MODEL_CARD.md` — 37 KB, and the document a reviewer is most likely to open after the
manuscript — references ~20 non-existent scripts (`BS_brain_predict.py`, `BS_fetch_endpoints.py`,
`BS_train_endpoints.py`, `BS_external_validation.py`, …), three deleted model directories
(`models_brain/`, `models_brain_reg/`, `models_genuine/`) and eight absent JSON reports. Same for
`docs/BS_LLM_benchmark_protocol.md`, `docs/BS_BENCHMARK_ANALYSIS.md`, `docs/EXPANSION_PLAN.md` and
`docs/decisions_log.md` (which point at `archive/legacy/`, since consolidated).

## BS-M-36 · No test suite and no continuous integration `[verified]`

```
ls .github                                   → No such file or directory
find . -name "test_*.py" -o -name conftest.py -o -type d -name tests   → no tracked results
ls pytest.ini tox.ini setup.py pyproject.toml Makefile noxfile.py      → all absent
```
The substitute, `app_health.py`, is a genuinely good 20-check smoke test, but it is invoked by hand and
nothing enforces it. **Fix:** one GitHub Actions workflow running `docker build` + `app_health.py`. Add
unit tests for `featurize`, the label thresholds and the split logic — the last would have caught
BS-C-06.

## BS-M-37 · `app_health.py` downgrades a package-version mismatch to a warning `[verified]`

`app_health.py:92-120` implements exactly the right check, with a correct rationale ("for pickled
estimators it can change predictions silently"), then puts a mismatch into `warns` at line 117 so
`main()` still returns 0. Under this repository's pinned `brainsafe_env` the models load cleanly and
**no mismatch exists today** — I verified this. The finding is that the gate would not stop a future
drift. **Fix:** promote a mismatch on `scikit-learn`, `rdkit` and `numpy` to a failure.

### Data layer

## BS-M-38 · Held-out measured inactives are eligible to be drawn as decoys `[reported]`

`train_binders_hybrid.py:100-105` excludes only the target's **actives** from decoy eligibility:
```python
aset = set(act)
idx = [i for i in np.where(cand)[0] if bg_ok[i] not in aset]
```
The target's measured **inactives** — half of which are held out as `ina_hold` (`:92`) and used both to
compute the headline AUROC *and* to set the deployed threshold — are not excluded. Measured
eligibility: A1 79 of 219 (36%), GBA1 123 of 262 (47%). Whether a given compound is actually drawn
depends on a seeded shuffle over ~159k candidates, so the realised overlap is not recoverable without
replaying training; the channel is open and unguarded. The comparison is also a **raw SMILES string**
test, so a salt form of an active is not excluded. **Fix:** exclude `ina_all` too, keyed on InChIKey,
then re-derive every binder threshold and AUROC.

## BS-M-39 · `year` aggregation is inconsistent across fetchers, biasing every expansion temporal result `[reported]`

`rebuild_endpoints.py:59,79` uses `year=("year","min")` for the 11 core targets; `fetch_batch2.py:92`,
`batch3:69`, `batch4:120`, `batch5:69` and `fetch_new_targets.py:75` use **`max`** for the ~48
expansion endpoints. The deleted original fetcher was explicit that `min` is correct ("per-compound
EARLIEST document_year kept for time-splits"). Using `max` dates an old, well-known compound by its
most recent re-measurement and pushes it into the "future" test set, so **every temporal result
reported for an expansion endpoint is optimistically biased.** Compounds themselves do not cross the
cutoff (BS-P-05); the bias is in which compounds land on which side.

## BS-M-40 · A superseded provenance table is still shipped as current, and the script that produced it is still runnable `[reported]`

The project correctly identified, measured and reverted a decoy/source-bias leak:
`results/tables/inactives_audit.csv` records GSK3B scaffold AUROC inflating from **0.937 to 0.989**
when 4,276 bulk PubChem HTS inactives at median Tanimoto 0.288 were added. That is exactly the failure
mode worth catching, and `merge_pubchem_inactives.py:43`'s `MIN_SIM = 0.40` hard-negative filter is the
right replacement. **But** `results/tables/inactives_merge_provenance.csv` still records the reverted
state as fact (`GSK3B,4617,4617,4276,0.931,0.5`), while live `GSK3B.csv` has 341 inactives, not 4,617.
And `src/brainsafe/data/rebuild_with_inactives.py` — which performs the unfiltered bulk merge with no
similarity filter (`:65`) — is still listed in `SCRIPT_INDEX.csv` under "Data acquisition" with no
superseded flag. Running it would silently re-introduce the leak. **Fix:** archive both.

## BS-M-41 · The applicability-domain and read-across indexes fingerprint salts `[reported]`

`build_ad_reference.py:27-41` and `build_readacross_index.py:66-75` use bare `Chem.MolFromSmiles` with
no fragment stripping, unlike the featuriser. Measured: **3,061 of 158,890** AD reference structures
and **4,780 of 153,038** read-across structures carry a counter-ion. For those, the deployed Tanimoto
similarity is computed on the salt while the query is salt-stripped, so the "nearest measured analogue"
and the in/out-of-domain flag are not a like-for-like comparison.

## BS-M-42 · Backup files are globbed into the deployed read-across index as phantom targets `[reported]`

`build_readacross_index.py:37-38` globs `data/endpoints/*.csv` and takes the filename stem as the
target name. Verified in the shipped pickle: `n_targets 72`, including **`OX1.chembl_only` and
`OX2.chembl_only`** — backup copies written by `merge_pubchem_inactives.py:66-68`. Users see
`OX1.chembl_only` presented as a target.

## BS-M-43 · ChEMBL is queried live with no version pin `[reported]`

Every fetcher hits `https://www.ebi.ac.uk/chembl/api/data/activity.json`, which serves the current
release; `grep -rn "CHEMBL_VERSION\|chembl_37" src/brainsafe/data` returns nothing. The claim "ChEMBL
version 37 (release 2026-05-01)" lives only in `measured_endpoints_SOURCE.md`. It happens to hold today
(a live `/status.json` check returned `ChEMBL_37`, 2026-05-01), but on the day ChEMBL 38 ships every
script returns different data and, with no cache committed and no checksum, nothing reveals the drift.
**Fix:** call `/status.json`, record the version into every output, and abort on mismatch.

## BS-M-44 · Provenance is absent for the ~48 expansion endpoints, and identifier verification is inconsistent `[reported]`

`data/raw/measured_endpoints_SOURCE.md`, `data/external/SOURCE.md` and `data/adme/SOURCE.md` are good:
release versions, access dates, URLs, target and UniProt IDs, licences. But `data/endpoints/` gained 48
files between 2026-07-30 and 2026-08-06 with provenance existing only as Python literals — no
SOURCE.md, no access date, no query record. No retrieval date is written into any data file. And while
`fetch_batch4.py:63-74` has an exemplary `verify()` that aborts if the live ChEMBL preferred name or
organism does not match (written after four of twenty-six recalled identifiers pointed at the wrong
protein), `fetch_batch2.py`, `fetch_batch3.py`, `fetch_new_targets.py` and `fetch_readacross_targets.py`
have **no `verify()` call** despite docstrings claiming each identifier was name-checked. No checksum
or manifest exists for any raw download.

## BS-M-45 · A single hard-coded activity threshold is applied to all 59 targets `[reported]`

`rebuild_endpoints.py:40-41` — `1 if pvalue >= 6.0 else (0 if pvalue < 5.0 else -1)` — repeated as a
literal, not imported, in six other fetchers, with a second cut (`ACTIVE_P = 7.0`) governing the 37
binder models and a third in `fetch_readacross_targets.py:58`. 1 µM is a defensible generic cut, but it
is applied unchanged to a kinase, an ion channel, a peptide-receptor antagonist class (CGRP, potencies
to pChEMBL 10.3) and a lysosomal hydrolase, with no per-target justification anywhere in `docs/`.
Related: no fetcher filters on `assay_type`, `confidence_score` or `target_organism`, so binding
constants and functional potencies from different technologies are median-pooled —
`fetch_new_targets.py:5` claims a "human single-protein binding/functional assay" filter that its code
does not implement.

## BS-M-46 · Silent page-limit truncation and exception swallowing in the data layer `[reported]`

`adme/fetch_pgp_substrate.py:48` caps at 6,000 activities (`for _ in range(6)`) and exits silently if
more remain; `fetch_pubchem_inactives.py:41-43` has undocumented caps (`MAX_ASSAYS = 120`,
`CAP_INACTIVE_CIDS = 30000`, `FETCH_SMILES_CAP = 12000`). Ten `except Exception: pass` / `continue`
blocks turn network and parse failures into missing data — the worst is
`build_ad_reference.py:29-33`, where a malformed endpoint CSV drops silently out of the applicability
domain, after which every compound from that endpoint would be reported out-of-domain with no trace.
Related: `requests.get(..., verify=False)` appears in eight data-layer fetchers, whereas
`fetch_bindingdb.py:48-53` and `fetch_pubchem_inactives.py:47-52` do it properly (certifi first).

---

# MINOR FINDINGS

**BS-N-01 · `featurize.py` documents a property folded fingerprints do not have `[verified]`.**
`featurize.py:6-8` claims the fingerprint is *"collision-free-by-construction … bit k always means the
same environment"*. Folding to 1024 bits is what creates collisions. Measured on the first 2,000 BBB
compounds: **11,986 distinct Morgan radius-2 environments** map into 1,024 bits, so at least **10,962**
must share a bit. Minor in code, **Major if it reaches the manuscript**. **Fix:** delete the clause.

**BS-N-02 · Scaffold groups are computed on raw SMILES while features are computed on the desalted
parent `[verified]`.** `train_rf.py:66` calls `Chem.MolFromSmiles` directly; `featurize.py:60-62` strips
salts first. Measured: 45 feature-identical clusters (133 rows) in BBB can be split across scaffold
folds; BACE1 4, AChE 3, D2 2, hERG 0. **This is the good news in the leakage story** — the scaffold
split contains the duplicates almost entirely (1.7% of BBB), so the scaffold AUROC of 0.9197 is **not**
materially contaminated. (Deduplicated scaffold CV gives 0.8826, but that reflects a changed evaluation
set, not leakage, and must not be reported as such.) **Fix:** route the scaffold computation through
`_mol_from_smiles`.

**BS-N-03 · The scaffold is atom-preserving, not generic, contradicting the docstring and the
manuscript `[verified]`.** `train_rf.py:62` says *"generic Bemis-Murcko scaffold"* and the manuscript
says *"generic scaffolds"*, but line 66 calls `MurckoScaffoldSmiles`, which retains atom types; the
generic framework needs `MakeScaffoldGeneric`. This matters in the unfavourable direction: more groups
means pyridine and benzene analogues of one framework land in different folds, so the split is
**weaker** than claimed. Also, acyclic compounds each get their **own** group
(`train_rf.py:69`, `f"_none_{len(mapping)}"`) — 311 in BBB, and 1,446 of BBB's 2,441 groups are
singletons — so for the acyclic fraction the scaffold split degenerates to a random split. The same
code is duplicated verbatim at `binder_cv_per_fold.py:81-90`. **Fix:** correct the wording (cheap) or
apply `MakeScaffoldGeneric` and rerun (what the claim implies); assign acyclic compounds one shared
group.

**BS-N-04 · `REPOSITORY_MAP.md` says the training data is not committed; it is `[verified]`.**
`REPOSITORY_MAP.md:8-11` states *"Neither is committed"* of `data/` and `models_rf/`. `git ls-files`
shows **165 tracked files under `data/`, including all 61 endpoint CSVs**, and 64 metadata JSONs under
`models_rf/`. Committing the training tables is **good**; the map understates what a reviewer gets.

**BS-N-05 · ECE off-by-one silently drops `p == 1.0` `[reported]`.** `calibrate.py:44-48` —
`np.digitize(1.0, linspace(0,1,11))` returns 11, so `idx = 10` falls outside `range(n_bins)`. Measured:
**1,023 / 40,040 (2.6%)** of OOF predictions excluded, up to 7.4% for GSK3B. Harmless today only
because those samples happen to have `y.mean() == 1.0` — the bug hides exactly the confident-and-wrong
case ECE exists to catch. **Fix:** `np.clip(np.digitize(p, bins) - 1, 0, n_bins - 1)`.

**BS-N-06 · Mean calibrated ECE rounds to 0.013, not 0.012 `[verified]`.** Recomputed: 0.01274.
Stated as 0.012 in `README.md:33`, the manuscript `:264` and three docs.

**BS-N-07 · Decoy selection uses 2048-bit fingerprints while the model uses 1024 `[reported]`.**
`train_binders_hybrid.py:45` — `GetMorganGenerator(radius=2, fpSize=2048)` against
`featurize.py:29` `MORGAN_BITS = 1024`. Probably deliberate; nothing documents it.

**BS-N-08 · Positional feature indices with no assertion `[reported]`.** `train_binders_hybrid.py:75`
and `train_receptor_binders.py:63-64` — `bg_mw, bg_lp = bgX[:, -12], bgX[:, -11]` with the comment
*"'mw' is first descriptor after 1024 fp bits"*. Silently wrong if `_DESCRIPTORS` order changes; no
check against `feature_names()`.

**BS-N-09 · Decoy:active ratio is not the documented 3:1 `[reported]`.**
`train_binders_hybrid.py:104` — `need = max(DECOY_RATIO*len(act) - len(ina_train), len(act))` gives
`3 − n_ina_train/n_act`, floored at 1:1 (KEAP1: 2.57:1). Class balance therefore varies per endpoint
while `class_weight="balanced"` absorbs it, so the calibrated probability scale is not comparable
across the panel — yet the app applies a single `BINDER_FLOOR = 0.40` (`app.py:501`).

**BS-N-10 · Per-endpoint hold-out splits are not reproducible in isolation `[reported]`.**
`scaffold_holdout_panel.py:47` uses one module-level `default_rng(2026)` advanced by every target in
sequence (`:104`), so a single-target run differs from a full-panel run. `models_rf/holdout/` is a
composite of several runs, and the pooled recall cannot be reproduced by any single invocation. The
correct pattern is already used at `binder_cv_per_fold.py:108`.

**BS-N-11 · Hold-out inactives are split randomly, not by scaffold `[reported]`.**
`scaffold_holdout_panel.py:114-116`, undercutting the docstring's claim that the threshold is never
informed by the compounds it scores.

**BS-N-12 · Recall reconstructed from a rounded rate `[reported]`.** `scaffold_holdout_report.py:44-46`
— `k = int(round(r * n))` where `r` was stored to 3 dp, so the Wilson intervals and the pooled
11,914/15,069 carry avoidable rounding error. Also `FLOOR = 0.06` (`:26`) detects a floor set at `0.05`
(`panel.py:152`), so thresholds in (0.05, 0.06] are misclassified as collapsed.

**BS-N-13 · H6's "never seen in training" index omits 11.1% of the corpus `[reported]`.**
`inv_clinical_indication.py:161-179` builds it from two derived artefacts rather than the training
tables; 17,997 of 161,814 endpoint SMILES are in neither, and **8 of the 162 drugs (4.9%)** in the
"never seen" stratum are in fact present.

**BS-N-14 · Permutation p-values reported at their floor as point estimates `[reported]`.**
`inv_disease_layer.py:144` with `N_PERM = 200` yields 0.004975; `REPORT.md:9` renders it as "p=0.005".
Should be `p < 0.005`, and 200 permutations is low for a headline claim.

**BS-N-15 · HTTP 204 sent with a body `[verified]`.** `api.py:142-143` calls `_send({}, 204)`, and
`_send` (`:131-140`) always writes a body and `Content-Length`. RFC 7230 §3.3.2 forbids this; with
keep-alive a conforming client reads `{}` as the next response. This is the CORS preflight, so **every
browser client hits it first**.

**BS-N-16 · Under-escaped URL path segment `[verified]`.** `app.py:2747` — `urllib.parse.quote(name)`
defaults to `safe="/"`. The host is fixed so cross-host SSRF is impossible, but
`x/../../../compound/cid/2244` can reach other PubChem REST paths. **Fix:** `quote(name, safe="")`.

**BS-N-17 · HTML escapers do not escape quotes and are used in an attribute `[verified]`.**
`app.py:2173-2174` and `:1199-1200` escape `&`, `<`, `>` only; the result is interpolated at `:2233`
into `alt="structure of {esc(name)}"`. Not reachable today (PubChem will not resolve a name containing
a quote), but the guard is incidental and the report is built to be shared. **Fix:**
`html.escape(s, quote=True)`.

**BS-N-18 · Unbounded input length `[verified]`.** No cap in `app.py:2762`, `:2778`, or `api.py:163`,
`:180-188`. A multi-megabyte "SMILES" goes straight into `Chem.MolFromSmiles`.

**BS-N-19 · Network failure reported as a spelling mistake `[verified]`.** `app.py:2757-2759` catches
everything and returns `None`; `:2774-2775` then advises *"Check the spelling"*. A PubChem outage, a
DNS failure and a typo are indistinguishable.

**BS-N-20 · Batch error rows put the unparsed input in the `smiles` column `[verified]`.**
`app.py:2404`. `api.py:187` gets this right by omitting the field.

**BS-N-21 · Base image and several dependencies unpinned `[verified]`.** `Dockerfile:7` —
`FROM python:3.13-slim`, a floating tag. `requirements.txt:17-23` leaves `requests`, `joblib`,
`streamlit`, `pypandoc_binary` unpinned and `xgboost>=3.3` open. Streamlit is the risk:
`app.py:897-899` and `:1020-1090` target private test IDs (`[data-testid="stMainBlockContainer"]`) and
BaseWeb attributes, so a minor release can silently break the entire visual identity. The header
comment at `requirements.txt:2-8` argues carefully for pinning scikit-learn and does not apply the
argument to the framework that renders everything.

**BS-N-22 · `.dockerignore` omits `models_rf/` `[verified]`.** A developer machine holding 0.78 GB of
models uploads all of it into the build context for `RUN python model_fetch.py` to no-op.
`.gcloudignore:6` gets this right.

**BS-N-23 · Dead code that still advertises a withdrawn endpoint `[reported]`.** `COVERAGE_YES`
(`app.py:1806-1821`) is superseded by `coverage_modelled()` (`:1848-1891`); the comment at `:1822-1828`
explains why hand-maintained lists were abandoned, and the abandoned list was left in place. It still
lists `Nav1.1`, documented as **withdrawn** at `:204-207`. The `tags` list at `:1130-1131` is similarly
stale.

**BS-N-24 · Personal operational detail in a public deployment document `[verified]`.**
`deploy/gcloud/deploy_cloud_run.md:41` embeds a real GCP project id
(`project-47d4c9c7-5bfb-4527-b67`); `:14` records the author's billing credit (`₹28,694`); `:80`
contains `cd D:\BRAINSAFE_AI`. This is also the **only** hard-coded absolute path anywhere — `app.py`,
`api.py`, `serve.py` and `model_fetch.py` all correctly derive `ROOT` from `Path(__file__)`.

**BS-N-25 · API cannot report missing models `[reported]`.** `api.py:51` calls `app.load_models()`,
whose failure path calls `st.stop()`; outside Streamlit that raises `StopException`, caught at `:169-171`
and returned as a bare `{"error": "internal error"}` 500.

**BS-N-26 · Model-fetch lock not released on every path `[reported]`.** `model_fetch.py:101-115` —
`lock.touch()` with no `try/finally`; the path-violation `return False` at `:147` exits without
unlinking. A stale lock blocks startup for 900 s.

**BS-N-27 · Silent exceptions that weaken exclusion guarantees `[reported]`.** 49 broad `except` blocks
across `src/`, `app.py`, `api.py`. Most are defensible RDKit guards, but six sit on paths where firing
*shrinks an exclusion set and inflates apparent specificity or novelty* —
`external_100.py:129-130,163-165`, `noncns_specificity.py:66-68,85-87,93-95`,
`noncns_specificity_fast.py:62-64,74-76` — and none increments a counter. That is the wrong direction
for a failure to be silent in.

**BS-N-28 · Truncated, non-deterministic canonical safety re-check `[reported]`.**
`noncns_specificity_fast.py:85` — `{c for c in (canon(s) for s in list(active_raw)[:40000]) if c}`.
`active_raw` is a `set`, so `list(set)[:40000]` is an arbitrary slice in an order not guaranteed stable
across processes, and the exclusion set is far larger than 40,000.

**BS-N-29 · `reliable == False` comparisons treat `None` inconsistently `[reported]`.**
`final_thresholds.py:111` reports silently-unevaluated endpoints as reliable-by-omission
(`None == False` is `False`); `calibrate_background_specificity.py:113` coerces the same case to 0 and
reports it as a failure — opposite conventions in sibling scripts.

**BS-N-30 · Docstring/code count mismatches `[reported]`.** `adme/train_adme.py:5` says "One classifier
and five regressors"; the code has **two** classifiers and **seven** regressors (`:35-37`).
`_train_small.py:1-8` carries `train_measured_label_holdout.py`'s docstring naming SIRT1 and Nav1.5
while training TAAR1, GluA2 and Nav1_1.

**BS-N-31 · `docs/METHODS.md:91` overstates the ADME set `[reported]`.** "~23,000 measured compounds"
against the true 21,694 (sum of `adme_cv_summary.csv`); `docs/ADME_RESULTS.md:9` says "About 21,700".

**BS-N-32 · "Knowledge graph of 52 targets" counts a non-target `[reported]`.** One of the 52 keys is
`NEURO`, the neuroprotection axis. The accurate statement is 51 targets plus an axis.

**BS-N-33 · Duplicate endpoint tables shipped beside their successors `[verified]`.**
`data/endpoints/` contains both `OX1.csv` / `OX1.chembl_only.csv` and `OX2.csv` /
`OX2.chembl_only.csv`, with nothing stating which is authoritative.

**BS-N-34 · "Supplementary training record" names no existing file `[reported]`.** Manuscript `:114`;
`supplementary/` contains `STable0`–`STable15`, none so titled.

**BS-N-35 · `provenance_audit.py:45-55` allowlist omits most quoted tables `[reported]`.** `PANEL_FILES`
lists 6 artefacts, omitting `external_bbb_validation.csv`, `rf_conformal.csv`, `rf_temporal.csv`,
`noncns_specificity_summary.csv`, `integrity_audit.csv` and `external_100_summary.csv`. The file's own
docstring argues for enumerating `results/tables/*.csv` rather than hand-maintaining a list.

**BS-N-36 · Streamlit's skip-to-content affordance is hidden `[reported]`.** `app.py:898-899` sets
`display:none !important` on `[data-testid="stHeader"]`, removing the skip link with the deploy menu.

**BS-N-37 · Colour-only encoding in the model-comparison table `[reported]`.** `app.py:2488-2489` marks
the best algorithm per endpoint **only** by green text — the one place the interface breaks its own
stated rule at `:276-277`.

**BS-N-38 · Ref 20's DOI is very likely non-resolving `[reported]`.** `10.1351/goldbook.11443` — IUPAC
Gold Book DOIs take the form `goldbook.XNNNNN` with a letter prefix. The entry also has no author, and
it stands in for the ECFP-4 primary method paper.

---

# TRIVIAL

- **BS-T-01** `packages.txt` duplicates `Dockerfile:17-19`; two places to update.
- **BS-T-02** `DEPLOY.md:117-118` instructs `curl https://<your-url>:8000/health`, contradicting
  `serve.py:29-31`, which correctly states a Space publishes only one port.
- **BS-T-03** `app.py:783-786` — `_b64` returns `""` on any exception; benign, with a documented
  fallback at `:1126-1129`.
- **BS-T-04** Last-digit disagreements between manuscript Table 1 and `docs/RF_CV_RESULTS.md`
  (MAO_B 0.955 vs 0.954, SERT R² 0.389 vs 0.388, D2 SD 0.051 vs 0.052) — all rounding of the same CSV.
- **BS-T-05** Reference typography: double space (ref 5), lower-case `10.1016/b978-` (refs 9–11), curly
  apostrophe (ref 10), `insomnia?.` (ref 11).
- **BS-T-06** Truncated URL at `docs/BS_BENCHMARK_ANALYSIS.md:105`.
- **BS-T-07** `inversion/H1H2H3.log` is from a superseded run recording `VERDICT H2: WEAKENED` against
  `REPORT.md`'s `REFUTED`; the two emitters disagree (`inv_disease_layer.py:194` vs `summarise.py:47`).
- **BS-T-08** Dead imports: `average_precision_score` and `featurize_one` in
  `train_receptor_binders.py`; `_scaffold_groups` in `feature_analysis.py:37-38`.
- **BS-T-09** `train_batch2.py:128-133` references `hard`, defined only in the `else` branch at `:90`;
  correct today purely by the guard.

---

# What is genuinely good

These are true, and a reviewer should be pointed at them.

- **BS-P-01 — Exact reproducibility.** Five headline CV numbers regenerated to four decimal places from
  a cold start. Seed 42 set and honoured. Rarer than it should be.
- **BS-P-02 — `model_fetch.py` is exemplary.** DOI-pinned Zenodo source, SHA-256 on the archive (`:137`)
  *and* all 195 extracted files (`:161-166`), fatal on mismatch, `tarfile` `filter="data"` plus an
  explicit path-containment check (`:143-152`), and a lock file written after a real concurrency
  incident. **The "unpickling unverified remote data" risk is not present.**
- **BS-P-03 — Calibration is done correctly.** `cross_val_predict` on out-of-fold probabilities
  (`calibrate.py:59`); deployment models use `CalibratedClassifierCV(cv=5)`, never `cv="prefit"`. The
  classic leak I went looking for is absent.
- **BS-P-04 — The conformal implementation is correct.** `rf_conformal_temporal.py:41-57`: `tr`, `cal`,
  `te` are disjoint slices of one permutation, the model is fitted on `tr` only, and the Mondrian
  quantile is the correct finite-sample `ceil((n+1)(1−ε))` form.
- **BS-P-05 — The temporal split is verifiably compound-disjoint.** Zero compounds cross the cutoff in
  any of the 11 temporally-tested endpoints, because the endpoint tables are aggregated one row per
  SMILES (`fetch_new_targets.py:74-75`). Worth stating in the manuscript.
- **BS-P-06 — The ADME layer has no leakage.** `train_adme.py` **imports** `RF_COMMON`, `SEED` and
  `N_SPLITS` from `train_rf.py` rather than copying them — the only training script that does — fits no
  calibrator, selects no thresholds, and has zero duplicate SMILES across all nine files.
- **BS-P-07 — Clean install works** from `requirements.txt` with every pin resolving.
- **BS-P-08 — No secrets, no hard-coded paths in source.** 679 tracked files, no credentials, no
  binaries; `brainsafe_env/` correctly ignored.
- **BS-P-09 — Training data is committed.** All 61 endpoint CSVs are in git, so a reviewer can rerun
  training without regenerating from live APIs.
- **BS-P-10 — Model artefacts verify**: all 195 manifest files present with matching sizes.
- **BS-P-11 — `app_health.py` is a real gate**: 20 substantive checks including directional
  pharmacology and export well-formedness, exiting non-zero on failure.
- **BS-P-12 — The manuscript's server description and limitations sections are strong** (`:780-824`,
  `:834-915`) and meet the NAR requirements they address.
- **BS-P-13 — Deliberate accessibility thinking** in the palette (`app.py:272-281`), and robustness
  engineering throughout: RDKit drawing libraries guarded at import with a stated rationale
  (`app.py:36-41`), invalid SMILES, empty input and over-limit batches all handled cleanly.
- **BS-P-14 — Statistical instrumentation is better than most.** Wilson intervals in seven independent
  and algebraically identical implementations, `n` reported beside every proportion, permutation nulls,
  Cohen's κ and Fisher exact tests where appropriate. The gap is clustering and multiplicity
  (BS-M-12, BS-M-13), not absence.
- **BS-P-15 — `compress_models.py` measures its own nondeterminism** (`:74-81`) rather than assuming it
  away, and `provenance_audit.py`'s graph-fingerprint mechanism (`:106-116`) is genuinely good work that
  should be the model for the rest.

---

# Ranked summary

| Rank | ID | Sev | Finding |
|---|---|---|---|
| 1 | BS-C-01 | Critical | Deployed thresholds set from measured inactives the model was trained on — 40/49 endpoints |
| 2 | BS-C-02 | Critical | Every `sensitivity_at_threshold` measured on training positives; five docstrings claim the opposite |
| 3 | BS-C-03 | Critical | Reported false-positive rates are in-sample and bounded by construction |
| 4 | BS-C-04 | Critical | The "independent" background sample is the decoy pool |
| 5 | BS-C-05 | Critical | Headline specificity 0.875 measured on the applicability-domain reference itself (all 1000 at T=1.000) |
| 6 | BS-C-06 | Critical | Feature-identical duplicates split across random folds; BBB 48.3%, AUROC 0.9605 → 0.8990 |
| 7 | BS-C-07 | Critical | "306 FDA drugs absent from training" is 21.2% memorised |
| 8 | BS-C-08 | Critical | 25.9% of inversion hold-out compounds are training actives of another target |
| 9 | BS-C-09 | Critical | Inversion duplicate check tests the wrong artefact; three of six checks cannot fail |
| 10 | BS-C-10 | Critical | `integrity_audit.py` computes calibration on training data and labels it "held-out" |
| 11 | BS-C-11 | Critical | Manuscript has zero in-text citations and an unfilled reference placeholder |
| 12 | BS-C-12 | Critical | The documented rebuild command silently empties the 11 core endpoint tables on a fresh clone |
| 13 | BS-C-13 | Critical | Core tables unreproducible: fetch scripts deleted in `fea5029`, raw caches gitignored |
| 14 | BS-C-14 | Critical | BindingDB censored measurements stripped of their relation; 1,148 records used as exact |
| 15 | BS-C-15 | Critical | The `pchembl_value` filter discards the measured negative class; 39/60 endpoints >90% active, 5 with <25 inactives |
| 16 | BS-C-16 | Critical | 48 of 59 endpoints deduplicate on raw SMILES — the root cause of BS-C-06 |
| 17–62 | BS-M-01…46 | Major | Numbers and provenance (01–04), calibration (05–07), validation framing (08–11), statistics (12–14), reproducibility (15–21), server and compliance (22–37), data layer (38–46) |
| 63–100 | BS-N-01…38 | Minor | See above |
| 101–109 | BS-T-01…09 | Trivial | See above |

**Where the Criticals cluster.** Eleven of the sixteen are one of two failures repeated across
components: **a number is measured on the data it was derived from** (C-01 to C-05, C-10) or
**compounds the model cannot distinguish are counted as independent** (C-06 to C-09, C-16). The
remaining five are reproducibility and submission blockers (C-11 to C-15). Fixing the two root causes —
persist and honour holdout sets; standardise to InChIKey before splitting — resolves most of the list.

---

# Reproducing this audit

```bash
brainsafe_env/Scripts/python.exe -m pip freeze
brainsafe_env/Scripts/python.exe src/brainsafe/evaluation/app_health.py
brainsafe_env/Scripts/python.exe audit/evidence/leak_external.py       # BS-C-07
brainsafe_env/Scripts/python.exe audit/evidence/leak_internal.py       # BS-C-06
brainsafe_env/Scripts/python.exe audit/evidence/repro_cv.py            # BS-C-06 impact
brainsafe_env/Scripts/python.exe audit/evidence/repro_scaffold.py      # BS-N-02
brainsafe_env/Scripts/python.exe audit/evidence/scaffold_mismatch.py   # BS-N-02
```
All scripts are read-only and write nothing outside `audit/evidence/`.

---

# Open questions and gaps in this audit

Stated explicitly rather than papered over.

1. **The Zenodo DOI `10.5281/zenodo.21858576` was not verified as live.** I made no network calls. If
   the record is not public, `docker build` fails at `Dockerfile:45`. **Verify before submission.**
2. **BS-C-08's exact magnitude is unsettled** — 25.9% by my computation, 31.4% by a component audit
   using a different active-set definition. Pin it down before quoting.
3. **The execution order of the six `binder_modes.json` writers is unrecoverable** from the artefacts.
   File mtimes suggest the *older* `final_thresholds.py` wrote the majority of live `threshold_basis`
   values, which would mean the later `deployed`-guard fix is not reflected in the shipped thresholds —
   but this cannot be proven read-only. A reviewer will ask. Add `produced_by` and `produced_at` fields.
4. **Censored data — now answered (BS-C-14), but the ChEMBL side rests on convention.** BindingDB
   qualifiers are stripped; ChEMBL is safe only because it assigns a pChEMBL value exclusively where the
   relation is `=`. That is an undocumented API behaviour the pipeline never asserts. Add
   `&standard_relation==` explicitly.
5. **Decoy / source bias — partially answered.** The project measured and reverted a bulk PubChem merge
   that inflated GSK3B scaffold AUROC from 0.937 to 0.989 (`inactives_audit.csv`), which is good work.
   What remains untested is whether the *retained* actives (ChEMBL/BindingDB) and inactives (PubChem)
   still differ systematically. The definitive test was not run: **train a classifier to predict which
   database a compound came from.** If it succeeds, the endpoint models are partly learning provenance
   rather than biology. This is the highest-value remaining check.
6. **`ad_reference.pkl` provenance was not verified** — described as a "74k measured background library"
   in three docstrings; its composition and overlap with the endpoint training sets beyond BS-C-04 is
   unknown.
7. **Whether the shipped `*_calibrated.joblib` files correspond to `calibration.csv`** was not
   confirmed; that requires loading and re-scoring them.
8. **Figures were not inspected** for legibility, resolution or caption self-containment.
9. **The two uncommitted manuscript files** were audited for placeholders and the one substantive hunk;
   their full diff was not reviewed.
10. **`app.py` is 2,892 lines in one module.** Not itself a defect, but it is why several findings here
    concern duplicated presentation logic, and it will make every fix above more expensive.

---

*Read-only diagnostic pass. No pre-existing repository file was modified. Every quantitative claim is
traceable to a `file:line` reference, or to a script under `audit/evidence/` with its captured output.
Findings are tagged `[verified]` where I re-executed them and `[reported]` where a component audit
supplied the evidence with code quoted but I did not re-run it; every Critical finding is `[verified]`.*
