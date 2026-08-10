# BrainSafe AI — Executive Summary

**Commit `5d338d6` · 2026-08-10 · read-only diagnostic · full detail in [`AUDIT_REPORT.md`](AUDIT_REPORT.md)**

## The headline

**Nothing is fabricated.** Every published figure I tested regenerated *exactly* — BBB random 10-fold
0.9605 against 0.9605, scaffold 0.9197 against 0.9197, four endpoints to four decimal places from a
cold start. The clean install works, the smoke test passes (20 checks, exit 0), all 195 model artefacts
verify against their manifest, and `model_fetch.py` is genuinely exemplary security engineering.

The problem is not the arithmetic. It is that **several headline numbers do not measure what they are
said to measure**, and the checks written to catch that cannot fail.

**16 Critical, 46 Major, 38 Minor, 9 Trivial.** Eleven of the sixteen Criticals are one of two root
causes repeated across components:

1. **A number is measured on the data it was derived from.**
2. **Compounds the model cannot tell apart are counted as independent.**

Fix those two and most of the list resolves.

---

## Critical (16) — must be fixed before submission

### Root cause 1: measured on its own data

| ID | Finding | Evidence |
|---|---|---|
| **C-01** | Deployed thresholds are set from measured inactives the model was **trained on**. `final_thresholds.py:69-73` re-reads all inactives, undoing the holdout `train_binders_hybrid.py` built. | A1: `n_train=109`, `n_holdout=110`, `n_for_threshold=**219**`. **40 of 49 endpoints** affected. |
| **C-02** | Every `sensitivity_at_threshold` is computed on **training positives**, while five docstrings state the opposite in so many words. | `train_binders_hybrid.py:115` puts `act` into training; `:139` evaluates on the same `act`. |
| **C-03** | Reported false-positive rates are in-sample and **bounded by construction** — threshold and rate come from the same sample. | `background_fpr_at_threshold` never exceeds 0.0500 = `BACKGROUND_FPR`; `screening_background_fpr` sits at exactly 0.01. |
| **C-04** | The "independent" background sample **is the decoy pool** — both are `ad_reference.pkl`, with no mutual exclusion. | `final_thresholds.py:8-9` claims independence; `:44-47` vs `train_binders_hybrid.py:69-70`. |
| **C-05** | The headline specificity **0.875** is measured on the applicability-domain reference set itself. | **All 1,000 compounds have max Tanimoto = 1.000.** The "in domain" stratum is numerically identical to the unstratified row. |
| **C-10** | `integrity_audit.py` computes calibration on training compounds under the header *"against held-out measured inactives"*. Its leakage check also silently examines only **8 of 46** targets. | `integrity_audit.py:83` vs `:92-103`; `:139` `[:8]`. |

### Root cause 2: indistinguishable compounds counted as independent

| ID | Finding | Evidence |
|---|---|---|
| **C-06** | No deduplication anywhere. Stereoisomers and salts give **byte-identical feature vectors** and split across random CV folds. | **BBB is 48.3% duplicate rows** (3,773 of 7,805); 13,846 duplicates and 288 contradictory-label groups panel-wide. **BBB AUROC 0.9605 → 0.8990 when collapsed.** |
| **C-16** | *Root cause of C-06:* **48 of 59 endpoints deduplicate on the raw SMILES string**, not InChIKey. No standardisation before the tables are written. | `fetch_batch2-5.py`, `fetch_new_targets.py` etc. group on `"smiles"`; 4,012 salted SMILES survive. |
| **C-07** | The "306 FDA drugs absent from training" external set is **21.2% memorised** — exclusion uses full InChIKey, features do not. | 65/306 at Tanimoto 1.0; 78/306 share an InChIKey skeleton. Warfarin, Prednisone, Quinidine, Naloxone… |
| **C-08** | **25.9%** of the inversion suite's "held-out" compounds are training actives of *another* panel target, while `summarise.py:37` states none was seen in training. | DAT 68.6%, D2 63.1%, SERT 53.9%. Affects exactly the homologous families the suite's own H8 flags. |
| **C-09** | The inversion duplicate check tests `compound_library.csv`, **not the tables trained on**, so it passes while C-06 holds. Two more of the six checks assert scikit-learn contracts and cannot fail; others rest on n=4 and n=1. | `inversion_validation.csv` says "0 duplicates" against my 13,846. |

### Reproducibility and submission blockers

| ID | Finding | Evidence |
|---|---|---|
| **C-12** | **The documented rebuild command silently destroys the training data on a fresh clone.** Missing cache returns an empty frame; `main()` writes it unconditionally. All 11 core CSVs become header-only, exit code 0. | `rebuild_endpoints.py:46-48` and `:115-116`. |
| **C-13** | The core tables **cannot be regenerated**: `BS_fetch_endpoints.py` was deleted in `fea5029`, and all raw caches are gitignored — yet `SOURCE.md` still cites the deleted script. | `git log --diff-filter=D fea5029`; `ls src/brainsafe/data \| grep BS_fetch` → nothing. |
| **C-14** | BindingDB **censored measurements are stripped of their relation** and used as exact values. `">10000"` becomes exactly 10,000 nM. | `fetch_bindingdb.py:87` `lstrip("><~=")`; **1,148 censored records** in the committed cache. |
| **C-15** | The `pchembl_value` filter **discards the measured negative class**, leaving most endpoints near-degenerate. | **39 of 60 endpoints >90% active; 21 >96%; 5 with fewer than 25 inactives** (GABAA_a5 has **4**). AChE alone discards 2,210 measured inactives. |
| **C-11** | The manuscript has **zero in-text citations**, a `[To be completed]` reference placeholder, and a closing list of uncited essentials including **B3DB** — the BBB training source. | `:1003-1005`, `:1035`; regex sweep of lines 1–1000 → 0 matches. |

---

## Major (46) — grouped

- **Numbers (M-01 to M-04).** Record count stated three ways — manuscript 67,984, README and the **live app** 64,474, true sum **67,982**. `BS_MODEL_CARD.md:215` calls its own figures authoritative; **six of eight AUROCs disagree** with the live table. `MASTER_validation_summary.csv` names CGRP and RIPK1 — the two flagship endpoints — as failures while the manuscript reports them as successes. Endpoint counts appear as 13, 40, 52, 62, 63, 67, 69 and 71 across the documentation.
- **Calibration (M-05 to M-07).** Calibration **degrades BBB** (ECE 0.0403→0.0424, Brier 0.0788→0.0947) and the degraded model ships; only the mean is reported. "Every classifier is isotonically calibrated" is false — the 47 binders use sigmoid. The deployed calibrated artefact is not the one whose ECE is published.
- **Validation framing (M-08 to M-11).** **BBB has no temporal validation at all** (no `year` column). Temporal regression R² collapses to 0.042 (D2) and 0.009 (DPPH) and is not surfaced. Conformal coverage describes a refit model on a random split, and no prediction set reaches the server. Out-of-domain AUROC (0.816) *exceeds* in-domain (0.772), contradicting the script's own stated expectation, silently.
- **Statistics (M-12 to M-14).** Wilson intervals throughout treat clustered observations as independent. No multiple-comparison control anywhere across 39-target and 47-endpoint families. Twelve verdict thresholds are hard-coded with no pre-registration; two pairs actively disagree.
- **Reproducibility (M-15 to M-21, M-38 to M-46).** Six deployed endpoints have **no recorded training command**. RF hyper-parameters are copy-pasted into **four incompatible settings**; `random_state=42` is written literally 14 times. `binder_modes.json` is accreted by six scripts with contradictory semantics and no run stamp. Held-out inactives are eligible to be drawn as decoys. `year` uses `max` for 48 endpoints and `min` for 11, biasing every expansion temporal result. ChEMBL is queried live with no version pin.
- **Server and compliance (M-22 to M-37).** Four `[URL]` placeholders; no contact address; no maintenance commitment; no interface screenshot; `CITATION.cff` author is `[To be completed]`. **The documented `docker run` command yields a dead server** (`serve.py` binds 7860, Dockerfile exposes and health-checks 8501). TLS verification is silently disabled on fallback. The API accepts an unbounded body with no rate limit. Text contrast 2.56:1 and the focus ring 1.3:1 fail WCAG AA. **No tests, no CI.**

---

## What is genuinely strong

Worth saying, and worth pointing a reviewer at:

- **Exact reproducibility** — five headline numbers to four decimal places, cold start.
- **`model_fetch.py`** — DOI-pinned, SHA-256 on the archive *and* all 195 files, fatal on mismatch, `tarfile` `filter="data"` plus path containment. The remote-unpickling risk is **not present**.
- **Calibration and conformal implementations are correct** — `cross_val_predict` on out-of-fold predictions, never `cv="prefit"`; disjoint conformal calibration slices with the correct finite-sample quantile.
- **The temporal split is verifiably compound-disjoint** — zero compounds cross the cutoff in 11 endpoints.
- **The ADME layer has no leakage** and is the only training script that *imports* its configuration rather than copying it.
- **The training data is committed** (all 61 endpoint CSVs), and `app_health.py` is a real 20-check gate.
- **The decoy/source-bias leak was found, measured and reverted** by the authors themselves (`inactives_audit.csv`: GSK3B 0.937 → 0.989 when contaminated). The instinct is right; it just was not applied everywhere.

---

## Recommended order of work

1. **Stop the bleeding.** Make `rebuild_endpoints.py` refuse to write an empty table (C-12); restore the three deleted fetch scripts and commit or deposit the caches (C-13).
2. **Fix the two root causes.** Standardise to InChIKey before splitting (C-16 → C-06, C-07); persist and honour holdout sets for positives, inactives and background (C-01 to C-05, C-10).
3. **Recompute and re-report.** BBB random-split AUROC, every binder sensitivity, every specificity and FPR, the external set at n=241 — which *raises* the number to ≈0.802 and is a free win.
4. **Repair the falsification suite** so it can fail (C-09, C-08), then re-run it.
5. **Fix the data labels.** Honour censoring; recover the negative class via `standard_relation=>` (C-14, C-15).
6. **Complete the manuscript.** Citations and references (C-11); reconcile every count; deploy and fill the four `[URL]` placeholders.
7. **Then the server work**: port contract, TLS fallback, API limits, accessibility, CI.

**Two things to note about tone.** First, three of these fixes make the paper *stronger*, not weaker —
the external validation improves, and the temporal and calibration limitations are the kind of honest
reporting reviewers reward. Second, the fastest way to lose a reviewer is the docstrings in
`train_binders_hybrid.py:15-16`, `train_measured_label_holdout.py:5-8` and `final_thresholds.py:8-9`,
which assert no-leakage while the code beneath them leaks. **Correct those sentences first**, whatever
else is deferred.

---

## Caveats on this audit

- **The Zenodo DOI was not verified as live** (no network calls made). If it is not public,
  `docker build` fails at `Dockerfile:45`.
- **C-08's magnitude is unsettled**: 25.9% by my computation, 31.4% by a component audit using a
  different active-set definition. Pin it before quoting.
- **Untested:** whether actives and inactives remain separable by *source database* (the definitive
  decoy-bias test); figure legibility; the full manuscript diff.
- Findings are tagged `[verified]` (I re-executed them) or `[reported]` (evidence supplied with code
  quoted, not re-run by me) in the full report. **Every Critical is `[verified]`.**
