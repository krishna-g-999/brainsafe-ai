# Regeneration pass

Applying the committed fixes to the shipped artefacts, once, in dependency order. Every step records
what it wrote and how the numbers moved. Started 2026-08-11.

**Recovery point.** Before anything was written: `models_rf/` (801 MB, all 71 estimators), `results/`,
`data/endpoints/` and `inversion/` were copied to `/d/brainsafe_regen_backup_20260811T151047`, outside
the repository. Tracked files are additionally recoverable from git; the model binaries are not in
git, which is why the copy exists. `rebuild_endpoints.py` also wrote its own timestamped copy to
`archive/legacy/endpoints_before_rebuild_2026-08-11T181846`.

---

## Step 0 — baseline captured

| Quantity | Before regeneration |
|---|---|
| Classifier AUROC, random 10-fold | 0.947–0.969, mean **0.9604** |
| Classifier AUROC, scaffold 10-fold | 0.868–0.956, mean **0.9186** |
| Regression R², random | mean 0.6381 |
| Regression R², scaffold | mean 0.4743 |
| Binder panel | 49 endpoints, mean sensitivity **0.8942**, median 0.9280 |
| External BBB | n=306, AUROC **0.7741** |
| Endpoint tables | 65,122 rows |

Baseline tables copied to the scratchpad for row-by-row comparison.

---

## Step 1 — BindingDB re-parsed, endpoint tables rebuilt · DONE

Regenerated the per-compound BindingDB tables from the committed raw responses under the fixed
parser, then pooled the endpoint tables again.

- 21,791 compounds from exact measurements; **1,141 censored records preserved** with their relation
- Endpoint tables **65,122 → 64,018 rows**, matching the figure measured in isolation beforehand
- GSK3B −381, A2A −364, HT2A −153, SERT −65, D2 −63, MAO_B −49, BACE1 −12, AChE −7, BChE −6,
  MAO_A −4; hERG and BBB unchanged
- No label moved. These are compounds whose only evidence was a bound.

---

## Step 2 — 13-endpoint random forest panel · IN PROGRESS

Retraining with deduplication before the split and before the deployed refit, the scaffold computed
on the desalted parent, and acyclic compounds sharing one group.

---

## BLOCKED — the measured negative class (BS-C-15) cannot be fetched here

The code fix is committed and correct. The **data** step cannot run in this environment.

Every request to `https://www.ebi.ac.uk` fails TLS verification. This is not a missing certificate
bundle: the handshake completes and a peer certificate is returned, but the chain does not verify.
Building a bundle from certifi plus all 53 certificates in the Windows ROOT and CA stores moves the
error from

```
unable to get local issuer certificate
```

to

```
Basic Constraints of CA cert not marked critical
```

which means the intercepting CA **is** present in the machine's trust store and OpenSSL is rejecting
it as malformed. The network between this machine and EBI is intercepting TLS with a certificate that
does not satisfy RFC 5280.

**I did not work around this.** Passing `verify=False`, or relaxing the OpenSSL verification flags,
would silently accept whatever the interceptor returned, and this pipeline would then be assigning
labels from data of unverified origin. That is the same defect recorded as **BS-M-25** in the audit,
where `app.py:2748` retries with verification disabled and can render a complete pharmacological
report for a molecule an attacker chose. Disabling verification to collect training labels is worse,
not better, than doing it to resolve a name.

**What is needed:** run the fetch from a network without TLS interception, or install a correctly
formed corporate root, then

```bash
python src/brainsafe/data/fetch_endpoints.py        # writes data/_chembl_cache/<target>_inactive.json
python src/brainsafe/data/rebuild_endpoints.py      # pools them in as label 0
```

Until then the panel keeps the class balance it has: 39 of 60 endpoints above 90 per cent active,
five with fewer than 25 inactives, GABAA_a5 with four. Everything else in this pass is unaffected,
because the negative-class recovery only ever adds compounds that are currently absent.

---

## Deliberately not regenerated

**B3DB / `data/endpoints/BBB.csv`.** `fetch_endpoints.main()` would re-download B3DB and rewrite it.
It was not run. The purpose of this pass is to apply the fixes to the data already held; pulling newer
upstream data at the same time would make every BBB comparison ambiguous, since a movement could be
the fix or could be a different B3DB. The shipped `BBB.csv` (7,807 rows) is carried through unchanged,
so the BBB numbers before and after are about the same compounds.

---

## Steps 2–6 · DONE

### The 13-endpoint panel

| | Before | After |
|---|---|---|
| Classifier AUROC, random 10-fold | 0.947–0.969, mean **0.9604** | 0.899–0.977, mean **0.9532** |
| Classifier AUROC, scaffold 10-fold | 0.868–0.956, mean **0.9186** | 0.864–0.965, mean **0.9122** |
| Regression R², random | 0.6381 | **0.6403** |
| Regression R², scaffold | 0.4743 | **0.4763** |

**BBB 0.9605 → 0.8990** random, **0.9197 → 0.8777** scaffold. Excluding BBB nothing moves by more
than 0.011; BACE1 rises in both splits, SERT's scaffold R² by 0.029. The README ranges survive for
the target panel and do not survive for BBB.

### The binder panel, all 49 endpoints

| | Before | After |
|---|---|---|
| Mean sensitivity | **0.8942** | **0.7513** |
| Median sensitivity | 0.9280 | 0.8210 |
| Mean AUROC vs measured inactives | **0.9563** | **0.9003** |
| Background FPR, median | 0.0273 | 0.0297 |
| Background FPR, **maximum** | **0.0500** | **0.0590** |
| Endpoints failing the reliability gate | 2 | **6** |

The maximum is the point. Before, the background false-positive rate could not exceed 0.0500, which
is exactly `BACKGROUND_FPR`, because the threshold was a quantile of the sample it was then scored
on. Measured on the disjoint evaluation pool it reaches 0.0590. **A number that cannot exceed its own
target is not a measurement**, and that it now can is the evidence the fix took.

Every sensitivity records `held_out_actives_by_scaffold`; 45 of 49 thresholds record
`held_out_measured_inactives_and_background`. COX2, GBA1, GluA2, Nav1_6 and TAAR1 join Nav1_1 below
the gate. Nine endpoints that had no recorded training command (BS-M-15) were regenerated with one.

### External validation

| Subset | n | AUROC |
|---|---|---|
| Not in B3DB by InChIKey | 306 | 0.7608 |
| **Also novel in feature space** | **241** | **0.7882** |
| Of which memorised despite the key | 65 | 0.7102 |

### Calibration, conformal, temporal, ADME

Mean ECE 0.0127 → 0.0160 after calibration. **BBB's ECE now improves (0.0657 → 0.0412)** where
calibration previously made it worse, which is what removing the ties its duplicates created was
expected to do; its Brier score still worsens slightly, so BS-M-05 is reduced rather than closed.
Conformal coverage 0.9026 → 0.8974 against a 0.90 target. Temporal essentially unchanged, classifier
mean 0.752 → 0.753, as expected since that split was never the contaminated one.

### Inversion suite: 5 of 6 pass

Reproducibility is exact (0.864 vs 0.864). The domain-flag check still fails, and more sharply: with
the external set restricted to compounds the model can distinguish, unseen drugs sit at median
similarity 0.44 to training while non-drug-like chemistry sits at 0.48.

### Two checks repaired rather than models

- The integrity audit's leakage section read `heldout_actives.json`, which goes stale the moment the
  panel is retrained; it reported three targets with shared scaffolds that were an artefact of
  comparing a current panel with a previous run's list. Pointed at the per-endpoint records, it
  checks 49 targets and finds none.
- The inversion reproducibility check retrained without deduplication, so it was reproducing a
  different experiment.

---

## Still to run

- `noncns_specificity_fast.py` — in progress. This is the **BS-C-05** headline: the 0.875 specificity
  was measured on the applicability-domain reference itself, where every compound sits at Tanimoto
  1.000 to the reference. It now draws from 8,418 DrugBank structures absent from that reference.
- `scaffold_holdout_panel.py` — regenerates `heldout_actives.json` and `scaffold_holdout_results.csv`,
  which feed manuscript Table 6 and `MASTER_validation_summary.csv`.
- The `inversion/` H1–H8 suite, including the BS-C-08 restriction.
- Manuscript tables and figures, then `app_health.py` as the release gate.

---

## Step 7 — specificity on independent chemistry · DONE

**BS-C-05.** The 0.875 figure was measured on a sample drawn from `ad_reference.pkl`, the library the
applicability domain measures against, so all 1,000 compounds sat at Tanimoto 1.000 to it and the
stratification meant to qualify the number was empty.

Drawn instead from the 8,418 DrugBank structures absent from that reference (8,207 eligible):

| | Before | After |
|---|---|---|
| False-positive rate | 0.125 | **0.0800** (95% CI 0.0647–0.0985) |
| Specificity | 0.875 | **0.9200** (95% CI 0.9015–0.9353) |

The stratification now has content, because the sample spans the domain instead of sitting on it
(similarity 0.059–1.000, median 0.514, only 8.8% at 1.000):

| Stratum | k / n | rate |
|---|---|---|
| in domain (T≥0.5) | 52 / 543 | 0.0958 |
| near domain (0.3–0.5) | 25 / 398 | 0.0628 |
| out of domain (T<0.3) | 3 / 59 | 0.0508 |

The headline improves. But the rate **falls** as chemistry gets more distant, which is the reverse of
what the applicability-domain argument predicts.

---

## The applicability domain: three independent measurements now disagree with the claim

1. Out-of-domain AUROC **0.816** exceeds in-domain **0.772** on the external set (BS-M-10).
2. The inversion domain check fails: non-drug-like chemistry sits **nearer** training (median 0.48)
   than genuinely novel drugs (0.44), and only 22% of it falls below the 0.30 cut.
3. Specificity is **better** out of domain (0.0508) than in it (0.0958).

Whatever `AD_THRESHOLD = 0.30` is selecting, it is not distance from the training chemistry, and it
qualifies every prediction the server returns. **This needs its own investigation before submission.**

---

## Step 8 — health gate · ONE FAILURE

`app_health.py` exits 1 on a single check:

```
[FAIL] directional pharmacology: levodopa top-3 does not include Parkinson's disease:
       ['Neuroprotection / oxidative stress', 'Epilepsy', "Alzheimer's disease"]
```

Before the regeneration levodopa read `Neuroprotection 0.34, Alzheimer's 0.14, Parkinson's 0.14` —
Parkinson's was third **on a tie at 0.14**. It is now displaced by Epilepsy at 0.17. The margin was
always inside noise, and levodopa is a prodrug rather than a binder, so a target-engagement panel was
never going to score it well. That explains the fragility; it does not excuse the failure. **The gate
fails and the release is not clean.** Either the check is testing something this panel cannot deliver
and should say so, or the disease layer needs work for prodrugs. That is a judgement for the authors.

Everything else passes: 52 targets, all probes score, distinct profiles, peripheral controls below
threshold, exports well formed, ADME units consistent.

---

## Step 9 — model manifest · DONE, deposit outstanding

Regeneration changed 93 of 195 files, so `model_fetch` judged them wrong and tried to re-download.
On this network that failed; on a working one it would have **overwritten the regenerated models with
the previous ones**. The manifest is rebuilt (244 files), version bumped to **1.1**, and `doi`/`urls`
emptied with the old DOI kept under `supersedes_doi`. A manifest pointing nowhere fails honestly; one
pointing at a record serving different bytes does not.

---

## Remaining

- `scaffold_holdout_panel.py` — Table 6 and `MASTER_validation_summary.csv`
- the `inversion/` H1–H8 suite, with the BS-C-08 restriction
- manuscript tables and figures
- deposit models v1.1 and the source-cache archive, then fill both manifests
- the blocked ChEMBL negative-class fetch
