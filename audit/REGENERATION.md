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
