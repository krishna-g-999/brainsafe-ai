# Ladder A2: is there more clean measured data? An empirical exhaustion check

The learning curves showed BACE1, MAO-A and the receptor regressors are still improving with data, so
A2 asked whether more clean measured data exists. Four independent public sources were checked
directly (not assumed). Date: 2026-07-21.

| Source | Quality | Finding |
|---|---|---|
| ChEMBL 37 | dose-response | Live target counts confirmed ~97% already captured; a few hundred more compounds at most. |
| BindingDB | dose-response | Already pooled in (+4,246 net-new; `results/tables/bindingdb_yield.csv`). Exhausted for the useful range. |
| PubChem BioAssay HTS | single-concentration | Inactives are easy negatives that inflate metrics (A1, `docs/INACTIVES_EXPERIMENT.md`); primary-screen actives carry high false-positive rates. Not clean. |
| IUPHAR / Guide to Pharmacology | expert-curated | Only 10-185 interactions per target with affinity (~500 total across six targets), largely overlapping ChEMBL. Negligible net-new. |

## Conclusion

The clean, dose-response public measured data for these well-studied CNS targets is **exhausted**.
ChEMBL and BindingDB together already hold essentially all of it, and the remaining sources are either
too noisy (PubChem HTS) or too small and redundant (GtoPdb) to help. The learning-curve headroom for
BACE1, MAO-A and the receptors is therefore **data-limited**: realising it requires **new experimental
measurements**, not further data mining.

This is a genuine, defensible result rather than a shortfall: we have demonstrably used the available
measured evidence, and we can say precisely why more of the same is not obtainable. The productive next
step is not more target-activity data (ladder A) but a new kind of measured endpoint that speaks to how
much drug actually reaches the brain (ladder B: ADME / exposure).
