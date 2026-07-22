# Ladder A1: adding measured inactives, a tested and (for the naive form) rejected approach

Date: 2026-07-21. This records an experiment that did **not** make it into the primary model, and why.
Keeping it on the record is deliberate: the negative result is as informative as a positive one.

## Motivation

Several classifiers are strongly active-skewed in the public data (GSK-3-beta ~93% active), because
ChEMBL and BindingDB report inhibitors. The practical symptom is an unrealistic base rate: the
GSK-3-beta model predicted **71.6%** of the 11,723 DrugBank drugs as active, which is nonsensical , 
most drugs are not kinase inhibitors. The proposed fix was to add genuine measured negatives from
PubChem high-throughput screens.

## What was done

`fetch_pubchem_inactives.py` pulled measured inactives for the protein-target classifiers, excluding
any compound ever active for the target. Large screens were captured for GSK-3-beta (11,966
inactives) and hERG (11,953; not used, hERG is already inactive-majority). `rebuild_with_inactives.py`
added 4,276 of the GSK-3-beta inactives (capped to balance to 50% active, never overriding a
dose-response active). The model was retrained and audited.

## Result

| GSK-3-beta | Before | After adding inactives |
|---|---|---|
| Active fraction | 93% | 50% |
| DrugBank predicted-active | 71.6% | **16.4%** |
| Scaffold AUROC | 0.937 | **0.989** |
| Scaffold MCC | 0.559 | 0.901 |
| Median Tanimoto of added inactives to actives | n/a | **0.288** (1.7% >= 0.5) |

Two things happened at once:
1. **The base rate was genuinely corrected** (71.6% -> 16.4% predicted active on DrugBank). This part
   is real and valuable.
2. **The discrimination metric was inflated** (AUROC 0.937 -> 0.989). The audit
   (`results/tables/inactives_audit.csv`) shows the added inactives are chemically **unlike** the
   inhibitors (median similarity 0.29). The model is separating "looks like a kinase inhibitor" from
   "random screening compound", which is trivial. This is the classic decoy / artificial-enrichment
   bias, not a real improvement in telling active inhibitors from *similar* inactive ones.

## Decision

The naive augmentation was **reverted from the primary model**. Reporting AUROC 0.989 would overstate
the model; the honest discrimination estimate remains **0.937** (scaffold), measured against
dose-response inactives that are chemically similar to the actives. The class imbalance is recorded as
a genuine limitation of the public data, and the base-rate issue is handled at prediction time by the
calibration layer and the applicability-domain flag rather than by inflating the training set.

## The correct way to do this (future work)

The negatives that fix the base rate (bulk HTS randoms) are the easy ones; genuine drug-like measured
"hard inactives" near the decision boundary are scarce. A defensible augmentation would keep only
**property- or similarity-matched hard negatives** (for example measured inactives with maximum
Tanimoto to an active above ~0.4), which do not inflate discrimination, but these are few, so they
correct the metric honestly without fully fixing the base rate. The base rate itself is best handled by
modelling prevalence explicitly and always reporting the calibrated probability with its domain flag.

Evidence: `results/tables/inactives_audit.csv`, `inactives_merge_provenance.csv`,
`data/_pubchem_cache/`.
