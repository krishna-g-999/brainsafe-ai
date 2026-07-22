# Ladder B: measured ADME / exposure layer — results

Ladder B adds the properties that decide whether an achievable dose puts *free* drug on a brain
target, i.e. what connects target binding to a real brain effect. Six measured endpoints, trained with
the same random-forest + 10-fold protocol as the target panel. Date: 2026-07-22.

## Data

18,699 measured compounds from established public benchmarks (TDC / Harvard Dataverse and
MoleculeNet), standardised identically to the rest of the project (`data/adme/`, `data/adme/SOURCE.md`).

## Cross-validation (10-fold)

| Endpoint | Task | n | Random | Scaffold | Confidence |
|---|---|---|---|---|---|
| P-gp inhibition | classification | 1,212 | AUROC 0.955 | **AUROC 0.937** | high |
| Aqueous solubility | regression | 9,573 | R2 0.804 | **R2 0.763** | high |
| Caco-2 permeability | regression | 897 | R2 0.734 | R2 0.593 | moderate |
| Lipophilicity (logD) | regression | 4,200 | R2 0.639 | R2 0.564 | moderate |
| Plasma-protein binding | regression | 1,797 | R2 0.434 | R2 0.374 | low |
| Hepatocyte clearance | regression | 1,020 | R2 0.230 | R2 0.193 | low |

Full detail: `results/tables/adme_cv_summary.csv`, `adme_cv_folds.csv`; models in `models_rf/adme/`;
fold predictions in `data/processed/cv_predictions/adme/`.

**Honest reading.** P-gp inhibition and solubility are strong; permeability and logD are usable;
plasma-protein binding is modest and hepatocyte clearance is weak (R2 ~0.19). The clearance result is
expected — structure-only metabolic-clearance prediction is a known hard problem — and it is reported
as low-confidence rather than presented as equal to the others.

## The combined CNS free-exposure readout (a K_p,uu proxy)

`cns_exposure.py` chains BBB penetration, passive permeability (Caco-2), free fraction (plasma-protein
binding) and a P-gp flag into a qualitative free-brain-exposure call. It is a heuristic proxy for
K_p,uu (unbound brain-to-plasma), not a measured K_p,uu — no large public K_p,uu dataset exists to
train on directly. On known drugs (`results/tables/cns_exposure_demo.csv`):

| Compound | BBB | Caco-2 | free frac. | Call | Correct? |
|---|---|---|---|---|---|
| Diazepam (CNS) | 0.98 | -4.40 | 0.106 | favourable | yes |
| Donepezil (CNS) | 0.96 | -4.96 | 0.386 | favourable | yes |
| Caffeine (CNS) | 0.84 | -4.48 | 0.731 | favourable | yes |
| Atenolol (peripheral) | 0.37 | -5.35 | 0.649 | limited | yes |
| Loperamide (peripheral) | 0.87 | -5.02 | 0.061 | favourable | **no** |

The readout correctly separates central drugs from a peripheral low-BBB drug (atenolol). It **fails on
loperamide**, which is not centrally active in reality because P-glycoprotein actively pumps it out of
the brain. Our P-gp model is trained on **inhibition**, not **substrate** status, so it cannot capture
this efflux. This single known-drug failure is the clearest possible motivation for the next endpoint:

## Next step (B, continued)

Add a measured **P-gp substrate** dataset (distinct from inhibition) so the efflux liability that keeps
loperamide out of the brain is modelled. That closes the most important remaining gap between "crosses
the barrier" and "achieves free brain exposure". Protein binding and clearance would also benefit from
larger measured sets. K_p,uu itself remains a proxy until a dedicated measured dataset is sourced.
