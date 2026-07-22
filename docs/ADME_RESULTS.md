# Ladder B: measured ADME / exposure layer: results

Ladder B adds the properties that decide whether an achievable dose puts *free* drug on a brain
target, i.e. what connects target binding to a real brain effect. Nine measured endpoints, trained with
the same random-forest + 10-fold protocol as the target panel. Date: 2026-07-22.

## Data

About 21,700 measured compounds from established public sources (TDC / Harvard Dataverse, MoleculeNet,
B3DB and ChEMBL), standardised identically to the rest of the project (`data/adme/`,
`data/adme/SOURCE.md`).

## Cross-validation (10-fold)

| Endpoint | Task | n | Random | Scaffold | Confidence |
|---|---|---|---|---|---|
| P-gp inhibition | classification | 1,212 | AUROC 0.955 | **AUROC 0.937** | high |
| P-gp substrate (efflux ratio >= 2) | classification | 1,371 | AUROC 0.858 | **AUROC 0.808** | good |
| Aqueous solubility | regression | 9,573 | R2 0.804 | **R2 0.763** | high |
| **K_p,uu (unbound brain/plasma)** | regression (log10) | 566 | R2 0.404 | **R2 0.352** | moderate (small n) |
| logBB (total brain/plasma) | regression | 1,058 | R2 0.577 | R2 0.455 | moderate |
| Caco-2 permeability | regression | 897 | R2 0.734 | R2 0.593 | moderate |
| Lipophilicity (logD) | regression | 4,200 | R2 0.639 | R2 0.564 | moderate |
| Plasma-protein binding | regression | 1,797 | R2 0.434 | R2 0.374 | low |
| Hepatocyte clearance | regression | 1,020 | R2 0.230 | R2 0.193 | low |

Full detail: `results/tables/adme_cv_summary.csv`, `adme_cv_folds.csv`; models in `models_rf/adme/`;
fold predictions in `data/processed/cv_predictions/adme/`.

**Honest reading.** P-gp inhibition and solubility are strong; permeability and logD are usable;
plasma-protein binding is modest and hepatocyte clearance is weak (R2 ~0.19). The clearance result is
expected, structure-only metabolic-clearance prediction is a known hard problem, and it is reported
as low-confidence rather than presented as equal to the others.

## The combined CNS free-exposure readout

`cns_exposure.py` combines the measured models into a single free-brain-exposure call. The primary
signal is the directly modelled **K_p,uu** (described in the next section); BBB penetration, logBB,
passive permeability (Caco-2), the P-gp substrate (efflux) model and free fraction (plasma-protein
binding) are reported alongside as interpretable supporting predictions. The P-gp substrate model is
what lets the readout distinguish a compound that crosses passively but is actively pumped back out.

## Data note

The P-gp substrate labels come from ChEMBL bidirectional-transport efflux ratios (target CHEMBL4302),
thresholded at ratio >= 2 (the standard substrate cut); 1,371 compounds, 840 substrate / 531 non-
substrate (`src/brainsafe/adme/fetch_pgp_substrate.py`). Efflux-ratio cut-offs are assay-dependent, so
this is a defensible but approximate label.

## K_p,uu: a directly measured model (not a proxy)

A dedicated measured **K_p,uu** (unbound brain-to-plasma) dataset was found in ChEMBL: activities with
the standard type `K(p,uu,brain)` in "unbound brain" assays. After standardisation and per-compound
median aggregation this gives **566 unique compounds with measured K_p,uu** (range 0.005-47, median
0.40); the values sanity-check against known pharmacology (erlotinib 0.05, a P-gp substrate with low
brain exposure; diazepam ~1). This is small but real, so it is modelled directly as a regression on
log10(K_p,uu) (`src/brainsafe/adme/fetch_pgp_substrate.py`-style fetch; `data/adme/kpuu.csv`). Scaffold
R2 is **0.352** - honest for a small, integrative endpoint, in line with published K_p,uu models.

The directly-modelled K_p,uu now drives the exposure call (K_p,uu >= 0.3 = meaningful free brain
exposure, the standard cut). On known drugs (`results/tables/cns_exposure_demo.csv`):

| Compound | predicted K_p,uu | Call | Correct? |
|---|---|---|---|
| Diazepam | 0.94 | favourable | yes (K_p,uu ~ 1) |
| Donepezil | 0.84 | favourable | yes (works in AD) |
| Caffeine | 0.16 | borderline | yes (moderate) |
| Atenolol | 0.07 | limited | yes (peripheral) |
| Loperamide | 0.04 | limited | yes (effluxed, non-central) |

The direct model resolves **both** hard cases the earlier heuristic could not: loperamide is correctly
limited (efflux) and donepezil is correctly favourable (works despite being a weak P-gp substrate),
because it is trained on the actual integrated endpoint rather than assembled from components. logBB and
the P-gp substrate model are retained as interpretable supporting predictions.

## Where B stands

The exposure layer now includes a **directly measured K_p,uu model** as the primary free-brain-exposure
readout, supported by BBB penetration, logBB, passive permeability (Caco-2), active efflux (P-gp
substrate) and free fraction (plasma-protein binding). Remaining honest gaps: K_p,uu would benefit from
more than 566 compounds (it is the smallest endpoint), and plasma-protein binding and hepatocyte
clearance remain the weakest regressors.
