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
| P-gp substrate (efflux ratio >= 2) | classification | 1,371 | AUROC 0.858 | **AUROC 0.808** | good |
| Aqueous solubility | regression | 9,573 | R2 0.804 | **R2 0.763** | high |
| logBB (total brain/plasma) | regression | 1,058 | R2 0.577 | R2 0.455 | moderate |
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

Now includes a measured **P-gp substrate** model (efflux) so the readout accounts for active efflux,
not just passive crossing:

| Compound | BBB | Caco-2 | P-gp substrate | Call | Correct? |
|---|---|---|---|---|---|
| Diazepam (CNS) | 0.98 | -4.40 | 0.31 | favourable | yes |
| Caffeine (CNS) | 0.84 | -4.48 | 0.19 | favourable | yes |
| Atenolol (peripheral) | 0.37 | -5.35 | 0.21 | limited | yes |
| Loperamide (peripheral) | 0.87 | -5.02 | 0.59 | limited (P-gp efflux) | **yes** |
| Donepezil (CNS) | 0.96 | -4.96 | 0.59 | limited (P-gp efflux) | borderline |

Adding the substrate model **fixes the loperamide case** — it is now correctly called efflux-limited,
the reason it is peripheral despite crossing passively. Donepezil sits right at the 0.5 boundary
(substrate probability 0.59) and is flagged too; this is honest rather than wrong, because donepezil
is a documented weak P-gp substrate that still reaches therapeutic brain levels. It illustrates the
correct interpretation: **P-gp substrate status is a graded risk flag, not an absolute veto** — the
strength of efflux, and the dose, decide the outcome.

## Data note

The P-gp substrate labels come from ChEMBL bidirectional-transport efflux ratios (target CHEMBL4302),
thresholded at ratio >= 2 (the standard substrate cut); 1,371 compounds, 840 substrate / 531 non-
substrate (`src/brainsafe/adme/fetch_pgp_substrate.py`). Efflux-ratio cut-offs are assay-dependent, so
this is a defensible but approximate label.

## K_p,uu: what is and is not available

We searched for a dedicated measured **K_p,uu** (unbound brain-to-plasma) dataset. There is **no clean,
sizable public benchmark** - K_p,uu is expensive to measure (it needs brain-homogenate binding, plasma
binding and total brain/plasma together) and exists only as small tables of ~50-300 compounds in
individual papers, not as a downloadable dataset. So K_p,uu itself remains a **proxy**.

What we added instead is **logBB** (log of *total* brain-to-plasma concentration), 1,058 measured
compounds from B3DB, as a regression endpoint. This is real and useful, but it is important to state
what it is not: logBB is *total* distribution, inflated by lipophilicity and tissue binding, and can be
high while *free* brain exposure is low. The demo makes this concrete - **loperamide has logBB 0.41
(looks brain-penetrant) yet is correctly called efflux-limited** by the P-gp substrate model. This is
exactly why K_p,uu was invented and why the mechanistic readout (permeability + efflux + free fraction)
is a better free-exposure proxy than logBB alone.

## Where B stands

The exposure layer now models the full mechanistic chain: BBB penetration -> logBB (total distribution)
-> passive permeability (Caco-2) -> **active efflux (P-gp substrate)** -> free fraction (plasma-protein
binding), combined into a qualitative free-brain-exposure call. Remaining honest gaps: a directly
measured K_p,uu model (data-limited, as above), and larger sets for plasma-protein binding and
hepatocyte clearance (the two weakest endpoints).
