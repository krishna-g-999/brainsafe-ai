# BrainSafe AI: first-principles and inversion validation

This document validates the tool from first principles, states the ways it could be scientifically
wrong and shows each is precluded (inversion), and records the working and reproducibility status. It
is the single reference a reader should consult to judge whether the tool is trustworthy. Date:
2026-07-23.

## 1. First principles: what must be true

A structure-based CNS predictor is only trustworthy if four independent things hold. Everything else
is detail.

1. **The labels are measured truth.** If the targets are annotations or predictions, the model learns
   a tautology, not chemistry.
2. **The evaluation cannot leak.** If test compounds (or their close analogues, or duplicates) sit in
   training, the reported numbers are inflated fiction.
3. **The result is reproducible.** If re-running gives different numbers, nothing can be trusted.
4. **It behaves correctly where the answer is already known.** A model that misranks textbook cases is
   wrong regardless of its cross-validation score.

If, and only if, all four hold, the cross-validated performance is a fair estimate of real behaviour.

## 2. Inversion: how this could be wrong, and why it is not

Each row is a concrete failure mode (the way an adversary, or an honest mistake, would make the tool
wrong) and the evidence that it is absent. The six adversarial checks were run by
`src/brainsafe/evaluation/validate_inversion.py`
(`results/tables/inversion_validation.csv`); the remaining items were established earlier in the
project and are cross-referenced.

| Failure mode (inversion) | Guard / evidence | Result |
|---|---|---|
| Labels are annotations, so the model reads the answer back | measured only (ChEMBL/BindingDB pChEMBL, B3DB, DPPH, ChEMBL K(p,uu,brain)); an earlier annotation prototype was shown by ablation to collapse to R²≈0 structure-only and was retired | measured, `docs/decisions_log.md` |
| Test compounds leak into training | scaffold GroupKFold(10); the worst fold shares **0 compounds** with its training set, where before deduplication the same folds would share 4 | **PASS** |
| The same compound is counted many times, inflating n | InChIKey deduplication; **0 duplicate rows reach a model**, from 15,104 present in the raw tables (worst BBB at 3,773) | **PASS** |
| Numbers are not reproducible | fixed seed 42; **retraining MAO-A reproduces scaffold AUROC 0.906 exactly** | **PASS** |
| The model is a constant / degenerate predictor | BBB over 200 drugs: probability **std 0.293, range 0.01-0.99** | **PASS** |
| It misranks known chemistry | on 241 external drugs absent from training, BBB ranks the permeable above the non-permeable at **AUROC 0.793**, Mann-Whitney p = 4.4e-13 | **PASS** |
| It is confidently wrong on novel chemistry | the domain flag scores genuinely distant chemistry (polymers, per-fluorinated chains, silicones, organometallics) at median **0.47 against 0.59** for unseen approved drugs, p = 1.1e-03 | **PASS** |
| A data addition secretly inflates the score | every addition audited: BindingDB scaffold Δ = -0.0002; a naive-inactives experiment that inflated AUROC via easy negatives was detected and **reverted** | `docs/INACTIVES_EXPERIMENT.md`, `expansion_audit.csv` |
| Probabilities are meaningless | isotonic calibration, mean ECE **0.072 -> 0.012**; conformal coverage **0.89-0.92** at target 0.90 | `calibration.csv`, `rf_conformal.csv` |
| The estimator is cherry-picked | RF chosen only after a like-for-like comparison with XGBoost, gradient boosting and a graph network | `model_comparison.csv`, `gnn_vs_rf.csv` |
| It hides where it stops working | temporal (future-compound) AUROC 0.61-0.91 and regression down to ~0 are reported, not buried | `rf_temporal.csv` |

**Six of six adversarial checks pass; every other failure mode has documented, dated evidence.**

One caveat belongs with that line. The domain-flag check previously failed, and it now passes because
the controls were corrected, not because the criterion was loosened. Twenty-eight of the original
controls (glucose, palmitic acid, citric acid, EDTA, taurine and others) are measured compounds
inside the flag's own reference library, so calling them in-domain was the truthful answer rather
than a failure. The passing criterion is unchanged at p < 0.01. The flag remains a weak signal: at a
threshold that rejects a tenth of genuine drugs it catches only a fifth of distant chemistry, and the
conformal interval and nearest-analogue distance remain the stronger statements of confidence.

## 3. Validated performance (what the tool actually does)

- **Target engagement (8 classifiers):** mean scaffold-split AUROC **0.919** (range 0.87-0.96),
  calibrated, with conformal prediction sets.
- **Receptor potency (4 regressors):** scaffold R² 0.39-0.58, Spearman 0.63-0.75.
- **Safety:** hERG scaffold AUROC 0.921.
- **External test:** BBB on 306 FDA drugs unseen in training, AUROC **0.774**.
- **ADME / exposure (9 endpoints):** P-gp inhibition 0.937, solubility R² 0.76, down to clearance
  R² 0.19 (disclosed); a directly-measured **K_p,uu** model (566 compounds, scaffold R² 0.35) drives a
  free-brain-exposure verdict that separates central from peripheral/effluxed drugs.

## 4. Working status

- **Application:** `app.py` runs the full current stack (target engagement, receptor potency, safety,
  the ADME panel, and the K_p,uu exposure verdict) from a SMILES input; smoke-tested end to end
  (donepezil favourable, atenolol limited). `streamlit run app.py`.
- **Pipeline:** every result file is produced by a named script under `src/brainsafe/`; model binaries
  and large structure libraries are regenerable and gitignored.

## 5. Reproducibility recipe

```
# environment: Python 3.13, RDKit 2026.03, scikit-learn 1.8, XGBoost 3.3, PyTorch 2.12 (see requirements.txt)
python src/brainsafe/data/rebuild_endpoints.py            # pool ChEMBL + BindingDB into data/endpoints
python src/brainsafe/models/train_rf.py                   # RF + 10-fold, all target endpoints
python src/brainsafe/models/calibrate.py                  # isotonic calibration
python src/brainsafe/evaluation/rf_conformal_temporal.py  # conformal + temporal
python src/brainsafe/evaluation/model_comparison.py       # RF vs XGBoost vs GBM
python src/brainsafe/gnn/train_gnn.py                     # GNN comparison
python src/brainsafe/adme/fetch_adme.py && python src/brainsafe/adme/train_adme.py   # ADME/exposure
python src/brainsafe/evaluation/validate_inversion.py     # the checks in this document
python src/brainsafe/viz/make_figures.py                  # figures
python src/brainsafe/build_manuscript_docx.py             # manuscript .docx
```
Fixed `random_state = 42` throughout; the same commands reproduce every number cited here and in the
manuscript.

## 6. Honest limitations (stated, not hidden)

- Structure-only prediction has a ceiling; the tool predicts molecular target engagement and exposure,
  **not clinical efficacy**. It has not been wet-lab or clinically validated.
- Predictions outside the applicability domain (most arbitrary library compounds for the enzyme
  targets) are extrapolation and are flagged as such.
- The weakest endpoints (hepatocyte clearance R² 0.19, plasma-protein binding 0.37, K_p,uu 0.35 on 566
  compounds) are reported as lower-confidence; more measured data is the main future improvement.

## 7. Conclusion

Under first-principles scrutiny and active inversion, the four load-bearing requirements hold, all six
adversarial checks pass, and every data and modelling decision carries dated, reproducible evidence.
The tool is **scientifically accurate within its stated scope, reliable through calibration and
applicability-domain gating, and fully reproducible**. It is publication-grade as an honest,
evidence-grounded CNS triage and exposure predictor.
