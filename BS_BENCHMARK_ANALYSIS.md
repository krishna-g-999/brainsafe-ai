# BrainSafe AI — Final Honest Benchmark & Competitive Analysis

*Fact-based. Our numbers are read from the saved reports; external numbers are from
the cited literature. No assumptions — where split methodology differs, it is stated.*

## 1. What exists (verified from files)
- 7 measured-data endpoint models + 7 conformal calibrations + antioxidant model +
  engine (`BS_brain_predict.py`) + app + model card.
- **33,146** measured compound–endpoint datapoints (ChEMBL pChEMBL + B3DB).
- Validation artifacts: `endpoints_report.json`, `BS_external_validation_report.json`,
  `BS_predictive_report.json`, `BS_validation_report.json`.

## 2. Per-endpoint performance vs published literature

| Endpoint | BrainSafe (scaffold-CV / **strict cluster-split**) | Published range | Honest position |
|---|---|---|---|
| BBB (B3DB) | 0.921 / **0.906** AUROC | 0.88–0.96 on B3DB (mostly **random** splits) | Competitive; ours is under a *stricter* split, so not directly comparable — random-split numbers inflate. |
| hERG | 0.907 / **0.900** | 0.86–0.93 (SVC ~0.93, DL 0.88–0.91) | Top of the typical range. |
| AChE | 0.915 / **0.912** | comparable ChEMBL RF/SVM studies (n~2.5k) | Strong; larger dataset (n=4,324), stricter validation. |
| BACE1 | 0.950 / **0.918** | multi-target Alzheimer's RF models exist | Strong; class-imbalanced (89% active) — AUROC/MCC reported, not accuracy. |
| MAO-B | 0.885 / **0.873** | 0.88–0.96 (small pharmacophore sets, n~126) | Solid on a far larger set (n=3,455); small-set literature numbers are optimistic. |
| MAO-A | 0.867 / **0.890** | — | Solid. |
| GSK-3β | 0.920 / **0.915** | — | Strong; imbalanced (MCC 0.47). |
| antioxidant | R²=0.27 (scaffold-CV) | — | **Weak** — honest; coarse curated labels cap it. |
| conformal coverage | **0.885–0.905 @ 0.90 target** | valid CP gives ~target coverage | Empirically valid. |

**Honest takeaway:** per endpoint we are **competitive, not state-of-the-art-beating**.
Our distinguishing strength is **validation rigor** (leak-free scaffold CV + strict
leave-cluster-out + verified conformal coverage), where many published tools report
random-split metrics that overstate real generalisation.

## 3. Comparison with existing online tools (facts)

| Capability | SwissADME | ADMETlab 3.0 | admetSAR/pkCSM | SwissTargetPrediction / PPB2 | **BrainSafe** |
|---|---|---|---|---|---|
| BBB penetration | ✔ (BOILED-EGG rule) | ✔ | ✔ | ✘ | ✔ (ML, B3DB) |
| hERG / safety | ✘ | ✔ | ✔ | ✘ | ✔ (hERG) |
| Druggability / CNS-MPO | ✔ | partial | partial | ✘ | ✔ |
| **CNS target activity (AChE/BACE1/MAO/GSK-3β)** | ✘ | ✘ | ✘ | ✔ (generic targets, similarity) | ✔ (measured QSAR) |
| **BBB-gated → disease-level brain effect** | ✘ | ✘ | ✘ | ✘ | ✔ (novel integration) |
| **Nearest measured-analog evidence per call** | ✘ | ✘ | ✘ | ✔ (similarity-based) | ✔ |
| **Conformal calibrated confidence** | ✘ | uncertainty est. | ✘ | ✘ | ✔ (coverage-verified) |
| Open / transparent / local | ✔/✘ | web | web | web | ✔ (fully inspectable) |

**What no single existing tool does (our integrative contribution):** unify *measured*
CNS-target polypharmacology **gated by BBB penetration** into **disease-level** brain-effect
scores, **plus** a safety anti-target, **plus** conformal confidence, **plus**
measured-analog evidence — in one transparent tool.

## 4. Honest novelty assessment
- **Not novel:** the individual methods. ECFP/descriptor + RF/ExtraTrees/HistGB ensembles,
  BBB ML, hERG ML, target QSAR, applicability domain, and conformal prediction
  (Norinder 2014) are all **standard, established** techniques. Multi-target Alzheimer's
  QSAR (AChE+MAO-B+BACE1) has been published.
- **Genuinely novel (incremental, integrative):** the *combination* — a BBB-gated,
  evidence-grounded, conformal-calibrated, safety-aware **CNS disease-effect profiler**
  driven entirely by measured public data, as an open transparent tool. This is an
  **application/integration contribution, not a methodological breakthrough.**

## 5. Verdicts (honest)
- **Scientifically correct?** **Yes.** Leak-free scaffold CV, strict leave-cluster-out
  (AUROC barely drops → real generalisation), isotonic calibration (Brier 0.04–0.14),
  conformal coverage empirically ≈ target, competitive AUROCs, documented limitations.
- **Approvable / publishable?** **Yes, as an application/resource paper** (e.g.,
  *J. Cheminformatics*, *Scientific Reports*, *Molecular Informatics* tier) **if framed
  honestly** as an integrative open CNS-profiling tool with rigorous validation — **not**
  as a SOTA-beating or novel-method paper. A top-tier predictor paper would additionally
  require temporal + prospective validation.
- **Useful to researchers?** **Yes.** Free, transparent, evidence-backed multi-endpoint
  CNS triage with calibrated confidence + safety + BBB gating is directly useful for
  hypothesis generation, compound prioritisation, and teaching — especially for
  natural products/flavonoids where such integrated CNS read-outs are scattered.

## 6. Honest remaining gaps (what a reviewer would flag)
1. **No temporal/prospective validation** — the strongest test of real-world use; needs a
   time-split ChEMBL pipeline + wet-lab follow-up (cannot be fabricated).
2. **Analog density** — ChEMBL target sets have median scaffold-split test→train Tanimoto
   0.55–0.71; cluster-split mitigates but a reviewer will want this stated (it is).
3. **Target activity ≠ clinical efficacy** — we predict molecular engagement, not outcomes.
4. **Antioxidant model is weak (R²0.27)** and disease coverage is a limited 5-target panel
   (no e.g. NMDA, D2, 5-HT receptors, BChE, Nrf2).
5. **Imbalanced endpoints** (BACE1, GSK-3β) — reported via AUROC/MCC, but PR-curves and
   decision-threshold sensitivity should be shown in a manuscript.

## 7. Bottom line
BrainSafe is a **scientifically valid, rigorously-validated, genuinely useful integrative
CNS profiler** — competitive per-endpoint performance with **above-average validation
honesty**, and an integration (BBB-gated disease effect + safety + conformal + evidence)
not offered as a unit by existing tools. It is **publishable as an application/resource**
with honest framing, and **not** a methodological breakthrough. The path to a flagship
predictor paper is temporal + prospective validation and a broader, balanced target panel.

## Sources
- BBB / B3DB benchmarks: [Sci Rep 2024](https://www.nature.com/articles/s41598-024-66897-y), [Chem Res Toxicol](https://pubs.acs.org/doi/10.1021/acs.chemrestox.0c00343), [DeePred-BBB](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112838/)
- hERG benchmarks: [JCIM benchmark study](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00744), [Front. Pharmacol. 2022](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2022.951083/full)
- ADMET tools: [ADMETlab 3.0 (NAR 2024)](https://academic.oup.com/nar/article/52/W1/W422/7640525), [ADMETlab 2.0 (NAR 2021)](https://academic.oup.com/nar/article/49/W1/W5/6249611), [admetSAR 3.0](https://lmmd.ecust.edu.cn/admetsar3/about/endpoint.php)
- Target prediction: [SwissTargetPrediction (NAR 2019)](https://academic.oup.com/nar/article/47/W1/W357/5491750), [Polypharmacology Browser PPB2](https://www.researchgate.net/publication/329741123_The_Polypharmacology_Browser_PPB2_Target_Prediction_Combining_Nearest_Neighbors_with_Machine_Learning)
- CNS target QSAR: [AChE ML screening (Front. Neurosci. 2022)](https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2022.1007389/full), [Dual AChE/MAO-B ML (Mol. Divers. 2024)](https://link.springer.com/article/10.1007/s11030-024-11061-x)
- Conformal prediction: [Norinder et al. 2014 / CP in drug discovery (J Pharm Sci 2020)](https://jpharmsci.org/article/S0022-3549(20)30589-X/fulltext), [Large-scale QSAR + CP (J Cheminform 2018)](https://link.springer.com/article/10.1186/s13321-018-0325-4)
