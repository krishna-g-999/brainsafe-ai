# BrainSafe AI — Benchmark & Competitive Analysis

*Our numbers are read straight from the saved validation reports; external numbers come from the
cited literature. Where a comparison rests on a different split methodology, we say so rather than
paper over it.*

## 1. What the system contains (verified from files)
Eight measured-data classification endpoints (BBB, AChE, BChE, BACE1, GSK-3β, MAO-A, MAO-B, hERG),
each with its own isotonic calibration and Mondrian conformal layer; four receptor potency
regressions (D2, A2A, 5-HT2A, SERT); a measured-DPPH antioxidant model; a deterministic
druggability/CNS-MPO layer; and a 504-compound clinical-precedent layer. All of it is trained on
**64,474 measured compound–endpoint records** (ChEMBL_37 pChEMBL + B3DB + ChEMBL DPPH). The inference
engine is `BS_brain_predict.py`; validation artifacts are the `*_report.json` files and the
per-endpoint `*_meta.json` files.

## 2. Per-endpoint performance vs the published literature

| Endpoint | BrainSafe (scaffold-CV / **cluster-split**) | Published range | Honest position |
|---|---|---|---|
| BBB (B3DB) | 0.921 / **0.906** AUROC | 0.88–0.96 on B3DB (mostly **random** splits) | Competitive; ours is under a stricter split, so a like-for-like random comparison would read higher. |
| hERG | 0.901 / **0.870** | 0.86–0.93 (SVC ~0.93, DL 0.88–0.91) | Within the typical range, tested more strictly. |
| AChE | 0.915 / **0.912** | comparable ChEMBL RF/SVM studies (n~2.5k) | Strong, on a larger set (n=4,324) and a stricter split. |
| BACE1 | 0.950 / **0.940** | multi-target Alzheimer's RF models exist | Strong; class-imbalanced (~89% active), so we report AUROC/MCC rather than accuracy. |
| MAO-B | 0.885 / **0.873** | 0.88–0.96 (small pharmacophore sets, n~126) | Solid on a much larger set (n=3,455); small-set literature numbers tend to be optimistic. |
| MAO-A | 0.867 / **0.890** | — | Solid. |
| GSK-3β | 0.920 / **0.915** | — | Strong; imbalanced (MCC 0.47). |
| antioxidant (measured DPPH) | R² = 0.43 (scaffold-CV) | — | Honest and moderate; it replaces an earlier text-derived proxy that fit far worse (R² ≈ 0.25). |
| conformal coverage | **0.885–0.905 @ 0.90 target** | valid CP gives ~target coverage | Empirically valid. |

**Takeaway.** On like-for-like random splits our AUROC (0.94–0.98) sits at or above the published
ranges; under the stricter scaffold and cluster splits it holds at 0.87–0.95. What sets the work apart
is less any single headline number than the **validation rigour** — leak-free scaffold CV, strict
leave-cluster-out, and empirically verified conformal coverage — next to a literature that often
reports only the friendlier random split.

## 3. Does the ensemble earn its keep? (comparative)
Under an identical scaffold-split protocol the deployed ensemble reaches a mean AUROC of **0.912**,
against **0.867** for a k-nearest-neighbour Tanimoto read-across and **0.808** for logistic regression,
and it is best on every one of the eight endpoints. Beating a pure nearest-neighbour baseline is the
check that the model has learned structure–activity relationships rather than simply memorising
look-alikes (Supplementary Table S9).

## 4. Comparison with existing online tools

| Capability | SwissADME | ADMETlab 3.0 | admetSAR/pkCSM | SwissTargetPrediction / PPB2 | **BrainSafe** |
|---|---|---|---|---|---|
| BBB penetration | ✔ (BOILED-EGG rule) | ✔ | ✔ | ✘ | ✔ (ML, B3DB) |
| hERG / safety | ✘ | ✔ | ✔ | ✘ | ✔ (hERG) |
| Druggability / CNS-MPO | ✔ | partial | partial | ✘ | ✔ |
| CNS target activity (AChE/BACE1/MAO/GSK-3β) | ✘ | ✘ | ✘ | ✔ (generic, similarity) | ✔ (measured QSAR) |
| BBB-gated → disease-level brain effect | ✘ | ✘ | ✘ | ✘ | ✔ (novel integration) |
| Nearest measured-analogue evidence per call | ✘ | ✘ | ✘ | ✔ (similarity-based) | ✔ |
| Conformal calibrated confidence | ✘ | uncertainty est. | ✘ | ✘ | ✔ (coverage-verified) |
| Open / transparent / local | ✔/✘ | web | web | web | ✔ (fully inspectable) |

No single existing tool does the thing BrainSafe is built around: unify *measured* CNS-target
polypharmacology, gate it by BBB penetration, roll it up into disease-level brain-effect scores, and
carry a safety anti-target, conformal confidence, and measured-analogue evidence alongside — all in one
transparent tool.

## 5. Why not just use a general-purpose LLM?
We ran a pre-registered head-to-head against four LLMs (Gemini Pro, ChatGPT/GPT-4o, Perplexity,
Claude), scored against a frozen measured-data key with live ChEMBL verification (Supplementary Table
S13). The result splits cleanly in two. On well-known drugs the LLMs are strong — three of four matched
or beat BrainSafe at BBB and hERG classification. But 45% of the ChEMBL identifiers they cited as
evidence were fabricated or pointed to the wrong molecule, and all four invented a target and potency
for an unpublished compound. BrainSafe fabricated nothing and reported honest uncertainty on the novel
structure. LLMs are useful for summarising known chemistry; they cannot be trusted for verifiable
provenance or for genuinely new molecules, which is where a grounded tool matters.

## 6. Where the novelty is, and isn't
The individual methods are standard: ECFP/descriptor features with RF/ExtraTrees/HistGB ensembles, BBB
and hERG ML, target QSAR, an applicability domain, and conformal prediction (Norinder et al., 2014).
Multi-target Alzheimer's QSAR has been published before. What we have not seen assembled elsewhere is
the *combination* — a BBB-gated, evidence-grounded, conformal-calibrated, safety-aware CNS
disease-effect profiler driven entirely by measured public data, released as an open tool. This is an
application and integration contribution, not a new algorithm, and we frame it that way.

## 7. Honest gaps a reviewer will (rightly) raise
1. **No prospective wet-lab validation** — the strongest test of real-world use; it needs experiments
   and cannot be manufactured on paper. (We do report a temporal split as a proxy for prospective use.)
2. **Analogue density** — ChEMBL target sets have moderate scaffold-split test-to-train similarity;
   the cluster split mitigates this, and we state it.
3. **Engagement ≠ efficacy** — we predict molecular engagement, not clinical outcome.
4. **Scope** — receptor endpoints are ranking-grade, and safety is represented by hERG alone; other
   liabilities are out of scope.
5. **Imbalanced endpoints** (BACE1, GSK-3β) are reported via AUROC/MCC, with threshold sensitivity in
   Supplementary Table S4.

## 8. Bottom line
BrainSafe is a scientifically valid, rigorously validated and genuinely useful integrative CNS
profiler. Its per-endpoint performance is competitive-to-strong, its validation is more honest than
most, and its integration — BBB-gated disease effect plus safety plus conformal confidence plus
measured-analogue evidence — is not offered as a unit by existing tools. It is well suited to an
application/resource publication with honest framing, and it is not a methodological breakthrough. The
route to a flagship predictor paper runs through prospective validation and a broader, more balanced
target panel.

## Sources
- BBB / B3DB benchmarks: [Sci Rep 2024](https://www.nature.com/articles/s41598-024-66897-y), [Chem Res Toxicol](https://pubs.acs.org/doi/10.1021/acs.chemrestox.0c00343), [DeePred-BBB](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112838/)
- hERG benchmarks: [JCIM benchmark study](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00744), [Front. Pharmacol. 2022](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2022.951083/full)
- ADMET tools: [ADMETlab 3.0 (NAR 2024)](https://academic.oup.com/nar/article/52/W1/W422/7640525), [ADMETlab 2.0 (NAR 2021)](https://academic.oup.com/nar/article/49/W1/W5/6249611), [admetSAR 3.0](https://lmmd.ecust.edu.cn/admetsar3/about/endpoint.php)
- Target prediction: [SwissTargetPrediction (NAR 2019)](https://academic.oup.com/nar/article/47/W1/W357/5491750), [Polypharmacology Browser PPB2](https://www.researchgate.net/publication/329741123_The_Polypharmacology_Browser_PPB2_Target_Prediction_Combining_Nearest_Neighbors_with_Machine_Learning)
- CNS target QSAR: [AChE ML screening (Front. Neurosci. 2022)](https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2022.1007389/full), [Dual AChE/MAO-B ML (Mol. Divers. 2024)](https://link.springer.com/article/10.1007/s11030-024-11061-x)
- Conformal prediction: [Norinder et al. 2014 / CP in drug discovery (J Pharm Sci 2020)](https://jpharmsci.org/article/S0022-3549(20)30589-X/fulltext), [Large-scale QSAR + CP (J Cheminform 2018)](https://link.springer.com/article/10.1186/s13321-018-0325-4)
