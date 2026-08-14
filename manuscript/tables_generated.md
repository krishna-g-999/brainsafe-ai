# Generated tables (all numbers computed from saved cross-validation predictions)


## Table 1. Core target panel (13 endpoints), 10-fold cross-validation

| Endpoint | Target | Task | Compounds | Scaffolds | Active fraction | Train/fold | Test/fold | Metric | Random 10-fold | Scaffold 10-fold | Why this endpoint |
|---|---|---|---|---|---|---|---|---|---|---|---|
| BBB | Blood-brain barrier penetration | classification | 7,805 | 2,439 | 0.64 | 7,024 | 780 | AUROC | 0.961 ± 0.007 | 0.920 ± 0.037 | Gate: no CNS effect without brain entry |
| AChE | Acetylcholinesterase | classification | 4,387 | 2,116 | 0.72 | 3,948 | 439 | AUROC | 0.963 ± 0.008 | 0.921 ± 0.021 | Symptomatic Alzheimer therapy (donepezil class) |
| BChE | Butyrylcholinesterase | classification | 2,621 | 1,345 | 0.70 | 2,359 | 262 | AUROC | 0.968 ± 0.012 | 0.937 ± 0.021 | Rises as Alzheimer progresses; selective-inhibitor target |
| BACE1 | Beta-secretase 1 | classification | 8,501 | 3,088 | 0.91 | 7,651 | 850 | AUROC | 0.967 ± 0.010 | 0.956 ± 0.021 | Rate-limiting step of amyloid-beta generation |
| GSK3B | Glycogen synthase kinase-3 beta | classification | 4,958 | 2,118 | 0.93 | 4,462 | 496 | AUROC | 0.969 ± 0.013 | 0.937 ± 0.030 | Tau hyperphosphorylation; neuroprotection |
| MAO_A | Monoamine oxidase A | classification | 2,228 | 827 | 0.41 | 2,005 | 223 | AUROC | 0.947 ± 0.017 | 0.868 ± 0.046 | Serotonin/noradrenaline catabolism; depression |
| MAO_B | Monoamine oxidase B | classification | 3,665 | 1,196 | 0.68 | 3,298 | 366 | AUROC | 0.955 ± 0.007 | 0.889 ± 0.033 | Dopamine catabolism; Parkinson therapy (selegiline) |
| hERG | hERG potassium channel | classification | 5,875 | 3,377 | 0.41 | 5,287 | 587 | AUROC | 0.954 ± 0.007 | 0.921 ± 0.035 | Cardiotoxicity liability; principal safety filter |
| D2 | Dopamine D2 receptor | regression | 7,734 | 3,329 | n/a | 6,961 | 773 | R2 | 0.601 ± 0.018 | 0.483 ± 0.051 | Antipsychotic efficacy; motor control |
| A2A | Adenosine A2A receptor | regression | 6,785 | 2,654 | n/a | 6,106 | 678 | R2 | 0.682 ± 0.023 | 0.576 ± 0.066 | Non-dopaminergic Parkinson target (istradefylline) |
| HT2A | Serotonin 5-HT2A receptor | regression | 5,989 | 2,254 | n/a | 5,390 | 599 | R2 | 0.636 ± 0.028 | 0.490 ± 0.054 | Atypical antipsychotics; psychedelics |
| SERT | Serotonin transporter | regression | 4,572 | 1,459 | n/a | 4,115 | 457 | R2 | 0.602 ± 0.038 | 0.389 ± 0.124 | SSRI antidepressant target |
| antioxidant_DPPH | Radical-scavenging capacity | regression | 2,862 | 919 | n/a | 2,576 | 286 | R2 | 0.669 ± 0.065 | 0.434 ± 0.100 | Oxidative stress in neurodegeneration |


## Table 2. Binder classifiers validated against held-out measured inactives

| Endpoint | Target | Measured binders | Training negatives | AUROC (held-out measured inactives) | Threshold | Sensitivity | Why this endpoint |
|---|---|---|---|---|---|---|---|
| D2 | Dopamine D2 receptor | 3,676 | 10,895 | 0.938 | 0.986 | 0.785 | Antipsychotic efficacy |
| A2A | Adenosine A2A receptor | 4,352 | 12,941 | 0.960 | 0.994 | 0.861 | Non-dopaminergic Parkinson target |
| HT2A | Serotonin 5-HT2A receptor | 3,951 | 11,813 | 0.953 | 0.983 | 0.850 | Atypical antipsychotic profile |
| SERT | Serotonin transporter | 3,028 | 9,012 | 0.983 | 0.768 | 0.970 | SSRI antidepressant target |
| HT1A | 5-HT1A receptor | 3,703 | 11,012 | 0.989 | 0.758 | 0.971 | Anxiety, depression (buspirone) |
| HT6 | 5-HT6 receptor | 2,741 | 8,195 | 0.977 | 0.947 | 0.942 | Cognition enhancement in Alzheimer |
| HT7 | 5-HT7 receptor | 1,478 | 4,377 | 0.960 | 0.921 | 0.913 | Mood, circadian rhythm, sleep |
| H3 | Histamine H3 receptor | 3,212 | 9,548 | 0.990 | 0.919 | 0.964 | Wakefulness (pitolisant), cognition |
| DAT | Dopamine transporter | 1,337 | 3,855 | 0.986 | 0.531 | 0.984 | ADHD, addiction, stimulant liability |
| NET | Noradrenaline transporter | 1,572 | 4,583 | 0.990 | 0.715 | 0.971 | Depression, ADHD (atomoxetine) |
| Sigma1 | Sigma-1 receptor | 1,928 | 5,762 | 0.940 | 0.955 | 0.856 | Neuroprotection, ER-stress chaperone |
| CB1 | Cannabinoid CB1 receptor | 2,680 | 7,936 | 0.965 | 0.988 | 0.878 | Pain, appetite, mood |
| OPRK1 | Kappa-opioid receptor | 3,059 | 9,022 | 0.968 | 0.991 | 0.903 | Analgesia, dysphoria, mood |
| OPRM1 | Mu-opioid receptor | 3,718 | 10,971 | 0.981 | 0.874 | 0.969 | Analgesia, addiction liability |
| D3 | Dopamine D3 receptor | 4,127 | 12,289 | 0.978 | 0.960 | 0.938 | Addiction, Parkinson motor complications |
| A1 | Adenosine A1 receptor | 1,943 | 5,720 | 0.954 | 0.978 | 0.847 | Neuroprotection, epilepsy, sedation |
| a7nAChR | Alpha-7 nicotinic receptor | 337 | 921 | 0.982 | 0.642 | 0.932 | Cognition, neuroinflammation in Alzheimer |
| LRRK2 | LRRK2 kinase | 1,173 | 3,510 | n/a | 0.400 | n/a | Most common genetic cause of Parkinson disease |


## Table 3. ADME / exposure layer (9 endpoints)

| Endpoint | Property | Task | Compounds | Metric | Random 10-fold | Scaffold 10-fold | Why this endpoint |
|---|---|---|---|---|---|---|---|
| kpuu | Unbound brain/plasma ratio (Kp,uu) | regression | 566 | R2 | 0.404 ± 0.099 | 0.352 ± 0.158 | Free drug available to CNS targets |
| logbb | Total brain/plasma ratio (logBB) | regression | 1,058 | R2 | 0.577 ± 0.081 | 0.423 ± 0.150 | Bulk brain distribution |
| caco2_permeability | Caco-2 permeability | regression | 897 | R2 | 0.734 ± 0.051 | 0.584 ± 0.158 | Passive membrane permeability |
| pgp_substrate | P-glycoprotein substrate | classification | 1,371 | AUROC | 0.858 ± 0.032 | 0.806 ± 0.038 | Active efflux out of the brain |
| pgp_inhibition | P-glycoprotein inhibition | classification | 1,212 | AUROC | 0.955 ± 0.018 | 0.935 ± 0.029 | Efflux-mediated drug interactions |
| solubility | Aqueous solubility (logS) | regression | 9,573 | R2 | 0.804 ± 0.017 | 0.728 ± 0.077 | Formulation and absorption |
| lipophilicity | Lipophilicity (logD) | regression | 4,200 | R2 | 0.639 ± 0.028 | 0.566 ± 0.042 | Permeability/promiscuity balance |
| plasma_protein_binding | Plasma protein binding | regression | 1,797 | R2 | 0.434 ± 0.090 | 0.364 ± 0.056 | Determines free fraction |
| clearance_hepatocyte | Hepatocyte clearance | regression | 1,020 | R2 | 0.230 ± 0.104 | 0.206 ± 0.064 | Metabolic stability, exposure duration |


## Table 4. Between-fold error-bar decomposition

Observed between-fold SD separated into sampling noise (finite test set) and genuine fold-to-fold heterogeneity, by within-fold bootstrap.


| Endpoint | Random SD | of which sampling | heterogeneity share | Scaffold SD | of which sampling | heterogeneity share |
|---|---|---|---|---|---|---|
| BBB | 0.0075 | 0.0062 | 32% | 0.0367 | 0.0102 | 92% |
| AChE | 0.0083 | 0.0094 | 0% | 0.0206 | 0.0137 | 56% |
| BChE | 0.0118 | 0.0096 | 34% | 0.0212 | 0.0145 | 53% |
| BACE1 | 0.0099 | 0.0106 | 0% | 0.0209 | 0.0118 | 68% |
| GSK3B | 0.0128 | 0.0112 | 23% | 0.0302 | 0.0204 | 54% |
| MAO_A | 0.0169 | 0.0147 | 24% | 0.0464 | 0.0251 | 71% |
| MAO_B | 0.0075 | 0.0103 | 0% | 0.0330 | 0.0185 | 69% |
| hERG | 0.0066 | 0.0082 | 0% | 0.0351 | 0.0109 | 90% |
| D2 | 0.0177 | 0.0255 | 0% | 0.0515 | 0.0299 | 66% |
| A2A | 0.0228 | 0.0251 | 0% | 0.0662 | 0.0289 | 81% |
| HT2A | 0.0276 | 0.0264 | 8% | 0.0537 | 0.0311 | 67% |
| SERT | 0.0379 | 0.0329 | 25% | 0.1239 | 0.0404 | 89% |
| antioxidant_DPPH | 0.0648 | 0.0457 | 50% | 0.0995 | 0.0592 | 65% |


## Table 5. Temporal (future-compound) validation

Models are trained only on compounds published before the cutoff year and tested on compounds published after it. This is the most demanding regime and the closest analogue of prospective use.


| Endpoint | Target | Cutoff year | Train | Test | Metric | Score |
|---|---|---|---|---|---|---|
| AChE | Acetylcholinesterase | 2020 | 4,141 | 1,028 | AUROC | 0.760 |
| BChE | Butyrylcholinesterase | 2021 | 2,696 | 617 | AUROC | 0.801 |
| BACE1 | Beta-secretase 1 | 2017 | 6,901 | 1,673 | AUROC | 0.910 |
| GSK3B | Glycogen synthase kinase-3 beta | 2021 | 4,205 | 546 | AUROC | 0.776 |
| MAO_A | Monoamine oxidase A | 2020 | 2,859 | 819 | AUROC | 0.727 |
| MAO_B | Monoamine oxidase B | 2020 | 3,719 | 923 | AUROC | 0.843 |
| hERG | hERG potassium channel | 2020 | 8,158 | 2,059 | AUROC | 0.715 |
| D2 | Dopamine D2 receptor | 2019 | 6,796 | 1,556 | R2 | 0.080 |
| A2A | Adenosine A2A receptor | 2020 | 4,600 | 1,516 | R2 | 0.226 |
| HT2A | Serotonin 5-HT2A receptor | 2020 | 4,741 | 986 | R2 | 0.157 |
| SERT | Serotonin transporter | 2015 | 3,647 | 1,147 | R2 | 0.111 |
| antioxidant_DPPH | Radical-scavenging capacity | 2016 | 2,340 | 522 | R2 | 0.009 |


Classifier endpoints: mean AUROC 0.790 (range 0.715 to 0.910). Regression endpoints: mean R2 0.117 (range 0.009 to 0.226).


## Table 6. Prospective validation under a scaffold hold-out

Twenty per cent of Bemis-Murcko scaffolds were withheld per target and every model retrained on the remainder, so no held-out compound shares a scaffold with anything its model saw. Thresholds were recalibrated on held-out negatives and an independent background sample. Targets marked excluded produced a threshold at the permitted floor, meaning no separation from background chemistry, and do not contribute to the pooled estimate.


| Target | Train actives | Held-out actives | Held-out scaffolds | Threshold | Recall | 95% CI | Note |
|---|---|---|---|---|---|---|---|
| OX2 | 3,228 | 624 | 233 | 0.050 | 0.990 | [0.979, 0.996] | excluded |
| LRRK2 | 752 | 421 | 90 | 0.050 | 0.957 | [0.933, 0.973] | excluded |
| PDE4B | 954 | 174 | 99 | 0.262 | 0.954 | [0.912, 0.977] |  |
| mTOR | 2,429 | 555 | 216 | 0.494 | 0.950 | [0.928, 0.965] |  |
| CSF1R | 1,510 | 274 | 152 | 0.190 | 0.927 | [0.890, 0.952] |  |
| H3 | 2,615 | 597 | 308 | 0.758 | 0.926 | [0.903, 0.945] |  |
| GABA_A | 166 | 35 | 16 | 0.098 | 0.914 | [0.776, 0.970] |  |
| DAT | 1,130 | 207 | 93 | 0.315 | 0.903 | [0.855, 0.937] |  |
| PDE10A | 3,141 | 943 | 326 | 0.979 | 0.902 | [0.882, 0.920] |  |
| SERT | 2,423 | 548 | 185 | 0.560 | 0.901 | [0.874, 0.924] |  |
| Nav1_5 | 194 | 48 | 27 | 0.367 | 0.875 | [0.753, 0.941] |  |
| NLRP3 | 176 | 46 | 27 | 0.953 | 0.870 | [0.743, 0.939] |  |
| OPRM1 | 2,931 | 787 | 302 | 0.842 | 0.859 | [0.833, 0.882] |  |
| HT6 | 2,023 | 718 | 167 | 0.961 | 0.858 | [0.830, 0.882] |  |
| a7nAChR | 266 | 71 | 35 | 0.494 | 0.845 | [0.743, 0.911] |  |
| Nav1_7 | 2,172 | 572 | 179 | 0.957 | 0.844 | [0.812, 0.872] |  |
| D3 | 3,283 | 844 | 351 | 0.965 | 0.842 | [0.816, 0.865] |  |
| HDAC1 | 2,421 | 634 | 292 | 0.905 | 0.841 | [0.810, 0.867] |  |
| OX1 | 2,736 | 576 | 178 | 0.983 | 0.833 | [0.801, 0.862] |  |
| KEAP1 | 105 | 17 | 15 | 0.430 | 0.824 | [0.590, 0.938] |  |
| NET | 1,307 | 265 | 90 | 0.646 | 0.819 | [0.768, 0.861] |  |
| HT1A | 2,902 | 801 | 327 | 0.980 | 0.814 | [0.786, 0.839] |  |
| A2A | 3,353 | 748 | 337 | 0.993 | 0.813 | [0.783, 0.839] |  |
| OPRK1 | 2,575 | 484 | 234 | 0.985 | 0.800 | [0.762, 0.833] |  |
| CB1 | 2,161 | 519 | 196 | 0.969 | 0.775 | [0.737, 0.808] |  |
| A1 | 1,535 | 408 | 177 | 0.980 | 0.762 | [0.719, 0.801] |  |
| HDAC6 | 2,928 | 787 | 335 | 0.995 | 0.729 | [0.697, 0.759] |  |
| HT7 | 1,166 | 312 | 133 | 0.954 | 0.721 | [0.669, 0.768] |  |
| Sigma1 | 1,571 | 357 | 186 | 0.961 | 0.675 | [0.625, 0.722] |  |
| MT1 | 571 | 126 | 48 | 0.971 | 0.667 | [0.581, 0.743] |  |
| SIRT1 | 141 | 21 | 17 | 0.766 | 0.667 | [0.454, 0.828] |  |
| HT2A | 3,102 | 715 | 306 | 0.996 | 0.639 | [0.603, 0.674] |  |
| GluN2B | 742 | 162 | 35 | 0.986 | 0.617 | [0.541, 0.689] |  |
| D2 | 2,989 | 668 | 340 | 0.992 | 0.554 | [0.516, 0.591] |  |
| TAAR1 | 78 | 15 | 9 | 0.397 | 0.533 | [0.301, 0.752] |  |
| P2X7 | 2,450 | 488 | 139 | 0.999 | 0.531 | [0.486, 0.575] |  |
| GBA1 | 101 | 19 | 10 | 0.453 | 0.526 | [0.317, 0.727] |  |
| COX2 | 945 | 237 | 67 | 0.951 | 0.371 | [0.312, 0.434] |  |
| mGluR5 | 986 | 214 | 106 | 0.999 | 0.308 | [0.250, 0.373] |  |
| GluA2 | 78 | 19 | 9 | 0.695 | 0.105 | [0.029, 0.314] |  |


Pooled recall 11,794/15,011 = 0.786; median per-target 0.817; 22 of 38 targets at or above 0.80.


## Table 7. Falsification analysis: each claim paired with a null model

Every hypothesis was stated so that it could fail. Where predictive power was at issue, scoring used scaffold hold-out models that never saw the compounds they scored. A refuted hypothesis is reported as prominently as a confirmed one; the purpose of the exercise was to find failures.


| Hypothesis | Verdict | Evidence |
|---|---|---|
| H1 the disease score is informative | **SUPPORTED** | top-3 accuracy 0.769 vs permutation null 0.152 (p=0.005) and frequency null 0.560 |
| H2 the curated edge weights add value | **REFUTED** | curated 0.7691, uniform 0.7678, permuted 0.7670 |
| H3 BBB gating discriminates between diseases | **REFUTED (by construction)** | the gate multiplies every disease equally and cannot change their order |
| H4 specificity transfers to novel chemistry | **SUPPORTED** | false-positive rate 0.033 on 61 distant compounds against 0.080 measured on library chemistry |
| H5 read-across beats a frequency baseline | **SUPPORTED** | recall 0.973 against 0.060 |
| H6 the disease scores match real clinical indications | **WEAKENED** | top-3 accuracy 0.352 on 162 drugs never seen in training, against permutation null 0.145 (p=0.001) and frequency null 0.654 |
| H7 some panel targets are non-discriminative and explain the silent antiepileptics | **REFUTED** | none of 35 targets ranks below AUROC 0.70; the cause is the operating point, with median deployed sensitivity 0.83 and 3 targets under 0.50 |
| H8 engaged targets are independent observations | **REFUTED** | 38 targets fire across approved drugs but span only 15 independent directions; 5 homologous pairs correlate above 0.5 |


## Table 8. Recovery of approved clinical indications, by condition

Ground truth is ChEMBL's drug_indication table restricted to phase 4, mapped to the panel through a keyword list fixed before any prediction was computed. 'Ranking only' removes the reporting threshold, separating the ability to rank conditions from the decision to stay silent. 'Silent' counts drugs for which no condition reached the reporting threshold.


| Indication | Drugs | Recovered in top 3 | Ranking only | Silent |
|---|---|---|---|---|
| Alzheimer's disease | 10 | 0.200 | 0.800 | 5 |
| Parkinson's disease | 39 | 0.154 | 0.538 | 16 |
| Depression / anxiety | 113 | 0.611 | 0.779 | 30 |
| Psychosis / schizophrenia | 73 | 0.644 | 0.740 | 14 |
| Addiction | 17 | 0.529 | 0.588 | 5 |
| ADHD | 31 | 0.226 | 0.258 | 16 |
| Chronic pain | 144 | 0.229 | 0.250 | 83 |
| Sleep / wakefulness | 27 | 0.296 | 0.296 | 11 |
| Epilepsy | 58 | 0.103 | 0.224 | 47 |


## Table 9. Targeted expansion: candidates, audit and outcome

Candidates were selected because the clinical-indication test identified epilepsy and chronic pain as the conditions served worst, and a systematic ChEMBL query found these mechanisms clear both a volume bar of 800 measured activities and a source-diversity bar. Each was then audited on its own fetched data before training, and again on deployed specificity afterwards.


| Candidate | Actives | Scaffolds | Measured inactives | AUROC vs inactives | Sensitivity | Outcome | Reason if rejected |
|---|---|---|---|---|---|---|---|
| a4b2nAChR | 796 | 278 | 146 | 0.937 | 0.899 | deployed |  |
| a3b4nAChR | 398 | 153 | 187 | 0.973 | 0.955 | deployed |  |
| Nav1_6 | 681 | 187 | 45 | 0.634 | 0.565 | deployed |  |
| Nav1_8 | 501 | 163 | 39 | 0.876 | 0.728 | deployed |  |
| Cav3_2 | 633 | 201 | 33 | 0.971 | 0.907 | withdrawn after training | active band compressed near zero, calibrated threshold 0.065, atenolol scores 0.084 |
| GABAA_a5 | 672 | 284 | 4 | n/a | n/a | not trained | only 4 measured inactives, cannot set a threshold honestly |
| CGRP | 761 | 333 | 26 | 0.985 | 0.997 | deployed | only 26 measured inactives, cannot set a threshold honestly |
| DHODH | 1421 | 360 | 145 | 0.964 | 0.967 | deployed |  |
| RIPK1 | 2349 | 826 | 719 | 0.983 | 0.969 | deployed |  |
