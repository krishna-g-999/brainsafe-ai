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
| kpuu | Unbound brain/plasma ratio (Kp,uu) | regression | 566 | R2 | 0.406 ± 0.096 | 0.352 ± 0.158 | Free drug available to CNS targets |
| logbb | Total brain/plasma ratio (logBB) | regression | 1,058 | R2 | 0.584 ± 0.079 | 0.413 ± 0.145 | Bulk brain distribution |
| caco2_permeability | Caco-2 permeability | regression | 897 | R2 | 0.736 ± 0.050 | 0.582 ± 0.174 | Passive membrane permeability |
| pgp_substrate | P-glycoprotein substrate | classification | 1,371 | AUROC | 0.856 ± 0.032 | 0.807 ± 0.041 | Active efflux out of the brain |
| pgp_inhibition | P-glycoprotein inhibition | classification | 1,212 | AUROC | 0.955 ± 0.018 | 0.935 ± 0.029 | Efflux-mediated drug interactions |
| solubility | Aqueous solubility (logS) | regression | 9,573 | R2 | 0.801 ± 0.013 | 0.726 ± 0.084 | Formulation and absorption |
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
| AChE | Acetylcholinesterase | 2020 | 4,141 | 1,028 | AUROC | 0.761 |
| BChE | Butyrylcholinesterase | 2021 | 2,696 | 617 | AUROC | 0.801 |
| BACE1 | Beta-secretase 1 | 2017 | 6,901 | 1,673 | AUROC | 0.910 |
| GSK3B | Glycogen synthase kinase-3 beta | 2021 | 4,205 | 546 | AUROC | 0.764 |
| MAO_A | Monoamine oxidase A | 2020 | 2,859 | 819 | AUROC | 0.727 |
| MAO_B | Monoamine oxidase B | 2020 | 3,719 | 923 | AUROC | 0.843 |
| hERG | hERG potassium channel | 2020 | 8,158 | 2,059 | AUROC | 0.720 |
| D2 | Dopamine D2 receptor | 2019 | 6,796 | 1,556 | R2 | 0.080 |
| A2A | Adenosine A2A receptor | 2020 | 4,600 | 1,516 | R2 | 0.229 |
| HT2A | Serotonin 5-HT2A receptor | 2020 | 4,741 | 986 | R2 | 0.155 |
| SERT | Serotonin transporter | 2015 | 3,647 | 1,147 | R2 | 0.111 |
| antioxidant_DPPH | Radical-scavenging capacity | 2016 | 2,340 | 522 | R2 | 0.007 |


Classifier endpoints: mean AUROC 0.789 (range 0.720 to 0.910). Regression endpoints: mean R2 0.116 (range 0.007 to 0.229).


## Table 6. Prospective validation under a scaffold hold-out

Twenty per cent of Bemis-Murcko scaffolds were withheld per target and every model retrained on the remainder, so no held-out compound shares a scaffold with anything its model saw. Thresholds were recalibrated on held-out negatives and an independent background sample. Targets marked excluded produced a threshold at the permitted floor, meaning no separation from background chemistry, and do not contribute to the pooled estimate.


| Target | Train actives | Held-out actives | Held-out scaffolds | Threshold | Recall | 95% CI | Note |
|---|---|---|---|---|---|---|---|
| HDAC6 | 3,052 | 663 | 335 | 0.710 | 0.970 | [0.954, 0.980] |  |
| mTOR | 2,354 | 630 | 216 | 0.050 | 0.968 | [0.951, 0.979] | excluded |
| DAT | 901 | 436 | 93 | 0.223 | 0.938 | [0.911, 0.957] |  |
| H3 | 2,578 | 634 | 308 | 0.437 | 0.931 | [0.908, 0.948] |  |
| CSF1R | 1,406 | 378 | 152 | 0.816 | 0.921 | [0.889, 0.944] |  |
| HDAC1 | 2,352 | 703 | 292 | 0.706 | 0.919 | [0.896, 0.937] |  |
| PDE10A | 3,305 | 779 | 326 | 0.649 | 0.917 | [0.895, 0.934] |  |
| OX2 | 3,136 | 716 | 233 | 0.716 | 0.911 | [0.887, 0.929] |  |
| NET | 1,244 | 328 | 90 | 0.288 | 0.909 | [0.872, 0.935] |  |
| PDE4B | 908 | 220 | 99 | 0.146 | 0.909 | [0.864, 0.940] |  |
| OX1 | 2,648 | 664 | 178 | 0.413 | 0.898 | [0.872, 0.918] |  |
| OPRM1 | 2,964 | 754 | 302 | 0.472 | 0.897 | [0.873, 0.916] |  |
| Nav1_5 | 196 | 46 | 27 | 0.050 | 0.891 | [0.770, 0.953] | excluded |
| Nav1_1 | 40 | 25 | 9 | 0.223 | 0.880 | [0.700, 0.958] |  |
| LRRK2 | 936 | 237 | 90 | 0.250 | 0.878 | [0.830, 0.913] |  |
| CB1 | 2,235 | 445 | 196 | 0.144 | 0.867 | [0.833, 0.896] |  |
| HT6 | 2,151 | 590 | 167 | 0.921 | 0.864 | [0.834, 0.890] |  |
| TAAR1 | 68 | 25 | 9 | 0.513 | 0.840 | [0.653, 0.936] |  |
| D3 | 3,272 | 855 | 351 | 0.932 | 0.830 | [0.804, 0.854] |  |
| SERT | 2,347 | 624 | 185 | 0.552 | 0.829 | [0.797, 0.856] |  |
| a7nAChR | 260 | 77 | 35 | 0.281 | 0.818 | [0.718, 0.888] |  |
| A2A | 3,303 | 798 | 337 | 0.917 | 0.815 | [0.786, 0.840] |  |
| HT1A | 2,902 | 801 | 327 | 0.937 | 0.814 | [0.786, 0.839] |  |
| Nav1_7 | 2,365 | 379 | 179 | 0.903 | 0.794 | [0.751, 0.832] |  |
| OPRK1 | 2,501 | 558 | 234 | 0.915 | 0.767 | [0.730, 0.800] |  |
| HT7 | 1,214 | 264 | 133 | 0.911 | 0.746 | [0.690, 0.795] |  |
| GBA1 | 102 | 18 | 10 | 0.294 | 0.722 | [0.491, 0.875] |  |
| NLRP3 | 170 | 52 | 26 | 0.409 | 0.692 | [0.557, 0.801] |  |
| Sigma1 | 1,518 | 410 | 186 | 0.972 | 0.683 | [0.636, 0.726] |  |
| D2 | 2,971 | 686 | 340 | 0.939 | 0.681 | [0.645, 0.715] |  |
| HT2A | 2,983 | 834 | 306 | 0.941 | 0.676 | [0.644, 0.707] |  |
| A1 | 1,486 | 457 | 177 | 0.951 | 0.656 | [0.612, 0.699] |  |
| mGluR5 | 990 | 210 | 106 | 0.910 | 0.652 | [0.586, 0.714] |  |
| P2X7 | 2,123 | 815 | 139 | 0.995 | 0.531 | [0.497, 0.565] |  |
| SIRT1 | 141 | 21 | 17 | 0.542 | 0.476 | [0.283, 0.676] |  |
| COX2 | 991 | 191 | 67 | 0.827 | 0.424 | [0.356, 0.495] |  |
| GABA_A | 162 | 39 | 16 | 0.878 | 0.410 | [0.271, 0.566] |  |
| MT1 | 570 | 127 | 48 | 0.974 | 0.315 | [0.241, 0.400] |  |
| KEAP1 | 85 | 37 | 15 | 0.841 | 0.297 | [0.175, 0.458] |  |
| GluN2B | 772 | 132 | 35 | 0.996 | 0.250 | [0.184, 0.330] |  |
| GluA2 | 65 | 32 | 9 | 0.835 | 0.000 | [0.000, 0.107] |  |


Pooled recall 12,858/16,014 = 0.803; median per-target 0.815; 21 of 39 targets at or above 0.80.


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
| a4b2nAChR | 796 | 278 | 146 | 0.923 | 0.863 | deployed |  |
| a3b4nAChR | 398 | 153 | 187 | 0.974 | 0.970 | deployed |  |
| Nav1_6 | 681 | 187 | 45 | 0.862 | 0.764 | deployed |  |
| Nav1_8 | 501 | 163 | 39 | 0.956 | 0.947 | deployed |  |
| Cav3_2 | 633 | 201 | 33 | 0.982 | 0.980 | deployed |  |
| GABAA_a5 | 672 | 284 | 4 | n/a | n/a | not trained | only 4 measured inactives, cannot set a threshold honestly |
| CGRP | 761 | 333 | 26 | 0.985 | 0.997 | deployed | only 26 measured inactives, cannot set a threshold honestly |
| DHODH | 1421 | 360 | 145 | 0.966 | 0.981 | deployed |  |
| RIPK1 | 2349 | 826 | 719 | 0.966 | 0.984 | deployed |  |
