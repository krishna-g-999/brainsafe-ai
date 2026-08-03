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
| logbb | Total brain/plasma ratio (logBB) | regression | 1,058 | R2 | 0.577 ± 0.081 | 0.455 ± 0.145 | Bulk brain distribution |
| caco2_permeability | Caco-2 permeability | regression | 897 | R2 | 0.734 ± 0.051 | 0.593 ± 0.126 | Passive membrane permeability |
| pgp_substrate | P-glycoprotein substrate | classification | 1,371 | AUROC | 0.858 ± 0.032 | 0.808 ± 0.054 | Active efflux out of the brain |
| pgp_inhibition | P-glycoprotein inhibition | classification | 1,212 | AUROC | 0.955 ± 0.018 | 0.937 ± 0.024 | Efflux-mediated drug interactions |
| solubility | Aqueous solubility (logS) | regression | 9,573 | R2 | 0.804 ± 0.017 | 0.763 ± 0.066 | Formulation and absorption |
| lipophilicity | Lipophilicity (logD) | regression | 4,200 | R2 | 0.639 ± 0.028 | 0.564 ± 0.054 | Permeability/promiscuity balance |
| plasma_protein_binding | Plasma protein binding | regression | 1,797 | R2 | 0.434 ± 0.090 | 0.374 ± 0.104 | Determines free fraction |
| clearance_hepatocyte | Hepatocyte clearance | regression | 1,020 | R2 | 0.230 ± 0.104 | 0.193 ± 0.048 | Metabolic stability, exposure duration |


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
| AChE | Acetylcholinesterase | 2020 | 3,392 | 839 | AUROC | 0.785 |
| BChE | Butyrylcholinesterase | 2021 | 2,015 | 527 | AUROC | 0.737 |
| BACE1 | Beta-secretase 1 | 2017 | 6,497 | 1,604 | AUROC | 0.908 |
| GSK3B | Glycogen synthase kinase-3 beta | 2021 | 3,270 | 417 | AUROC | 0.657 |
| MAO_A | Monoamine oxidase A | 2020 | 1,667 | 446 | AUROC | 0.611 |
| MAO_B | Monoamine oxidase B | 2020 | 2,697 | 755 | AUROC | 0.781 |
| hERG | hERG potassium channel | 2019 | 4,570 | 1,248 | AUROC | 0.785 |
| D2 | Dopamine D2 receptor | 2019 | 6,081 | 1,396 | R2 | 0.042 |
| A2A | Adenosine A2A receptor | 2021 | 4,479 | 1,057 | R2 | 0.338 |
| HT2A | Serotonin 5-HT2A receptor | 2020 | 4,276 | 895 | R2 | 0.182 |
| SERT | Serotonin transporter | 2015 | 3,380 | 999 | R2 | 0.100 |
| antioxidant_DPPH | Radical-scavenging capacity | 2016 | 2,340 | 522 | R2 | 0.009 |


Classifier endpoints: mean AUROC 0.752 (range 0.611 to 0.908). Regression endpoints: mean R2 0.134 (range 0.009 to 0.338).


## Table 6. Prospective validation under a scaffold hold-out

Twenty per cent of Bemis-Murcko scaffolds were withheld per target and every model retrained on the remainder, so no held-out compound shares a scaffold with anything its model saw. Thresholds were recalibrated on held-out negatives and an independent background sample. Targets marked excluded produced a threshold at the permitted floor, meaning no separation from background chemistry, and do not contribute to the pooled estimate.


| Target | Train actives | Held-out actives | Held-out scaffolds | Threshold | Recall | 95% CI | Note |
|---|---|---|---|---|---|---|---|
| OX2 | 2,890 | 962 | 233 | 0.050 | 0.992 | [0.984, 0.996] | excluded |
| LRRK2 | 906 | 267 | 90 | 0.050 | 0.985 | [0.962, 0.994] | excluded |
| mTOR | 2,421 | 563 | 216 | 0.377 | 0.977 | [0.961, 0.986] |  |
| H3 | 2,615 | 597 | 308 | 0.758 | 0.926 | [0.903, 0.945] |  |
| PDE4B | 919 | 209 | 99 | 0.148 | 0.923 | [0.879, 0.952] |  |
| DAT | 1,130 | 207 | 93 | 0.315 | 0.903 | [0.855, 0.937] |  |
| OX1 | 2,601 | 711 | 178 | 0.955 | 0.895 | [0.870, 0.915] |  |
| OPRM1 | 3,004 | 714 | 303 | 0.838 | 0.895 | [0.870, 0.915] |  |
| NLRP3 | 186 | 36 | 27 | 0.053 | 0.889 | [0.747, 0.956] | excluded |
| PDE10A | 3,212 | 872 | 326 | 0.980 | 0.877 | [0.854, 0.897] |  |
| SERT | 2,293 | 735 | 190 | 0.573 | 0.861 | [0.834, 0.884] |  |
| HT6 | 2,023 | 718 | 167 | 0.961 | 0.858 | [0.830, 0.882] |  |
| CSF1R | 1,543 | 241 | 152 | 0.726 | 0.846 | [0.796, 0.887] |  |
| D3 | 3,337 | 790 | 351 | 0.980 | 0.846 | [0.819, 0.869] |  |
| Nav1_1 | 40 | 25 | 9 | 0.325 | 0.840 | [0.653, 0.936] |  |
| OPRK1 | 2,497 | 562 | 235 | 0.948 | 0.838 | [0.805, 0.866] |  |
| HDAC1 | 2,414 | 641 | 292 | 0.900 | 0.836 | [0.806, 0.863] |  |
| CB1 | 2,201 | 479 | 200 | 0.809 | 0.827 | [0.790, 0.858] |  |
| NET | 1,307 | 265 | 90 | 0.646 | 0.819 | [0.768, 0.861] |  |
| Nav1_7 | 2,250 | 494 | 179 | 0.941 | 0.818 | [0.781, 0.849] |  |
| HT1A | 2,902 | 801 | 327 | 0.980 | 0.814 | [0.786, 0.839] |  |
| P2X7 | 2,432 | 506 | 139 | 0.979 | 0.800 | [0.763, 0.833] |  |
| A2A | 3,428 | 924 | 364 | 0.994 | 0.777 | [0.749, 0.803] |  |
| A1 | 1,514 | 429 | 177 | 0.972 | 0.769 | [0.727, 0.807] |  |
| a7nAChR | 251 | 86 | 35 | 0.648 | 0.767 | [0.668, 0.844] |  |
| GABA_A | 184 | 17 | 16 | 0.089 | 0.765 | [0.527, 0.904] |  |
| HT7 | 1,166 | 312 | 133 | 0.954 | 0.721 | [0.669, 0.768] |  |
| HDAC6 | 2,960 | 755 | 335 | 0.990 | 0.699 | [0.666, 0.731] |  |
| Nav1_5 | 198 | 44 | 27 | 0.197 | 0.682 | [0.534, 0.800] |  |
| Sigma1 | 1,571 | 357 | 186 | 0.961 | 0.675 | [0.625, 0.722] |  |
| GluN2B | 652 | 252 | 35 | 0.999 | 0.619 | [0.558, 0.677] |  |
| HT2A | 3,161 | 790 | 317 | 0.987 | 0.618 | [0.583, 0.651] |  |
| COX2 | 918 | 264 | 67 | 0.964 | 0.587 | [0.527, 0.645] |  |
| D2 | 2,880 | 796 | 341 | 0.989 | 0.580 | [0.546, 0.614] |  |
| SIRT1 | 134 | 28 | 20 | 0.529 | 0.464 | [0.295, 0.642] |  |
| mGluR5 | 975 | 225 | 106 | 0.995 | 0.444 | [0.381, 0.510] |  |
| MT1 | 545 | 152 | 48 | 0.996 | 0.283 | [0.217, 0.359] |  |
| KEAP1 | 97 | 25 | 15 | 0.978 | 0.240 | [0.115, 0.434] |  |
| GluA2 | 74 | 23 | 9 | 0.816 | 0.217 | [0.097, 0.419] |  |


Pooled recall 12,325/15,609 = 0.790; median per-target 0.807; 19 of 36 targets at or above 0.80.
