# Model inventory, 2026-08-19

Commit `6dd47bd`. One row per deployed estimator, taken from the estimators and their metadata on disk. Four receptors carry both a potency regression and a binder classifier, so the model count exceeds the endpoint count.

**75 estimators, 70 deployed.** Fitted between 2026-08-04 10:01 and 2026-08-17 21:50.

| family | estimators | deployed |
|---|---|---|
| binder | 52 | 47 |
| exposure | 10 | 10 |
| safety | 1 | 1 |
| target | 12 | 12 |

## Every estimator

| model | predicts | task | n train | metric | random | scaffold | calibration | deployed | fitted |
|---|---|---|---|---|---|---|---|---|---|
| A2A | potency, pChEMBL | regression | 6756 | R2 | 0.7243 | 0.6118 | none | yes | 2026-08-13 19:05 |
| AChE | probability of activity | classification | 5125 | AUROC | 0.9655 | 0.9227 | isotonic, out-of-fold | yes | 2026-08-13 18:56 |
| antioxidant_DPPH | potency, pChEMBL | regression | 2782 | R2 | 0.6592 | 0.4161 | none | yes | 2026-08-13 19:09 |
| BACE1 | probability of activity | classification | 8213 | AUROC | 0.9775 | 0.9652 | isotonic, out-of-fold | yes | 2026-08-13 18:57 |
| BBB | probability of activity | classification | 3901 | AUROC | 0.899 | 0.8777 | isotonic, out-of-fold | yes | 2026-08-13 18:56 |
| BChE | probability of activity | classification | 3278 | AUROC | 0.9724 | 0.9451 | isotonic, out-of-fold | yes | 2026-08-13 18:57 |
| D2 | potency, pChEMBL | regression | 7907 | R2 | 0.6436 | 0.5237 | none | yes | 2026-08-13 19:02 |
| GSK3B | probability of activity | classification | 5443 | AUROC | 0.965 | 0.9389 | isotonic, out-of-fold | yes | 2026-08-13 18:58 |
| hERG | probability of activity | classification | 9934 | AUROC | 0.9572 | 0.9323 | isotonic, out-of-fold | yes | 2026-08-13 18:59 |
| HT2A | potency, pChEMBL | regression | 6076 | R2 | 0.6988 | 0.5626 | none | yes | 2026-08-13 19:07 |
| MAO_A | probability of activity | classification | 3595 | AUROC | 0.9638 | 0.899 | isotonic, out-of-fold | yes | 2026-08-13 18:58 |
| MAO_B | probability of activity | classification | 4548 | AUROC | 0.9638 | 0.9189 | isotonic, out-of-fold | yes | 2026-08-13 18:58 |
| pka_basic | pKa | regression | 6384 | R2 |  |  | none | yes | 2026-08-04 10:01 |
| SERT | potency, pChEMBL | regression | 4479 | R2 | 0.6897 | 0.4612 | none | yes | 2026-08-13 19:09 |
| A1_binder | probability this compound binds this target | classification | 7384 | AUROC vs measured non-binders |  | 0.913 | sigmoid, prefit | yes | 2026-08-13 20:44 |
| A2A_binder | probability this compound binds this target | classification | 15682 | AUROC vs measured non-binders |  | 0.943 | sigmoid, prefit | yes | 2026-08-13 20:56 |
| CB1_binder | probability this compound binds this target | classification | 10244 | AUROC vs measured non-binders |  | 0.948 | sigmoid, prefit | yes | 2026-08-13 20:41 |
| CGRP_binder | probability this compound binds this target | classification | 2251 | AUROC vs measured non-binders |  | 0.985 | sigmoid, prefit | yes | 2026-08-13 21:02 |
| COX2_binder | probability this compound binds this target | classification | 4437 | AUROC vs measured non-binders |  | 0.788 | sigmoid, prefit | yes | 2026-08-13 20:46 |
| CSF1R_binder | probability this compound binds this target | classification | 6778 | AUROC vs measured non-binders |  | 0.961 | sigmoid, prefit | yes | 2026-08-13 20:47 |
| Cav3_2_binder | probability this compound binds this target | classification | 1318 | AUROC vs measured non-binders |  | 0.971 | sigmoid, prefit | WITHDRAWN | 2026-08-13 21:02 |
| D2_binder | probability this compound binds this target | classification | 13887 | AUROC vs measured non-binders |  | 0.862 | sigmoid, prefit | yes | 2026-08-13 20:55 |
| D3_binder | probability this compound binds this target | classification | 15723 | AUROC vs measured non-binders |  | 0.95 | sigmoid, prefit | yes | 2026-08-13 20:44 |
| DAT_binder | probability this compound binds this target | classification | 5090 | AUROC vs measured non-binders |  | 0.962 | sigmoid, prefit | yes | 2026-08-13 20:39 |
| DHODH_binder | probability this compound binds this target | classification | 3559 | AUROC vs measured non-binders |  | 0.964 | sigmoid, prefit | yes | 2026-08-13 21:02 |
| GABA_A_binder | probability this compound binds this target | classification | 762 | AUROC vs measured non-binders |  | 0.835 | sigmoid, prefit | yes | 2026-08-13 20:50 |
| GBA1_binder | probability this compound binds this target | classification | 446 | AUROC vs measured non-binders |  | 0.845 | sigmoid, prefit | yes | 2026-08-13 20:53 |
| GluA2_binder | probability this compound binds this target | classification | 68 | AUROC vs measured non-binders |  | 0.61 | sigmoid, prefit | yes | 2026-08-13 21:04 |
| GluN2B_binder | probability this compound binds this target | classification | 3418 | AUROC vs measured non-binders |  | 0.764 | sigmoid, prefit | yes | 2026-08-13 20:50 |
| H3_binder | probability this compound binds this target | classification | 12233 | AUROC vs measured non-binders |  | 0.966 | sigmoid, prefit | yes | 2026-08-13 20:38 |
| HDAC1_binder | probability this compound binds this target | classification | 11663 | AUROC vs measured non-binders |  | 0.961 | sigmoid, prefit | yes | 2026-08-13 20:48 |
| HDAC6_binder | probability this compound binds this target | classification | 14177 | AUROC vs measured non-binders |  | 0.972 | sigmoid, prefit | yes | 2026-08-13 20:49 |
| HT1A_binder | probability this compound binds this target | classification | 14179 | AUROC vs measured non-binders |  | 0.936 | sigmoid, prefit | yes | 2026-08-13 20:36 |
| HT2A_binder | probability this compound binds this target | classification | 14421 | AUROC vs measured non-binders |  | 0.93 | sigmoid, prefit | yes | 2026-08-13 20:57 |
| HT6_binder | probability this compound binds this target | classification | 10548 | AUROC vs measured non-binders |  | 0.959 | sigmoid, prefit | yes | 2026-08-13 20:37 |
| HT7_binder | probability this compound binds this target | classification | 5614 | AUROC vs measured non-binders |  | 0.907 | sigmoid, prefit | yes | 2026-08-13 20:38 |
| KEAP1_binder | probability this compound binds this target | classification | 471 | AUROC vs measured non-binders |  | 0.904 | sigmoid, prefit | yes | 2026-08-13 20:53 |
| LRRK2_binder | probability this compound binds this target | classification | 4357 | AUROC vs measured non-binders |  | 0.985 | sigmoid, prefit | yes | 2026-08-13 20:45 |
| MT1_binder | probability this compound binds this target | classification | 2611 | AUROC vs measured non-binders |  | 0.867 | sigmoid, prefit | yes | 2026-08-13 20:52 |
| NET_binder | probability this compound binds this target | classification | 6012 | AUROC vs measured non-binders |  | 0.954 | sigmoid, prefit | yes | 2026-08-13 20:39 |
| NFKB1_binder | probability this compound binds this target | classification | 102 | AUROC vs measured non-binders |  | 0.392 | sigmoid, prefit | WITHDRAWN | 2026-08-17 21:50 |
| NLRP3_binder | probability this compound binds this target | classification | 849 | AUROC vs measured non-binders |  | 0.872 | sigmoid, prefit | yes | 2026-08-13 20:45 |
| NR3C1_binder | probability this compound binds this target | classification | 37 | AUROC vs measured non-binders |  | 0.479 | sigmoid, prefit | WITHDRAWN | 2026-08-17 21:50 |
| NRF2_binder | probability this compound binds this target | classification | 285 | AUROC vs measured non-binders |  | 0.539 | sigmoid, prefit | WITHDRAWN | 2026-08-17 21:50 |
| Nav1_1_binder | probability this compound binds this target | classification | 227 | AUROC vs measured non-binders |  | 0.918 | sigmoid, prefit | WITHDRAWN | 2026-08-13 21:04 |
| Nav1_5_binder | probability this compound binds this target | classification | 870 | AUROC vs measured non-binders |  | 0.921 | sigmoid, prefit | yes | 2026-08-13 21:04 |
| Nav1_6_binder | probability this compound binds this target | classification | 1013 | AUROC vs measured non-binders |  | 0.634 | sigmoid, prefit | yes | 2026-08-13 21:02 |
| Nav1_7_binder | probability this compound binds this target | classification | 10496 | AUROC vs measured non-binders |  | 0.96 | sigmoid, prefit | yes | 2026-08-13 21:03 |
| Nav1_8_binder | probability this compound binds this target | classification | 1070 | AUROC vs measured non-binders |  | 0.876 | sigmoid, prefit | yes | 2026-08-13 21:03 |
| OPRK1_binder | probability this compound binds this target | classification | 11767 | AUROC vs measured non-binders |  | 0.934 | sigmoid, prefit | yes | 2026-08-13 20:41 |
| OPRM1_binder | probability this compound binds this target | classification | 14134 | AUROC vs measured non-binders |  | 0.96 | sigmoid, prefit | yes | 2026-08-13 20:43 |
| OX1_binder | probability this compound binds this target | classification | 12655 | AUROC vs measured non-binders |  | 0.955 | sigmoid, prefit | yes | 2026-08-13 20:51 |
| OX2_binder | probability this compound binds this target | classification | 14823 | AUROC vs measured non-binders |  | 0.954 | sigmoid, prefit | yes | 2026-08-13 20:52 |
| P2X7_binder | probability this compound binds this target | classification | 11191 | AUROC vs measured non-binders |  | 0.76 | sigmoid, prefit | yes | 2026-08-13 20:46 |
| PDE10A_binder | probability this compound binds this target | classification | 15570 | AUROC vs measured non-binders |  | 0.967 | sigmoid, prefit | yes | 2026-08-13 20:48 |
| PDE4B_binder | probability this compound binds this target | classification | 4308 | AUROC vs measured non-binders |  | 0.973 | sigmoid, prefit | yes | 2026-08-13 20:54 |
| RIPK1_binder | probability this compound binds this target | classification | 5626 | AUROC vs measured non-binders |  | 0.983 | sigmoid, prefit | yes | 2026-08-13 21:03 |
| SERT_binder | probability this compound binds this target | classification | 11164 | AUROC vs measured non-binders |  | 0.936 | sigmoid, prefit | yes | 2026-08-13 20:58 |
| SIRT1_binder | probability this compound binds this target | classification | 276 | AUROC vs measured non-binders |  | 0.795 | sigmoid, prefit | yes | 2026-08-13 21:03 |
| Sigma1_binder | probability this compound binds this target | classification | 7298 | AUROC vs measured non-binders |  | 0.81 | sigmoid, prefit | yes | 2026-08-13 20:40 |
| TAAR1_binder | probability this compound binds this target | classification | 82 | AUROC vs measured non-binders |  | 0.725 | sigmoid, prefit | yes | 2026-08-13 21:04 |
| a3b4nAChR_binder | probability this compound binds this target | classification | 750 | AUROC vs measured non-binders |  | 0.973 | sigmoid, prefit | yes | 2026-08-13 21:03 |
| a4b2nAChR_binder | probability this compound binds this target | classification | 1817 | AUROC vs measured non-binders |  | 0.937 | sigmoid, prefit | yes | 2026-08-13 21:03 |
| a7nAChR_binder | probability this compound binds this target | classification | 1277 | AUROC vs measured non-binders |  | 0.896 | sigmoid, prefit | yes | 2026-08-13 20:44 |
| mGluR5_binder | probability this compound binds this target | classification | 4554 | AUROC vs measured non-binders |  | 0.924 | sigmoid, prefit | yes | 2026-08-13 20:50 |
| mTOR_binder | probability this compound binds this target | classification | 11382 | AUROC vs measured non-binders |  | 0.992 | sigmoid, prefit | yes | 2026-08-13 20:53 |
| adme_caco2_permeability | measured value | regression | 897 | R2 | 0.7342 | 0.5837 | none | yes | 2026-08-13 22:12 |
| adme_clearance_hepatocyte | measured value | regression | 1020 | R2 | 0.2302 | 0.2062 | none | yes | 2026-08-13 22:13 |
| adme_kpuu | measured value | regression | 566 | R2 | 0.4044 | 0.3525 | none | yes | 2026-08-13 22:14 |
| adme_lipophilicity | measured value | regression | 4200 | R2 | 0.6389 | 0.5659 | none | yes | 2026-08-13 22:12 |
| adme_logbb | measured value | regression | 1058 | R2 | 0.5774 | 0.4235 | none | yes | 2026-08-13 22:14 |
| adme_pgp_inhibition | probability | classification | 1212 | AUROC | 0.9549 | 0.9346 | none | yes | 2026-08-13 22:06 |
| adme_pgp_substrate | probability | classification | 1371 | AUROC | 0.8576 | 0.8059 | none | yes | 2026-08-13 22:07 |
| adme_plasma_protein_binding | measured value | regression | 1797 | R2 | 0.4336 | 0.3645 | none | yes | 2026-08-13 22:13 |
| adme_solubility | measured value | regression | 9573 | R2 | 0.8036 | 0.7279 | none | yes | 2026-08-13 22:10 |