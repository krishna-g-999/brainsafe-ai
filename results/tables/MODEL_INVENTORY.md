# Model inventory, 2026-08-19

Commit `120f2c7`. One row per deployed estimator, taken from the estimators and their metadata on disk. Four receptors carry both a potency regression and a binder classifier, so the model count exceeds the endpoint count.

**75 estimators, 70 deployed.** Fitted between 2026-08-04 10:01 and 2026-08-19 17:54.

| family | estimators | deployed |
|---|---|---|
| binder | 52 | 47 |
| exposure | 10 | 10 |
| safety | 1 | 1 |
| target | 12 | 12 |

## Every estimator

| model | predicts | task | n train | metric | random | scaffold | calibration | deployed | fitted |
|---|---|---|---|---|---|---|---|---|---|
| A2A | potency, pChEMBL | regression | 6743 | R2 | 0.7231 | 0.6205 | none | yes | 2026-08-19 13:56 |
| AChE | probability of activity | classification | 5125 | AUROC | 0.9659 | 0.9212 | isotonic, out-of-fold | yes | 2026-08-19 13:48 |
| antioxidant_DPPH | potency, pChEMBL | regression | 2782 | R2 | 0.6589 | 0.4153 | none | yes | 2026-08-19 14:01 |
| BACE1 | probability of activity | classification | 8207 | AUROC | 0.9764 | 0.9648 | isotonic, out-of-fold | yes | 2026-08-19 13:49 |
| BBB | probability of activity | classification | 3901 | AUROC | 0.899 | 0.8777 | isotonic, out-of-fold | yes | 2026-08-19 13:47 |
| BChE | probability of activity | classification | 3278 | AUROC | 0.9724 | 0.9451 | isotonic, out-of-fold | yes | 2026-08-19 13:48 |
| D2 | potency, pChEMBL | regression | 7905 | R2 | 0.6403 | 0.5311 | none | yes | 2026-08-19 13:54 |
| GSK3B | probability of activity | classification | 5439 | AUROC | 0.9649 | 0.9425 | isotonic, out-of-fold | yes | 2026-08-19 13:49 |
| hERG | probability of activity | classification | 9933 | AUROC | 0.9565 | 0.927 | isotonic, out-of-fold | yes | 2026-08-19 13:51 |
| HT2A | potency, pChEMBL | regression | 6075 | R2 | 0.6996 | 0.556 | none | yes | 2026-08-19 13:58 |
| MAO_A | probability of activity | classification | 3585 | AUROC | 0.9619 | 0.9059 | isotonic, out-of-fold | yes | 2026-08-19 13:49 |
| MAO_B | probability of activity | classification | 4534 | AUROC | 0.963 | 0.917 | isotonic, out-of-fold | yes | 2026-08-19 13:50 |
| pka_basic | pKa | regression | 6384 | R2 |  |  | none | yes | 2026-08-04 10:01 |
| SERT | potency, pChEMBL | regression | 4479 | R2 | 0.6897 | 0.4612 | none | yes | 2026-08-19 14:00 |
| A1_binder | probability this compound binds this target | classification | 7352 | AUROC vs measured non-binders |  | 0.914 | sigmoid, prefit | yes | 2026-08-19 14:22 |
| A2A_binder | probability this compound binds this target | classification | 15636 | AUROC vs measured non-binders |  | 0.949 | sigmoid, prefit | yes | 2026-08-19 14:31 |
| CB1_binder | probability this compound binds this target | classification | 10198 | AUROC vs measured non-binders |  | 0.949 | sigmoid, prefit | yes | 2026-08-19 14:19 |
| CGRP_binder | probability this compound binds this target | classification | 2251 | AUROC vs measured non-binders |  | 0.985 | sigmoid, prefit | yes | 2026-08-19 17:52 |
| COX2_binder | probability this compound binds this target | classification | 4429 | AUROC vs measured non-binders |  | 0.783 | sigmoid, prefit | yes | 2026-08-19 14:23 |
| CSF1R_binder | probability this compound binds this target | classification | 6772 | AUROC vs measured non-binders |  | 0.95 | sigmoid, prefit | yes | 2026-08-19 14:24 |
| Cav3_2_binder | probability this compound binds this target | classification | 1325 | AUROC vs measured non-binders |  | 0.982 | sigmoid, prefit | yes | 2026-08-19 17:52 |
| D2_binder | probability this compound binds this target | classification | 13861 | AUROC vs measured non-binders |  | 0.872 | sigmoid, prefit | yes | 2026-08-19 14:30 |
| D3_binder | probability this compound binds this target | classification | 15831 | AUROC vs measured non-binders |  | 0.956 | sigmoid, prefit | yes | 2026-08-19 14:21 |
| DAT_binder | probability this compound binds this target | classification | 5113 | AUROC vs measured non-binders |  | 0.959 | sigmoid, prefit | yes | 2026-08-19 14:17 |
| DHODH_binder | probability this compound binds this target | classification | 3619 | AUROC vs measured non-binders |  | 0.966 | sigmoid, prefit | yes | 2026-08-19 17:52 |
| GABA_A_binder | probability this compound binds this target | classification | 771 | AUROC vs measured non-binders |  | 0.719 | sigmoid, prefit | yes | 2026-08-19 14:27 |
| GBA1_binder | probability this compound binds this target | classification | 441 | AUROC vs measured non-binders |  | 0.932 | sigmoid, prefit | yes | 2026-08-19 14:29 |
| GluA2_binder | probability this compound binds this target | classification | 68 | AUROC vs measured non-binders |  | 0.696 | sigmoid, prefit | WITHDRAWN | 2026-08-19 17:54 |
| GluN2B_binder | probability this compound binds this target | classification | 3461 | AUROC vs measured non-binders |  | 0.761 | sigmoid, prefit | yes | 2026-08-19 14:26 |
| H3_binder | probability this compound binds this target | classification | 12295 | AUROC vs measured non-binders |  | 0.978 | sigmoid, prefit | yes | 2026-08-19 14:17 |
| HDAC1_binder | probability this compound binds this target | classification | 11728 | AUROC vs measured non-binders |  | 0.96 | sigmoid, prefit | yes | 2026-08-19 14:25 |
| HDAC6_binder | probability this compound binds this target | classification | 14155 | AUROC vs measured non-binders |  | 0.975 | sigmoid, prefit | yes | 2026-08-19 14:26 |
| HT1A_binder | probability this compound binds this target | classification | 14179 | AUROC vs measured non-binders |  | 0.936 | sigmoid, prefit | yes | 2026-08-19 14:15 |
| HT2A_binder | probability this compound binds this target | classification | 14369 | AUROC vs measured non-binders |  | 0.925 | sigmoid, prefit | yes | 2026-08-19 14:32 |
| HT6_binder | probability this compound binds this target | classification | 10532 | AUROC vs measured non-binders |  | 0.96 | sigmoid, prefit | yes | 2026-08-19 14:16 |
| HT7_binder | probability this compound binds this target | classification | 5623 | AUROC vs measured non-binders |  | 0.947 | sigmoid, prefit | yes | 2026-08-19 14:16 |
| KEAP1_binder | probability this compound binds this target | classification | 463 | AUROC vs measured non-binders |  | 0.88 | sigmoid, prefit | yes | 2026-08-19 14:29 |
| LRRK2_binder | probability this compound binds this target | classification | 4470 | AUROC vs measured non-binders |  | 0.968 | sigmoid, prefit | yes | 2026-08-19 14:22 |
| MT1_binder | probability this compound binds this target | classification | 2670 | AUROC vs measured non-binders |  | 0.896 | sigmoid, prefit | yes | 2026-08-19 14:28 |
| NET_binder | probability this compound binds this target | classification | 5933 | AUROC vs measured non-binders |  | 0.92 | sigmoid, prefit | yes | 2026-08-19 14:18 |
| NFKB1_binder | probability this compound binds this target | classification | 102 | AUROC vs measured non-binders |  | 0.459 | sigmoid, prefit | WITHDRAWN | 2026-08-19 17:54 |
| NLRP3_binder | probability this compound binds this target | classification | 843 | AUROC vs measured non-binders |  | 0.93 | sigmoid, prefit | yes | 2026-08-19 14:22 |
| NR3C1_binder | probability this compound binds this target | classification | 37 | AUROC vs measured non-binders |  | 0.41 | sigmoid, prefit | WITHDRAWN | 2026-08-19 17:54 |
| NRF2_binder | probability this compound binds this target | classification | 285 | AUROC vs measured non-binders |  | 0.789 | sigmoid, prefit | WITHDRAWN | 2026-08-19 17:54 |
| Nav1_1_binder | probability this compound binds this target | classification | 227 | AUROC vs measured non-binders |  | 0.952 | sigmoid, prefit | WITHDRAWN | 2026-08-19 17:54 |
| Nav1_5_binder | probability this compound binds this target | classification | 870 | AUROC vs measured non-binders |  | 0.921 | sigmoid, prefit | yes | 2026-08-19 14:33 |
| Nav1_6_binder | probability this compound binds this target | classification | 1027 | AUROC vs measured non-binders |  | 0.862 | sigmoid, prefit | yes | 2026-08-19 17:52 |
| Nav1_7_binder | probability this compound binds this target | classification | 10397 | AUROC vs measured non-binders |  | 0.957 | sigmoid, prefit | yes | 2026-08-19 17:53 |
| Nav1_8_binder | probability this compound binds this target | classification | 1030 | AUROC vs measured non-binders |  | 0.956 | sigmoid, prefit | yes | 2026-08-19 17:53 |
| OPRK1_binder | probability this compound binds this target | classification | 11642 | AUROC vs measured non-binders |  | 0.945 | sigmoid, prefit | yes | 2026-08-19 14:20 |
| OPRM1_binder | probability this compound binds this target | classification | 14058 | AUROC vs measured non-binders |  | 0.954 | sigmoid, prefit | yes | 2026-08-19 14:20 |
| OX1_binder | probability this compound binds this target | classification | 12724 | AUROC vs measured non-binders |  | 0.964 | sigmoid, prefit | yes | 2026-08-19 14:27 |
| OX2_binder | probability this compound binds this target | classification | 14679 | AUROC vs measured non-binders |  | 0.964 | sigmoid, prefit | yes | 2026-08-19 14:28 |
| P2X7_binder | probability this compound binds this target | classification | 11122 | AUROC vs measured non-binders |  | 0.813 | sigmoid, prefit | yes | 2026-08-19 14:23 |
| PDE10A_binder | probability this compound binds this target | classification | 15463 | AUROC vs measured non-binders |  | 0.962 | sigmoid, prefit | yes | 2026-08-19 14:24 |
| PDE4B_binder | probability this compound binds this target | classification | 4332 | AUROC vs measured non-binders |  | 0.966 | sigmoid, prefit | yes | 2026-08-19 14:30 |
| RIPK1_binder | probability this compound binds this target | classification | 5619 | AUROC vs measured non-binders |  | 0.966 | sigmoid, prefit | yes | 2026-08-19 17:54 |
| SERT_binder | probability this compound binds this target | classification | 11429 | AUROC vs measured non-binders |  | 0.945 | sigmoid, prefit | yes | 2026-08-19 14:33 |
| SIRT1_binder | probability this compound binds this target | classification | 276 | AUROC vs measured non-binders |  | 0.792 | sigmoid, prefit | yes | 2026-08-19 14:33 |
| Sigma1_binder | probability this compound binds this target | classification | 7427 | AUROC vs measured non-binders |  | 0.881 | sigmoid, prefit | yes | 2026-08-19 14:18 |
| TAAR1_binder | probability this compound binds this target | classification | 82 | AUROC vs measured non-binders |  | 0.78 | sigmoid, prefit | yes | 2026-08-19 17:54 |
| a3b4nAChR_binder | probability this compound binds this target | classification | 760 | AUROC vs measured non-binders |  | 0.974 | sigmoid, prefit | yes | 2026-08-19 17:54 |
| a4b2nAChR_binder | probability this compound binds this target | classification | 1835 | AUROC vs measured non-binders |  | 0.923 | sigmoid, prefit | yes | 2026-08-19 17:54 |
| a7nAChR_binder | probability this compound binds this target | classification | 1285 | AUROC vs measured non-binders |  | 0.763 | sigmoid, prefit | yes | 2026-08-19 14:22 |
| mGluR5_binder | probability this compound binds this target | classification | 4511 | AUROC vs measured non-binders |  | 0.893 | sigmoid, prefit | yes | 2026-08-19 14:26 |
| mTOR_binder | probability this compound binds this target | classification | 11510 | AUROC vs measured non-binders |  | 0.983 | sigmoid, prefit | yes | 2026-08-19 14:29 |
| adme_caco2_permeability | measured value | regression | 897 | R2 | 0.7359 | 0.5818 | none | yes | 2026-08-19 14:38 |
| adme_clearance_hepatocyte | measured value | regression | 1020 | R2 | 0.2302 | 0.2062 | none | yes | 2026-08-19 14:39 |
| adme_kpuu | measured value | regression | 566 | R2 | 0.4056 | 0.3523 | none | yes | 2026-08-19 14:39 |
| adme_lipophilicity | measured value | regression | 4200 | R2 | 0.6389 | 0.5659 | none | yes | 2026-08-19 14:37 |
| adme_logbb | measured value | regression | 1058 | R2 | 0.5836 | 0.4131 | none | yes | 2026-08-19 14:39 |
| adme_pgp_inhibition | probability | classification | 1212 | AUROC | 0.9549 | 0.9346 | none | yes | 2026-08-19 14:33 |
| adme_pgp_substrate | probability | classification | 1371 | AUROC | 0.8561 | 0.807 | none | yes | 2026-08-19 14:34 |
| adme_plasma_protein_binding | measured value | regression | 1797 | R2 | 0.4336 | 0.3645 | none | yes | 2026-08-19 14:38 |
| adme_solubility | measured value | regression | 9573 | R2 | 0.8008 | 0.7263 | none | yes | 2026-08-19 14:36 |