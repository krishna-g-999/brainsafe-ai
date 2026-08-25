# BrainSafe AI: exposure-gated prediction of small-molecule target engagement in the brain

**Authors.** Krishnasalini Gunanathan, Raghunatha Sarma, Sai Shyam, Ramya E. M.,
Venketesh Sivaramakrishnan. Sri Sathya Sai Institute of Higher Learning (SSSIHL), Puttaparthi,
Andhra Pradesh 515134, India. Corresponding author: Prof. Venketesh Sivaramakrishnan,
svenketesh@sssihl.edu.in.

**Server.** https://huggingface.co/spaces/Krishnag999/brainsafe-ai (live, no registration, MIT licence)

**Keywords.** target prediction; blood-brain barrier; ADMET; applicability domain

---

**Input.** A single SMILES string, or a compound name resolved to structure through PubChem. A batch
mode accepts up to 300 structures. No registration, no login, and no email address is requested.
Submitted structures are processed in memory and are not stored.

**Output.** For one submission the server returns, in a few seconds: engagement of 54 molecular
targets, each as a calibrated probability, with conformal prediction intervals on the eight core
classifiers; predicted blood-brain barrier penetration; nine ADME and exposure endpoints, among them
a directly modelled unbound brain-to-plasma ratio and a cardiac hERG liability; and the distance to
the nearest measured analogue, with that analogue shown.
Engaged targets are traced through a curated pathway graph anchored to KEGG, Reactome and IUPHAR to
the conditions those mechanisms touch. Output includes an interactive mechanism map, per-endpoint
tables, and downloadable CSV, JSON and a self-contained HTML report.

**Method.** 75 estimators were trained and 70 deployed, on 228,200 measured compound-endpoint records
from ChEMBL, BindingDB, B3DB, Therapeutics Data Commons and MoleculeNet, covering 169,341 compounds
keyed by InChIKey. Each endpoint is fitted on its own measured set alone, of median 3,789 rows.
Compounds are represented as a 1,024-bit ECFP-4 fingerprint with twelve physicochemical descriptors,
and a random forest is fitted per endpoint. That choice was made on a like-for-like comparison over
the thirteen core endpoints against XGBoost, histogram gradient boosting, L2 logistic regression and
a nearest-neighbour read-across, and against a graph neural network on four of them, which the
random forest won on all four. On the scaffold split the forest is best on seven of the eight
classification endpoints and on none of the five regressions, where gradient boosting scores
higher; it was deployed for stability under hyperparameters, for not extrapolating at the edge
of the applicability domain, and because TreeSHAP is exact for it.

Three choices distinguish the method. The negative class is recovered from censored measurements
rather than simulated with decoys, which returns experimentally tested non-binders to 57 endpoints.
Decision thresholds are measured on a background pool disjoint from the one that set them, so the
reported false-positive rate is an observation rather than the quantile that defined it. Target
scores are admitted only in proportion to predicted barrier penetration, so potency at a target a
compound cannot reach contributes nothing.

**Validation.** Scaffold-grouped 10-fold cross-validation gives a mean AUROC of 0.925 (0.878 to
0.965) across the measured-label classifiers. The binder panel, validated against compounds tested at
the same target and found inactive rather than against decoys, reaches a mean AUROC of 0.917 (0.719
to 0.985). Expected calibration error falls from 0.080 to 0.015 after isotonic calibration. On 306
FDA-curated approved drugs absent from the training source the barrier model reaches AUROC 0.764,
and 0.793 on the 241 also distinguishable from training in feature space. On 1,000 compounds with no
recorded activity at any modelled target the server stays silent 94.9% of the time; those compounds
are presumed rather than proven inactive, so that figure is a lower bound. Nine falsification
hypotheses were tested and four were refuted; all are reported, including that the curated pathway
edge weights carry no predictive information beyond the graph topology.

**Scope.** The server is organ-directed rather than disease-directed: its unifying constraint is the
blood-brain barrier, a pharmacokinetic property, and its targets span neurodegeneration, psychiatry,
analgesia, sleep and neuroinflammation rather than one disease group. Its core outputs, target
engagement and brain exposure, are the same class of quantity as ADMETlab (PMID 38572755) and
SwissTargetPrediction (PMID 31106366), both published in this issue. We would welcome the editor's
view on whether this places it within scope.

**Previous publications.** None. The method and the server have not been published, in NAR or
elsewhere, and no PubMed ID exists for either.

**This is not an update** of an application published previously in the Web Server Issue or in any
other venue.

**This is not a resubmission.** No proposal or manuscript for this server has been submitted to, or
rejected by, the Web Server Issue in any previous year.
