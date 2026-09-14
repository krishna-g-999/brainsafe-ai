# BrainSafe AI: a calibrated, exposure-gated web server for multi-endpoint prediction of small-molecule action in the human brain

**Authors:** Krishnasalini Gunanathan¹, Raghunatha Sarma¹, Sai Shyam¹, Ramya E. M.¹,
Venketesh Sivaramakrishnan¹

¹Sri Sathya Sai Institute of Higher Learning (SSSIHL), Puttaparthi, Andhra Pradesh 515134, India

**Correspondence:** Prof. Venketesh Sivaramakrishnan, svenketesh@sssihl.edu.in

**Manuscript type:** NAR Web Server Issue.

**Graphical abstract:** submitted as a separate file, `manuscript/figures/GraphicalAbstract.png`
(and `.pdf`).

**Key Points**

- BrainSafe AI profiles small-molecule mechanism in the human brain from structure alone, gating
  every target score by predicted exposure so potency at a target a compound cannot reach
  contributes nothing to the output.
- Its 74 cross-validated estimators are each validated under both random and scaffold-grouped
  splits and carry a calibrated probability, the core classifiers additionally reporting a
  conformal interval.
- The negative class is recovered from compounds measured and found inactive rather than simulated
  with decoys, and every validation is reported whichever way it falls, which led to withdrawing
  two endpoints that could not separate a real ligand from an unrelated metabolite.

## Abstract

BrainSafe AI predicts, from structure alone, how a small molecule may act on the human brain. For a
submitted SMILES string or compound name it returns engagement of 54 molecular targets spanning the
principal neurodegenerative, psychiatric, neuroinflammatory, analgesic and sleep-related mechanisms,
predicted blood-brain barrier penetration, nine ADME and exposure endpoints including a directly
modelled unbound brain-to-plasma ratio, and two cardiac safety liabilities. Every target score is
admitted only in proportion to predicted brain exposure, so potency at a target a compound cannot
reach contributes nothing, and engaged targets are traced through a curated pathway graph to the
conditions those mechanisms touch. The server is built on 75 estimators, 70 deployed, trained on
228,200 measured compound-endpoint records from ChEMBL (1), BindingDB (2) and B3DB
(3). Under scaffold-grouped 10-fold cross-validation the measured-label classifiers reach a mean
AUROC of 0.925 (0.958 under a random split), and expected calibration error falls from 0.0801 to
0.0147 after isotonic calibration (4). The binder panel is validated against compounds
tested at the same target and found inactive rather than against decoys, reaching a mean AUROC of
0.917. Every prediction carries a calibrated probability, a conformal interval, and an
applicability-domain distance to the nearest measured analogue, and the server reports silence
rather than a guess for compounds outside its competence: on non-CNS chemistry its specificity is
0.925 (95% CI 0.907 to 0.940). Disease-level scores are presented as a route from a mechanism to the
conditions it touches, not as an indication prediction. BrainSafe AI is freely available without
registration at https://huggingface.co/spaces/Krishnag999/brainsafe-ai, with source code, trained
models and every validation artefact at https://github.com/krishna-g-999/brainsafe-ai.

## Introduction

Central nervous system drug discovery fails more often, and later, than any other therapeutic area.
Two reasons dominate: a compound potent at its target may never reach the brain, and a compound that
reaches the brain may engage more than the target it was designed for. Existing servers address these
questions separately. ADMET platforms predict barrier penetration without saying what the compound
would do on arrival; target-prediction servers rank probable targets without asking whether the
compound can reach them.

BrainSafe AI answers both in one pass, and couples them: a target score is admitted only in
proportion to predicted exposure (Figure 1). It reports what it cannot do as prominently as what it
can. Five endpoints were trained, tested and withdrawn because no threshold separated real ligands
from trivial chemistry; every validation is written so that it could fail and is reported whichever
way it falls; and the recall achieved on chemistry distant from the training set is reported beside
every result rather than left to be discovered.

## Materials and Methods

**Training data.** Labels are measured experimental values only, never qualitative annotation.
Potency data are ChEMBL pChEMBL values augmented with BindingDB affinities pooled at compound level;
blood-brain barrier labels come from B3DB augmented with FDA-curated approved drugs; the nine ADME
endpoints use measured sets from Therapeutics Data Commons (5), MoleculeNet (6), B3DB
and ChEMBL. The panel
holds 228,200 measured compound-endpoint records over 169,341 unique compounds keyed by the InChIKey
of the desalted parent. Each endpoint is trained on its own measured set alone; across the deployed
panel those sets hold a median of 3,789 rows and span 387 compounds (KEAP1) to 10,276 (hERG).
That set is every deployed model owning a table in data/endpoints, the barrier model included;
excluding it, since it is an exposure rather than a target endpoint, gives 54 tables and a
median of 3,587.

**Recovery of the negative class.** A compound assayed and found inactive is often deposited only as
a censored bound, and the conventional query discards exactly those rows, leaving a training set of
actives against synthetic decoys. A bound settles a label when the whole interval it defines falls on
one side of the activity cut, and is discarded as undecidable when it spans both. Recovering these
returned experimentally tested non-binders to 57 endpoints.

**Representation.** Each compound is a 1,036-column vector: a 1,024-bit folded ECFP-4 fingerprint
(7) and twelve physicochemical descriptors. Structures are reduced to the largest organic fragment and
neutralised. Neutralisation is part of the representation rather than a detail of it, because a drug
and its salt must give the same answer: removing a counter-ion without it leaves the parent carrying
the salt's charge, and haloperidol hydrochloride then scored a barrier probability of 0.613 against
0.993 for the free base. A permanent charge is retained, because a quaternary ammonium's charge is
precisely what prevents it crossing. Chirality is excluded, so two enantiomers give identical rows;
rows identical in feature space are collapsed before any split. Of 228,198 training structures, 40.8
per cent carry a stereocentre, but where one skeleton appears as several stereoisomers measured at
the same endpoint the labels agree in 94.6 per cent of 8,013 cases, so the share of the panel where
chirality could change a class call is 0.19 per cent.

**Models.** A random forest (8) is fitted per endpoint. That choice was made on a
like-for-like comparison over thirteen core endpoints against XGBoost (9), histogram gradient
boosting, L2 logistic regression and a nearest-neighbour read-across, and against a graph neural
network on four of them, which the forest won on all four. Under the scaffold split the forest is
best on seven of eight classification endpoints, losing AChE to histogram gradient boosting, and on
none of the five regressions, where boosting scores higher. It was deployed for its stability under
hyperparameters, for not extrapolating beyond the training range, and because TreeSHAP (10)
is exact for it rather than approximate. Classifiers are isotonically calibrated (4) on
out-of-fold predictions, so no compound contributes to the calibrator that scores it. Binder models
use Platt scaling (11), the withheld set for one target often being too small to fit a step
function.

**Thresholds.** The background library is partitioned into three disjoint pools by a stable hash of
the structure: one supplies property-matched decoys (12), one sets thresholds, one measures the
false-positive rate. Choosing a threshold as a quantile of a sample and then measuring the rate on
that same sample restates the target rather than testing it.

**Exposure gating and the disease layer.** A target score is admitted in proportion to predicted
barrier penetration. Engaged targets are traced through a curated graph anchored to KEGG (13),
Reactome (14) and IUPHAR (15) to the conditions they touch. The graph's edge weights were
ablated and carry no predictive information beyond the topology (curated 0.7901, uniform 0.7897,
permuted 0.7874), so they are reported as structure rather than as tuned parameters.

## Results

**Cross-validation.** The panel and its per-endpoint performance are shown in Figure 2. Under random
10-fold cross-validation the measured-label classifiers reach a mean AUROC of 0.958 (0.899 to
0.976); under a scaffold-grouped split that withholds entire structural classes, 0.925 (0.878 to
0.965). Expected calibration error falls from 0.0801 to 0.0147 after isotonic calibration, and
conformal prediction (16) on the eight core classifiers, on the deduplicated matrix the
classifiers are trained on, achieves empirical coverage of 0.876 to 0.933 against a 0.90 target,
with mean set size from 0.956 to 1.215 on a two-class problem where 1.0 is a confident single label
(Figure 3A).

**The binder panel.** The 52 binder classifiers are validated not against the decoys used to train
them but against compounds experimentally tested at the same target and found inactive. Across the 47
deployed they reach a mean AUROC of 0.917 and a mean sensitivity of 0.764, both on actives withheld
by scaffold. Both are means over 47 endpoints and the spread is wide: AUROC ranges from 0.719 at
GABA-A to 0.985 at CGRP, sensitivity from 0.303 at GABA-A to 0.993 at CGRP with a median of 0.835.

That sensitivity figure carries a correction we report rather than leave for a reader to find. Four
scripts write it in sequence, and the last scored every active in the endpoint table, so roughly four
fifths of the scoring set were compounds the model had been fitted on, while the field's own basis
annotation went on reading that the measurement was held out. The panel mean was published as 0.898
against a held-out 0.764, an inflation of 0.134. The registry has been recomputed from the
scaffold-held-out partition, it agrees with the threshold table on all 47 deployed endpoints, and a
regression test fails if any endpoint reports sensitivity on an undeclared basis. Six endpoints fall
below the reliability gate on the corrected figure, against one before; each holds its background
false-positive rate at or below target, so each is weak rather than misleading and stays deployed
with a low-power marker on any negative call. Five further endpoints are withdrawn outright for
firing on trivial metabolites at every usable threshold.

**Leakage and null models.** On the deduplicated matrix the pipeline fits, no InChIKey, no feature
vector and no scaffold appears on both sides of any fold. With labels permuted the same pipeline
returns a mean AUROC of 0.496 (random) and 0.503 (scaffold) over the eight core classifiers, every
one of the sixteen values within 0.02 of chance and a smallest margin over its own null 19-fold
larger than the largest such departure.

**External and prospective validation.** The barrier model was tested on FDA-curated approved drugs
absent from B3DB by InChIKey: AUROC 0.764 on all 306, and 0.767 on the 227 also distinguishable from
training in feature space, once compounds the featuriser cannot tell apart from a training row are
also excluded. For the target panel no external set of comparable size exists, because
for most of these targets the public measured chemistry is the training set. Two kinds of
independence were therefore constructed, with every model refitted rather than scored. By date: each
endpoint refitted on its pre-cutoff rows with its decision threshold also frozen before the cutoff,
and tested on compounds first published afterwards; 39 of 47 deployed endpoints qualify, giving
45,244 test compounds. By curator: compounds deposited in BindingDB and absent from ChEMBL, withheld
entirely (Figure 5).

Read in aggregate the time split suggests prospective decay, with mean AUROC 0.823 against 0.951 for
a size-matched random control. It is not decay. The false-positive rate on background chemistry is
unchanged, and the gap closes once test compounds are stratified by maximum Tanimoto similarity to
the training actives. A random split of medicinal-chemistry data holds out 83 per cent close
analogues of its own training set, because the published record is series; a time split holds out 28
per cent chemistry below Tanimoto 0.40. Three test sets built by unrelated rules trace one recall
curve: 0.16, 0.55, 0.74 and 0.86 by date across four novelty bands, against 0.12, 0.52, 0.77 and 0.93
at random and 0.05, 0.46, 0.83 and 0.90 by curator. Recall is therefore a function of chemical
distance rather than of publication date, and the expected recall for a submitted compound is
reported at query time.

**Specificity.** One thousand compounds with no recorded activity at any modelled target were scored
through the deployed pipeline; 925 returned no actionable disease signal, a specificity of 0.925 (95%
CI 0.907 to 0.940). These compounds are presumed inactive because nothing is recorded about them, not
proven inactive, so this is a lower bound.

**Adversarial checks and falsification.** Six checks were written so that each could fail, and all
six pass. The domain-flag check initially failed and passes only after its control set was corrected:
28 of the original controls are measured compounds inside the flag's own reference library, where
calling them in domain is truthful rather than a failure. The passing criterion was not moved. A
further falsification analysis found that the curated edge weights add nothing beyond topology
(curated 0.7901, uniform 0.7897, permuted 0.7874); that silence at a target reflects the operating
point rather than a non-discriminative model, every target separating its own actives at AUROC 0.91
or better; and that engaged targets are not independent observations, 37 firing targets spanning
only 16 independent directions. It also found a deployed endpoint calling common metabolites
binders at its calibrated threshold, and that endpoint is withdrawn.

**Use case.** A worked profile and the silence behaviour are shown in Figure 4. For donepezil the
server returns AChE engagement at 1.00 against a training base rate of 0.596, barrier penetration
0.991, and Alzheimer's disease as the top condition at 0.991 with AChE named as the driver. It also
returns a hERG probability of 0.734, an enrichment of 0.652 over a 0.236 base rate, so the compound
that scores highest on the efficacy axis carries a liability a medicinal chemist would need to see.
For atenolol, a beta-blocker optimised not to enter the brain, every target probability sits below
its base rate and the top condition scores exactly 0, well below the reporting threshold: the server
is correctly silent.

## Discussion

BrainSafe AI couples exposure and engagement in one pass, and reports the confidence of each. Its
distinguishing choices are that the negative class is recovered from measurement rather than
simulated, that thresholds are measured on a pool disjoint from the one that set them, and that
target scores are gated by predicted exposure.

Five limitations bound its use. First, the applicability-domain flag is a weak discriminator of
non-drug-like chemistry: in the adversarial check it scores genuinely absent chemistry at a median
maximum similarity of 0.47 against 0.57 for unseen approved drugs (n = 25, one-sided Mann-Whitney
p = 1.8e-03), but at a threshold rejecting a tenth of genuine drugs it catches only a fifth of
distant chemistry. What it
does predict well is sensitivity, the distance it measures being the variable recall tracks. Recall
on chemistry beyond Tanimoto 0.40 is also near 0.16, so a negative result on a novel scaffold is close
to uninformative, and the server reports the expected recall beside it. The specificity estimate
rests on compounds presumed rather than proven inactive. Terpenoid and steroidal natural products,
in turn, are largely outside the training library, whose median fraction-sp3 is 0.34, and such
compounds are flagged as outside the domain. Last, the disease layer does not predict indication:
27 of the 51 targets in the pathway graph drive more than one condition, and what selects among them,
dose, regimen and patient population, is not present in a structure.

The server does not distinguish an agonist from an antagonist. The training label is a potency value
measuring affinity, which an agonist and an antagonist at the same receptor can share, and ChEMBL's
action_type field appears nowhere in this project's data. The honest description of what the panel
predicts is engagement, not modulation.

## Data availability

BrainSafe AI is freely available without registration or login at
https://huggingface.co/spaces/Krishnag999/brainsafe-ai, served over HTTPS; no submitted structure is
retained beyond the lifetime of its request. Source code, trained models, every training table and
every validation artefact are at https://github.com/krishna-g-999/brainsafe-ai under the MIT
licence; underlying data retain their own sources' licences. Trained estimators and the raw API
responses are deposited separately with a manifest recording the SHA-256 of every file. **[TO BE
SUPPLIED BEFORE SUBMISSION]** the archive deposit's own DOI.

## Funding

[TO BE SUPPLIED]

## Conflict of interest

None declared.

## References

Every entry was resolved by a live query against CrossRef or Europe PMC and accepted only on a title match, or, where the identity is known and the registered title is a short form, by resolving the DOI and confirming the title and first author. The requested title, the matched title and the score are recorded in `manuscript/references_verified.json`, so the list can be re-checked mechanically. None is written from memory. A PubMed Central or PubMed abstract link is given where NCBI indexes the work (`manuscript/references_links.json`); neither exists for a work outside PubMed's coverage, which is expected for some conference proceedings and for the software citations.

1. Zdrazil B, Felix E, Hunter F et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Research. 2024. doi:10.1093/nar/gkad1004. https://pmc.ncbi.nlm.nih.gov/articles/PMC10767899/ https://pubmed.ncbi.nlm.nih.gov/37933841/
2. Gilson M, Liu T, Baitaluk M et al. BindingDB in 2015: A public database for medicinal chemistry, computational chemistry and systems pharmacology. Nucleic Acids Research. 2016. doi:10.1093/nar/gkv1072. https://pmc.ncbi.nlm.nih.gov/articles/PMC4702793/ https://pubmed.ncbi.nlm.nih.gov/26481362/
3. Meng F, Xi Y, Huang J et al. A curated diverse molecular database of blood-brain barrier permeability with chemical descriptors. Scientific Data. 2021. doi:10.1038/s41597-021-01069-5. https://pmc.ncbi.nlm.nih.gov/articles/PMC8556334/ https://pubmed.ncbi.nlm.nih.gov/34716354/
4. Niculescu-Mizil A, Caruana R. Predicting good probabilities with supervised learning. Proceedings of the 22nd international conference on Machine learning  - ICML '05. 2005. doi:10.1145/1102351.1102430.
5. Huang K, Fu T, Gao W, Zhao Y, Roohani Y, Leskovec J, Coley CW, Xiao C, Sun J, Zitnik M. Artificial intelligence foundation for therapeutic science. Nature chemical biology. 2022. doi:10.1038/s41589-022-01131-2. https://pmc.ncbi.nlm.nih.gov/articles/PMC9529840/ https://pubmed.ncbi.nlm.nih.gov/36131149/
6. Wu Z, Ramsundar B, Feinberg E et al. MoleculeNet: a benchmark for molecular machine learning. Chemical Science. 2018. doi:10.1039/c7sc02664a. https://pmc.ncbi.nlm.nih.gov/articles/PMC5868307/ https://pubmed.ncbi.nlm.nih.gov/29629118/
7. Rogers D, Hahn M. Extended-Connectivity Fingerprints. Journal of Chemical Information and Modeling. 2010. doi:10.1021/ci100050t. https://pubmed.ncbi.nlm.nih.gov/20426451/
8. Breiman L. Random Forests. Machine Learning. 2001. doi:10.1023/a:1010933404324.
9. Chen T, Guestrin C. XGBoost: A Scalable Tree Boosting System. Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. 2016. doi:10.1145/2939672.2939785.
10. Lundberg S, Erion G, Chen H et al. From local explanations to global understanding with explainable AI for trees. Nature Machine Intelligence. 2020. doi:10.1038/s42256-019-0138-9. https://pmc.ncbi.nlm.nih.gov/articles/PMC7326367/ https://pubmed.ncbi.nlm.nih.gov/32607472/
11. Lin H, Lin C, Weng R. A note on Platt’s probabilistic outputs for support vector machines. Machine Learning. 2007. doi:10.1007/s10994-007-5018-6.
12. Mysinger M, Carchia M, Irwin J et al. Directory of Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking. Journal of Medicinal Chemistry. 2012. doi:10.1021/jm300687e. https://pmc.ncbi.nlm.nih.gov/articles/PMC3405771/ https://pubmed.ncbi.nlm.nih.gov/22716043/
13. Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research. 2000. doi:10.1093/nar/28.1.27. https://pmc.ncbi.nlm.nih.gov/articles/PMC102409/ https://pubmed.ncbi.nlm.nih.gov/10592173/
14. Milacic M, Beavers D, Conley P, Gong C, Gillespie M, Griss J, Haw R, Jassal B, Matthews L, May B, Petryszak R, Ragueneau E, Rothfels K, Sevilla C, Shamovsky V, Stephan R, Tiwari K, Varusai T, Weiser J, Wright A, Wu G, Stein L, Hermjakob H, D'Eustachio P. The Reactome Pathway Knowledgebase 2024. Nucleic acids research. 2024. doi:10.1093/nar/gkad1025. https://pmc.ncbi.nlm.nih.gov/articles/PMC10767911/ https://pubmed.ncbi.nlm.nih.gov/37941124/
15. Harding SD, Armstrong JF, Faccenda E, Southan C, Alexander SPH, Davenport AP, Spedding M, Davies JA. The IUPHAR/BPS Guide to PHARMACOLOGY in 2024. Nucleic acids research. 2024. doi:10.1093/nar/gkad944. https://pmc.ncbi.nlm.nih.gov/articles/PMC10767925/ https://pubmed.ncbi.nlm.nih.gov/37897341/
16. Norinder U, Carlsson L, Boyer S et al. Introducing Conformal Prediction in Predictive Modeling. A Transparent and Flexible Alternative to Applicability Domain Determination. Journal of Chemical Information and Modeling. 2014. doi:10.1021/ci5001168. https://pubmed.ncbi.nlm.nih.gov/24797111/
17. Wilson E. Probable Inference, the Law of Succession, and Statistical Inference. Journal of the American Statistical Association. 1927. doi:10.1080/01621459.1927.10502953.



## Figure legends

![Figure 1](figures/Figure1_architecture.png)

**Figure 1.** How a query is answered. (**A**) A submitted structure is standardised and represented
as one fixed 1,036-column vector, scored by four model families, as MODEL_INVENTORY.csv defines
them: the 52-endpoint binder panel of which 47 are deployed, twelve target potency and activity
endpoints, the ten-endpoint exposure and ADME layer including the barrier model itself, and the
single hERG safety classifier. Every target score is admitted only in proportion to the predicted
probability that the compound reaches the brain, and surviving scores are ranked by enrichment over
each endpoint's base rate rather than by raw probability. Every reported value carries a calibrated
probability and an applicability-domain distance, and the eight core classifiers additionally report
a conformal interval. (**B**) The counts in (**A**) are trained
estimators, 75 in total, of which 70 are deployed. Each was preceded by twenty fits that never serve
a prediction and exist only to measure how the twenty-first behaves on withheld compounds, 1,480
across the panel.

![Figure 2](figures/Figure9_model_atlas.png)

**Figure 2.** The panel, one mark per estimator, so that no claim rests on a mean a reader cannot
check. (**A**) Every estimator, deployed or withdrawn, placed by the number of compounds it was
trained on
and by the performance claimed for it, coloured by model family. Marker shape carries the metric:
AUROC and R² both run to 1.0 and are not the same quantity, since 0.5 is chance for one and a
respectable fit for the other, so they are distinguished rather than averaged. The five estimators
withdrawn after specificity testing are drawn in outline, because a panel showing only what survived
is a selection rather than an inventory. Training sets span two orders of magnitude, from 37 to
15,831 rows; binder training sets include property-matched decoys while the others are measured
compounds only, which the axis states. (**B**) The same population by family, with the family median
marked. The spread is the point: the binder classifiers have a median of 0.945 while the exposure and
ADME family, which mixes AUROC and R² endpoints, has a median of 0.574, and a single panel average
would describe neither. The complete inventory, with training-set composition, validation scheme,
calibration and fitting date for every estimator, is Supplementary Table S1.

![Figure 3](figures/Figure6_validation.png)

**Figure 3.** Four validations that a cross-validated score cannot replace. (**A**) Expected
calibration error before and after isotonic regression fitted on out-of-fold predictions, so no
compound contributes to the calibrator that scores it. (**B**) Recall on whole scaffold classes
withheld before training, with 95 per cent Wilson intervals (17) and marker area proportional
to the number of withheld actives, so an interval that is wide because the evidence is thin looks
thin. (**C**) Specificity on chemistry the server should stay quiet about, and external
discrimination on approved drugs absent from the training source. (**D**) The adversarial suite, in
which each check was written so that it could fail. All six pass, one of them, the
applicability-domain flag, only after its controls were corrected rather than its criterion
retuned; every check is drawn at the same size whatever its verdict, so a future failure would be
exactly this visible.

![Figure 4](figures/Figure8_use_case.png)

**Figure 4.** A worked profile and the silence behaviour. For four approved CNS drugs the server
recovers the pharmacologically correct driving mechanism and the corresponding condition; for four
peripherally acting compounds no disease score reaches the reporting threshold. Bars are disease
scores after exposure gating, and the driving target is named beside each. The two behaviours are the
same design decision seen from opposite sides: a target score is admitted only in proportion to
predicted barrier penetration, so a compound that does not arrive cannot generate a disease call.

Supplementary figures, each named by the file it is generated into so that the number and the
artefact cannot come apart:

![Figure 5](figures/Figure11_external_validation.png)

**Figure 5.** External validation, and an apparent temporal decay that is a composition effect.
(**A**) Per endpoint, AUROC under a size-matched random split against AUROC under a time split that
withholds the most recent quarter of the data and freezes the decision threshold before the cutoff.
The size match matters: a time split trains on less data as well as none of the future, so without a
control at the same n a drop cannot be attributed to either. (**B**) The same comparison for
sensitivity at the frozen operating point, where the gap is larger. (**C**) Why the gap exists. A
random split of medicinal-chemistry data holds out mostly close analogues of its own training set,
because the published record is series; a time split does not. The two are not testing comparable
populations. (**D**) The resolution. Recall against maximum Tanimoto similarity to the training
actives, for three test sets built by unrelated rules: withheld by publication date, withheld at
random, and withheld by curator, the last being compounds deposited in BindingDB and absent from
ChEMBL. They trace one curve, so recall is a function of chemical distance rather than of how the set
was held out, and the expected sensitivity for a submitted compound is knowable at query time.
