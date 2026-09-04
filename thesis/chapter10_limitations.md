# Chapter 10. Limitations, and what should be done next

> **Provenance.** Every figure in this chapter was read from a file in this repository during the
> session in which the chapter was written. Sources are named inline. Where a number in a document
> disagrees with the artefact it cites, the artefact is taken as correct and the disagreement is
> reported rather than reconciled silently.

## 10.1 Three different things called limitations

A closing chapter usually mixes three kinds of statement and is weaker for it. They are separated
here.

**What the system cannot do**, which is a property of the science and would survive any amount of
tidying. Sections 10.2 to 10.6.

**What the record gets wrong**, which is a property of this project's housekeeping and is fixable.
Section 10.7 identifies the one structural defect behind most of it, and section 10.8 is the full
ledger of 41 items raised across the preceding nine chapters, ten of which are now closed.

**What should be done next**, ordered by what it would change rather than by how interesting it is.
Section 10.9.

One item of that future work is done in this chapter rather than proposed, because Chapter 9 named
it as the falsification suite's largest gap and a gap named but not filled is only half an
observation. Section 10.3 is H10.

## 10.2 The limitation that dominates every other: accuracy is a function of chemical distance

If a reader takes one limitation from this thesis it should be this one, because it bounds every
number in Chapters 3 to 8 and because the system's own users will meet it first.

The panel's recall depends almost entirely on how close a query is to chemistry the panel has already
been shown. Measured on the prospective time split, over 39 endpoints
(`results/tables/external_novelty_strata.csv`):

| Nearest-training Tanimoto | Actives | AUROC | Recall at the deployed threshold |
|---|---:|---:|---:|
| Below 0.40, a different chemotype | 4,366 | 0.6711 | **0.1612** |
| 0.40 to 0.55, a related series | 5,228 | 0.7645 | 0.5532 |
| 0.55 to 0.70, the same series | 3,674 | 0.7871 | 0.7403 |
| 0.70 and above, a close analogue | 2,586 | 0.8617 | 0.8616 |

Recall runs from 0.16 to 0.86 across that range, a factor of five. On the cross-provenance set,
built by a rule that knows nothing about dates, the same curve is steeper still: 0.0494 in the
lowest band against 0.9011 in the highest.

**On genuinely novel chemistry the tool finds roughly one true activity in six, and across a database
boundary roughly one in twenty.** No analysis in this thesis improves that figure. What Chapter 8
establishes is that it is *predictable*, not that it is good, and the server now reports the expected
recall beside a negative result for exactly that reason. A silence at Tanimoto 0.30 is close to
uninformative, and a user who reads it as evidence of inactivity is making an error the interface
must keep working to prevent.

This is not a defect of the random forest. It is what a similarity-based featuriser does, and every
model family benchmarked in Chapter 3 shares it. It is the reason the composition finding of section
8.7 matters: the apparent decay between a random and a temporal split is mostly this curve, seen
from a different angle.

## 10.3 The barrier model, tested rather than assumed

Chapter 9 observed that the falsification suite has nine hypotheses and none of them tests the
blood-brain barrier classifier, which gates every disease score and names the architecture. Its only
evidence was an external validation, which is a confirmation and not an attempt at refutation.

**H10** states the null the project would least like to be true:

> Permeability is a bulk physicochemical property. The twelve descriptors already in the feature
> vector encode it, and a rule or a small model over those descriptors alone ranks compounds as well
> as a 1,036-column random forest does.

This is not a straw man. Polar surface area, molecular weight and hydrogen-bond donor count are the
terms of every published CNS-permeability heuristic, and barrier permeability is the endpoint in this
project with the strongest prior claim to being predictable without a fingerprint at all. If the null
holds, the barrier model should be replaced by a rule a reader can apply by hand.

`inversion/inv_barrier_necessity.py` scores every comparator on two populations: the same
Bemis-Murcko folds the deployed cross-validation used, read from the stored out-of-fold predictions
so the assignment is identical rather than merely similar, and the 241 FDA-curated approved drugs
that are absent from B3DB and also distinguishable from training chemistry in feature space.

| | Scaffold hold-out, n = 3,901 | External approved, n = 241 |
|---|---:|---:|
| **Deployed forest, 1,036 features** | **0.8779** | **0.7934** |
| Random forest, 12 descriptors only | 0.8441 | 0.7701 |
| Logistic regression, 12 descriptors | 0.8105 | 0.7420 |
| TPSA alone, no fitting of any kind | 0.7687 | 0.7440 |
| Hydrogen-bond donors alone | 0.7272 | 0.7174 |
| QED alone | 0.7138 | 0.6332 |
| CNS heuristic: TPSA ≤ 90 and MW ≤ 400 | 0.6751 | 0.6369 |
| Molecular weight alone | 0.6643 | 0.5992 |
| cLogP alone | 0.6122 | 0.6058 |

Read as point estimates the fingerprint wins on both populations, by 0.0338 and 0.0233. **Read with
its interval, the answer changes on the population that matters.** Every comparator carries a paired
bootstrap of 2,000 resamples over compounds, so both models are always judged on the same draw:

| Comparison | Margin | 95% interval | Bootstrap p |
|---|---:|---|---:|
| Descriptor forest, scaffold hold-out | 0.0338 | −0.0409 to −0.0263 | < 0.001 |
| Descriptor forest, **external approved** | 0.0233 | **−0.0553 to +0.0106** | **0.083** |

**The verdict is WEAKENED, not SUPPORTED.** On the training distribution the fingerprint demonstrably
earns its place. On 241 approved drugs the model has never seen, its advantage over twelve
descriptors is not established: the interval crosses zero.

The verdict rule is stated that way deliberately. Chapter 9 criticised H4 for deciding SUPPORTED on a
point-estimate comparison whose interval contains its own comparator, and a hypothesis written in
response to that criticism must not repeat it. Under the naive rule H10 would read SUPPORTED. Under a
rule that reads the interval it does not, and the difference is the whole value of having written the
criticism down.

Two consequences, and neither is that the barrier model should be removed.

**What the project is entitled to claim shrinks.** The fingerprint helps on chemistry resembling
training data and its contribution to unseen approved drugs is not demonstrated. Any wording implying
that the barrier model substantially outperforms physicochemical prediction should be qualified. What
it does clearly beat, on both populations and by a wide margin, is the published heuristic: 0.7934
against 0.6369 externally.

**The result converges with section 10.2 and with the calibration finding below.** The barrier model
is the component whose advantage vanishes soonest as chemistry becomes unfamiliar, and it is also the
worst-calibrated of the eight core classifiers. It is the single weakest link in the architecture and
it sits at the gate.

H10 is now in `inversion/results/H10_barrier_necessity.csv` and in `VERDICTS.csv`. It reads its
verdict from its own artefact rather than having it recomputed in `summarise.py`, which is the
one-rule arrangement Chapter 9 recommended; H1 to H9 have not been migrated.

## 10.4 What the panel does not cover

**16 conditions and 51 graph targets** (`inversion/results/GRAPH_FINGERPRINT.json`), of
52 endpoints with 47 deployed and **six of those 47 failing the reliability gate**: COX-2, GABA-A,
GluN2B, P2X7, SIRT1 and TAAR1. Any mechanism outside that set is invisible, and invisibility is
reported as silence, which is the failure mode section 10.2 describes.

**Natural-product chemistry is the largest and best-characterised gap**, and the reason is visible in
the training library before any model is scored. Across the **170,617** structures in the endpoint
tables that parse, the median fraction of sp3-hybridised carbon is **0.3448**, and only **9.23 per
cent** are both sp3-rich, at a fraction of 0.55 or above, and carry at most one aromatic ring. The
quartiles are 0.2308 and 0.4615, so the distribution is narrow as well as low
(`results/tables/library_sp3_coverage.csv`). Terpenoid and steroidal chemistry is therefore largely
outside the library, which is a statement about coverage rather than about the models.

Assembling an external test set from NPASS began with 3,547 candidate compounds across the panel and
produced three scorable endpoints
(`results/tables/external_natural_products_summary.csv`):

| Endpoint | n | Actives | AUROC | Median max Tanimoto |
|---|---:|---:|---:|---:|
| GBA1 | 89 | 2 | 0.7414 | 0.270 |
| AChE | 41 | 4 | 0.4628 | 0.500 |
| hERG | 80 | 5 | 0.4613 | 0.309 |

Two things about this table matter more than the AUROCs. **1,385 of the candidates had to be removed
as contaminated**, being training compounds written as a different tautomer or salt, which is the
same stereo-blind featuriser problem that Chapter 8 found in the external barrier set. And **23 of
the surveyed targets have no deployed estimator at all**, so most of the natural-product question
cannot be asked of this system.

With two to five actives per endpoint, those three intervals are far too wide to establish failure.
They are equally far from supporting a claim of natural-product competence, and the honest statement
is that the system is untested on this chemistry rather than that it fails on it.

## 10.5 Uncertainty is complete for eight estimators and largely absent for forty-seven

The uncertainty stack described in Chapter 4, isotonic calibration on out-of-fold predictions plus
Mondrian conformal prediction plus an applicability band, is complete for the **eight** core
classifiers and not for the target panel. Counted per layer rather than asserted, since the layers
do not stop at the same place:

| Layer | 8 core classifiers | 47 deployed binder endpoints | 11 ADME and auxiliary |
|---|---:|---:|---:|
| Isotonic calibration | 8 of 8 | none | none |
| Conformal set | 8 of 8 | none | none |
| Applicability band | 8 of 8 | 39 of 47 | 1 of 11 |
| Per-endpoint calibration error | 8 of 8 | 38 of 47 | none |

Two counts in that table deserve stating rather than leaving to be read off. **Nine of the 47
deployed binder endpoints carry no measured calibration error at all**: CGRP, GABA-A, GBA1, GluN2B,
KEAP1, Nav1.5, Nav1.8, SIRT1 and α3β4 nAChR. **Eight carry no per-endpoint applicability
reference**, so a query against them falls back to the global domain call. Neither absence is
reported anywhere in the interface.

The core figures are good and the headline is honest about the mean but not about the spread.
Expected calibration error falls from a mean of 0.0801 raw to **0.0147** after isotonic calibration,
over a range of **0.0049 at BACE1 to 0.0412 at BBB** (`results/tables/calibration.csv`).

**The barrier model is the worst-calibrated of the eight, at 0.0412, nearly three times the mean.**
That has not been stated anywhere in the project, and it matters more than the rank alone suggests,
because the barrier probability is not merely reported: it multiplies every disease score. A
mis-calibrated gate propagates into every condition at once.

Across the 38 Platt-scaled binder endpoints that carry one, the expected calibration error is
**mean 0.0762, median 0.0675, range 0.0300 to 0.1780**
(`results/tables/integrity_calibration_per_target.csv`), about five times the core figure. The
binder endpoints also have **no conformal coverage statement of any kind**. So the uncertainty
machinery is most complete for the eight endpoints a user is least likely to query and thinnest for
the 47 they are most likely to.

## 10.6 What the specificity figure is, and what it is not

Specificity on the non-CNS library is **0.925**, 925 of 1,000, interval 0.907 to 0.9397
(`results/tables/noncns_specificity_summary.csv`). The artefact labels it correctly and the label
should travel with the number: it is an **upper bound**, because the compounds are presumed inactive
on the grounds that nothing is recorded about them, not proven inactive.

Two further qualifications, both established elsewhere in this thesis:

- The compounds are drawn from within the applicability-domain reference, so the figure does not
  bound behaviour on genuinely distant chemistry. H4 addresses that and, as Chapter 9 shows, does so
  on 61 compounds carrying one false positive, which is not enough to separate the strata.
- The applicability flag itself is a weak instrument. Against chemistry genuinely absent from the
  reference it separates at a median maximum similarity of 0.59 for unseen drugs against 0.47 for
  non-drug-like controls (n = 25, Mann-Whitney p = 1.11 × 10⁻³), but only **20 per cent** of
  non-drug-like structures fall below the deployed threshold of 0.30
  (`results/tables/inversion_validation.csv`). The conformal interval and the nearest-analogue
  distance are the stronger statements and the interface presents them as such.

## 10.7 One structural defect explains most of the record's errors

Across nine chapters this thesis has reported roughly two dozen places where a document disagrees
with an artefact. They are not independent mistakes. Almost all of them are instances of one defect:

> **A number that lives in a generator rather than in an artefact cannot be checked, cannot be cited,
> and drifts silently.**

The instances, in the order they were found:

| Chapter | The number | Where it lived |
|---|---|---|
| 3 | The censored-bound recovery count, 21,994 against a recount of 29,751 | hard-coded in the report generator |
| 8 | The label-permutation null, 0.4938 and 0.4921 | hard-coded in a literal prose block |
| 8 | The independent-reproduction figure, 26 values to 4.7 × 10⁻⁵ | hard-coded in the same block |
| 9 | The H2 weight ablation, 0.7917 / 0.7911 / 0.7899 over 15,609 compounds | hard-coded in `app.py` |
| 9 | Five co-firing correlations and three conditional probabilities | hard-coded in `app.py`, user-visible, since derived from the artefact |
| 10 | The natural-product coverage figures | computed inside the report generator, held by no file until this chapter |

The last row is the mildest and the most instructive. Those two figures are genuinely *measured* at
report time, from the endpoint tables, and the generator's own comment records that they had
previously drifted as prose: the median was stated as 0.36 against a measured 0.34, and the sp3-rich
share as 3.3 per cent against a measured 9.2. Measuring them was the right response. Not writing them
to a file was not, because no other document could cite them and nothing checked them. They are now
in `results/tables/library_sp3_coverage.csv`, declared in the freshness graph, and they reproduce the
generator's values exactly at 0.3448 and 9.23 per cent. What remains is for the generator to read the
artefact instead of recomputing it, which is outstanding item 1 below.

The remedy is uniform and cheap: **every number a document states should be read from a file at build
time, and every file should be declared in the freshness graph.** This project already does that for
most of its figures, which is why the exceptions are visible at all.

Two related gaps in the checking machinery were found alongside them, and both are recorded in
section 10.8:

- **The freshness graph's inputs are data and models, never code.** A generator can be corrected
  without anything marking its output stale, which is why `background_specificity.csv` still holds
  pre-correction output from a script fixed on 27 August.
- **The falsification verdicts are pinned by no test.** A rerun that flipped one would be caught only
  by a person reading the report.

## 10.8 The ledger

41 items were raised across Chapters 1 to 9. Ten are now closed. They are listed as closed rather
than removed, because a thesis that silently drops its own outstanding items is doing the thing this
chapter is about.

**Closed while writing.** Six were closed in the course of writing Chapters 5 to 10.

| Item | Raised | Closed by |
|---|---|---|
| The published binder sensitivity is measured on training compounds | 5 | Commit `ab88039`; registry now holds the held-out figure on all 47 endpoints, reconciliation gap +0.000 |
| `calibrate_background_specificity.py` should score held-out actives | 5 | Fixed at the source, so the overwrite cannot recur |
| The specificity change from 0.949 to 0.925 was untraced | 7 | Traced to stale base rates corrected by `fd7aaa0`; not a degradation |
| The label-permutation null has no artefact | 8 | `results/tables/permutation_null.csv`, built for Chapter 9 |
| The falsification suite has no hypothesis for the barrier model | 9 | H10, section 10.3 |
| The natural-product coverage figures have no artefact | 10 | `results/tables/library_sp3_coverage.csv`, built for this chapter |

**Closed after writing.** Four more were repaired on 2026-09-02 in commit `1f91925`, which audited
Chapters 8 to 10 against the files they cite. They are recorded here because the ledger is the
chapter's own claim about the state of the project, and a ledger that goes stale is the defect
section 10.7 describes.

| Item | Raised | Closed by |
|---|---|---|
| The interface quotes co-firing correlations that no longer match the artefact | 9 | `FAMILY_COFIRE` is now derived from the artefact by a function. The repair found a worse defect than Chapter 9 had: the nicotinic entry asserted a family-level conclusion resting on **one** joint approved drug, which the file could not distinguish from a certainty because every rate in it is a multiple of 1/400. The artefact now carries `n_drugs`, `n_a`, `n_b` and `n_joint`, and suppresses any pair supported by fewer than ten drugs |
| `app.py:102-104` quotes the withdrawn H2 triple | 9 | Corrected to 0.7901 / 0.7897 / 0.7874 over 7,008, with a note recording what it used to say. Audit item BS-M-04 now closed in the source |
| `background_specificity.csv` holds pre-correction output | 9 | Refreshed by updating only its two reported columns from the registry. **No threshold and no background false-positive rate moved on any of the 47 endpoints**, and the shipped copy in the submission package was updated with it. Two tests pin both |
| The technical report's cross-provenance section is stale | 8 | The 0.868 was not hard-coded; the report computes it and had never been rebuilt, because `TECHNICAL_REPORT.md` was the one generated document outside the freshness graph. It is now declared against all 23 of its inputs and rebuilt at 0.716 |

The third of those deserves a sentence, because Chapter 9 posed it as a dilemma and the resolution is
narrower than either horn. Four of that file's six columns were already correct; only the two
reported ones were stale, so refreshing those two from the registry avoids re-running the sequence
that would rewrite every operating threshold to correct a reported one. And the version a reviewer
holds was corrected at the same time as the local one, which is the part that mattered: **the
disagreement was visible to reviewers before it was visible to us.**

**Open, ordered by consequence.** The first two change something a reader would quote; the rest are
housekeeping.

1. **The sensitivity correction has reached the registry and not all of the prose.** The manuscript
   and `EVIDENCE_MAP.md` still quote 0.8983 under the qualifier "on actives withheld by scaffold".
   The correct figure is 0.7638. The technical report has been rebuilt; these two have not. (5, 6)
2. **The headline calibration figure should be qualified everywhere it appears**, since 0.0801 to
   0.0147 covers eight of the deployed estimators and the 38 binder endpoints sit at a mean of
   0.0762. Section 10.5 adds that the barrier model is the worst of the eight. (4)
3. Give each hypothesis one verdict rule, as H10 now has, and pin the verdicts with a test. The
   suite has grown to 73 tests and none of them reads `inversion/`. (9)
4. Make the freshness graph aware of code, so a corrected generator marks its output stale. The
   repairs above were all found by reading rather than by a check. (9)
5. Two citation gaps remain in Chapter 1: a CNS attrition series, and verified citations for the
   named comparison servers. (1)
6. The unique-compound count of 169,341 keyed by InChIKey is quoted in the manuscript and the
   evidence map and has not been re-derived. Counting distinct SMILES gives 170,619, which is
   consistent with the two being different quantities but verifies nothing. (1, 2)
7. Per-endpoint AUROC and sensitivity carry no intervals anywhere, and the two extremes of the
   reported 0.719 to 0.985 range rest on 37 and 23 compounds. Wilson or bootstrap intervals are
   cheap and would settle how much of that spread is real. (4, 6)
8. TAAR1 should be reviewed for withdrawal or an explicit warning; it is the least trustworthy
   deployed endpoint on four independent measures. (6)
9. Nav1.5, SIRT1 and TAAR1 define "active" differently in two scripts, and a test should assert
   that `sensitivity_basis` matches the population the stored sensitivity was computed on. (5)
10. The H3 refutation has an unstated exception: gating is rank-invariant across the 14
    non-peripheral conditions but not between those and migraine and multiple sclerosis. (7)
11. The haloperidol figures in the manuscript describe a fixed defect without saying so, and the
    learning-curve conclusion in the technical report is not supported by its own table. (3)
12. Two evidence pointers in `docs/decisions_log.md` cite paths that do not exist and describe the
    superseded ensemble. (3)
13. `in_domain` in the artefacts means T ≥ 0.30, which is the interface's "in or near domain". The
    naming should be reconciled. (4)
14. The eight endpoints excluded from the prospective arm should be listed with their reasons. (8)
15. The three 4 August logs in `inversion/` record a superseded verdict and should be refreshed or
    deleted. (9)
16. The alternative median in the manuscript should read 3,587.5, not 3,587. (2)

## 10.9 What should be done next, ordered by what it would change

**First, the things that would change a user's decision.**

*Extend the uncertainty stack to the target panel.* Isotonic calibration for those binder endpoints
with enough held-out positives to support it, per-endpoint calibration error carried with every
prediction rather than a headline inherited from elsewhere, and a conformal coverage statement for
the 38 endpoints that have none. This is the largest single improvement available and needs no new
data.

*Report intervals on every per-endpoint figure.* An AUROC of 0.719 on 37 compounds and one of 0.985
on 23 are not comparable numbers, and presenting them in one column invites a reader to rank
endpoints on noise.

*Decide about the six endpoints that fail the reliability gate.* This thesis records the decision to
deploy them with a flag rather than making it. Section 6 of Chapter 5 sets out why the flag is
defensible; whether TAAR1 in particular should be served at all is a question the project has not
answered.

**Second, the things that would change what the project may claim.**

*Migrate H1 to H9 to the single-verdict-rule arrangement and pin the verdicts with a test.* The suite
is the strongest evidence in this thesis and it is the least defended part of the codebase.

*Write the missing hypotheses.* H10 filled the largest gap; others remain. The most valuable is a
null for the aggregation rule itself, asking whether the maximum over engaged targets beats a simple
count of them, since Chapter 7 justifies the maximum by argument and by H8's redundancy finding
rather than by a direct comparison.

*Make the freshness graph aware of code*, and add a check that every number in a built document is
traceable to a declared artefact. The second is ambitious; the first is an afternoon and would have
caught the `background_specificity.csv` case.

**Third, the things that need new data rather than new analysis.**

*Natural-product chemistry.* Three endpoints with two to five actives each cannot answer the
question. This needs a curated external set assembled with the feature-vector contamination check
from the start, not the InChIKey check.

*A prospective test in the proper sense.* Every result in Chapter 8 is retrospective: the dates are
real but no compound was predicted before it was measured. A pre-registered set of predictions,
lodged before the corresponding assays are run, is the only evidence that would settle what the
temporal simulation can only approximate.

*Measured brain exposure rather than a barrier classification.* The gate multiplies by a probability
of being CNS-positive, which is a coarser quantity than unbound brain-to-plasma ratio. H10 suggests
the classifier is close to what descriptors already give on unseen chemistry, so the improvement
here is likely to come from a better target variable rather than a better model of the current one.

## 10.10 What this thesis establishes

Stated as plainly as the evidence allows, with the spread rather than the mean wherever one exists.

- **The component models rank compounds well and the ranking is not an artefact of leakage.** Random
  10-fold AUROC runs 0.899 to 0.9764 across the eight core classifiers and scaffold-grouped
  10-fold runs 0.8777 to 0.9648. Under permuted labels the same pipeline on the same folds returns 0.4959 and
  0.5026, every one of the sixteen values within 0.020 of chance.
- **Accuracy is a function of chemical distance, and the function is now measured.** Recall runs from
  0.1612 to 0.8616 across four novelty bands, and three independently constructed test sets trace the
  same curve.
- **The disease layer carries real information about mechanism and does not predict indication.**
  0.7901 top-3 accuracy against the project's own map and a permutation null of 0.1625; a mean
  per-indication AUROC of 0.6158 against external clinical indications, ranging from 0.7981 for
  depression to 0.4898 for epilepsy.
- **Exposure gating is a filter and not a discriminator**, and the curated edge weights carry no
  measurable information.
- **The panel's sensitivity at its deployed operating point is 0.7638 on held-out actives**, median
  0.835, from 0.303 at GABA-A to 0.993 at CGRP, with six of 47 endpoints firing for fewer than half
  their own actives. Silence is not evidence of inactivity, and the number that makes a silence
  interpretable now travels with it.
- **Specificity is 0.925 on presumed-inactive non-CNS chemistry**, which is an upper bound.
- **Four of ten falsification hypotheses were refuted and two weakened.** The refutations cost the
  project a claim about the knowledge graph, a claim about the gate, a change to the interface, and
  four withdrawn endpoints that were calling glucose and urea binders.

What the thesis does not establish is anything prospective, anything about natural-product chemistry,
and any advantage for the barrier model over twelve descriptors on chemistry it has not seen. Those
three are the honest boundary of the work, and section 10.9 says what would move them.

---

## Outstanding items for this chapter

1. **`build_technical_report.py` should read `library_sp3_coverage.csv`** rather than recomputing
   the two figures inline. The artefact now exists and agrees with the generator exactly, so this is
   a substitution rather than a correction, but until it is made the number in the report is still
   one that nothing checks.
2. **H10 covers the barrier model's ranking, not its role as a gate.** Whether replacing the barrier
   probability with a descriptor rule would change which compounds receive a disease call, at the
   0.30 reporting threshold, is a separate and more directly user-relevant question that this
   chapter does not answer.
3. **The ledger in section 10.8 should be maintained as a file rather than as prose**, for exactly
   the reason section 10.7 gives.
