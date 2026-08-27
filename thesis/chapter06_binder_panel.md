# Chapter 6. The binder panel and the measured-inactive validation

> Every quantitative statement in this chapter was computed during this session from the artefact
> named beside it. The sensitivity figures follow the correction established in Chapter 5: where a
> published figure and a held-out figure differ, both are given and the held-out one is the claim.

---

## 6.1 What the binder panel is

The binder panel is the largest part of the system and the part a user is most likely to query. It is
**52 fitted classifiers, 47 of them deployed**, each answering one question about one target: does
this compound bind here at pChEMBL 7 or better, which is about 100 nM.

It is fitted in two modes (`submission_package/07_MODELS/binder_panel_registry.json`):

| Mode | Endpoints | What the negative class is |
|---|---:|---|
| `hybrid_decoys_plus_measured_inactives` | 44 | recovered censored bounds where they exist, topped up with property-matched decoys |
| `measured_labels_holdout` | 8 | measured labels only, no decoys needed |

Across the deployed endpoints the median training set holds 1,189 positives against 4,855 decoys, a
median ratio of 3.55 decoys per active. Those decoys are drawn from the decoy pool of Chapter 5, and
never from the pools that set or measure the threshold.

## 6.2 Why the validation had to change

Chapter 2 argued that a decoy is an assumption wearing the clothes of a data point. This chapter
measures what that assumption is worth, and the measurement is the reason the panel is reported the
way it is.

A model trained against property-matched decoys and then **validated** against property-matched
decoys is being asked to separate real ligands from chemistry selected to be dissimilar to them. It
is a genuine task, and it is far easier than the one a user needs, which is to separate compounds a
medicinal chemist thought worth testing at that target from those that turned out not to work.

Both numbers exist for every deployed endpoint, so the difference can be stated rather than argued:

| | Scaffold-grouped CV against decoys | Against measured inactives |
|---|---:|---:|
| Mean AUROC | **0.9784** | **0.9174** |
| Median | 0.9920 | 0.9470 |
| Minimum | 0.763 | 0.719 |
| Maximum | **1.000** | 0.985 |

The decoy-validated figure is higher on **45 of 47** endpoints, by a mean of **+0.0611**. Its median
of 0.992 and maximum of 1.000 are the signature of a saturated problem: on more than half the panel,
separating actives from decoys is very nearly solved, which tells a reader almost nothing about
whether the model is useful.

The endpoints where the two diverge most are the informative ones:

| Endpoint | Against decoys | Against measured inactives | Change |
|---|---:|---:|---:|
| GABA-A | 0.965 | 0.719 | **-0.246** |
| GluN2B | 0.993 | 0.761 | -0.232 |
| α7-nAChR | 0.977 | 0.763 | -0.214 |
| P2X7 | 0.998 | 0.813 | -0.185 |
| COX-2 | 0.935 | 0.783 | -0.152 |
| D2 | 0.989 | 0.872 | -0.117 |

GABA-A validated against decoys reads 0.965, which would place it comfortably in the middle of the
panel. Validated against compounds actually tested at GABA-A and found inactive, it reads 0.719, the
weakest endpoint deployed. Nothing about the model changed between those two numbers. Only the
question did.

**The panel therefore reports the second number everywhere**, and the 0.0611 it costs is the price of
a figure that means what a reader will take it to mean.

## 6.3 The measured-inactive hold-out, and how thin it is

The honesty of section 6.2 depends on there being enough measured inactives to validate against, and
for several endpoints there are not many.

Across the 44 deployed endpoints that have one, the measured-inactive hold-out has a **median of 228
compounds**, from **23 to 1,616**. Eleven endpoints have fewer than 100 and five have fewer than 50:
the CGRP receptor at 23, Nav1.8 at 33, GluN2B at 34, GABA-A at 37 and MT1 at 41.

That matters for how the headline range should be read. The panel's best AUROC against measured
inactives, 0.985, belongs to the CGRP receptor, and it is computed against 23 compounds. Its worst,
0.719 at GABA-A, is computed against 37. **The two ends of the reported range are the two thinnest
measurements in it**, and neither carries an interval in any current artefact. Adding Wilson or
bootstrap intervals to the per-endpoint AUROCs is the single cheapest improvement available to this
chapter and has not been done.

The active hold-outs are more comfortable, with a median of 299 and a range of 25 to 899.

## 6.4 Sensitivity, corrected

Chapter 5 established that the sensitivity figure stored in the panel registry, and quoted in every
document downstream of it, is computed over every active in the endpoint table rather than over the
held-out actives its own label claims. The correction belongs here too, because this is the chapter
where the panel's operating behaviour is reported.

| | As published | Held out only |
|---|---:|---:|
| Mean sensitivity | 0.8983 | **0.7638** |
| Median | 0.9320 | 0.8350 |
| Minimum | 0.639, COX-2 | **0.303, GABA-A** |
| Maximum | 0.997, CGRP | 0.993, CGRP |
| Endpoints firing for under half their own actives | 0 of 47 | **6 of 47** |

The claim this thesis makes is the right-hand column. The panel fires for about **0.76** of the
actives it should, not 0.90, and six deployed endpoints fire for fewer than half: GABA-A at 0.303,
COX-2 at 0.338, TAAR1 at 0.355, GluN2B at 0.374, SIRT1 at 0.378 and P2X7 at 0.465.

Note what those six have in common with section 6.2. **GABA-A, GluN2B, COX-2 and P2X7 appear in both
lists**: they are the endpoints that lose most when validated against measured inactives, and they
are the endpoints that fire for fewest of their own held-out actives. That is not a coincidence and
it is not two problems. A target whose measured inactives are hard to separate from its actives is a
target where any threshold must sit high to control false positives, and a high threshold is what
costs sensitivity. Chapter 5's operating-point argument and Chapter 6's discrimination argument are
the same fact seen from two sides.

## 6.5 The withdrawal set

Five endpoints were trained, tested and withheld. They are reported as part of the inventory rather
than removed from it, because a panel showing only what survived is a selection.

| Endpoint | AUROC vs measured inactives | Sensitivity | Positives | Why |
|---|---:|---:|---:|---|
| GluA2 | 0.696 | 0.103 | 68 | fires on glucose and atenolol at its calibrated threshold of 0.629 |
| Nav1.1 | 0.952 | 0.120 | 40 | fires on glucose, urea, glycine, lactate and atenolol at 0.571 |
| NRF2 | 0.789 | 0.545 | 70 | background false-positive rate 0.057, above the 5 per cent the panel holds to |
| NFKB1 | 0.459 | 0.000 | 39 | recovers no active while calling five trivial metabolites binders |
| NR3C1 | 0.410 | 0.167 | 60 | AUROC below chance against its own held-out measured inactives |

Two of these deserve comment because they show the withdrawal gate doing something a single metric
would not.

**Nav1.1 has an AUROC of 0.952**, which is better than eleven deployed endpoints. It was withdrawn
anyway. Its registry entry states the reasoning exactly: holding 5 per cent on random chemistry would
need a threshold of 0.594, its sensitivity at the calibrated cut is already 0.120, and no threshold
separates trivial metabolites from real ligands. "Its AUROC of 0.952 measures ranking, which is not
the quantity a deployed cut needs." That is the distinction Chapter 4 drew between ranking and recall,
arriving independently at the deployment gate.

**NR3C1 fails in the opposite direction.** It passes the specificity audit, firing on no trivial
molecule, and fails on discrimination instead, at an AUROC of 0.410, which is below chance. A panel
that only checked for firing on glucose would have deployed it.

### The natural-product endpoints failed for a diagnosable reason

NRF2, NFKB1 and NR3C1 were added specifically to extend coverage into natural-product chemistry, and
all three failed. The registry records why, and the diagnosis is about the labels rather than the
models (`results/tables/np_endpoint_assay_composition.csv`):

| Endpoint | Labelled records | Recorded as "Potency" | Direct binding constant (Ki or Kd) |
|---|---:|---:|---:|
| NRF2 | 1,029 | 1,019 | **0.0%** |
| NFKB1 | 362 | 344 | **0.3%** |
| NR3C1 | 461 | 160 | 14.8% |

For NRF2 and NFKB1 essentially every label is a pooled functional readout rather than a binding
constant. A functional potency does not define a binding class that a ligand fingerprint can
separate: two compounds can share a potency in a cell-based assay through entirely different
mechanisms, and the structure carries no signal for the label. NR3C1, the only one of the three with
real binding data, failed for the other available reason: 140 compounds after deduplication is too
few to fit.

This is the most useful negative result in the panel. It says that extending to natural-product
targets is not blocked by the chemistry being unusual, which was the assumed obstacle, but by the
assays in the public record for those targets not measuring binding.

## 6.6 Reproduction

The whole binder cross-validation was re-run from the endpoint tables with separately written scoring
code and compared against the recorded values (`results/tables/binder_cv_summary.csv`). Across the 47
deployed endpoints the re-run agrees with the recorded scaffold AUROC to a **mean absolute deviation
of 0.0041**, with a **maximum of 0.0386 at TAAR1**, which is also the panel's smallest binder training
set at 82 rows and its worst-calibrated endpoint from Chapter 4. The re-run mean is 0.9808 against a
recorded 0.9784.

A deviation of 0.039 on one endpoint is not a reproduction failure, but it is not nothing either, and
the endpoint where it occurs is the one where every other measurement in this thesis is also least
stable. TAAR1 should be read as the panel's least trustworthy deployed endpoint on all four counts:
smallest training set, worst calibration at 0.178 expected calibration error, held-out sensitivity of
0.355, and the largest reproduction deviation.

## 6.7 What the panel is, honestly summarised

Putting the corrected figures together, the deployed binder panel:

- discriminates against compounds measured and found inactive at the same target at a **mean AUROC of
  0.9174**, median 0.9470, from **0.719 at GABA-A to 0.985 at the CGRP receptor**, with the two
  extremes resting on 37 and 23 measured inactives respectively;
- fires for a **mean of 0.7638** of held-out actives at the deployed threshold, median 0.8350, with
  **six of 47 endpoints below 0.50**;
- would report 0.9784 and 0.8983 for those two quantities if validated the conventional way, and
  neither figure would mean what a reader takes it to mean;
- carries five withdrawn endpoints in its inventory, three of which failed because the public assay
  record for those targets does not measure binding.

---

## Outstanding items for this chapter

1. **Per-endpoint AUROC against measured inactives carries no interval.** The reported range, 0.719
   to 0.985, has its two extremes computed on 37 and 23 compounds. Wilson or bootstrap intervals
   would settle how much of that spread is real, and are cheap.
2. **The sensitivity correction from Chapter 5 must propagate.** The manuscript, the evidence map and
   the technical report all quote 0.8983 with the qualifier "on actives withheld by scaffold".
3. **TAAR1 should be reviewed for withdrawal or for a warning flag.** It is the least trustworthy
   deployed endpoint on four independent measures, and nothing in the interface says so beyond the
   generic reliability flag.
4. The four endpoints appearing in both the low-AUROC and low-sensitivity lists, GABA-A, GluN2B,
   COX-2 and P2X7, are already flagged unreliable by the threshold sequence. Whether they should be
   deployed at all is a decision this thesis records rather than makes.
