# Chapter 2. Data

> Every quantitative statement in this chapter was computed during this session from the artefact
> named beside it, principally the 63 endpoint tables in
> `submission_package/06_TRAINING_DATA/endpoints/`. Citation numbers refer to
> `manuscript/references.md`.

---

## 2.1 What counts as a label

The single decision that shapes everything in this chapter is that a label must be a measurement.

This is not the obvious choice, and it is expensive. ChEMBL [30] carries curator-assigned activity
comments alongside its numeric potency values, and using them would multiply the available training
data several times over. The reason for refusing them is recorded in `docs/decisions_log.md` and it
is a lesson learned rather than a principle asserted. An earlier prototype of this system did train
on curated annotation, and a feature ablation showed it was reading the answer back out of
disease-association features derived from the same body of knowledge that produced the annotation.
With those features removed, structure-only performance collapsed to near chance. The model had
learned to restate a curator's opinion, and it had learned nothing about chemistry.

A measured label cannot fail in that particular way. It can fail in others, and section 2.9 sets out
which. But the circularity is closed, and closing it is what makes the rest of the validation in this
thesis meaningful rather than self-confirming.

Concretely, an endpoint table holds five columns, `smiles`, `label`, `pchembl`, `year` and `source`,
and every row is one compound measured at one endpoint. The classification cut is pChEMBL 7, which is
100 nM: a compound is active at or above it and inactive at or below pChEMBL 5, which is 10 µM. The
interval between the two is discarded rather than assigned, because a measurement that lands there
does not settle the question the label is asking.

## 2.2 The sources, and what each contributes

Across the 63 endpoint tables the panel holds **228,200 measured compound-endpoint records**. Counted
by the `source` field, they divide as follows.

| Source | Rows | Share |
|---|---:|---:|
| ChEMBL, exact potency values | 167,856 | 73.6% |
| ChEMBL, censored bounds recovered as inactives | 29,751 | 13.0% |
| BindingDB and ChEMBL, the same compound in both | 18,552 | 8.1% |
| The barrier set, B3DB with FDA-curated approved drugs | 7,807 | 3.4% |
| BindingDB alone | 3,143 | 1.4% |
| NPASS 3.0 | 1,054 | 0.5% |
| PubChem high-throughput screening | 37 | 0.02% |

Three observations follow from that table, and none is visible in the aggregate figure of 228,200.

**The panel is a ChEMBL panel.** Taking the exact values, the recovered bounds and the shared rows
together, 94.7 per cent of the data passes through one curation pipeline. Section 2.6 describes what
was done to test whether that dependence matters, and Chapter 8 reports the result.

**The second-largest single contribution is the recovered negative class**, at 29,751 rows. That is
larger than every non-ChEMBL source combined. Section 2.3 is therefore the methodological centre of
this chapter rather than a detail within it.

**Every row that NPASS 3.0 contributed went into an endpoint that was later withdrawn.** NPASS
supplied NFKB1 (263 rows), NR3C1 (140) and NRF2 (651), which sum exactly to the 1,054 in the table,
and all three failed the deployment gate described in Chapter 6. The natural-product source was added
to extend coverage into chemistry the panel handles badly, and the honest summary is that the attempt
was made, measured and abandoned. It is reported here rather than quietly dropped, because a source
table that lists only the sources that worked would misrepresent the effort.

The PubChem contribution of 37 rows is similarly a residue rather than a policy, and section 2.5
explains why it is not larger.

## 2.3 Recovering the negative class from censored bounds

A bioactivity database is a record of what people found when they looked. It is not a balanced
sample, and the imbalance is not random.

The conventional query for building a QSAR training set filters on the pChEMBL field, because that is
where the potency lives. A compound that was assayed and found **inactive** frequently has no
pChEMBL value at all. It is deposited as a censored bound: `standard_relation` of `>` with a
concentration, meaning the assay ran to its highest tested concentration without reaching an effect.
Filtering on pChEMBL discards exactly those rows. The result is a training set whose positives came
from measurement and whose negatives have to be invented, usually as property-matched decoys [12].

The recovery rule is a piece of interval arithmetic, and its whole content is knowing when to refuse
to answer. A censored bound defines an interval of possible true potencies. It settles a label if and
only if that entire interval falls on one side of the activity cut.

- `IC50 > 10 µM` places the true potency strictly below pChEMBL 5.0. The whole interval is inactive,
  so the compound is a **measured non-binder**.
- `IC50 > 100 nM` places the true potency below pChEMBL 7.0, which spans both classes. The compound
  might be a weak binder or might be inert, and the record does not say. It is **discarded as
  undecidable** rather than guessed at.

Applying that rule recovers **29,751 measured non-binders across 57 of the 63 endpoints**. These are
compounds somebody put in an assay and observed not to bind, which is a categorically different
object from a compound assumed inactive because it resembles nothing active.

The distinction matters most where it is easiest to overlook. A model validated against decoys is
being asked to separate real ligands from chemistry that was selected to be dissimilar to them, and
it will do so easily. A model validated against measured non-binders is being asked to separate
compounds that a medicinal chemist thought worth testing at that target from those that turned out
not to work, which is the discrimination a user actually needs. Chapter 6 reports the binder panel
against measured inactives for exactly this reason, and the resulting mean AUROC of 0.9174 should be
read as a harder number than a decoy-validated figure of the same size.

Six endpoints have no recovered inactives. They are the ones whose deposited record contains no
usable bounds, and their negative class remains what section 2.4 describes.

## 2.4 Where the negative class is still an assumption

Recovery is partial, and the chapter would be dishonest if it stopped at the previous section.

Where an endpoint has enough recovered inactives, its negatives are measurements. Where it does not,
the shortfall is made up with property-matched decoys drawn from the background library: compounds
matched on the physicochemical properties of the actives but with a Tanimoto similarity below 0.35 to
any of them. A decoy is an assumption wearing the clothes of a data point. It asserts that a compound
nobody has tested at this target would be inactive if they did, and that assertion is usually but not
always right.

This has two consequences that are carried forward rather than resolved.

First, the negative class means different things at different endpoints, which is one of the four
reasons Chapter 3 gives for fitting each endpoint independently rather than as one multi-task model.
Pooling a measured non-binder and an assumed one into a single loss silently averages a measurement
with a guess.

Second, any metric computed against decoys is easier than the same metric computed against measured
inactives. The panel registry records both quantities per endpoint, and this thesis quotes the
measured-inactive figure wherever one exists.

## 2.5 An experiment that was run and reverted

The obvious way to enlarge the negative class further is to add bulk high-throughput screening
inactives from PubChem. It was tried, on GSK-3β, the worst-skewed endpoint in the panel, and the
result is instructive enough to belong in the main text rather than an appendix.

Adding 4,276 measured PubChem inactives did two things at once
(`results/tables/inactives_audit.csv`).

It **corrected a real problem**. The GSK-3β model had been predicting 71.6 per cent of the 11,723
DrugBank drugs as active, which is not a possible state of the world. After the addition, that fell
to 16.4 per cent (`docs/INACTIVES_EXPERIMENT.md`).

It also **manufactured an apparent improvement that was not one**. Scaffold-split AUROC rose from
0.9369 to 0.9891, and the Matthews correlation from 0.5593 to 0.9013. Those are enormous gains, and
they are an artefact. The added negatives have a median Tanimoto similarity to the actives of 0.288,
and only 1.7 per cent of them reach 0.5. They are chemically unlike inhibitors, so separating them
from inhibitors is an easy problem the model was not previously being asked to solve. The audit
labelled the result *interpret with care, some easy negatives*, and it was correct to.

The decision was to revert. The honest discrimination estimate for GSK-3β remains 0.9369, and the
base-rate skew is handled instead by calibration and by the applicability-domain flag, and disclosed
as a data limitation. The correct version of this work is similarity-matched hard negatives, which
would fix the base rate without inflating the metric, and it has not been done.

That episode is the clearest illustration in this thesis of why a rising number is not evidence. Had
the experiment been kept, the panel would report a better AUROC at GSK-3β and would be worse at the
job it exists to do.

## 2.6 Pooling a second curator

For eleven protein targets a second independent source, BindingDB [19], was pooled with ChEMBL at
compound level. The reason was that the review asked for more data and ChEMBL had been verified as
near-complete for these targets, so genuine growth required a different curator rather than a deeper
query of the same one.

Because BindingDB's affinity export returns actives, the addition could only move the class balance
in one direction, and it was therefore audited against the ChEMBL-only baseline before being adopted
(`results/tables/expansion_audit.csv`). Over thirteen endpoints under the scaffold-grouped split the
headline metric moved by a **mean of -0.0002, ranging from -0.0151 at 5-HT2A to +0.0129 at A2A**. The
addition neither inflates nor degrades performance. A2A gained the most, +0.0129 in R², on the
largest single addition of 1,238 compounds; BBB and the antioxidant endpoint gained nothing because
nothing was added to them.

A change that moves the mean by two ten-thousandths is a change worth making only if the reason is
something other than performance, and it is: pooling a second curator is what makes the
cross-provenance validation of Chapter 8 possible at all. Compounds deposited in BindingDB and absent
from ChEMBL become a test set curated by different people from different papers.

## 2.7 Deduplication, and the leak it exists to prevent

The featuriser described in Chapter 3 excludes stereochemistry, so two enantiomers produce
byte-identical feature vectors. This is a limitation, bounded in Chapter 10. It is also a hazard, and
the hazard is more immediate than the limitation.

If two rows are identical in feature space and land on opposite sides of a cross-validation fold, the
model is tested on a compound it was trained on. The measured AUROC then includes a component of
memorisation, and no amount of scaffold grouping removes it, because the scaffold split groups by
Bemis-Murcko skeleton [2] and these rows share a skeleton by construction.

Rows identical in the 1,036-column representation are therefore collapsed before any split is drawn,
and a collapsed group whose labels disagree is dropped entirely rather than resolved by majority
vote, because the featuriser cannot tell its members apart and a vote would be inventing an answer.

The scale is not trivial. **15,104 duplicate rows exist in the tables before deduplication**, the
worst single endpoint being the barrier set at 3,773, and **zero reach a model**
(`results/tables/inversion_validation.csv`). Those 15,104 rows are correct chemistry: stereoisomers
are genuinely distinct compounds, deposited separately and measured separately. They are duplicates
only from the point of view of a stereo-blind representation, which is precisely why the collapse
happens at the featurisation boundary and not in the source tables.

Three checks confirm the boundary holds. On the deduplicated matrix the pipeline actually fits, no
InChIKey, no feature vector and no scaffold appears on both sides of any fold. On the raw table the
feature-vector overlap reaches 544. And across the panel, 52 targets were checked for scaffold
overlap between training and test and **none shares a scaffold**
(`results/tables/integrity_audit.csv`).

## 2.8 The panel endpoint by endpoint

An aggregate of 228,200 records describes no endpoint. The quantity that governs what any single
prediction is worth is the size of the set behind it, and those sizes span more than an order of
magnitude.

Across the **55 deployed endpoints that own a table**, the median is **3,789 rows**, from **387 for
KEAP1** to **10,276 for hERG**.

That set requires naming, because the phrase "the deployed panel" admits two readings that differ.
The set above is every deployed model owning a table in the endpoints directory, the barrier model
included. Excluding the barrier model, on the grounds that it is an exposure rather than a target
endpoint, gives 54 tables and a median of **3,587.5**. The difference is one endpoint: BBB holds
7,807 rows, sits above the median, and moves it. Neither reading is wrong and the figure quoted
throughout this thesis is 3,789, but a reader cannot check it without knowing which set is meant, so
the set is named wherever the number appears.

The practical consequence of the range is that a prediction at KEAP1 and a prediction at hERG are not
comparable statements of confidence, even when both return the same probability. KEAP1 has 387
measured compounds behind it. The uncertainty machinery of Chapter 4 exists so that this difference
reaches the user rather than being averaged away, and the model atlas figure plots every estimator
against its training-set size for the same reason.

## 2.9 Class balance, and what the deposition process does to it

Of the 228,200 rows, **176,123 are active and 52,077 inactive**, an overall active fraction of 0.772.
Across the 55 deployed endpoints the median active fraction is **0.825**, running from **0.122 at
Nav1.5** to **0.962 at P2X7**. Fourteen of the 55 are above 90 per cent active and 32 are above 80
per cent; only five sit below 50 per cent, the most balanced being MAO-B at 0.505, BChE at 0.543 and
KEAP1 at 0.558.

This skew is a property of the deposition process and not of the chemistry. People publish what
binds. A target whose set is 96 per cent active is not a target at which most compounds are active;
it is a target whose non-binders were never written down, or were written down as the censored bounds
section 2.3 recovers.

Three consequences run through the rest of the thesis.

**A high probability may say more about the endpoint than about the compound.** This is why the
engagement signal of Chapter 7 is a base-rate enrichment rather than a raw probability, and why a
prediction below an endpoint's own base rate is treated as evidence of inactivity rather than as weak
evidence of activity.

**Accuracy is uninformative here** and is not reported. A model that answers "active" at P2X7 scores
0.962. Every classifier is fitted with `class_weight="balanced"` for the same reason.

**Four receptor endpoints are not classification problems at all.** D2, A2A, 5-HT2A and SERT have
ChEMBL sets that are 96 to 98 per cent active, because only binders are reported, and a binary task
on such a set is ill-posed: Matthews correlations of 0.21 to 0.44 failed the deployment quality gate.
They are served as potency regressions on pChEMBL instead. That is a case where the right response to
a data limitation was to change the question rather than the model.

## 2.10 The temporal structure

Every row carries the year of the document it came from, and 212,227 of the 228,200 are dated. They
span **1976 to 2025** with a median of **2016**.

This field does no work at training time. It exists so that Chapter 8 can freeze the panel at a
cutoff year, refit on what was known then, and test on compounds first published afterwards. That is
the closest this thesis comes to a prospective test, and section 8.6 states plainly why it is not
one: the dates are real, but the analysis is retrospective, because no compound here was predicted
before it was measured.

The median of 2016 is worth noting for a different reason. It means roughly half the panel's evidence
postdates the deep-learning turn in cheminformatics and roughly half precedes it, and the
distribution is not uniform across endpoints. An endpoint whose chemistry appeared in a single burst
of publications has nothing on the far side of any wall that could be drawn, which is why 8 of the 47
deployed endpoints are excluded from the prospective arm.

## 2.11 What this data cannot support

Four claims are unavailable from this data whatever is done downstream, and they are stated here
because a modelling chapter that begins with an honest inventory of its inputs cannot be surprised
later.

**Direction of effect.** The label is a potency value, IC50, Ki, Kd or EC50. Those measure affinity,
which an agonist and an antagonist at the same receptor can share. ChEMBL's `action_type` field
carries direction, is sparsely populated, and was never pulled: it appears nowhere in these tables,
which retain only structure, label, potency, year and source. No modelling choice recovers what was
not collected.

**Stereochemistry.** Held out of the representation, for reasons given in Chapter 3 and bounded in
Chapter 10.

**Absolute potency below the cut.** A compound discarded as undecidable is genuinely unknown, not
weakly active. The panel has no opinion about it.

**Coverage of natural-product chemistry.** The training library has a median fraction-sp3 of 0.34,
and only 9.2 per cent of it is both sp3-rich and free of aromatic rings. The three endpoints added
specifically to address this, drawn from NPASS, all failed. This is a gap in the data that a better
model cannot close.

---

## Outstanding items for this chapter

1. The manuscript now states the alternative median as 3,587. The exact value is **3,587.5**, the
   median of an even-numbered set of integers. It should read 3,587.5, or 3,588 if a whole number is
   wanted, but not 3,587, which is a truncation.
2. The figure of 169,341 unique compounds, keyed by InChIKey of the desalted parent, is quoted in the
   manuscript and the evidence map. Counting distinct SMILES across the endpoint tables gives
   170,619, which is consistent with the two being different quantities but does not verify the
   InChIKey figure. It is not used in this chapter and should be re-derived before it is used.
3. A third figure for deployed sensitivity now exists. The panel registry gives a mean of 0.898 over
   47 endpoints, `inversion/results/H7_target_discrimination.csv` a median of 0.792 over 37, and
   `results/tables/integrity_audit.csv` a mean of 0.846 over its own set. All three are the deployed
   operating point measured on different held-out actives. Chapter 6 must say which is which, and no
   chapter should quote one without naming its hold-out.
