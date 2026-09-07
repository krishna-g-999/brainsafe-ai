# Chapter 9. The falsification suite, and what it cost

> **Provenance.** Every figure in this chapter was read from a file in this repository during the
> session in which the chapter was written. Sources are named inline. Where a number in a document
> disagrees with the artefact it cites, the artefact is taken as correct and the disagreement is
> reported rather than reconciled silently.

## 9.1 The premise

The suite's own plan states the problem it was built for, and states it better than a paraphrase
would:

> Every result so far has been produced by people who wanted the tool to work, including me. The
> validation to date establishes that the target models rank compounds correctly. It does not
> establish that the layer a user actually reads, the per-disease brain-relevance score, means
> anything at all.
>
> `inversion/PLAN.md`

That is the correct diagnosis of what Chapters 3 to 8 leave open. Cross-validation, conformal
coverage, applicability banding and the prospective simulation are all measurements of the target
models. The disease layer sits on top of them, is assembled from a curated graph and a gating rule,
and until August 2026 nothing in the project had tested it against an outcome at all. A system can
have well-calibrated components and an output layer that means nothing, and no amount of component
validation will detect that.

The suite is therefore not more validation. It is an attempt to make the central claims fail, run by
the same person who wrote them, which is a weaker instrument than independent replication but a much
stronger one than another favourable measurement.

## 9.2 What makes a test capable of failing

Four constraints separate this exercise from the rest of the evidence, and they are worth stating
because they are the reason its results are usable.

**Every hypothesis is paired with a null that could produce the same apparent success by accident.**
This is the constraint that does the work. A disease layer that recovers the right condition 79 per
cent of the time sounds decisive until one asks what always naming the three commonest conditions
would score, and the answer is 0.5508 (`inversion/results/H1_disease_layer.csv`). Most of the suite's
value comes from having asked that question before seeing the result rather than after.

**The analysis is read-only.** No trained model, no curated dataset and no threshold was changed in
order to obtain any number in it. Acting on a finding is a separate step recorded in the git history,
so a wording change made in response to a result cannot be mistaken for part of the evidence for it.

**Results are written under `inversion/`, not `results/`.** A falsification cannot be filed where a
validation would be found and quoted as one.

**Where predictive power is at issue, scoring uses `models_rf/holdout/`**, the scaffold-hold-out
twins, which never saw the compounds they score.

The reporting rule is stated in the plan as well: a refuted hypothesis is the most valuable outcome
available, and is to be reported as prominently as a confirmation. Section 9.4 is the test of whether
that was honoured, and it largely was.

## 9.3 The nine hypotheses

Four hypotheses were refuted, four supported and one weakened
(`inversion/results/VERDICTS.csv`). The graph under test is fingerprinted, at SHA-256 `25cae8c8...`
over 16 conditions and 51 targets (`inversion/results/GRAPH_FINGERPRINT.json`), so a verdict cannot
silently come to describe a different graph.

| | Hypothesis | Verdict | Argued in |
|---|---|---|---|
| H1 | The disease score is informative | SUPPORTED | 7.7 |
| H2 | The curated edge weights add value | **REFUTED** | 7.6 |
| H3 | BBB gating discriminates between diseases | **REFUTED by construction** | 7.5 |
| H4 | Specificity transfers to novel chemistry | SUPPORTED | **9.6** |
| H5 | Read-across beats a frequency baseline | SUPPORTED | **9.6** |
| H6 | The disease scores match real clinical indications | WEAKENED | 7.7 |
| H7 | Some panel targets are non-discriminative | **REFUTED** | 5.9 |
| H8 | Engaged targets are independent observations | **REFUTED** | 7.4, 7.8 |
| H9 | The disease layer discriminates between compounds | SUPPORTED | 7.7 |

Seven of the nine are argued in full elsewhere in this thesis and are not re-derived here. This
chapter takes the four questions the individual hypotheses cannot answer: what the refutations cost,
what the suite found that nothing else did, whether its two untreated hypotheses are as strong as
their verdicts read, and what the suite does not cover.

## 9.4 What the refutations cost

A falsification suite is worth what the project paid for it. On this test the record is good, and the
payments are traceable in the git history rather than asserted.

**H2 changed how the knowledge graph is described.** The curated edge weights score 0.7901 against
0.7897 for uniform weights and 0.7874 for randomly permuted ones, a spread of 0.0027
(`inversion/results/H2_weight_ablation.csv`). Commit `f40b5d7`, "Act on the inversion findings:
describe the weights and the BBB term accurately", relabelled them in the source, the interface and
the manuscript as a mechanistic prior rather than as tuned parameters. The weights were kept, which is
defensible because they still express which link is the more direct, but no predictive claim is made
for them.

**H3 turned a design claim into a filter claim.** Because the gate multiplies the 14 non-peripheral
conditions by the same factor it is rank-invariant among them and cannot sharpen a disease call. The
same commit corrected the wording wherever it appeared.

**H7 cost an endpoint and a disclosure.** Commit `8bb402b`, "Expand the panel where it was measured
to be weakest, and withdraw an endpoint that was firing on sugar", followed from the H7 line of
enquiry. A direct audit of the served models, `results/deployed_specificity_audit.csv`, then found
four failures, three of them by firing on molecules that are not drugs at all:

| Endpoint | FPR on random chemistry | Trivial molecules called binders |
|---|---:|---|
| Nav1.1 | 0.0800 | glucose, urea, glycine, lactate, atenolol |
| GluA2 | 0.0717 | glucose, atenolol |
| NRF2 | 0.0567 | none, but above the 5 per cent target |
| NF-kappa-B1 | 0.0267 | glucose, urea, acetate, glycine, lactate |

All four were withdrawn, together with NR3C1, which passed this audit and was withdrawn on other
grounds. The panel stands at 52 endpoints with 47 deployed.

**H8 changed the interface.** Across 400 approved drugs, 37 targets fire at least once but span only
16 independent directions by Kaiser's rule on the correlation matrix of the firing indicators
(`inversion/results/H8_panel_independence.csv`), and five homologous pairs co-fire above phi = 0.5
(`inversion/results/H8_family_correlation.csv`):

| Pair | phi | P(second engaged given first) |
|---|---:|---:|
| mu-opioid and kappa-opioid | 0.798 | 0.657 |
| SERT and NET | 0.732 | 0.671 |
| D2 and D3 | 0.720 | 0.717 |
| 5-HT2A and 5-HT7 | 0.671 | 0.659 |
| 5-HT1A and 5-HT7 | 0.583 | 0.629 |

The interface now groups engaged targets by homology family and quotes the measured correlation
rather than leaving a reader to assume independence. Section 9.8 reports that the numbers it quotes
are no longer the ones in the artefact.

**Two costs were paid on the suite itself rather than on the tool**, and they matter more than they
look. Commit `7fa7203`, "Make the inversion checks capable of failing, and let one of them fail",
repaired checks that could not have returned a negative result. Commits `da82d80` and `a3d184f`
restricted H1 to compounds that are training actives of no other panel target, because a compound
active at two targets is memorised by one of the models scoring it, and shuffling the
target-to-disease map does not control for that: the memorisation survives the shuffle and is merely
sent to the wrong condition. That correction removed roughly three fifths of the evaluation set, from
16,690 held-out entries across 41 endpoints (`models_rf/holdout/heldout_actives.json`) down to the
7,008 H1 now scores, and it lowered the headline. A suite that tightens its own null after the null
has already been passed is behaving correctly.

## 9.5 The result the suite bought: H7 found a defect nobody was looking for

The strongest argument for the whole exercise is not any of the nine verdicts. It is that H7, asking
an unrelated question, measured the quantity Chapter 5 would later find had been published wrongly
for months.

H6 had left 47 of 58 approved antiepileptics scoring exactly zero. H7 asked whether the responsible
models could rank at all, scoring each target's held-out actives against 600 random PubChem
structures using hold-out models. The answer is that they can. Across 37 testable targets the AUROC
runs from 0.9174 at NLRP3 to 0.9996 at HDAC1, median 0.9888, mean 0.9847
(`inversion/results/H7_target_discrimination.csv`). Not one falls near chance, so the
non-discriminative-model explanation is refuted.

What the same table records incidentally is deployed sensitivity, and it runs from 0.000 at P2X7 to
0.990 at mTOR, median 0.7922, mean 0.7331. Six targets fire for fewer than half their own held-out
actives: COX-2, GABA-A, KEAP1, P2X7, SIRT1 and TAAR1.

The panel registry at that time reported a mean sensitivity of 0.8983 with no deployed endpoint below
0.50, under a basis field reading `held_out_actives_by_scaffold`. Chapter 5 traced that to a
sequencing defect and established the held-out mean as 0.7638 with six endpoints below 0.50, five of
them the same six H7 names.

**H7 and the Chapter 5 audit share no construction.** H7 applies the deployed threshold to hold-out
models scored against random PubChem chemistry; the audit compares two scripts' active-selection
rules against the endpoint tables. They landed 0.031 apart on the panel mean, and the number they
agreed on was not the published one. The suite designed to embarrass the tool had recorded the
correct figure in a committed artefact for three weeks before anybody read it as a contradiction.

The lesson is about reading rather than about measurement. The evidence was present, in a committed
file, in a report the project had written itself. Nothing compared it against the registry, because
nothing had been built to.

## 9.6 The two hypotheses no other chapter treats

H4 and H5 are the suite's two untreated supported verdicts, and both are weaker than they read.

### H4: specificity outside the training neighbourhood

The null is precise. The measured specificity was an artefact of drawing negatives from the same
library the models were trained on, where every compound is in domain. The test scores 600 structures
drawn by random PubChem identifier, independent of every set used to build the tool, and stratifies
them by maximum Tanimoto to training chemistry
(`inversion/results/H4_distant_specificity.csv`):

| Stratum | n | False-positive rate | 95% interval | Median max Tanimoto |
|---|---:|---:|---|---:|
| In domain, T >= 0.50 | 75 | 0.1200 | 0.0644 to 0.2126 | 0.550 |
| Near domain, 0.30 to 0.50 | 464 | 0.0711 | 0.0511 to 0.0982 | 0.387 |
| Distant, T < 0.30 | 61 | **0.0164** | 0.0029 to 0.0872 | 0.275 |
| All | 600 | 0.0717 | 0.0536 to 0.0951 | 0.389 |

The comparator is the library-chemistry rate of 0.0750, 75 of 1,000, interval 0.0603 to 0.0930
(`results/tables/noncns_specificity_summary.csv`). It is read from that artefact at run time rather
than written into the script as a literal, so it cannot go stale, and that is good practice worth
noting.

The point estimate is a fourfold improvement on the comparator and the ordering is monotone in
distance. **But the verdict rests on a rule that ignores its own interval.** The script decides
`SUPPORTED if f <= ref_fpr * 1.5` (`inversion/inv_distant_specificity.py:162`), comparing two point
estimates. The distant stratum contains **one** false positive in 61 compounds, and its 95 per cent
interval, 0.0029 to 0.0872, **contains the comparator it is being compared against**. Tested rather
than compared:

| Comparison | Fisher exact, two-sided |
|---|---:|
| Distant against library chemistry | p = 0.119 |
| Distant against near domain | p = 0.161 |
| Distant against in domain | **p = 0.023** |
| Homogeneity across the three strata | chi-square = 5.44, 2 d.f., p = 0.066 |

So the defensible claim is narrower than the verdict. **Specificity does not degrade outside the
training neighbourhood**, which is what the null predicted it would do, and that is a real negative
result worth having. The stronger reading, that specificity improves with distance, rests on a single
count and does not separate from the comparator at conventional significance. The honest form is that
H4's null is not supported and H4's alternative is not established, on 61 compounds.

The in-domain row deserves a sentence of its own. At 0.1200 it is the **worst** of the three strata,
and it is higher than the 0.0750 measured on the non-CNS library, whose compounds are also in domain.
Random PubChem structures resembling training chemistry are harder for this system than curated
non-CNS drugs are. The intervals overlap, so this is a direction rather than a finding, but it points
the opposite way from the reassuring reading of the table.

### H5: read-across against a frequency baseline

Read-across recovers the true target in the top three of five similarity-weighted nearest neighbours
for **0.9726** of 4,092 held-out compounds across 41 targets, interval 0.9672 to 0.9772, against
**0.0587** for always answering with the three commonest targets in the index
(`inversion/results/H5_readacross_value.csv`). Per target the recall runs from 0.826 at Nav1.5 to
1.000, which nine of the 41 targets reach, with a median of 0.983
(`inversion/results/H5_readacross_per_target.csv`).

Two things constrain what this can be quoted for, and the second is not stated anywhere in the
project.

**The frequency null is bounded by the design of the evaluation set, not by the difficulty of the
task.** The set caps each target at 120 compounds, so with 41 targets contributing 4,092 compounds
between them, any constant answer naming three targets can be right at most about 360 times, roughly
0.088. The measured 0.0587 is exactly 240 of 4,092: the two commonest index targets that appear in
the evaluation set, D2 and A2A, each contributing its full block of 120, and nothing else
contributing at all. A null with a ceiling near 0.09 cannot distinguish a good retrieval method from
a mediocre one. Beating it by 0.91 is not evidence of quality; it is evidence that the null was the
wrong instrument.

**The index is not filtered by scaffold.** The query compound and any structure at Tanimoto >= 0.999
are excluded, and nothing else (`inversion/inv_readacross_value.py:76-79`). The compounds are held
out from the *models*, but the read-across index is the full deployed index, so a held-out compound's
own scaffold siblings are among the neighbours available to answer for it. The script's docstring is
candid about this, saying the scaffold is absent "only in the sense that the models never trained on
it", but neither the verdict line nor the report carries the qualification.

So H5 measures read-across **in its deployment regime**, where a query's chemical family is already
represented in the index, and 0.9726 is the right number for that regime. It says nothing about a
novel scaffold or an unrepresented target class and must not be quoted as if it did. The report
already carries the target-class half of this caveat. The scaffold half is missing and belongs beside
it.

## 9.7 The claim the suite never made: a permutation null

Chapter 8 closed by flagging that the technical report's null models have no artefact. Section 6.9
states that with labels permuted the same pipeline on the same folds returns a mean AUROC of 0.4938
random and 0.4921 scaffold, worst single endpoint 0.5174. Those three numbers are hard-coded inside a
literal prose block at `build_technical_report.py:1326`. No file in the repository held them and no
script computed them, which made this the load-bearing claim of the whole validation programme and
also the only one a reader could not check: it is the evidence that the cross-validated figures are
not inflated by leakage.

It has now been computed. `src/brainsafe/evaluation/permutation_null.py` imports `train_rf` rather
than reimplementing it, so the featuriser, the deduplication, the Bemis-Murcko grouping, the fold
objects and the forest hyper-parameters are identical by construction rather than by resemblance. One
permutation is drawn per endpoint from a seed derived from the endpoint name, preserving the class
balance, and the same ten-fold random and scaffold cross-validation is then run on the permuted
labels.

What a permuted label destroys, and what it deliberately does not, is the point of the test:

| | |
|---|---|
| Destroyed | any association between a compound's features and its class |
| Preserved | the class balance, the fold sizes, the scaffold grouping, and the tendency of whole scaffold classes to be sampled together |

The second row is the hypothesis under test. A scaffold-grouped fold could report a respectable AUROC
without any chemistry being learned, if scaffold classes differed enough in class frequency for the
grouping alone to be exploitable. Permuting the labels breaks the chemistry and keeps the grouping,
so a scaffold figure that stays at chance under permutation cannot be an artefact of the grouping.

The result is `results/tables/permutation_null.csv`, over the eight core classification endpoints:

| Endpoint | n | Random null | Scaffold null |
|---|---:|---:|---:|
| BBB | 3,901 | 0.4898 | 0.4917 |
| AChE | 5,125 | 0.4941 | 0.5071 |
| BChE | 3,278 | 0.4982 | **0.5199** |
| BACE1 | 8,207 | 0.4877 | **0.4830** |
| GSK-3-beta | 5,439 | **0.5200** | 0.5175 |
| MAO-A | 3,585 | 0.4892 | 0.5109 |
| MAO-B | 4,534 | 0.5015 | 0.4918 |
| hERG | 9,933 | **0.4868** | 0.4988 |
| **Mean** | | **0.4959** | **0.5026** |

**Every one of the sixteen values sits within 0.020 of chance.** The random mean is 0.4959 over a
range of 0.4868 to 0.5200; the scaffold mean is 0.5026 over a range of 0.4830 to 0.5199. Quoting the
means alone would hide that spread, and the spread is what a single permutation per endpoint costs:
the per-fold standard deviations run from 0.0197 to 0.0457, so an endpoint's null carries roughly
plus or minus 0.02 of sampling noise and the deviations observed are of that size and of both signs.
The panel mean is the reliable statistic here, not any individual row.

Two conclusions follow, and the second is the one the test was built for.

**The cross-validated figures are not obtainable without the labels.** Against each endpoint's own
permuted null, the margin runs from 0.4092 at BBB to 0.4887 at BACE1 under the random split, mean
0.4616, and from **0.3860 at BBB** to 0.4818 at BACE1 under the scaffold split, mean 0.4226. The
largest deviation any null shows from 0.500 is 0.0200, so the smallest margin in the panel is more
than nineteen times the largest departure from chance the null itself produces.

**Scaffold grouping alone confers nothing.** This is the specific leakage route the test exists to
close. Paired by endpoint, the scaffold null exceeds the random null by a mean of only 0.0067
(Wilcoxon p = 0.25, paired t-test p = 0.16), and the sign is negative for three of the eight
endpoints. Whole scaffold classes do not carry enough class-frequency information for a label-free
model to exploit the grouping, so the scaffold-split figures in Chapter 3 are not inflated by that
route.

One qualification, and it matters for how this section should be read. **The report's specific
figures are corroborated rather than reproduced.** Its 0.4938 random and 0.4921 scaffold sit 0.0021
and 0.0105 from the values measured here, and its worst single endpoint of 0.5174 sits 0.0026 from
the 0.5200 measured at GSK-3-beta. Those gaps are the size of the sampling noise above, so the two
runs agree, but the report records no seed and a permutation null is a random quantity, so its exact
triple cannot be recovered by anyone. The right correction is for section 6.9 to cite the artefact
rather than restate a number, which is the same fix Chapter 3 recommended for the censored-bound
count.

## 9.8 Where the suite is weaker than it reads

Six defects, found while writing this chapter. None changes a verdict. Together they describe a suite
that was built carefully and maintained less carefully than the pipeline it audits.

**1. The verdict rule for three hypotheses lives in two places, and for one of them the two
disagree.** Each test script prints its own verdict, and `inversion/summarise.py` recomputes every
verdict from the CSVs when it writes `REPORT.md` and `VERDICTS.csv`. For H1, H2 and H5 the two rules
use different label sets: the test script has a WEAKENED band that the summary drops. For H1 and H5
the outcome is the same either way. For H2 it is not. `inv_disease_layer.py:228` prints **WEAKENED**,
"curation adds little over uniform"; `summarise.py:60` publishes **REFUTED**. The published verdict is
the harsher of the two, which is the right direction to err in, but a verdict should not depend on
which of two copies of a rule a reader happens to run. This is the defect class that commit `d8da016`
fixed for the reliability gate, where the gate "lived in six places and disagreed with itself in
three".

**2. No test pins any verdict.** The repository's five test files contain no reference to
`inversion/`. Every other load-bearing quantity in this project is pinned by an assertion somewhere.
The falsification results are not, so a rerun that silently flipped a verdict would be caught only by
a person reading the report.

> **Resolved, items 3, 4 and 5.** All three were repaired on 2026-09-02 in commit `1f91925`, after
> this chapter was written. They are retained as found rather than deleted, because the record of a
> defect and the record of its repair are both evidence, and because the repair of item 4 turned up
> something worse than this chapter had noticed. Each item carries a note below saying what changed.

**3. `app.py` still quotes a withdrawn version of H2.** The comment at `app.py:102-104` reads: an
ablation "over 15,609 scaffold-held-out compounds found that curated, uniform and randomly permuted
weights give top-3 disease accuracy of 0.7917, 0.7911 and 0.7899". The artefact it cites now reads
0.7901, 0.7897 and 0.7874 on 7,008 compounds. The triple is the 4 August run and the artefact has
been regenerated six times since. The population is worse than stale: 15,609 never matched this
experiment at all. On the day the comment was written the artefact said 7,524, and
`audit/AUDIT_REPORT.md:824` identifies 15,609 as the denominator of a pooled-recall calculation, a
different quantity from a different table. The project's own audit recorded this as item BS-M-04,
"weight ablation, two triples for one experiment", and named `app.py` explicitly. It was corrected in
the manuscript and not in the source.

*Fixed.* The comment now states the current triple over 7,008 compounds and records what it used to
say and why the wrong denominator survived: 15,609 agrees with the H2 accuracy to three decimals, so
a wrong number sat beside a right-looking one.

**4. The interface quotes co-firing correlations that no longer match the artefact, and this one is
user-visible.** `FAMILY_COFIRE` at `app.py:529` supplies both a badge reading "correlated, r = 0.81"
and a sentence of explanation, rendered on a result page whenever two members of a homology family
both fire (`app.py:1563-1570`). Against `H8_family_correlation.csv` as it now stands:

| Family | Shown to the user | In the artefact |
|---|---|---|
| Dopamine D2-like | r = 0.81; D2 and D3 fire together for 78% of compounds engaging D2 | phi = 0.720; 0.717 |
| Opioid | r = 0.79; kappa fires for 65% of compounds engaging mu | phi = 0.798; 0.657 |
| Monoamine transporters | r = 0.64; NET fires for 84% of compounds engaging DAT | strongest pair is SERT and NET at phi = 0.732; the DAT-to-NET conditional is 0.618 |
| Serotonin receptors | r = 0.70; 5-HT7 fires for 78% of compounds engaging 5-HT2A | phi = 0.671; 0.659 |
| Neuronal nicotinic | r = 0.35; the least correlated family measured | phi = 0.498; still the least correlated of the five estimable families |
| Voltage-gated sodium | not estimable | correct: no pair has an estimable phi |

The accompanying comment says 36 targets fire spanning 14 independent directions, "roughly 2.6-fold
redundant", where the artefact says 37 and 16, a factor of 2.31. Every discrepancy is small and none
reverses a conclusion, but these are numbers presented to a scientist deciding how much weight two
engaged targets deserve, and the monoamine row names the wrong pair as the family's strongest.

*Fixed, and the fix found a worse defect than this chapter had.* `FAMILY_COFIRE` is now derived from
the artefact by a function rather than restated, so it cannot drift again. Deriving it required the
artefact to carry denominators, because every rate in it is a multiple of 1/400 and a conditional
probability of 1.000 computed from a single joint compound is indistinguishable in the file from a
well-supported certainty. That is exactly what the nicotinic row was: the interface asserted
"r = 0.35, the least correlated family measured", a family-level conclusion, on a pair whose joint
engagement rests on **one** approved drug. `inv_panel_independence.py` now records `n_drugs`, `n_a`,
`n_b` and `n_joint`, and a pair is reported only when at least ten approved drugs engage both.
Regenerating reproduced every pre-existing column exactly and left `VERDICTS.csv` unchanged.

**5. `results/tables/background_specificity.csv` is the output of a script that has since been fixed
and never re-run.** It still reports a mean sensitivity of 0.8983 across 47 endpoints and marks 46 of
them reliable. The registry, `final_thresholds.csv` and the reviewer package all say 0.7638 and 41 of
47. This is not a second occurrence of the Chapter 5 defect. The source was repaired on 27 August and
the registry patched directly on 30 August, deliberately, because re-running the threshold sequence
had regressed the project once before. The artefact was simply left holding the old output. The
freshness graph does not catch it because the graph's inputs are data and models, never code: a
generator can be corrected without anything marking its output stale. That is a real gap in a tool
this project relies on. It is stated here rather than fixed, because regenerating this file means
re-running the sequence the correction was designed to avoid.

*Fixed, without re-running the sequence.* The dilemma above was real and the resolution is narrower
than either horn of it. Four of the file's six columns were already correct; only the two reported
ones were stale. `refresh_background_specificity.py` updates those two from the registry and leaves
every operating threshold untouched, which is verifiable: after the repair no threshold and no
background false-positive rate moved on any of the 47 endpoints, and the file now agrees with the
registry on all 47 with zero disagreements. A byte-identical copy shipped in the submission package
and was updated with it, which matters more than the local file: **the disagreement was visible to
reviewers before it was visible to us.** Two tests now pin both, since the freshness graph cannot
catch this class for the reason recorded at `check_freshness.py:110`.

**6. Three pieces of explanatory prose have drifted from their artefacts.** The H7 script's docstring
says deployed sensitivity "ranges from 0.26 to 0.98"; the artefact says 0.000 to 0.990. The H4
comparator docstring says the distant-chemistry rate is 0.080; the artefact says 0.075. And
`inversion/H1H2H3.log`, `H4.log` and `H5.log` are the 4 August run, retained unchanged through four
regenerations, so the stored log records "VERDICT H2: WEAKENED" beside a published verdict of
REFUTED. None is load-bearing, and all three are the kind of drift that makes a reader distrust the
figures that are.

## 9.9 What a falsification suite cannot do

Three limits, stated because the suite is the strongest evidence in this thesis and overreading it
would be easy.

**It tests what its author thought to doubt.** All nine hypotheses were written by the person who
built the system. The suite is an excellent instrument against self-deception about things one has
noticed and no instrument at all against things one has not. Its own history illustrates the point:
H8 and H9 do not appear in `PLAN.md`, which stops at H7, and H9 was written only after H6 returned a
result its author judged unfairly scored. Three of the nine hypotheses exist because earlier ones
produced surprises, which is the right way for a suite to grow and also an admission that its initial
coverage was set by intuition.

**Nothing in it tests the barrier model.** H1, H2, H3, H6 and H9 test the disease layer; H5 tests
read-across; H7 and H8 test the target panel; H4 tests system-level specificity. The BBB classifier,
which gates every disease score and is the component the whole architecture is named for, has no
hypothesis. Its evidence is the external set in Chapter 8, AUROC 0.7666 on 227 unseen drugs, which is
a validation and not a falsification. The obvious missing test is a null asking whether the barrier
model's contribution to the gated score could be replaced by a molecular-weight or cLogP rule without
loss. That test does not exist and should.

**A verdict is a summary of one measurement, not a property of the system.** H4 is SUPPORTED on 61
compounds; H5 is SUPPORTED against a null with a ceiling near 0.09; H9 is SUPPORTED on a mean
per-indication AUROC of 0.6158 that conceals a range from 0.7981 for depression down to 0.4898 for
epilepsy. The verdict column is the least informative part of `VERDICTS.csv`, and any use of this
suite that reads the column rather than the tables will overstate what it found.

## 9.10 What this chapter establishes

- **The suite is the most valuable evidence in this thesis, and it is valuable in proportion to what
  it refuted.** Four of nine central claims did not survive. Two changed how the system is described,
  one changed the interface, and one contributed to withdrawing four endpoints that were firing on
  glucose, urea and lactate.
- **Its best result was incidental.** H7's sensitivity column recorded the correct panel figure three
  weeks before anything compared it with the published one, and it agreed with the eventual audit to
  within 0.031 by an entirely unrelated construction.
- **Two of its four supported verdicts are weaker than the word suggests.** H4 rests on one false
  positive in 61 compounds, with an interval containing its own comparator; the defensible claim is
  that specificity does not degrade with distance, not that it improves. H5 measures retrieval in a
  regime where the query's chemical family is present in the index, against a null the design of the
  evaluation set caps near 0.09.
- **The permutation null now exists.** The one claim in the validation programme that no artefact
  supported has been computed from the same pipeline and written to a file.
- **The suite is not maintained to the standard of the pipeline it audits.** A verdict rule that
  disagrees with itself, no test pinning any verdict, two stale citations in the served application,
  one of them user-visible, and an artefact holding the output of a script that has since been fixed.

---

## Outstanding items for this chapter

1. ~~Correct `app.py:102-104` to the current H2 triple and population.~~ **Done, 2026-09-02**
   (`1f91925`). Audit item BS-M-04 is now closed in the source as well as the manuscript.
2. ~~Derive `FAMILY_COFIRE` from `H8_family_correlation.csv` rather than restating it.~~
   **Done, 2026-09-02.** The artefact now carries denominators and suppresses any pair supported by
   fewer than ten approved drugs, which removed a family-level claim resting on one compound.
3. **Give each hypothesis one verdict rule.** H10 reads its verdict from its own artefact, which is
   the arrangement to adopt; H1 to H9 still derive theirs twice, and for H2 the two copies disagree.
4. **Pin the verdicts with a test**, so a rerun that flips one fails loudly. Still open: the test
   suite has grown to 73 and none of it reads `inversion/`.
5. ~~Add a hypothesis for the barrier model.~~ **Done.** H10, argued in section 10.3.
6. **Make the freshness graph aware of code.** An artefact whose generator has changed since the
   artefact was written is stale, and nothing currently says so.
7. ~~Regenerate `background_specificity.csv`.~~ **Done, 2026-09-02**, by refreshing the two reported
   columns from the registry rather than by re-running the threshold sequence. No threshold moved.
8. **Refresh or delete the three 4 August logs** in `inversion/`, which now record a superseded
   verdict.
