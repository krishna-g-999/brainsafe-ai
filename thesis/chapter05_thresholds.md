# Chapter 5. Thresholds and the three disjoint background pools

> Every quantitative statement in this chapter was computed during this session from the artefact
> named beside it. Section 5.8 reports a defect found while writing it, and the artefact that now
> quantifies that defect, `results/tables/sensitivity_reconciliation.csv`, was written during this
> session for that purpose.

---

## 5.1 A threshold is a separate problem from a probability

Chapter 4 ended with a calibrated probability and a coverage statement. Neither makes a call. Between
a probability and a reported engagement sits a threshold, and the threshold is where a panel of
independent models becomes a system that says yes or no.

It is also where the easiest self-deception in the whole project lives, and the design of this
chapter is organised around refusing it.

## 5.2 The circularity a single pool creates

Suppose a threshold is chosen as the 90th percentile of the scores a model gives to a sample of
background chemistry, so that ten per cent of that sample scores above it. Now measure the
false-positive rate of that threshold on the same sample. It is ten per cent. It was always going to
be ten per cent. The measurement restates the quantile and tests nothing about the model.

The failure is not hypothetical carelessness; it is the natural thing to do when one background
library is available and three jobs need doing. The binder pipeline asks that library for three
different things:

- **decoys**, trained on as presumed negatives where measured inactives are scarce,
- **a threshold sample**, whose score distribution sets the cut,
- **an evaluation sample**, on which the achieved false-positive rate is measured.

Drawing the threshold and the evaluation samples from one pool makes the reported rate an arithmetic
consequence of the quantile. Drawing decoys from the same pool as either is worse, because the model
was explicitly trained to score those compounds as negative and is then congratulated for doing so.

## 5.3 One partition, three roles

The background library of **158,890 compounds** is therefore partitioned once, into three pools that
no compound can belong to twice. Recomputing the assignment directly from
`src/brainsafe/models/pools.py` gives:

| Pool | Compounds | Share | Target share | Purpose |
|---|---:|---:|---:|---|
| Decoy | **95,515** | 60.11% | 60% | supply property-matched negatives during training |
| Threshold | **31,694** | 19.95% | 20% | set the decision cut |
| Evaluation | **31,681** | 19.94% | 20% | measure the false-positive rate the cut achieves |

Each threshold and each rate is computed on a sample of 3,000 drawn from the relevant pool.

The consequence is one sentence and it is the point of the chapter: **the measured rate can now
disagree with its target**, and that it can disagree is the evidence that it is a measurement rather
than a restatement. Section 5.6 reports by how much it does.

## 5.4 Why a hash and not a shuffle

Pool membership is assigned by a stable hash of the canonical structure rather than by shuffling.
A compound's role is `blake2b(salt + smiles) mod 100`, banded 0 to 59, 60 to 79, 80 to 99, with the
salt fixed at `brainsafe-background-v1`.

Two properties follow that a random shuffle does not give, and both matter more than they look.

**The assignment does not depend on the order the library happens to be in**, so it is reproducible
across machines and across runs without storing a seed alongside the data.

**Adding compounds later leaves every existing assignment untouched.** A threshold set today remains
comparable with a false-positive rate measured next year, because the library growing does not
reshuffle what came before. Under a shuffle, every expansion of the background library would silently
invalidate every previously measured rate.

The software suite pins both properties. `TestBackgroundPools` holds five subtests: that assignment
is a pure function of the structure and independent of call order, that every structure gets exactly
one role, that the bands cover exactly one hundred, and that the hash band is stable and in range.
These are not tests that the code runs. They are tests that a threshold set on one machine means the
same thing on another.

## 5.5 Two constraints, and which one binds

A threshold must satisfy two conditions at once, and the operating point is the stricter of them.

**On the target's own measured inactives**, no more than `TARGET_FPR = 0.10` may be called a binder.
The cut is the corresponding quantile of the scores the model gives to the half of the endpoint's
measured inactives that were withheld from fitting.

**On random library chemistry**, no more than `BACKGROUND_FPR = 0.05` may be called a binder.

The second condition exists because of a specific failure that per-target validation missed. Nav1.1
was trained on its own measured actives and inactives, both of which are sodium-channel screening
chemistry, and scored an AUROC of 0.979 against that narrow negative set. Shown arbitrary chemistry
it gave glucose a binder probability of 0.806. A model validated only against its own assay's
negatives can be confidently wrong about everything else. Nav1.1 was subsequently withdrawn, as
Chapter 6 records, but the constraint it motivated applies to every endpoint.

The 5 per cent level was itself tested. Tightening it to 2 per cent was tried and rejected, because
it cost too much sensitivity: it pushed the μ-opioid threshold to 0.999 and caused morphine to be
missed. That is the right kind of reason to reject a tightening, and it is recorded rather than
asserted.

Across the 47 deployed binder endpoints (`results/tables/final_thresholds.csv`) the measured-inactive
constraint binds on **39** and the background constraint on **8**. The measured-inactive cut is the
higher of the two on 39 of 47, with a median of 0.7353 against 0.0903 for the background cut. That
ordering is the same fact Chapter 2 reported from the other end: compounds a chemist thought worth
testing at a target are much harder to separate from real ligands than compounds drawn at random from
a library.

The threshold is finally clipped to [0.05, 0.999], so that no endpoint can be given a cut that admits
nearly everything or nothing at all.

## 5.6 The rate that is measured, against the rate that was targeted

This is what the partition buys.

| | On the threshold pool, in sample | On the disjoint evaluation pool |
|---|---:|---:|
| Mean | 0.1057 | **0.0228** |
| Median | 0.1020 | 0.0259 |
| Minimum | 0.0730 | 0.0001 |
| Maximum | 0.1880 | 0.0621 |
| Above 0.05 | most | 4 of 44 |
| Above 0.10 | most | **0 of 44** |

The in-sample rate sits near its 0.10 target by construction, though not exactly on it: it is exactly
0.10 on only 11 of 47 endpoints, because the clip and the discreteness of a quantile over a finite
sample move it. The rate measured on the disjoint pool is **lower on all 44 of 44** endpoints where
both are recorded, by a mean of 0.083.

The direction is systematic and it is not an error. It is the same asymmetry as section 5.5: for 39
of 47 endpoints the binding constraint is the measured-inactive quantile, which is a much stricter cut
than a background quantile would be, so the rate that cut achieves on ordinary library chemistry is
far below ten per cent. The two numbers answer different questions. The in-sample figure says *the
cut is where we put it*. The evaluation figure says *the panel fires on about 2 per cent of unrelated
chemistry*, and only the second is evidence.

The highest measured background rates are at the aminergic receptors, 5-HT2A at 0.0621, 5-HT1A at
0.0606, D2 at 0.0564 and D3 at 0.0544, which is the expected place for them: those targets bind a
large, chemically ordinary region of drug-like space, so more of a random library genuinely resembles
their ligands. The lowest are at Cav3.2, 0.0001, Nav1.8, 0.0005, and the CGRP receptor, 0.0008.

## 5.7 A second operating point

The deployed threshold is a triage cut. A second, stricter cut is derived and stored for every
endpoint (`results/tables/screening_thresholds.csv`), holding the background false-positive rate at
0.01 rather than 0.05, for the different task of screening a large library where the cost of a false
lead is multiplied by the size of the library.

The price is stated per endpoint. At A1 the triage cut of 0.8083 gives a sensitivity of 0.736 and the
screening cut of 0.9818 gives 0.523; at A2A, 0.8888 and 0.824 against 0.9979 and 0.344. Roughly a
third of the sensitivity is spent to reduce the background rate fivefold. Having both, and reporting
what each costs, is more honest than choosing one and calling it the operating point.

## 5.8 One quantity, four numbers, and a mislabelled set

Writing this chapter surfaced the most serious documentation defect found in this thesis so far. It
concerns the sensitivity reported at the deployed threshold, which is the number that tells a user
what a silence is worth, and it is stated here in full because the correction changes a published
figure.

Four artefacts report sensitivity at the deployed threshold for the binder panel, and they disagree:

| Source | n | Mean | Median |
|---|---:|---:|---:|
| `models_rf/binder_modes.json`, the panel registry, and everything quoted from it | 47 | **0.8983** | 0.9320 |
| `results/tables/final_thresholds.csv` | 47 | **0.7638** | 0.8350 |
| `inversion/results/H7_target_discrimination.csv` | 37 | 0.7330 | 0.7922 |
| `results/tables/integrity_audit.csv`, panel mean | – | 0.846 | – |

The registry's figure is the one the manuscript, the evidence map, the technical report and Chapter 1
of this thesis all quote, and it is **the highest of the four**.

The cause is a sequencing defect, not a calculation error. The threshold sequence runs
`final_thresholds.py` and then `calibrate_background_specificity.py`, and **both write
`sensitivity_at_threshold` into the registry**. The second overwrites the first. But the two select
their actives differently:

- `final_thresholds.py` scores the **held-out actives** listed in `models_rf/holdout/<ep>_binder_holdout.json`. For A1 that is 420 compounds.
- `calibrate_background_specificity.py` scores **every active in the endpoint table** at pChEMBL ≥ 7. For A1 that is 1,943 compounds, of which **1,523 were used to fit the model**.

The registry keeps the second script's number and the first script's label. Every deployed entry
still carries `sensitivity_basis: "held_out_actives_by_scaffold"`, and on all 47 of 47 the number
beside it is not that.

`results/tables/sensitivity_reconciliation.csv` quantifies the gap per endpoint:

- The published figure is higher on **47 of 47** endpoints, by a mean of **+0.1344**.
- The set it is measured on is **80.5 per cent training compounds at the median** over the 44
  endpoints where both scripts agree on what an active is.
- The largest gaps are SIRT1, held out 0.378 against a published 0.967, a gap of +0.589; GABA-A,
  0.303 against 0.826; TAAR1, 0.355 against 0.763; GluN2B, 0.374 against 0.737.

A second inconsistency surfaced in the same pass. Three endpoints, **Nav1.5, SIRT1 and TAAR1**, define
"active" differently in the two scripts: the binder pipeline selects training actives by the `label`
column while the specificity script selects by pChEMBL ≥ 7. Nav1.5 has 242 actives by label and 40 at
pChEMBL ≥ 7, so its published sensitivity is measured on a set built by a different rule from the one
the model was fitted with.

## 5.9 What the honest figure is, and what it changes

The defect described in section 5.8 has been corrected at the source: `models_rf/binder_modes.json`
now stores the held-out figure on all 47 deployed endpoints, and `sensitivity_reconciliation.csv`
reports a gap of exactly zero across the panel. The fix was made in the script that computes the
figure rather than by patching the stored value, because the same correction had been applied once
before by hand, reaching a mean of 0.7513, and was then silently undone the next time the calibration
stage ran. Section 5.8 is retained as the record of how the defect was found; the figures below are
the current, published ones.

Measured on held-out actives only, across the 47 deployed binder endpoints:

| | Published | Held out |
|---|---:|---:|
| Mean sensitivity | 0.8983 | **0.7638** |
| Median | 0.9320 | 0.8350 |
| Minimum | 0.639, COX-2 | **0.303, GABA-A** |
| Maximum | 0.997, CGRP | 0.993, CGRP |
| Endpoints below 0.50 | **0 of 47** | **6 of 47** |

The last row is the one that matters, and it is the reason this is a correction rather than a
quibble. On the published figure no deployed endpoint fires for fewer than half its own actives. On
the held-out figure six do: GABA-A at 0.303, COX-2 at 0.338, TAAR1 at 0.355, GluN2B at 0.374, SIRT1
at 0.378 and P2X7 at 0.465.

**The falsification suite already found this, from an entirely different direction.** H7 asked whether
the silent antiepileptics were explained by non-discriminative targets, concluded they were not, and
reported that the cause was the operating point with a median deployed sensitivity of 0.79 and six
targets under 0.50. Its six are COX-2, GABA-A, KEAP1, P2X7, SIRT1 and TAAR1, five of which are the
same six found here. H7 scored held-out actives against random PubChem chemistry, an unrelated
construction, and landed 0.031 from the held-out mean. The published figure was the outlier all
along, and the suite designed to embarrass the tool had already said so.

The consequences are three, and none of them is that the system is worse than believed.

**The panel's behaviour is unchanged.** No threshold moves and no prediction changes. What changes is
the number that describes them.

**The manuscript's qualifier is false and must be corrected.** It reads "a mean sensitivity of 0.898
on actives withheld by scaffold". The mean of 0.898 is not on actives withheld by scaffold. The
figure on actives withheld by scaffold is 0.764.

**Chapter 1 of this thesis and Chapter 6 must both be revised**, and the claim that silence at a
close analogue is reasonably strong evidence of inactivity has to be re-examined against 0.764 rather
than 0.898.

## 5.10 The endpoints that are marked unreliable

The threshold sequence marks an endpoint unreliable when its sensitivity at the operating point falls
below `MIN_SENS = 0.50` or its AUROC against measured inactives falls below 0.75. Six of the 47
deployed endpoints carry that flag: COX-2, GABA-A, GluN2B, P2X7, SIRT1 and TAAR1.

Those are, with one substitution, the same six that section 5.9 identifies as firing for fewer than
half their own held-out actives. So the flag is doing its job, and the panel already knows which
endpoints are weak. What the panel does not do is stop deploying them, and that is a defensible
choice only because the flag travels with the prediction. An unreliable endpoint reporting nothing is
not evidence of inactivity, and the interface has to say so at exactly those six.

## 5.11 What the tests pin, and what they missed

`TestThresholdSequenceIsAtomic` holds two subtests: that the four threshold steps stay one unit, and
that no member of the sequence depends on a file the sequence itself rewrites. Both pass. They pin
the *ordering* of the sequence, which is the property that stops a threshold being set from a file
that a later step will change.

They do not pin what section 5.8 found, which is that two steps in a correctly ordered sequence write
the same field from different populations, and that the label written by the first survives the value
written by the second. A test for that is straightforward and is the specific piece of work this
chapter recommends: assert that `sensitivity_basis` describes the set the stored
`sensitivity_at_threshold` was actually computed on, and that the two scripts agree on what an active
is. Four of the 47 thresholds also differ between the registry and `final_thresholds.csv`, at D3,
5-HT1A, Nav1.5 and SIRT1, the registry being higher in each case because it holds the later value;
that is correct behaviour, but nothing asserts it.

---

## Outstanding items for this chapter

1. ~~The published binder sensitivity is measured on a set that is mostly training compounds.~~
   **Done, 2026-08-30, and propagated.** The registry now stores the held-out figure on all 47
   endpoints and the reconciliation reports a gap of +0.000. Both the manuscript and
   `docs/TECHNICAL_REPORT.md` now state the corrected mean of 0.764 directly and name 0.898 as the
   superseded published figure; `EVIDENCE_MAP.md` was not found at the top level of the repository
   during this pass and could not be checked.
2. ~~`calibrate_background_specificity.py` should score held-out actives.~~ **Done at the source**,
   so the overwrite cannot recur.
3. Nav1.5, SIRT1 and TAAR1 define "active" differently in the two scripts. One definition should win.
4. A test should assert that `sensitivity_basis` matches the population the stored sensitivity was
   computed on.
5. `results/tables/sensitivity_reconciliation.csv` is new and is not yet declared in
   `tools/check_freshness.py`.
