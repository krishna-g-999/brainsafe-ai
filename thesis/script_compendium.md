# The scripts: what each one does, and why it exists

> A guide to every script in the live BrainSafe AI codebase, written for a reader who needs to know
> what produced a number and what would have to be re-run to change it. Counts and descriptions were
> taken from the files themselves on 7 September 2026, at commit `d25245b`. Where a script carries a
> caveat that a reader needs in order to use its output correctly, the caveat is given here rather
> than left in the source.

## 0. How to read this

The repository holds **197 live Python files**. That is a large number for one project, and the
reason is worth stating at the outset, because it is the organising principle of the whole codebase:
**almost nothing here is a general-purpose library.** Each script produces one named artefact, states
in its own docstring what that artefact means and what it does not mean, and writes it to
`results/tables/` or `models_rf/`. Documents then read the artefact rather than restating it. The
scripts that check the documents against the artefacts are themselves in this list.

Four counts, all verified during the audit that accompanies this document:

- **52 endpoints** in the registry, **47 of them deployed**.
- **252 model files** in the manifest, all checksums verified.
- **229 quantitative claims** across the thesis are pinned to an artefact and machine-checked.
- **81 tests**, plus 57 subtests, pass on the current tree.

A further 75 files sit under `_ARCHIVE_2026-08-10_Brainsafe/`. They are the superseded engine and are
not part of this system. Nothing in the live tree imports them.

## 1. The pipeline, in the order it runs

```
  fetch_*.py            measured activities from ChEMBL, BindingDB, PubChem, B3DB, NPASS
        |
  rebuild_endpoints.py  pool the sources, resolve replicates, apply the potency cuts
        |
  featurize.py          desalt, neutralise, ECFP-4 (1,024 bits) + 12 descriptors = 1,036 columns
        |
  train_rf.py           8 classifiers + 5 regressors, 10-fold random AND scaffold-grouped
  train_binders_*.py    47 binder endpoints, hybrid decoys plus measured inactives
        |
  calibrate.py          isotonic calibration on out-of-fold predictions
  final_thresholds.py   operating thresholds from held-out inactives and a disjoint background pool
        |
  evaluation/*.py       external, prospective, temporal, specificity, conformal, null
  inversion/*.py        H1 to H10: ten hypotheses stated so that they could fail
        |
  app.py                the served predictor
        |
  check_freshness.py    refuse to commit if any artefact is older than its inputs
```

Reproduction is one command: `tools/reproduce.py`.

---

## 2. Data acquisition and curation (`src/brainsafe/data`, 33 files)

This is the largest single group, which is the correct proportion for a project whose central claim is
that its labels are measured rather than annotated.

### The fetchers

| Script | What it retrieves |
|---|---|
| `fetch_endpoints.py` | ChEMBL activities for the core target panel, and the B3DB barrier set |
| `fetch_bindingdb.py` | Binding affinities from BindingDB, reporting the net-new yield per target |
| `fetch_batch2.py` … `fetch_batch5.py` | Four target expansions: ALS, Huntington, neuroinflammation, epilepsy, sleep, pain, migraine, multiple sclerosis |
| `fetch_new_targets.py` | Activities for the new brain-target endpoints |
| `fetch_antioxidant.py` | The DPPH radical-scavenging regression set |
| `fetch_pka.py` | Measured pKa, separating basic from acidic ionisation |
| `fetch_pubchem_inactives.py`, `fetch_pubchem_expansion.py` | Measured **inactives** for targets whose negative class is thin |
| `fetch_pgp_substrate.py` | P-glycoprotein substrate data, kept distinct from inhibition |
| `fetch_natural_products.py`, `ingest_npass.py` | Natural-product activity against panel targets, from NPASS 3.0 |
| `fetch_clinical.py` | ATC class N molecules that reached a clinical phase, as the CNS precedent reference |
| `fetch_readacross_targets.py` | Additional targets to widen read-across beyond the modelled panel |

`find_targets.py` deserves separate mention. It resolves ChEMBL target identifiers **by name search,
never by a memorised identifier**. That is a guard against the single most common silent error in
cheminformatics data collection, which is confidently fetching the wrong protein.

`_tls.py` handles verified HTTPS on a network whose interception certificate is malformed but trusted.
It is infrastructure, not science, but it is why the fetches are reproducible on the author's network.

### Assembly and curation

- **`rebuild_endpoints.py`** is the centre of this group. It pools ChEMBL and BindingDB evidence,
  resolves replicate measurements, applies the pChEMBL cuts, and writes the endpoint tables that
  everything downstream reads. If one script had to be understood before the others, it is this one.
- **`build_compound_library.py`** builds the library from measured sources and owns `standardise()`,
  the InChIKey routine the leakage checks use.
- **`merge_pubchem_inactives.py`**, **`rebuild_with_inactives.py`**, **`recover_expansion_inactives.py`**
  merge the measured negative class into the endpoint tables. The first is explicitly
  *similarity-matched*, which matters: an earlier bulk merge raised GSK-3β scaffold AUROC from 0.9369
  to 0.9891 and was reverted, because the added negatives sat at a median Tanimoto of 0.288 to the
  actives and the model was learning to recognise a chemotype rather than an activity.
- **`integrate_external.py`** builds the external test set and computes the novelty flags.
  **This script currently needs re-running**: its flags were computed before `parent_mol` began
  neutralising, so 14 compounds it marks novel are in fact feature-identical to a training row. See
  Finding A of the audit report.
- **`build_np_endpoints.py`** builds endpoint tables for the targets added to close the
  natural-product gap.
- **`package_caches.py`** packages the raw API caches as one citable archive with a verifying manifest.

### The decisions that were made by measurement

Four scripts exist only to answer a question that could have been answered by assertion, and are worth
knowing about for that reason alone:

- **`validate_extraction.py`** measures the error rate of literature extraction *before any of it is
  trusted*.
- **`verify_pathways.py`** checks every target-to-pathway assertion in the knowledge graph against KEGG
  and Reactome, so the disease layer is not built on remembered biology.
- **`assess_unmodelled_targets.py`**, **`survey_np_targets.py`**, **`assess_state_dependence.py`**,
  **`assess_remaining_gaps.py`** ask how much measured data actually exists for mechanisms the panel
  does not cover. They are the evidence behind which targets became endpoints and which did not.
- **`audit_batch4.py`** audits six new endpoints *before* any of them is trained.

---

## 3. Representation (`src/brainsafe/features`, 2 files)

**`featurize.py`** is the component every model depends on. It converts a SMILES string into 1,036
numbers: a 1,024-bit ECFP-4 (Morgan, radius 2) fingerprint plus twelve interpretable descriptors
(molecular weight, cLogP, TPSA, H-bond donors and acceptors, rotatable bonds, aromatic rings, fraction
sp3, ring count, heavy atoms, formal charge, QED).

Two properties of it govern how its outputs must be read.

**`parent_mol` desalts, neutralises and sanitises.** This is not cosmetic. Before neutralisation,
haloperidol hydrochloride written the way PubChem serves it scored BBB 0.613 against 0.993 for the free
base, and hERG 0.295 against 0.914, so a user pasting the salt form lost a cardiac liability flag on a
compound that has one. Permanent charges are deliberately left alone, so choline and neostigmine keep
their quaternary nitrogen, which is exactly the property that stops them crossing the barrier. Of
170,617 unique structures, 198 change representation and 1,155 charged ones are correctly untouched.

**The fingerprint is folded, therefore lossy.** Folding is a modulo, so distinct atomic environments
collide. On 20,000 structures from the project's own tables, **52,882 distinct environments map onto
1,024 bits, and every one of the 1,024 carries more than one**, a median of 52 and a maximum of 73
(`fingerprint_collisions.csv`). The consequence is for interpretation: a bit a tree ranks important
names a *set* of environments, so a SHAP attribution to a single bit cannot be read as a substructure.
The representation is also stereo-blind, which is what makes deduplication necessary before any split.

`encodings.py` holds reversible integer encodings for non-numeric metadata. That metadata never enters
the feature matrix.

---

## 4. Model fitting (`src/brainsafe/models`, 18 files; `src/brainsafe/gnn`, 3; `src/brainsafe/adme`, 4)

### The core panel

**`train_rf.py`** trains the eight classification endpoints (BBB, AChE, BChE, BACE1, GSK3B, MAO_A,
MAO_B, hERG) and five regression endpoints (D2, A2A, HT2A, SERT, antioxidant DPPH). Every model is
300 trees, `min_samples_leaf` 2, `class_weight="balanced"`, `random_state=42`, and every one is
cross-validated **twice**: 10-fold random (stratified) and 10-fold scaffold-grouped.

Three functions inside it do the scientific work and are reused across the codebase:

- **`_scaffold_groups`** assigns a Bemis-Murcko scaffold group per compound, computed on the same
  desalted parent the featuriser uses so that a salt and its free base cannot land in different folds.
  Acyclic compounds share one group rather than getting one each, because giving each its own group
  quietly turns the scaffold split back into a random split for that part of the set.
- **`_dedup_features`** collapses rows whose feature vectors are byte-identical, which is necessary
  because the featuriser is stereo-blind. Classification groups whose members disagree on the label are
  dropped rather than voted on. BBB goes from 7,807 rows to 3,901 under this step, so it is not a
  marginal correction.
- **`_cv`** runs both schemes and saves every compound's out-of-fold prediction to
  `data/processed/cv_predictions/`, which is what lets calibration and the null be computed later
  without refitting.

### The binder panel

**`train_binders_hybrid.py`** trains the 47 deployed binder endpoints on a hybrid negative class:
property-matched decoys drawn from a background pool at a similarity ceiling, **plus** experimentally
measured inactives. Actives are held out **by scaffold**, a fifth of the scaffold groups never entering
training, so reported sensitivity measures recall of structurally distinct compounds rather than of
memorised ones. Both halves of every split are written to `models_rf/holdout/` so that the later
threshold step cannot silently reach past the holdout by re-reading the endpoint table.

*Known asymmetry:* the measured inactives are split in half **at random**, not by scaffold, which is
Finding B of the audit report. The memorisation this permits is measurable in 18 of 40 testable
endpoints, but it produces no demonstrated bias in the delivered thresholds (median quantile shift
-0.0315, upward on only 16 of 40).

Companions: **`train_receptor_binders.py`** (D2, A2A, HT2A, SERT), **`train_new_binders.py`** (14 new
brain targets), **`train_batch2.py`**, **`train_measured_label_holdout.py`**, and the two small
single-target scripts **`_train_nav17.py`** and **`_train_ox.py`**.

**`train_pka.py`** trains the most-basic-pKa regressor used by the CNS MPO desirability score.

### Calibration, pools and thresholds

- **`calibrate.py`** fits isotonic calibration per classification endpoint and reports the honest
  improvement, measured on already-saved out-of-fold predictions under an inner 5-fold so nothing is
  scored on data it was fit on. Mean expected calibration error falls from 0.0801 raw to 0.0151.
  Deployment models wrap the forest in `CalibratedClassifierCV`.
- **`pools.py`** partitions the background chemistry into three disjoint pools by
  `blake2b(salt + smiles) mod 100`, split 60/20/20, **so that no two roles can draw the same compound**.
  One pool supplies training decoys, one sets thresholds, one measures the false-positive rate. This is
  the single most important structural guard in the project: without it, the reported specificity would
  be a restatement of the threshold that produced it.
- **`final_thresholds.py`** sets every binder threshold from two independent requirements and takes the
  stricter: at most 10 per cent of the target's **held-out** measured inactives may be called binders,
  and at most 5 per cent of the threshold pool may be. The rate is then reported on the *evaluation*
  pool, not the threshold pool, because a quantile of the threshold pool cannot exceed its own target
  and reporting it would say nothing. Thresholds are recomputed from scratch each run so repeated
  calibration cannot ratchet them upward.
- **`calibrate_binder_thresholds.py`**, **`calibrate_background_specificity.py`**,
  **`screening_thresholds.py`** (a second, high-precision set for low-prevalence screening),
  **`apply_specificity_decisions.py`** act on the specificity audit.
- **`package_models.py`** and **`compress_models.py`** package and recompress the deployed models,
  the latter verifying that **not one prediction changes**.

### The graph network and the ADME layer

**`gnn/`** holds a Graph Isomorphism Network in pure PyTorch (`gin_model.py`, `graph_features.py`,
`train_gnn.py`). `train_gnn.py` compares it to the random forest **on the same split**, which is what
makes the comparison mean anything.

**`adme/`** trains the exposure layer with the same protocol as the target panel (`train_adme.py`),
and `cns_exposure.py` combines target, barrier and ADME models into the interpretable free-exposure
readout that produces the reported Kp,uu.

---

## 5. Evaluation and validation (`src/brainsafe/evaluation`, 33 files)

The design principle is that a cross-validated score is not evidence of generalisation, so each script
here tests a different way the panel could be wrong.

### Independence from training

| Script | Question |
|---|---|
| `external_validation.py` | How does the panel do on FDA-curated drugs held entirely out of training? |
| `external_100.py`, `external_100_stats.py` | Prospective validation on compounds absent from every training set |
| `external_prospective.py`, `external_prospective_core.py` | What would the panel have said about chemistry that did not yet exist? |
| `external_cross_source.py` | Does the panel survive a change of curator? |
| `external_natural_products.py` | How does it behave on chemistry it has not seen a class of? |
| `external_novelty_strata.py` | Does accuracy survive when the test compound does not resemble training? |
| `scaffold_holdout_panel.py`, `scaffold_holdout_report.py` | Retrain the whole binder panel with scaffolds withheld |

`external_validation.py` reports the external result under **both** novelty criteria and says which is
which, on the explicit reasoning that excluding overlap by InChIKey is not enough to call a compound
unseen, because the InChIKey separates stereoisomers and salts while the model's features do not. That
reasoning is right; the flags it relies on are currently stale (Finding A).

### Uncertainty and domain

- **`rf_conformal_temporal.py`** produces the Mondrian conformal prediction sets at 90 per cent and the
  temporal validation. It reports `frac_ambiguous`, `frac_empty` and `frac_singleton` separately,
  because on two classes the mean set size is `1 + P(ambiguous) - P(empty)` and the excess over 1.0
  cannot be read as ambiguity unless no set is empty. It has a scaffold arm in
  `rf_conformal_scaffold.csv`.
- **`applicability_domain.py`**, **`applicability_measures.py`** measure distance from training
  chemistry. *Caveat:* the domain does not separate more accurate predictions from less accurate ones
  on the external set (in-domain AUROC 0.7597 against out-of-domain 0.8059, difference not significant).
  It does correlate with the size of the probability error (Spearman -0.1422, p = 0.0128). Read it as a
  statement about calibration, not discrimination.

### Nulls, baselines and comparisons

- **`permutation_null.py`** is the label-permutation null: what the pipeline scores when there is
  nothing to learn. The worst deviation from chance across the panel is **0.0200** of AUROC, so no
  part of the pipeline manufactures signal from a shuffled label.
- **`model_comparison.py`** puts five families on the same features and the same folds: random forest,
  XGBoost, histogram gradient boosting, kNN read-across, and scaled logistic regression.
- **`model_family_significance.py`** applies Wilcoxon signed-rank to that comparison **and reports the
  minimum attainable p at each sample size**, so a unanimous result at n = 5 is labelled "underpowered"
  rather than "not distinguishable". This distinction cuts both ways: it stops the project claiming a
  null against XGBoost, which is higher on 5 of 5 regression endpoints.
- **`feature_analysis.py`** is the feature-block ablation. Its result is against the project's own
  interest and is reported anyway: over 13 endpoints the twelve descriptors add a mean of **+0.0014**
  over the fingerprint alone. They are kept because they cost nothing and the exposure layer needs
  them, not because they were shown to help. The converse, +0.1049, is what justifies the fingerprint.
- **`learning_curve.py`** asks, scaffold-honestly, what more data would buy.

### Specificity and integrity

- **`noncns_specificity.py`** and **`noncns_specificity_fast.py`** test roughly 1,000 compounds with no
  known engagement of any modelled target. This exists because a prospective test found Nav1.1 giving
  glucose a binder probability of 0.806 while scoring 0.979 against its own assay's negatives: a model
  judged only on its own chemistry can be confidently wrong about everything else.
- **`integrity_audit.py`**, **`provenance_audit.py`** ask which reported results are still true of the
  model that is actually deployed.
- **`sensitivity_reconciliation.py`** reconciles the four different figures this project reports for
  deployed sensitivity, which is exactly the kind of thing that otherwise becomes a viva ambush.
- **`stereochemistry_audit.py`** measures what ignoring chirality actually costs.
- **`fingerprint_collisions.py`** measures the folding collision rate.
- **`audit_inactives.py`**, **`audit_expansion.py`** test whether a data change was a genuine
  improvement or an easy-negative artefact.
- **`binder_cv_per_fold.py`**, **`binder_vs_measured_inactives.py`**,
  **`refresh_background_specificity.py`**, **`deployed_specificity_audit.py`**, **`app_health.py`**,
  **`library_sp3_coverage.py`**, **`validate_inversion.py`** complete the group.

---

## 6. The falsification suite (`inversion`, 10 files)

Ten hypotheses stated so that they could fail, each in its own script, each writing its own verdict
into its own artefact so that `summarise.py` reads rather than re-derives it. **Four are refuted and
two weakened**, and those verdicts are reported rather than buried.

| Script | Hypothesis |
|---|---|
| `inv_disease_layer.py` | H1 to H3: does the disease layer carry information, and are its design choices justified? |
| `inv_distant_specificity.py` | H4: does the reported specificity survive outside the training neighbourhood? |
| `inv_readacross_value.py` | H5: does read-across beat a frequency baseline? |
| `inv_clinical_indication.py` | H6: do the disease scores predict real clinical indications? |
| `inv_target_discrimination.py` | H7: what does silence mean? |
| `inv_panel_independence.py` | H8: does a larger panel add information, or only more chances to fire? |
| `inv_disease_discrimination.py` | H9: does the disease layer discriminate, or only echo which indications are common? |
| `inv_barrier_necessity.py` | H10: does the barrier model earn its place over a descriptor rule? |

`inv_barrier_necessity.py` is worth studying as a model of how these should be written. Its verdict
rule reads the **bootstrap interval**, not the point estimate, explicitly because Chapter 9 criticises
H4 for deciding SUPPORTED on a point comparison whose interval contained its comparator. A hypothesis
added in response to that criticism must not repeat it. Applying its own standard returned WEAKENED
rather than SUPPORTED.

`summarise.py` collects the verdicts into one table and a written report, and is idempotent:
regenerating it from unchanged artefacts produces a byte-identical file. `audit_mesh_map.py` audits the
MeSH-to-panel mapping H6 depends on.

---

## 7. Independent reproduction (`validation/repro`, 7 files; `audit/evidence`, 5 files)

These re-derive the results **without using the pipeline's own code**, which is the only way to catch a
bug that is shared between a result and its check.

- **`r00_environment.py`** captures everything needed to say what produced a number, before producing
  any.
- **`r01_leakage.py`** checks split integrity independently of the pipeline's claims.
- **`r02_recompute_cv.py`** re-runs the core cross-validation from the endpoint tables and scores it
  independently.
- **`r03_ledger.py`** builds the reproduction ledger: every reported number against an independently
  produced one.
- **`r04_null_models.py`**, **`r05_calibration_importance.py`**, **`r06_shap.py`** complete the ladder.

Under `audit/evidence`: **`leak_external.py`** (is the 306-compound external set really disjoint?),
**`leak_internal.py`** (within-endpoint duplicates), **`repro_cv.py`** and **`repro_scaffold.py`**
(reproduce as published, then repeat with duplicates collapsed), and **`scaffold_mismatch.py`** (do
feature-identical compounds always land in the same scaffold group?).

*Note:* `leak_external.py` computes its Tanimoto overlap on raw SMILES rather than on the desalted,
neutralised parent the model sees, which is why it reports 65 rather than 79. Its own
InChIKey-skeleton criterion, printed directly above, gives 78 and is nearly right.

---

## 8. Serving (`app.py`, `api.py`, `serve.py`, `model_fetch.py`, `panel.py`)

**`app.py`** (3,680 lines) is the interactive predictor and the canonical entry point. It loads the
**calibrated** models where they exist. The scoring transforms a reader is most likely to misread are:

- **`enrichment(ep, p)`** converts a calibrated probability into signed enrichment over the endpoint
  base rate, piecewise linear with a kink at the base rate: `(p - b)/(1 - b)` above it and `(p - b)/b`
  below. It is continuous at `p = b`, returns +1 at certainty and -1 at zero, and is strictly
  increasing, so it cannot reorder compounds within an endpoint. It exists because a probability below
  prevalence is evidence of inactivity, not of weak engagement.
- **`disease_scores(r)`** takes the **strongest** engaged mechanism per condition and multiplies by
  predicted barrier penetration. Because the barrier term multiplies every gated disease identically it
  **cannot change their relative order**; it decides whether anything is reported at all. Conditions in
  `PERIPHERAL_MECHANISM_DISEASES` are exempt, because for them the assumption the gate encodes is false.
  Describing the gate as sharpening the disease call would overstate what it does.
- **`target_signal`** routes each target to the right readout: enrichment for measured-label
  classifiers, a threshold-referenced score for decoy-aware binders, and a percentile for the
  antioxidant model.

**`panel.py`** is the single place that says what the server is made of, and `verify()` checks that the
served panel, the registry and the knowledge graph agree. **`api.py`** exposes the same models
programmatically, **`serve.py`** runs both in one container, and **`model_fetch.py`** ensures the model
files are present before anything tries to load them.

---

## 9. Reporting, figures and thesis machinery

**`src/brainsafe/analysis`** (27 files) builds the technical report, the reviewer workbooks, the audit
package, the manuscript and its tables. Two are guards rather than builders:
**`check_manuscript_numbers.py`** compares every headline number in the manuscript against the artefact
that produces it, and **`verify_references.py`** with **`verify_references_strict.py`** resolve every
citation against Europe PMC and CrossRef by exact title, **so that no reference is written from
memory**.

**`src/brainsafe/figures`** (16 files) renders eleven manuscript figures and four thesis figures
(T1 the scoring pipeline, T2 enrichment and the exposure gate, T3 the uncertainty stack, T4 the
falsification suite), all to the NAR house style defined in `style.py`: Arial, a 6.5 pt floor, 89 mm
and 183 mm column widths, 400 dpi, PNG and PDF. Every figure reads its values from an artefact.

**`thesis/`** holds `build_docx.py` and **`verify_chapter_numbers.py`**, the latter being the script
that checks all 229 quantitative claims in the chapters against the artefacts they cite. It is the
reason a stale number in a chapter fails a build rather than reaching a committee.

**`thesis/presentations`** (12 files) builds the ten defence decks and the master deck, with
`validate_deck.py` checking each for layout defects a build script cannot see.

---

## 10. The guards (`tools`, 8 files; `tests`, 5 files)

- **`check_freshness.py`** (585 lines) refuses to let any derived artefact be older than something it
  was derived from, and `install_hooks.py` makes it unskippable via a git pre-commit hook. It currently
  verifies 252 model entries with zero checksum mismatches. *Coverage gap:*
  `data/external/processed/external_bbb_test.csv` is declared only as an input, never as a derived
  artefact, which is how Finding A survived.
- **`build_provenance.py`**, **`build_evidence_map.py`**, **`build_start_here.py`** generate the
  reviewer's entry point: every claim, its number, and the file that produced it.
- **`reproduce.py`** regenerates the results in one command on any platform.
- **`tests/`** holds five suites, and the one to know about is
  **`test_pipeline_invariants.py`**, described in its own docstring as "invariants the scientific
  results depend on, **each of which has been violated at least once**". That is the honest reason a
  test suite exists. `test_panel_app_consistency.py` (407 lines) pins that the served panel is the
  registered panel and that the graph names only models that exist.

---

## 11. What to say if asked "why so many scripts?"

Because each one is an answer to a question that could otherwise have been answered by assertion. The
count is high for the same reason the claim count is high: **the project chose to spend its effort on
the correctness of the evaluation rather than the performance of the estimator.**

That choice is defensible and should be stated first rather than conceded under questioning. It has a
cost, which is that no hyper-parameter search exists anywhere in this repository and an untuned XGBoost
already beats the deployed forest on all five regression endpoints, so no claim of optimality is
available. It has a corresponding benefit, which is structural: there was no tuning set, so there can
be no tuning-set leakage, and selecting hyper-parameters on the same folds a number is reported from is
the commonest way a figure in this field is quietly inflated.
