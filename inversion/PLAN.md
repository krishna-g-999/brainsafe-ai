# Inversion analysis: an attempt to falsify BrainSafe AI

## Premise

Every result so far has been produced by people who wanted the tool to work, including me. The
validation to date establishes that the target models rank compounds correctly. It does not
establish that the layer a user actually reads, the per-disease brain-relevance score, means
anything at all. That layer rests on a curated knowledge graph and a gating rule, neither of which
was learned from data and neither of which has ever been tested against an outcome.

This analysis therefore assumes each headline claim is **false** and asks what evidence would prove
it. A test that cannot fail is not evidence, so every hypothesis below is paired with a null model
that would produce the same apparent success by accident.

## Constraints

- **Read-only.** Nothing in `models_rf/`, `data/` or `app.py` is modified. Trained models and
  curated data are inputs to this analysis, never outputs of it.
- **All results are written under `inversion/`**, separate from `results/`, so a falsification can
  never be mistaken for a validation.
- **Prospective by construction.** Where a hypothesis concerns predictive power, it is tested with
  the scaffold-hold-out models in `models_rf/holdout/`, which never saw the compounds they score.

## Hypotheses under test

### H1. The disease score is informative
*Null: the disease layer carries no information. Apparent success comes from the frequency of common
diseases in the panel, or from target memorisation.*

Using hold-out models only, take each held-out active of target T and ask whether the disease that T
maps to appears among the top predicted diseases. Compare against two nulls: a permutation null in
which the target-to-disease map is shuffled, and a frequency null that always predicts the most
common disease. If the tool does not clearly beat both, the disease layer is decoration.

### H2. The curated knowledge-graph weights add value
*Null: the hand-assigned weights are arbitrary and uniform weights do as well.*

Recompute the same predictions with curated weights, with all weights set to one, and with weights
randomly permuted across edges. If curated and uniform are indistinguishable, the curation should be
described as a structure rather than as a tuned model. If curated and permuted are
indistinguishable, the weights are noise.

### H3. Multiplying by BBB penetration is the right operation
*Null: gating by BBB neither helps nor beats the alternatives.*

Compare ranking quality under the raw mechanism score, the score multiplied by BBB, and the score
multiplied by predicted unbound brain exposure. The current design is only justified if it wins.

### H4. Specificity survives outside the training neighbourhood
*Null: the measured specificity of 0.875 was an artefact of drawing negatives from the same library
the models were trained on, where every compound is in-domain.*

Sample compounds with low maximum Tanimoto to training chemistry and measure the false-positive rate
there. If it rises sharply, the headline specificity is not transferable and must be restated.

### H5. Read-across adds information beyond the models
*Null: the targets of a compound's nearest neighbours are no more informative than targets drawn at
random.*

For held-out compounds with a known target, ask how often read-across recovers that target, against
a random-target baseline matched to target frequency.

## Reporting rule

Each test writes its own CSV and a one-line verdict of SUPPORTED, WEAKENED or REFUTED. A refuted
hypothesis is the most valuable outcome available here and is reported as prominently as a
confirmation.
