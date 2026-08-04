# Inversion analysis: results

Each hypothesis below was stated so that it could fail, and paired with a null model capable of producing the same apparent success by accident. All scoring used the scaffold hold-out models where predictive power was at issue. Nothing in `models_rf/`, `data/` or `app.py` was modified by this analysis.

## Verdicts

| Hypothesis | Verdict | Evidence |
|---|---|---|
| H1 the disease score is informative | **SUPPORTED** | top-3 accuracy 0.792 vs permutation null 0.205 (p=0.005) and frequency null 0.558 |
| H2 the curated edge weights add value | **REFUTED** | curated 0.7917, uniform 0.7911, permuted 0.7899 |
| H3 BBB gating discriminates between diseases | **REFUTED (by construction)** | the gate multiplies every disease equally and cannot change their order |
| H4 specificity transfers to novel chemistry | **SUPPORTED** | false-positive rate 0.051 on 59 distant compounds against 0.125 measured on library chemistry |
| H5 read-across beats a frequency baseline | **SUPPORTED** | recall 0.970 against 0.060 |

## What each result means

**H1.** Scored with hold-out models only, so no compound was seen in training. The disease layer reaches 79.2% top-3 accuracy where shuffling the target-to-disease map gives 20.5% and always answering with the three commonest diseases gives 55.8%. The layer carries real information.

**H2.** Curated weights score 0.7917, uniform weights 0.7911 and randomly permuted weights 0.7899. The spread is 0.0018. The information lies in which target connects to which disease, not in how strongly. The weights should be described as structure rather than as tuned parameters, and the graph would be simpler and no less accurate without them.

**H3.** Multiplying every disease score by the same BBB probability leaves their ranking untouched. Gating therefore cannot improve or damage which disease is chosen; it changes only the absolute value and hence what crosses the reporting threshold. It is an exposure filter, and the manuscript should call it one rather than implying it sharpens the disease call.

**H4.** Structures drawn by random PubChem identifier, independent of every set used to build this tool. On compounds distant from training chemistry the false-positive rate is 0.051, against 0.125 measured on library compounds.

**H5.** Read-across recovers the true target for 97.0% of held-out compounds against 6.0% for always answering with the commonest targets. The query and any identical structure were excluded. This measures read-across in its intended regime, where the query's target family is represented in the index; it does not show that read-across works for a target class the index does not contain, and the figure should not be quoted as if it did.

## Consequences for the tool

1. The disease layer is validated for the first time, prospectively, and it works.
2. The curated edge weights are not doing measurable work and should be presented as graph structure rather than as tuned parameters.
3. BBB gating is an exposure filter, not a discriminator, and the manuscript wording should say so.
4. Read-across is validated only where the target family is already represented.
