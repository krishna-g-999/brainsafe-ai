# Inversion analysis: results

Each hypothesis below was stated so that it could fail, and paired with a null model capable of producing the same apparent success by accident. All scoring used the scaffold hold-out models where predictive power was at issue. The analysis itself is read-only: no trained model, no curated dataset and no scoring rule was changed in order to obtain any number reported here. Acting on the findings afterwards is a separate step and is recorded in the git history, so a wording change made in response to a result can never be mistaken for part of the evidence for it.

## Verdicts

| Hypothesis | Verdict | Evidence |
|---|---|---|
| H1 the disease score is informative | **SUPPORTED** | top-3 accuracy 0.786 vs permutation null 0.154 (p=0.005) and frequency null 0.548 |
| H2 the curated edge weights add value | **REFUTED** | curated 0.7865, uniform 0.7861, permuted 0.7844 |
| H3 BBB gating discriminates between diseases | **REFUTED (by construction)** | the gate multiplies every disease equally and cannot change their order |
| H4 specificity transfers to novel chemistry | **SUPPORTED** | false-positive rate 0.016 on 61 distant compounds against 0.051 measured on library chemistry |
| H5 read-across beats a frequency baseline | **SUPPORTED** | recall 0.973 against 0.059 |
| H6 the disease scores match real clinical indications | **WEAKENED** | top-3 accuracy 0.352 on 162 drugs never seen in training, against permutation null 0.145 (p=0.001) and frequency null 0.654 |
| H7 some panel targets are non-discriminative and explain the silent antiepileptics | **REFUTED** | none of 38 targets ranks below AUROC 0.70; the cause is the operating point, with median deployed sensitivity 0.79 and 7 targets under 0.50 |
| H8 engaged targets are independent observations | **REFUTED** | 37 targets fire across approved drugs but span only 17 independent directions; 5 homologous pairs correlate above 0.5 |

## What each result means

**H1.** Scored with hold-out models only, and restricted to compounds that are training actives of no other panel target, since scaffolds were withheld per target and a compound active at two targets would otherwise be memorised by one of the models scoring it. The disease layer reaches 78.6% top-3 accuracy where shuffling the target-to-disease map gives 15.4% and always answering with the three commonest diseases gives 54.8%. The layer carries real information.

**H2.** Curated weights score 0.7865, uniform weights 0.7861 and randomly permuted weights 0.7844. The spread is 0.0021. The information lies in which target connects to which disease, not in how strongly. The weights should be described as structure rather than as tuned parameters, and the graph would be simpler and no less accurate without them.

**H3.** Multiplying every disease score by the same BBB probability leaves their ranking untouched. Gating therefore cannot improve or damage which disease is chosen; it changes only the absolute value and hence what crosses the reporting threshold. It is an exposure filter, and the manuscript should call it one rather than implying it sharpens the disease call.

**H4.** Structures drawn by random PubChem identifier, independent of every set used to build this tool. On compounds distant from training chemistry the false-positive rate is 0.016, against 0.051 measured on library compounds.

**H5.** Read-across recovers the true target for 97.3% of held-out compounds against 5.9% for always answering with the commonest targets. The query and any identical structure were excluded. This measures read-across in its intended regime, where the query's target family is represented in the index; it does not show that read-across works for a target class the index does not contain, and the figure should not be quoted as if it did.

**H6.** H1 asked whether the disease layer recovers the disease that a compound's target maps to, using this project's own map. H6 asks the harder question: whether the score matches the condition the drug is actually approved to treat. Ground truth is ChEMBL's drug_indication table restricted to phase 4, mapped to the panel through a keyword list fixed before any prediction was made. On the 162 drugs whose exact structure appears nowhere in the training chemistry, top-3 accuracy is 35.2% against 14.5% for shuffling which drug carries which indication and 65.4% for always answering with the commonest indications. Two readings follow and both are true. The tool beats the permutation null decisively, so its output does depend on the compound and is not memorisation: drugs whose structure appears nowhere in training score the same as drugs that do. It does not beat the frequency null, because pain, depression and psychosis account for most approved CNS indications and a constant answer naming those three is right about 65% of the time. That constant answer carries no information about any particular compound, but it is a real bar and the tool does not clear it on this metric. Removing the reporting threshold and judging the ranking alone raises accuracy to 49.7%, which locates much of the gap in the decision to stay silent rather than in the ranking. Per indication the spread is wide: Psychosis / schizophrenia is recovered for 64% of its drugs, Epilepsy for 10%, which H7 explains.

**H7.** H6 left most approved antiepileptics scoring exactly zero, which looked like a broken model. It is not. Against 600 random PubChem structures, scored with models that never saw the compounds, every one of the 38 testable targets separates its own held-out actives from random chemistry at AUROC 0.91 or better. The cause is the operating point: thresholds hold the false-positive rate near 0.2 percent and the price is sensitivity, which ranges from 0.00 to 0.99 with a median of 0.79. 7 endpoints fire for under half their own actives (COX2, GABA_A, GluA2, KEAP1, P2X7, SIRT1, TAAR1). Old antiepileptics are small, simple and low-affinity, which is precisely the chemistry that sits under a strict cut. The user-facing consequence is that silence is not evidence of inactivity, and the tool should say so wherever it reports nothing.

**H8.** The interface reports how many targets a compound engages and draws an edge for each, which invites reading three engaged targets as three independent findings. For homologues that is wrong. Measured across approved drugs, SERT and NET r=0.73; OPRM1 and OPRK1 r=0.80; D2 and D3 r=0.72. Of the panel, 37 targets fire at least once but resolve into only 17 independent directions, so a raw count overstates the evidence by roughly a factor of two and a half. Co-firing is not itself an error, since a promiscuous ligand should engage both homologues; presenting it as corroboration is. Engaged targets are now grouped by homology family and the measured correlation is quoted wherever two members fire. Separately, a compound engaging nothing receives a reported finding 4.7% of the time on random chemistry, which bounds what a further expansion of the panel would cost in false leads.

## Consequences for the tool

1. The disease layer is validated for the first time, prospectively, and it works against the project's own target-to-disease map (H1) and, more weakly but far above chance, against external clinical indications (H6).
2. The curated edge weights are not doing measurable work and should be presented as graph structure rather than as tuned parameters.
3. BBB gating is an exposure filter, not a discriminator, and the manuscript wording should say so.
4. Read-across is validated only where the target family is already represented.
5. Agreement with clinical indications is real but does not beat a constant answer naming the three commonest CNS indications. The tool should be described as ranking mechanisms, not as predicting indications, and the per-indication table should be published alongside the aggregate so that the weak conditions are visible.
6. Silence is not evidence of inactivity. Deployed sensitivity has a median of 0.77 and falls to 0.26 at the strictest endpoints, so a null result must be reported with that number attached rather than as an absence of effect.
