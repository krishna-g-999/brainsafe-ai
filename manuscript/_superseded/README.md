# Superseded manuscript material

Nothing here is current. It is kept rather than deleted because a project that claims its results are
reproducible should be able to show what it previously believed, and because a deleted file cannot be
compared against the one that replaced it. Nothing in this directory is referenced by any live
document, and nothing here ships in the submission package.

## Manuscripts

| File | Date | Why it is here |
|---|---|---|
| `BS_MANUSCRIPT_FINAL.md` | 17 Aug 2026 | An earlier full manuscript, superseded by `NAR_WebServer_BrainSafe_draft.md`. Predates the 20 August retrain, the GluA2 withdrawal from the pathway graph, and the base-rate correction of 25 August. |
| `BrainSafe_AI_Manuscript.docx` | 23 Jul 2026 | Describes the earlier ensemble engine, before the binder panel existed. |
| `NAR_WebServer_BrainSafe_draft.docx` | 31 Jul 2026 | A build of the draft from before the retrain; the current build is `NAR_WebServer_BrainSafe.docx`. |

## Figures, July naming scheme

`figures_july_scheme/` holds eleven figures from a naming scheme that has since been replaced. They
are archived for two reasons rather than one.

The first is staleness: every one predates the 20 August retrain, so each shows numbers the panel no
longer produces.

The second matters more. Their names **collide** with the current scheme. There were two files called
`Figure1_*`, two called `Figure2_*`, and so on through `Figure11_*`. A document referring to
"Figure 3" had two candidates in one directory, one current and one two months old, distinguished
only by a suffix. That is how a stale image ends up in a submission, and it is why these were moved
rather than left beside their replacements.

| Archived (July scheme) | Current equivalent, if any |
|---|---|
| `Figure1_endpoint_rationale` | `Figure10_endpoint_selection` |
| `Figure2_mechanism` | `Figure1_architecture` |
| `Figure3_model_selection` | none; the comparison is now a table in the technical report |
| `Figure4_validation` | `Figure6_validation` |
| `Figure5_cv_design_and_errorbars` | `Figure3_cv_design` |
| `Figure6_all_endpoints` | `Figure9_model_atlas` |
| `Figure7_temporal_by_domain` | `Figure11_external_validation` |
| `Figure8_scaffold_holdout` | folded into `Figure6_validation` |
| `Figure9_specificity` | folded into `Figure6_validation` |
| `Figure10_performance` | `Figure9_model_atlas` |
| `Figure11_decision_analysis` | none |

## The live set

Eleven figures remain in `manuscript/figures/`, all regenerated on or after 20 August and all now
declared in `tools/check_freshness.py`, which fails if any falls behind the artefacts it draws from.
Five of them were outside that graph until this archive was made, and one of those five,
`Figure8_use_case`, had been drawn with the training base rates that were corrected on 25 August. It
has been regenerated.
