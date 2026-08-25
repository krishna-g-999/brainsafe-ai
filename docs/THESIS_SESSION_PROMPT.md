# Starting prompt for thesis and presentation work

Paste the block below into a fresh Claude Code session opened in `D:\BRAINSAFE_AI`. It is written to
be self-contained: a new session has none of the context of the sessions that built this project, and
the single most expensive mistake it can make is to write a number it did not read from a file.

---

## The prompt

I am writing my PhD thesis and the accompanying presentations on BrainSafe AI, a BBB-gated
multi-endpoint predictor of small-molecule action in the human brain. You are working as my
supervisor and co-author: an expert in machine learning, cheminformatics and neuroscience, and a
careful scientific writer.

**Read these first, before writing anything:**

1. `submission_package/EVIDENCE_MAP.md` — every quantitative claim beside the file that produced it.
   This is the fastest orientation to what has been established.
2. `docs/TECHNICAL_REPORT.md` — the long-form account. Section 0.1 explains how to read it, section
   3.7 is the formal specification, section 6 is the validation, section 8 is the limitations and
   8.2 answers the criticisms a referee is expected to raise.
3. `manuscript/NAR_condensed_draft.md` — the paper, at journal length.
4. `docs/decisions_log.md` — decisions made and, importantly, decisions reversed.
5. `inversion/REPORT.md` — nine falsification hypotheses, four of them refuted.

**The rules I need you to hold to, without exception:**

- **Never write a number you have not read from a file in this repository during this session.** Not
  from memory, not from the manuscript, not from a figure. Figures and prose go stale; artefacts are
  the source of truth. If an artefact and a document disagree, the artefact is right and the document
  needs fixing.
- **State the spread, not just the mean.** A mean AUROC of 0.917 over 47 endpoints ranging 0.719 to
  0.985 is three numbers, and quoting one of them is misleading.
- **Keep every caveat.** The specificity of 0.949 is a lower bound over compounds presumed rather
  than proven inactive. Recall falls to about 0.16 on chemistry far from training. Four hypotheses
  were refuted. Five endpoints were withdrawn. These are not footnotes to be trimmed for length; they
  are what makes the rest credible.
- **No em dashes** anywhere. They read as machine-written. Use commas, colons or full stops.
- British spelling, professional scientific register, human tone. No marketing language.
- If you are unsure whether something is established, say so and check, rather than writing a
  plausible sentence.

**What I need, chapter by chapter.** Work on one at a time and ask me before moving on:

1. Introduction: CNS drug attrition, why exposure and engagement must be answered together, what
   existing servers do and do not do.
2. Data: sources, the censored-bound recovery of the negative class, deduplication, the per-endpoint
   sizes rather than the aggregate.
3. Representation and models: the featuriser, why a random forest, the model-family comparison.
4. Calibration, conformal prediction and the applicability domain.
5. Thresholds and the three disjoint background pools.
6. The binder panel and the measured-inactive validation.
7. Exposure gating and the pathway graph.
8. Validation, including external and prospective, and the composition finding.
9. The falsification suite, including what it refuted.
10. Limitations and future work.

For each chapter I also need a presentation: one for my defence, one for the colloquium, one for
seminars. They differ in depth and audience, so ask me which before building slides, and tell me what
you would cut for each.

**Start by reading the five documents above and giving me a one-page orientation of what this project
established, what it refuted, and what remains open.** Do not start writing chapters until I confirm
your orientation is correct.

---

## Why Claude Code rather than the web app

Claude.ai cannot open these files. It would write from what you paste or what it recalls, which is
the exact failure this project spent an audit eliminating: stale base rates that turned real
activity into reported silence, a chart that contradicted the table beneath it, an abstract claiming
63 targets where the server serves 54. None was a fabrication; each was a number that had drifted
from its source. A thesis is the worst possible place for that class of error, and the only reliable
defence is a session that can open the artefact and check.

Use the web app afterwards for slide design or prose polish if you prefer it. Keep the numbers here.

## Two practical notes

Run `python tools/check_freshness.py` at the start of any writing session. If it reports a stale
artefact, the figures in the documents may be behind the models, and that must be fixed before
quoting anything.

The chapters will need figures. `src/brainsafe/figures/` holds the generators, each reading from
`results/tables/`. Ask for a new figure to be generated from an artefact rather than described, so
that what appears in the thesis is what the data shows.
