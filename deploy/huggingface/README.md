---
title: BrainSafe AI
emoji: 🧠
colorFrom: blue
colorTo: yellow
sdk: docker
app_port: 7860
pinned: false
license: mit
short_description: Multi-endpoint prediction of small-molecule effects on the human brain
---

# BrainSafe AI

Predicts, from chemical structure alone, whether a small molecule is likely to act on the brain:
blood-brain-barrier penetration, engagement of 52 disease-relevant targets, free brain exposure,
cardiac safety, and per-condition relevance across 16 brain conditions, traced through a curated
target-to-pathway-to-disease knowledge graph.

Every endpoint is trained on measured public bioactivity data (ChEMBL, BindingDB, B3DB; 64,474
records over 61,317 unique compounds), probability-calibrated, and reported with an
applicability-domain flag and the nearest measured analogue.

**Research decision-support for prioritisation and hypothesis generation. It predicts molecular
target engagement and physicochemical properties, not clinical efficacy, and has not undergone
wet-lab or clinical validation. Not for medical, diagnostic or treatment decisions.**

## How to read the output

Three quantities are easy to confuse and are not on the same scale:

- **Binder probability** is calibrated per target and is *not* comparable between targets, because
  each carries its own threshold and training base rate.
- **Engagement signal** is the distance above a target's own threshold. This is the only engagement
  quantity comparable across targets, and it is what the disease layer consumes.
- **Independent mechanisms** groups engaged targets by homology. Dopamine D2 and D3 are engaged
  together by the same ligands, so two engaged homologues are close to one observation. Weigh this
  number, not the raw target count.

**A silent result is not evidence of inactivity.** Thresholds are set for precision, holding the
false-positive rate near 0.2% on random chemistry, and the measured cost is a median sensitivity of
0.77 across the panel, falling to 0.26 at the strictest endpoints.

## Validation

| Test | Result |
|---|---|
| Prospective sensitivity, scaffold-held-out | 0.791 (11,914 / 15,069) |
| Specificity, non-CNS compounds | 0.875 (lower bound) |
| Binder AUROC vs measured inactives | 0.955 across 43 targets |
| Deployed false-positive rate, random chemistry | median 0.0017 |

Beyond conventional validation the server was subjected to a systematic falsification analysis, in
which each central claim was paired with a null model able to reproduce the same apparent success by
accident. Three came back against the tool and were acted on: the curated knowledge-graph edge
weights carry no measurable predictive content beyond the graph topology, blood-brain-barrier gating
is an exposure filter that cannot discriminate between conditions, and engaged targets are not
independent observations. The same analysis identified a deployed endpoint that assigned binder
probabilities between 0.801 and 0.816 to glucose, urea and acetate; it was withdrawn.

Source, models, validation artefacts and the scripts that regenerate every figure:
https://github.com/krishna-g-999/brainsafe-ai
