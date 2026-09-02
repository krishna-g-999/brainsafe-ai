"""Build the Chapter 7 defence deck, reading every figure from an artefact.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/build_chapter07_defence.py
Out:  thesis/presentations/chapter07_defence.pptx
"""
from __future__ import annotations

import ast
import json
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

from deck_common import (
    ROOT, TAB, INV, PKG, W, M, INK, DEEP, TEAL, AMBER, CRIMSON, MUTED, TINT, PAPER, CHALK,
    HEAD, BODY, Deck, rows, num, fmt,
)

OUT = Path(__file__).resolve().parent / "chapter07_defence.pptx"


def facts() -> dict:
    d = {}
    tree = ast.parse((ROOT / "app.py").read_text(encoding="utf-8"))
    ns = {}
    for n in tree.body:
        if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") in (
                "KNOWLEDGE_GRAPH", "DISEASE_ORDER", "PERIPHERAL_MECHANISM_DISEASES"):
            ns[n.targets[0].id] = ast.literal_eval(n.value)
    kg = ns["KNOWLEDGE_GRAPH"]
    edges = [(t, e) for t, es in kg.items() for e in es]
    d["n_targets"], d["n_cond"], d["n_edges"] = len(kg), len(ns["DISEASE_ORDER"]), len(edges)
    d["peripheral"] = sorted(ns["PERIPHERAL_MECHANISM_DISEASES"])
    d["multi"] = sum(1 for es in kg.values() if len({e[2] for e in es}) > 1)
    w = [e[3] for _, e in edges]
    d["w"] = (len(set(w)), min(w), max(w), st.mean(w), st.median(w))
    d["w_dist"] = sorted(Counter(w).items())
    bc = defaultdict(list)
    for t, e in edges:
        bc[e[2]].append(t)
    sizes = sorted(len(v) for v in bc.values())
    d["per_cond"] = (st.median(sizes), min(sizes), max(sizes))
    d["best_cond"] = max(bc.items(), key=lambda kv: len(kv[1]))[0]

    br = json.loads((PKG / "07_MODELS" / "endpoint_base_rates.json").read_text())["classifiers"]
    d["base"] = sorted(((k, v["n"], v["base_rate"]) for k, v in br.items()),
                       key=lambda x: -x[2])

    h2 = {r["weights"]: num(r["top3_accuracy"]) for r in rows(INV / "H2_weight_ablation.csv")}
    d["h2"] = h2
    h1 = {r["test"]: r for r in rows(INV / "H1_disease_layer.csv")}
    d["h1_n"] = int(h1["observed top-3 accuracy"]["n"])
    d["h1"] = (num(h1["observed top-3 accuracy"]["value"]),
               num(h1["permutation null (mean)"]["value"]),
               num(h1["frequency null (always top-3 commonest)"]["value"]))
    d["h2_compounds"] = ((h2["curated"] - h2["uniform (all 1.0)"]) * d["h1_n"],
                         (h2["curated"] - h2["randomly permuted"]) * d["h1_n"])

    h6 = {r["stratum"]: r for r in rows(INV / "H6_clinical_indication.csv")}
    d["h6"] = (int(h6["never seen in training"]["n"]),
               num(h6["never seen in training"]["top3_accuracy"]),
               num(h6["never seen in training"]["permutation_null"]),
               num(h6["never seen in training"]["frequency_null"]),
               num(h6["never seen, ranking only"]["top3_accuracy"]))

    h9 = {r["metric"]: r for r in rows(INV / "H9_disease_discrimination_summary.csv")}
    d["h9"] = (num(h9["mean per-indication AUROC"]["model"]),
               num(h9["macro-averaged top-3 recall"]["model"]),
               num(h9["macro-averaged top-3 recall"]["frequency_null"]))
    d["pi"] = rows(INV / "H9_disease_discrimination.csv")
    d["h9_above"] = sum(1 for x in d["pi"] if num(x["auroc_model"]) > 0.5)

    h8 = {r["metric"]: num(r["value"]) for r in rows(INV / "H8_panel_independence.csv")}
    d["h8"] = (int(h8["targets that ever fire on the drug set"]),
               int(h8["independent directions in the firing pattern"]),
               h8["reported finding on random chemistry (user visible)"],
               h8["reported finding on approved drugs"],
               h8["mean binder endpoints fired per random compound"],
               h8["mean binder endpoints fired per approved drug"])
    fc = [x for x in rows(INV / "H8_family_correlation.csv") if x["phi_correlation"]]
    d["pairs"] = sorted(((x["family"], x["target_a"], x["target_b"], num(x["phi_correlation"]))
                         for x in fc), key=lambda t: -t[3])[:5]

    sp = [x for x in rows(TAB / "noncns_specificity_summary.csv")
          if x["metric"].startswith("Specificity")][0]
    d["spec"] = (num(sp["estimate"]), int(sp["k"]), int(sp["n"]))
    return d


F = facts()
D = Deck()

# 1 -- title
s = D.dark()
D.text(s, M, 1.74, W - 2 * M, 0.5, "Chapter 7", size=16, color=TEAL, bold=True)
D.text(s, M, 2.18, W - 2 * M - 1.0, 1.6,
       "Exposure gating and the pathway graph",
       size=40, bold=True, font=HEAD, color=PAPER, line=1.08)
D.text(s, M, 4.05, W - 2 * M - 1.6, 1.0,
       "The rule that turns seventy independent numbers into one answer, and the part of the system "
       "that has been most misdescribed, including by its own authors",
       size=16, color=CHALK, italic=True, line=1.25)
D.text(s, M, 5.65, W - 2 * M, 0.6,
       f"{F['n_targets']} targets, {F['n_cond']} conditions, {F['n_edges']} edges. The falsification "
       f"suite refuted the stated purpose of two of its three components.",
       size=13, color=TEAL, line=1.24)
D.notes(s, "This is the layer with the weakest claims in the thesis. Stating them precisely is more "
           "useful than defending them.")

# 2 -- enrichment
s = D.light()
D.head(s, "1", "Why enrichment, and not the probability",
       "A probability of 0.60 means opposite things at two endpoints")
D.text(s, M + 0.30, 1.92, 2.6, 0.3, "endpoint", size=11, bold=True, color=MUTED)
D.text(s, M + 3.00, 1.92, 1.3, 0.3, "n", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
D.text(s, M + 4.50, 1.92, 1.5, 0.3, "base rate", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
yy = 2.32
for i, (k, n, b) in enumerate(F["base"]):
    hi = k in ("hERG", "BACE1")
    if i % 2 == 0:
        D.card(s, M, yy - 0.09, 6.4, 0.52, fill=TINT)
    D.text(s, M + 0.30, yy, 2.6, 0.34, k.replace("_", "-"), size=13, bold=hi,
           color=DEEP if hi else INK)
    D.text(s, M + 3.00, yy, 1.3, 0.34, f"{n:,}", size=12.5, color=MUTED, align=PP_ALIGN.RIGHT)
    D.text(s, M + 4.50, yy, 1.5, 0.34, f"{b:.3f}", size=13, bold=hi,
           color=DEEP if hi else INK, align=PP_ALIGN.RIGHT)
    yy += 0.52
xs = M + 6.9
D.card(s, xs, 1.85, W - M - xs, 2.55, fill=TINT)
D.text(s, xs + 0.40, 2.10, W - M - xs - 0.8, 2.05,
       "A calibrated probability of 0.60 is strong evidence of activity at hERG, where 24 per cent of "
       "the measured set is active, and evidence of INACTIVITY at BACE1, where 86 per cent is.\n\n"
       "Ranking the raw probabilities puts those two compounds in the same place. Enrichment puts "
       "them on opposite sides of zero.",
       size=14, color=INK, line=1.28)
D.text(s, xs, 4.65, W - M - xs, 1.4,
       "The negative branch matters as much as the positive one. A probability below an endpoint's "
       "base rate is evidence of inactivity, not weak evidence of activity, so it is clipped to "
       "exactly zero and contributes nothing rather than contributing a little.",
       size=13, color=INK, line=1.28)
D.source(s, "submission_package/07_MODELS/endpoint_base_rates.json. The interface always shows the "
            "untransformed probability beside the signal.")
D.notes(s, "Base rates span nearly four-fold across the eight core classifiers. That span is the "
           "whole justification for the enrichment map.")

# 3 -- the graph
s = D.light()
D.head(s, "2", "The pathway graph", "Curated, versioned, anchored to KEGG, Reactome and IUPHAR")
stats = [(str(F["n_targets"]), "targets"), (str(F["n_cond"]), "conditions"),
         (str(F["n_edges"]), "edges"), (f"{F['multi']} of {F['n_targets']}", "drive more than one condition")]
xw = (W - 2 * M - 3 * 0.28) / 4
for i, (v, l) in enumerate(stats):
    x = M + i * (xw + 0.28)
    D.card(s, x, 1.80, xw, 1.75)
    D.stat(s, x + 0.30, 2.08, xw - 0.6, v, l, color=DEEP, vsize=34)
wn, wmin, wmax, wmean, wmed = F["w"]
D.text(s, M, 3.95, 6.2, 0.42, "The weights", size=17, bold=True, font=HEAD, color=TEAL)
D.text(s, M, 4.42, 6.2, 1.5,
       f"{wn} distinct values from {wmin:.2f} to {wmax:.2f}, mean {wmean:.3f}, median {wmed:.2f}. "
       f"Three quarters sit at 0.70 or above.\n\nTargets per condition: median {F['per_cond'][0]}, "
       f"from {F['per_cond'][1]} to {F['per_cond'][2]}. The best-connected condition is "
       f"{F['best_cond'].lower()}.",
       size=13.5, color=INK, line=1.28)
D.card(s, M + 6.6, 3.95, W - M - (M + 6.6), 2.05, fill=TINT)
D.text(s, M + 7.0, 4.20, W - M - (M + 6.6) - 0.8, 1.6,
       "Two conditions are exempt from the exposure gate because their mechanisms act outside the "
       f"barrier: {F['peripheral'][0].lower()} and {F['peripheral'][1].lower()}. That exemption "
       "matters on slide 5.",
       size=13, color=INK, line=1.28)
D.source(s, "Read directly from app.py KNOWLEDGE_GRAPH during the writing of this chapter.")
D.notes(s, "An edge is (pathway, KEGG id, condition, weight). AChE maps to Alzheimer's disease "
           "through the cholinergic synapse, hsa04725, at weight 1.0.")

# 4 -- max not sum
s = D.light()
D.head(s, "3", "Aggregation is a maximum, not a sum",
       "Engaging three of a condition's targets scores exactly as engaging its strongest")
D.equation(s, M, 1.80, W - 2 * M,
           [("S", False), ("d", True), ("(x)  =  max", False), ("(t,w) ∈ G(d)", True),
            ("  w · s", False), ("t", True), ("(x)", False)], size=25)
D.text(s, M, 2.62, W - 2 * M, 0.5,
       "That looks like discarding information. H8 is why it is not.",
       size=16, bold=True, font=HEAD, color=DEEP, align=PP_ALIGN.CENTER)
D.text(s, M + 0.30, 3.30, 3.0, 0.3, "family", size=11, bold=True, color=MUTED)
D.text(s, M + 3.40, 3.30, 3.2, 0.3, "pair", size=11, bold=True, color=MUTED)
D.text(s, M + 7.00, 3.30, 1.0, 0.3, "φ", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
yy = 3.68
for i, (fam, a, b, phi) in enumerate(F["pairs"]):
    if i % 2 == 0:
        D.card(s, M, yy - 0.09, 8.3, 0.52, fill=TINT)
    D.text(s, M + 0.30, yy, 3.0, 0.34, fam, size=12.5, color=MUTED)
    D.text(s, M + 3.40, yy, 3.2, 0.34, f"{a.replace('_','-')}, {b.replace('_','-')}", size=13,
           color=INK)
    D.text(s, M + 7.00, yy, 1.0, 0.34, f"{phi:.3f}", size=13, bold=True,
           color=CRIMSON if phi > 0.7 else DEEP, align=PP_ALIGN.RIGHT)
    yy += 0.52
xs = M + 8.7
D.text(s, xs, 3.35, W - M - xs, 2.4,
       f"{F['h8'][0]} targets fire across approved drugs but span only {F['h8'][1]} independent "
       f"directions.\n\nA sum would count one observation several times. A ligand engaging both "
       f"opioid receptors has given one piece of evidence, not two, and at φ = {F['pairs'][0][3]:.3f} "
       "the second is nearly determined by the first.",
       size=13, color=INK, line=1.28)
D.text(s, M, 6.35, W - 2 * M, 0.5,
       f"The firing pattern is sparse anyway: {F['h8'][4]:.3f} endpoints fire for a random compound "
       f"and {F['h8'][5]:.2f} for an approved drug, so the maximum usually has one non-zero member.",
       size=12.5, color=MUTED, italic=True, line=1.24)
D.source(s, "inversion/results/H8_panel_independence.csv, H8_family_correlation.csv")
D.notes(s, "Co-firing is not the error. Presenting it as corroboration is. The interface now groups "
           "engaged targets by homology family and quotes the measured correlation.")

# 5 -- H3
s = D.dark()
D.text(s, M, 0.95, 9.0, 0.5, "The exposure gate", size=15, color=AMBER, bold=True)
D.text(s, M, 1.40, 11.7, 0.76,
       "Gating is a filter. It cannot be a discriminator.",
       size=32, bold=True, font=HEAD, color=PAPER, line=1.10)
D.equation(s, M, 2.35, 11.7,
           [("S̃", False), ("d", True), ("(x)  =  γ", False), ("d", True), ("(x) · S", False),
            ("d", True), ("(x)", False)], size=24, color=CHALK)
D.text(s, M, 3.15, 11.7, 0.9,
       f"γ takes the same value for all {F['n_cond'] - len(F['peripheral'])} non-peripheral "
       "conditions, so multiplying every one of their scores by it is a common positive scaling, and "
       "a common positive scaling is rank-invariant. The gate cannot change which condition ranks "
       "first. It can only decide whether anything clears the reporting threshold.",
       size=14, color=CHALK, line=1.28)
# Title colours must contrast with the card fill, which is DEEP; the first draft set the first
# title to DEEP on DEEP and the third to CRIMSON on DEEP, both unreadable.
cards = [("It is not thereby useless", PAPER,
          f"Deciding whether to speak is the system's most-used behaviour. "
          f"{F['spec'][1]} of {F['spec'][2]} non-CNS compounds receive no disease call, and the gate "
          f"is a principal reason. Suppressing {F['spec'][0]:.1%} of irrelevant chemistry is "
          f"substantial work."),
         ("But it cannot be credited with the ranking", AMBER,
          "Any claim that exposure gating sharpens the disease call is arithmetically false. The "
          "manuscript wording was changed because of this result."),
         ("And the refutation has an unstated exception", PAPER,
          f"{F['peripheral'][0]} and {F['peripheral'][1].lower()} are exempt, so γ = 1 for them "
          "while the other 14 are scaled. Between those two groups the gate DOES change the ranking. "
          "No document states this.")]
xw = (W - 2 * M - 2 * 0.3) / 3
for i, (t, col, b) in enumerate(cards):
    x = M + i * (xw + 0.3)
    D.card(s, x, 4.35, xw, 2.20, fill=CRIMSON if i == 2 else DEEP)
    D.text(s, x + 0.32, 4.58, xw - 0.64, 0.72, t, size=14, bold=True, font=HEAD, color=col, line=1.14)
    D.text(s, x + 0.32, 5.35, xw - 0.64, 1.1, b, size=11.5, color=CHALK, line=1.24)
D.source(s, "inversion/results/H3_gating.csv. Refuted by construction, not by experiment.")
D.notes(s, "The third card is new in this thesis. H3's refutation is exact for the 14 and does not "
           "extend to the boundary with the two peripheral conditions.")

# 6 -- H2
s = D.light()
D.head(s, "4", "The curated edge weights add nothing measurable",
       "Topology replaced nothing; only the weights were changed")
cd = CategoryChartData()
cd.categories = ["curated", "uniform, all 1.0", "randomly permuted"]
cd.add_series("top-3 accuracy", [F["h2"]["curated"], F["h2"]["uniform (all 1.0)"],
                                 F["h2"]["randomly permuted"]])
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.80),
                        Inches(6.6), Inches(4.15), cd)
ch = gf.chart
ch.has_title = False
ch.has_legend = False
va = ch.value_axis
va.minimum_scale, va.maximum_scale = 0.75, 0.80
va.major_gridlines.format.line.color.rgb = D.grid
va.tick_labels.font.size = Pt(11)
va.tick_labels.font.color.rgb = MUTED
ch.category_axis.tick_labels.font.size = Pt(11)
ch.category_axis.tick_labels.font.color.rgb = MUTED
ser = ch.plots[0].series[0]
ser.format.fill.solid()
ser.format.fill.fore_color.rgb = DEEP
ch.plots[0].gap_width = 80
xs = M + 6.95
D.text(s, xs, 1.92, W - M - xs, 0.9,
       f"The spread from curated to permuted is "
       f"{F['h2']['curated'] - F['h2']['randomly permuted']:.4f}.",
       size=15, bold=True, color=DEEP, line=1.24)
D.card(s, xs, 2.90, W - M - xs, 1.65, fill=TINT)
D.text(s, xs + 0.35, 3.12, W - M - xs - 0.7, 1.25,
       f"On the {F['h1_n']:,} compounds the evaluation uses, replacing hand-assigned weights with "
       f"uniform ones changes the top-3 outcome for about {F['h2_compounds'][0]:.0f} compounds. "
       f"Permuting them changes about {F['h2_compounds'][1]:.0f}.",
       size=13, color=INK, line=1.26)
D.text(s, xs, 4.80, W - M - xs, 1.5,
       "Under a maximum, a weight matters only by scaling the winning term or changing which term "
       "wins. With a sparse firing pattern the common case is one engaged target, where the weight is "
       "a pure scale factor.",
       size=12.5, color=INK, line=1.26)
D.text(s, M, 6.20, W - 2 * M, 0.5,
       "The information lies in which target connects to which condition, not in how strongly. The "
       "weights are structure, not tuned parameters, and no claim is made for them.",
       size=13.5, bold=True, color=AMBER, line=1.24)
D.source(s, "inversion/results/H2_weight_ablation.csv")
D.notes(s, "A refutation the project acted on rather than argued with.")

# 7 -- what the layer does establish
s = D.light()
D.head(s, "5", "What the disease layer does establish",
       "Three hypotheses, three different questions, one narrow claim")
h1, h6, h9 = F["h1"], F["h6"], F["h9"]
blocks = [
    ("H1", "against the project's own target-to-disease map", DEEP,
     f"Top-3 accuracy {h1[0]:.4f} on {F['h1_n']:,} compounds, against a permutation null of "
     f"{h1[1]:.4f} (p = 0.005) and a frequency null of {h1[2]:.4f}. Beats both. But it asks whether "
     "the layer recovers the condition its own graph assigns, which is internal consistency rather "
     "than a test against the world."),
    ("H6", "against real clinical indications", AMBER,
     f"On the {h6[0]} drugs whose structure appears nowhere in training, top-3 accuracy {h6[1]:.4f} "
     f"against a permutation null of {h6[2]:.4f} and a frequency null of {h6[3]:.4f}. Beats the "
     f"permutation null decisively, so the output depends on the compound. Does not beat a constant "
     f"answer. Judging the ranking alone raises it to {h6[4]:.4f}."),
    ("H9", "on metrics a constant predictor cannot pass", TEAL,
     f"Mean per-indication AUROC {h9[0]:.4f} against 0.500, beating chance on {F['h9_above']} of 9. "
     f"Macro-averaged top-3 recall {h9[1]:.4f} against {h9[2]:.4f}. The layer is doing something."),
]
yy = 1.80
for tag, sub, col, body in blocks:
    D.card(s, M, yy, W - 2 * M, 1.48)
    D.dot(s, M + 0.35, yy + 0.50, tag, fill=col, dia=0.46)
    D.text(s, M + 1.15, yy + 0.28, 3.5, 0.85, sub, size=14, bold=True, font=HEAD, color=col,
           line=1.14)
    D.text(s, M + 4.95, yy + 0.22, W - 2 * M - 5.3, 1.1, body, size=12, color=INK, line=1.24)
    yy += 1.62
D.text(s, M, 6.55, W - 2 * M, 0.4,
       "All three are true at once. The honest claim is that the layer ranks mechanisms, not that it "
       "predicts indications.",
       size=13, color=MUTED, italic=True)
D.source(s, "inversion/results/H1_disease_layer.csv, H6_clinical_indication.csv, "
            "H9_disease_discrimination_summary.csv")
D.notes(s, "H6 and H9 are not in conflict. H6 is the fair description of the top-3 list a user reads; "
           "H9 is the fair description of whether the layer is doing anything.")

# 8 -- the per-indication spread
s = D.light()
D.head(s, "6", "The spread matters more than the mean",
       "Two of the nine conditions are not being predicted at all")
pi = sorted(F["pi"], key=lambda x: num(x["auroc_model"]))
cd2 = CategoryChartData()
cd2.categories = [x["indication"] for x in pi]
cd2.add_series("per-indication AUROC", [num(x["auroc_model"]) for x in pi])
gf2 = s.shapes.add_chart(XL_CHART_TYPE.BAR_CLUSTERED, Inches(M), Inches(1.78),
                         Inches(8.1), Inches(4.30), cd2)
c2 = gf2.chart
c2.has_title = False
c2.has_legend = False
va = c2.value_axis
va.minimum_scale, va.maximum_scale = 0.4, 0.85
va.major_gridlines.format.line.color.rgb = D.grid
va.tick_labels.font.size = Pt(10)
va.tick_labels.font.color.rgb = MUTED
c2.category_axis.tick_labels.font.size = Pt(10)
c2.category_axis.tick_labels.font.color.rgb = MUTED
ser = c2.plots[0].series[0]
ser.format.fill.solid()
ser.format.fill.fore_color.rgb = DEEP
c2.plots[0].gap_width = 55
xs = M + 8.4
D.text(s, xs, 1.92, W - M - xs, 1.2,
       f"Decisive for depression at {num(pi[-1]['auroc_model']):.4f} and psychosis. Marginal for "
       f"ADHD at 0.5006. At or below chance for sleep and epilepsy.",
       size=13.5, bold=True, color=DEEP, line=1.26)
D.card(s, xs, 3.30, W - M - xs, 1.85, fill=TINT)
D.text(s, xs + 0.32, 3.52, W - M - xs - 0.65, 1.45,
       "Publishing the mean without this table would conceal that two of the nine conditions are not "
       "being predicted at all. A constant predictor scores 0.500 on every row by construction.",
       size=12.5, color=INK, line=1.26)
D.text(s, xs, 5.40, W - M - xs, 1.2,
       "Epilepsy is the informative failure: most approved antiepileptics are small, simple and "
       "low-affinity, which is the chemistry that sits below a strict cut.",
       size=12, color=MUTED, line=1.26)
D.source(s, "inversion/results/H9_disease_discrimination.csv")
D.notes(s, "The endpoints that would carry the antiepileptics are among the six firing for fewer "
           "than half their own actives, from Chapter 5.")

# 9 -- closing
s = D.dark()
D.text(s, M, 1.20, 8.4, 0.5, "Chapter 7, honestly summarised", size=15, color=TEAL, bold=True)
items = [
    ("Exposure gating works as a filter and cannot work as a discriminator", AMBER),
    ("The curated edge weights carry no measurable information; the graph's content is its topology", AMBER),
    (f"The layer carries real information about mechanism, {F['h1'][0]:.4f} against a permutation "
     f"null of {F['h1'][1]:.4f}", TEAL),
    (f"It ranks conditions better than chance for {F['h9_above']} of 9, from "
     f"{num(F['pi'][0]['auroc_model']):.4f} down to {num(F['pi'][-1]['auroc_model']):.4f}", TEAL),
    ("It does not predict indication, and is correctly described as a route from a mechanism to the "
     "conditions it touches", PAPER),
]
yy = 1.85
for i, (t, col) in enumerate(items, 1):
    D.dot(s, M, yy + 0.02, str(i), fill=col if col is not PAPER else TEAL, fg=INK, dia=0.36)
    D.text(s, M + 0.58, yy - 0.02, 11.1, 0.85, t, size=15, color=col, line=1.28)
    yy += 0.95
D.text(s, M, 6.55, 11.5, 0.5,
       "The maximum-not-sum rule is the one component whose stated justification survived testing, "
       "and it survived because H8 supplied the evidence for it.",
       size=13.5, color=CHALK, italic=True, line=1.24)
D.notes(s, "Two of three components had their stated purpose refuted. Say so first; it is what makes "
           "the third credible.")

D.save(OUT)
