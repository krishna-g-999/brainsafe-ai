"""Build the Chapter 5 defence deck, reading every figure from an artefact.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/build_chapter05_defence.py
Out:  thesis/presentations/chapter05_defence.pptx
"""
from __future__ import annotations

import json
import pickle
import statistics as st
import sys
from collections import Counter
from pathlib import Path

from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

from deck_common import (
    ROOT, TAB, PKG, W, M, INK, DEEP, TEAL, AMBER, CRIMSON, MUTED, TINT, PAPER, CHALK,
    HEAD, BODY, Deck, rows, num, fmt,
)

sys.path.insert(0, str(ROOT / "src"))
from brainsafe.models.pools import role_of, SHARES, AD_REFERENCE  # noqa: E402

OUT = Path(__file__).resolve().parent / "chapter05_defence.pptx"


def facts() -> dict:
    d = {}
    smi = pickle.load(open(AD_REFERENCE, "rb"))[0]
    c = Counter(role_of(s) for s in smi)
    d["lib"] = sum(c.values())
    d["pools"] = [(k, c[k], 100 * c[k] / d["lib"], SHARES[k]) for k in SHARES]

    reg = json.loads((PKG / "07_MODELS" / "binder_panel_registry.json").read_text())
    dep = {k: v for k, v in reg.items() if v.get("deployed")}
    d["n_dep"] = len(dep)
    ins = [v["fpr_in_sample_on_threshold_set"] for v in dep.values()
           if v.get("fpr_in_sample_on_threshold_set") is not None]
    ev = [v["background_fpr_held_out"] for v in dep.values()
          if v.get("background_fpr_held_out") is not None]
    d["ins"] = (st.mean(ins), st.median(ins), min(ins), max(ins), len(ins),
                sum(1 for v in ins if abs(v - 0.10) < 0.001))
    d["ev"] = (st.mean(ev), st.median(ev), min(ev), max(ev), len(ev),
               sum(1 for v in ev if v > 0.05), sum(1 for v in ev if v > 0.10))
    pair = [(k, v["fpr_in_sample_on_threshold_set"], v["background_fpr_held_out"])
            for k, v in dep.items()
            if v.get("fpr_in_sample_on_threshold_set") is not None
            and v.get("background_fpr_held_out") is not None]
    d["lower_on"] = sum(1 for _, a, b in pair if b < a)
    d["n_pair"] = len(pair)
    d["gap_fpr"] = st.mean(b - a for _, a, b in pair)
    d["hi_bg"] = sorted(pair, key=lambda x: -x[2])[:4]
    d["lo_bg"] = sorted(pair, key=lambda x: x[2])[:3]

    ft = [r for r in rows(TAB / "final_thresholds.csv") if r["target"] in dep]
    d["bind_mi"] = sum(1 for r in ft if r["binding_constraint"] == "measured inactives")
    d["bind_bg"] = sum(1 for r in ft if r["binding_constraint"] == "background")
    d["med_mi"] = st.median(num(r["from_measured_inactives"]) for r in ft)
    d["med_bg"] = st.median(num(r["from_background"]) for r in ft)
    d["unreliable"] = sorted(r["target"] for r in ft if r["reliable"] != "True")

    rec = rows(TAB / "sensitivity_reconciliation.csv")
    h = [num(r["sensitivity_heldout"]) for r in rec if r["sensitivity_heldout"]]
    p = [num(r["sensitivity_published"]) for r in rec if r["sensitivity_published"]]
    d["held"] = (st.mean(h), st.median(h), min(h), max(h), sum(1 for v in h if v < 0.5))
    d["pub"] = (st.mean(p), st.median(p), min(p), max(p), sum(1 for v in p if v < 0.5))
    d["n_sens"] = len(h)
    hm = {r["target"]: num(r["sensitivity_heldout"]) for r in rec if r["sensitivity_heldout"]}
    pm = {r["target"]: num(r["sensitivity_published"]) for r in rec if r["sensitivity_published"]}
    d["held_min_ep"] = min(hm, key=hm.get)
    d["pub_min_ep"] = min(pm, key=pm.get)
    dv = [num(r["published_minus_heldout"]) for r in rec if r["published_minus_heldout"]]
    d["higher_on"] = sum(1 for v in dv if v > 0)
    d["mean_gap"] = st.mean(dv)
    d["worst"] = sorted(((r["target"], num(r["sensitivity_heldout"]),
                          num(r["sensitivity_published"]), num(r["published_minus_heldout"]))
                         for r in rec if r["published_minus_heldout"]),
                        key=lambda x: -x[3])[:5]
    tf = [num(r["train_fraction_of_published_set"]) for r in rec
          if r["train_fraction_of_published_set"] and num(r["train_fraction_of_published_set"]) <= 1]
    d["train_frac"] = (st.median(tf), len(tf))
    d["low6"] = sorted(((k, v) for k, v in hm.items() if v < 0.5), key=lambda x: x[1])

    h7 = {r["target"] for r in rows(ROOT / "inversion" / "results" /
                                    "H7_target_discrimination.csv")
          if num(r["deployed_sensitivity"]) < 0.5}
    d["h7_low"] = sorted(h7)
    d["common"] = sorted(h7 & {k for k, _ in d["low6"]})
    return d


F = facts()
D = Deck()

# 1 -- title
s = D.dark()
D.text(s, M, 1.74, W - 2 * M, 0.5, "Chapter 5", size=16, color=TEAL, bold=True)
D.text(s, M, 2.18, W - 2 * M - 1.0, 1.6,
       "Thresholds, and the three disjoint background pools",
       size=38, bold=True, font=HEAD, color=PAPER, line=1.08)
D.text(s, M, 4.10, W - 2 * M - 1.6, 1.0,
       "Where a panel of independent models becomes a system that says yes or no, and where the "
       "easiest self-deception in the project lives",
       size=16, color=CHALK, italic=True, line=1.25)
D.text(s, M, 5.65, W - 2 * M, 0.6,
       f"One library of {F['lib']:,} compounds, partitioned once so that no compound can hold two "
       f"roles",
       size=13, color=TEAL)
D.notes(s, "This chapter contains a correction to a published number. Do not save it for the end of "
           "the viva; it is better volunteered.")

# 2 -- the circularity
s = D.light()
D.head(s, "1", "The circularity a single pool creates", None)
D.card(s, M, 1.80, W - 2 * M, 1.45, fill=TINT)
D.text(s, M + 0.45, 2.05, W - 2 * M - 0.9, 1.0,
       "Choose a threshold as the 90th percentile of the scores a model gives to a sample of "
       "background chemistry. Now measure the false-positive rate of that threshold on the same "
       "sample. It is ten per cent. It was always going to be ten per cent.",
       size=16, color=INK, line=1.30)
D.text(s, M, 3.50, W - 2 * M, 0.5,
       "The measurement restates the quantile and tests nothing about the model.",
       size=18, bold=True, font=HEAD, color=CRIMSON)
D.text(s, M, 4.20, W - 2 * M, 0.5,
       "It is not hypothetical carelessness. It is the natural thing to do when one library is "
       "available and three jobs need doing:",
       size=14, color=INK)
jobs = [("decoys", "trained on as presumed negatives where measured inactives are scarce"),
        ("a threshold sample", "whose score distribution sets the cut"),
        ("an evaluation sample", "on which the achieved false-positive rate is measured")]
yy = 4.90
for i, (a, b) in enumerate(jobs, 1):
    D.dot(s, M + 0.10, yy, str(i), fill=DEEP, dia=0.32)
    D.text(s, M + 0.60, yy - 0.02, 2.7, 0.34, a, size=14, bold=True, color=DEEP)
    D.text(s, M + 3.45, yy - 0.02, W - 2 * M - 3.55, 0.34, b, size=13, color=INK)
    yy += 0.48
D.text(s, M, 6.42, W - 2 * M, 0.4,
       "Sharing a pool between the second and third restates the target. Sharing with the first is "
       "worse: the model was trained to score those compounds as negative, and is then congratulated "
       "for doing so.",
       size=12.5, color=MUTED, italic=True, line=1.24)
D.notes(s, "Say the last line slowly. It is the argument for the whole partition.")

# 3 -- the partition
s = D.light()
D.head(s, "2", "One partition, three roles",
       f"Recomputed from src/brainsafe/models/pools.py during the writing of this chapter")
cd = CategoryChartData()
cd.categories = [k for k, _, _, _ in F["pools"]]
cd.add_series("compounds", [n for _, n, _, _ in F["pools"]])
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.80),
                        Inches(6.3), Inches(3.15), cd)
ch = gf.chart
ch.has_title = False
ch.has_legend = False
ch.value_axis.major_gridlines.format.line.color.rgb = D.grid
ch.value_axis.tick_labels.font.size = Pt(10)
ch.value_axis.tick_labels.font.color.rgb = MUTED
ch.category_axis.tick_labels.font.size = Pt(11)
ch.category_axis.tick_labels.font.color.rgb = MUTED
ser = ch.plots[0].series[0]
ser.format.fill.solid()
ser.format.fill.fore_color.rgb = DEEP
ch.plots[0].gap_width = 70

xs = M + 6.6
D.text(s, xs + 0.10, 1.88, 2.4, 0.3, "pool", size=11, bold=True, color=MUTED)
D.text(s, xs + 2.55, 1.88, 1.5, 0.3, "compounds", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
D.text(s, xs + 4.20, 1.88, 1.1, 0.3, "share", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
yy = 2.28
for k, n, pct, tgt in F["pools"]:
    D.card(s, xs, yy - 0.09, W - M - xs, 0.54, fill=TINT)
    D.text(s, xs + 0.10, yy, 2.4, 0.34, k, size=13.5, bold=True, color=INK)
    D.text(s, xs + 2.55, yy, 1.5, 0.34, f"{n:,}", size=13.5, color=DEEP, bold=True, align=PP_ALIGN.RIGHT)
    D.text(s, xs + 4.20, yy, 1.1, 0.34, f"{pct:.2f}%", size=13, color=MUTED, align=PP_ALIGN.RIGHT)
    yy += 0.58
D.text(s, xs, 4.10, W - M - xs, 0.85,
       f"{F['lib']:,} compounds, assigned once. Each threshold and each rate uses a sample of 3,000 "
       "from the relevant pool.",
       size=12.5, color=INK, line=1.24)

D.card(s, M, 5.20, W - 2 * M, 1.25, fill=TINT)
D.text(s, M + 0.45, 5.45, W - 2 * M - 0.9, 0.85,
       "The consequence is one sentence: the measured rate can now disagree with its target, and "
       "that it can disagree is the evidence that it is a measurement rather than a restatement.",
       size=15, bold=True, color=DEEP, line=1.28)
D.source(s, "Pool membership is blake2b(salt + smiles) mod 100, banded 60/20/20. Five subtests pin it.")
D.notes(s, "A hash rather than a shuffle: assignment does not depend on library order, and adding "
           "compounds later leaves every existing assignment untouched, so a threshold set today "
           "stays comparable with a rate measured next year.")

# 4 -- two constraints
s = D.light()
D.head(s, "3", "Two constraints, and which one binds",
       "The operating point is the stricter of them, clipped to [0.05, 0.999]")
D.card(s, M, 1.80, 5.85, 2.15)
D.text(s, M + 0.42, 2.05, 5.0, 0.4, "On the target's own measured inactives",
       size=16, bold=True, font=HEAD, color=DEEP)
D.text(s, M + 0.42, 2.55, 5.0, 1.25,
       "No more than 10 per cent may be called a binder. The cut is the quantile of the scores given "
       "to the half of that endpoint's measured inactives withheld from fitting.",
       size=13, color=INK, line=1.26)
D.card(s, M + 6.25, 1.80, W - 2 * M - 6.25, 2.15)
D.text(s, M + 6.67, 2.05, 5.0, 0.4, "On random library chemistry",
       size=16, bold=True, font=HEAD, color=TEAL)
D.text(s, M + 6.67, 2.55, 5.0, 1.25,
       "No more than 5 per cent may be called a binder. Nav1.1 scored AUROC 0.979 against its own "
       "narrow negative set and gave glucose a binder probability of 0.806. This constraint exists "
       "because of that.",
       size=13, color=INK, line=1.26)
D.text(s, M, 4.20, W - 2 * M, 0.45,
       f"The measured-inactive constraint binds on {F['bind_mi']} of {F['n_dep']}; background on "
       f"{F['bind_bg']}.",
       size=17, bold=True, font=HEAD, color=INK)
D.text(s, M, 4.78, W - 2 * M, 0.85,
       f"Median cut from measured inactives {fmt(F['med_mi'], 4)}, from background "
       f"{fmt(F['med_bg'], 4)}. That ordering is Chapter 2's finding from the other end: compounds a "
       "chemist thought worth testing are much harder to separate from real ligands than compounds "
       "drawn at random.",
       size=13.5, color=INK, line=1.28)
D.card(s, M, 5.80, W - 2 * M, 0.85, fill=TINT)
D.text(s, M + 0.42, 6.00, W - 2 * M - 0.85, 0.5,
       "The 5 per cent level was itself tested. Tightening it to 2 per cent was tried and rejected: "
       "it pushed the mu-opioid threshold to 0.999 and caused morphine to be missed.",
       size=12.5, color=INK, italic=True, line=1.24)
D.notes(s, "Nav1.1 was later withdrawn, but the constraint it motivated applies to every endpoint.")

# 5 -- the rate measured
s = D.light()
D.head(s, "4", "The rate that is measured, against the rate that was targeted",
       "This is what the partition buys")
i, e = F["ins"], F["ev"]
D.text(s, M + 0.30, 1.90, 3.4, 0.3, "", size=11)
D.text(s, M + 3.90, 1.90, 2.6, 0.32, "on the threshold pool,\nin sample", size=11, bold=True,
       color=MUTED, align=PP_ALIGN.RIGHT, line=1.1)
D.text(s, M + 7.00, 1.90, 2.6, 0.32, "on the DISJOINT\nevaluation pool", size=11, bold=True,
       color=DEEP, align=PP_ALIGN.RIGHT, line=1.1)
lines = [("mean", fmt(i[0], 4), fmt(e[0], 4)), ("median", fmt(i[1], 4), fmt(e[1], 4)),
         ("minimum", fmt(i[2], 4), fmt(e[2], 4)), ("maximum", fmt(i[3], 4), fmt(e[3], 4)),
         ("above 0.05", "most", f"{e[5]} of {e[4]}"),
         ("above 0.10", "most", f"{e[6]} of {e[4]}")]
yy = 2.60
for k, (lab, a, b) in enumerate(lines):
    if k % 2 == 0:
        D.card(s, M, yy - 0.09, 9.9, 0.54, fill=TINT)
    D.text(s, M + 0.30, yy, 3.4, 0.34, lab, size=13, color=INK)
    D.text(s, M + 3.90, yy, 2.6, 0.34, a, size=13, color=MUTED, align=PP_ALIGN.RIGHT)
    D.text(s, M + 7.00, yy, 2.6, 0.34, b, size=13, bold=True, color=DEEP, align=PP_ALIGN.RIGHT)
    yy += 0.54
D.text(s, M, 6.00, 9.9, 0.85,
       f"The measured rate is lower on all {F['lower_on']} of {F['n_pair']} endpoints, by a mean of "
       f"{abs(F['gap_fpr']):.3f}. The in-sample figure says the cut is where we put it. The "
       f"evaluation figure says the panel fires on about 2 per cent of unrelated chemistry, and only "
       "the second is evidence.",
       size=13.5, color=INK, line=1.28)
xs = M + 10.15
D.text(s, xs, 2.55, W - M - xs, 1.5,
       "Highest measured rates sit at the aminergic receptors, which bind a large and chemically "
       "ordinary region of drug-like space.",
       size=12, color=MUTED, line=1.24)
yy = 4.00
for k, a, b in F["hi_bg"][:3]:
    D.text(s, xs, yy, 1.4, 0.3, k.replace("_", "-"), size=11.5, color=INK)
    D.text(s, xs + 1.35, yy, 1.0, 0.3, fmt(b, 4), size=11.5, bold=True, color=AMBER,
           align=PP_ALIGN.RIGHT)
    yy += 0.34
D.source(s, "submission_package/07_MODELS/binder_panel_registry.json")
D.notes(s, "The systematic direction is not an error: for 39 of 47 endpoints the binding constraint "
           "is the measured-inactive quantile, a much stricter cut than a background quantile.")

# 6 -- the finding
s = D.dark()
D.text(s, M, 0.92, 9.0, 0.5, "One quantity, four numbers", size=15, color=AMBER, bold=True)
D.text(s, M, 1.34, 11.9, 0.76,
       "The published sensitivity is measured on the training set.",
       size=28, bold=True, font=HEAD, color=PAPER, line=1.10)
D.text(s, M, 2.24, 11.4, 0.85,
       "The threshold sequence runs final_thresholds.py and then calibrate_background_specificity.py. "
       "Both write sensitivity_at_threshold into the registry. The second wins, and the two select "
       "their actives differently.",
       size=13.5, color=CHALK, line=1.26)
D.card(s, M, 3.30, 5.6, 1.55, fill=DEEP)
D.text(s, M + 0.35, 3.52, 5.0, 0.35, "final_thresholds.py", size=13.5, bold=True, color=PAPER)
D.text(s, M + 0.35, 3.92, 5.0, 0.85,
       "scores the held-out actives listed in models_rf/holdout/. For A1 that is 420 compounds.",
       size=12.5, color=CHALK, line=1.24)
D.card(s, M + 6.0, 3.30, 5.7, 1.55, fill=DEEP)
D.text(s, M + 6.35, 3.52, 5.1, 0.35, "calibrate_background_specificity.py", size=13.5, bold=True,
       color=AMBER)
D.text(s, M + 6.35, 3.92, 5.1, 0.85,
       "scores every active in the endpoint table. For A1 that is 1,943 compounds, of which 1,523 "
       "were used to fit the model.",
       size=12.5, color=CHALK, line=1.24)
D.text(s, M, 5.15, 11.6, 0.55,
       f"The registry keeps the second script's number and the first script's label. All "
       f"{F['n_sens']} deployed entries still read sensitivity_basis: held_out_actives_by_scaffold.",
       size=14, bold=True, color=PAPER, line=1.24)
facts_line = [(f"published higher on {F['higher_on']} of {F['n_sens']}", ""),
              (f"by a mean of {F['mean_gap']:+.4f}", ""),
              (f"the scored set is {F['train_frac'][0]:.1%} training compounds at the median", "")]
yy = 5.90
for t, _ in facts_line:
    D.text(s, M + 0.15, yy, 11.4, 0.32, "•  " + t, size=13, color=AMBER)
    yy += 0.36
D.source(s, "results/tables/sensitivity_reconciliation.csv, written for this chapter")
D.notes(s, "This is a correction to a published number. The panel's behaviour is unchanged: no "
           "threshold moves and no prediction changes. What changes is the number describing them.")

# 7 -- the honest figure
s = D.light()
D.head(s, "5", "What the honest figure is, and what it changes",
       f"Measured on held-out actives only, across the {F['n_sens']} deployed binder endpoints")
p, h = F["pub"], F["held"]
D.text(s, M + 4.55, 1.88, 2.3, 0.3, "as published", size=11.5, bold=True, color=MUTED,
       align=PP_ALIGN.RIGHT)
D.text(s, M + 7.30, 1.88, 2.3, 0.3, "held out only", size=11.5, bold=True, color=DEEP,
       align=PP_ALIGN.RIGHT)
lines = [("mean sensitivity", fmt(p[0], 4), fmt(h[0], 4), False),
         ("median", fmt(p[1], 4), fmt(h[1], 4), False),
         ("minimum", f"{fmt(p[2])}, {F['pub_min_ep']}", f"{fmt(h[2])}, {F['held_min_ep'].replace('_','-')}", False),
         ("maximum", fmt(p[3]), fmt(h[3]), False),
         ("endpoints below 0.50", f"{p[4]} of {F['n_sens']}", f"{h[4]} of {F['n_sens']}", True)]
yy = 2.32
for k, (lab, a, b, hi) in enumerate(lines):
    if k % 2 == 0:
        D.card(s, M, yy - 0.10, 9.8, 0.58, fill=TINT)
    D.text(s, M + 0.30, yy, 4.1, 0.36, lab, size=13.5, bold=hi, color=CRIMSON if hi else INK)
    D.text(s, M + 4.55, yy, 2.3, 0.36, a, size=13.5, color=MUTED, bold=hi, align=PP_ALIGN.RIGHT)
    D.text(s, M + 7.30, yy, 2.3, 0.36, b, size=13.5, bold=True,
           color=CRIMSON if hi else DEEP, align=PP_ALIGN.RIGHT)
    yy += 0.58
D.text(s, M, 5.35, 9.8, 0.5,
       "The last row is why this is a correction and not a quibble.",
       size=15, bold=True, font=HEAD, color=CRIMSON)
low = ", ".join(f"{k.replace('_','-')} {v:.3f}" for k, v in F["low6"])
D.text(s, M, 5.90, 9.8, 0.85,
       f"On the published figure no deployed endpoint fires for fewer than half its own actives. On "
       f"the held-out figure six do: {low}.",
       size=13, color=INK, line=1.26)
xs = M + 10.05
D.text(s, xs, 2.30, W - M - xs, 3.0,
       "The falsification suite already found this, from an unrelated direction.\n\nH7 scored "
       "held-out actives against random PubChem chemistry and reported six targets under 0.50. Five "
       "of its six are the same six.\n\nThe published figure was the outlier all along, and the suite "
       "designed to embarrass the tool had already said so.",
       size=12, color=INK, line=1.26)
D.source(s, "results/tables/sensitivity_reconciliation.csv; inversion/results/H7_target_discrimination.csv")
D.notes(s, "The convergence with H7 is the strongest evidence that the held-out figure is the right "
           "one. Two unrelated constructions land within 0.03 of each other.")

# 8 -- unreliable endpoints
s = D.light()
D.head(s, "6", "The endpoints the panel already marks unreliable",
       "Sensitivity below 0.50 at the operating point, or AUROC below 0.75 against measured inactives")
D.text(s, M, 1.85, 5.7, 1.5,
       f"Six of the {F['n_dep']} deployed endpoints carry the flag: "
       f"{', '.join(x.replace('_','-') for x in F['unreliable'])}.",
       size=15, bold=True, color=DEEP, line=1.26)
D.text(s, M, 3.15, 5.7, 2.1,
       "With one substitution these are the same six that fire for fewer than half their own held-out "
       "actives. So the flag is doing its job, and the panel already knows which endpoints are weak.\n\n"
       "What it does not do is stop deploying them.",
       size=13.5, color=INK, line=1.28)
xs = M + 6.1
D.card(s, xs, 1.82, W - M - xs, 4.15, fill=TINT)
D.text(s, xs + 0.42, 2.06, 5.2, 0.42, "That is defensible only because the flag travels",
       size=16, bold=True, font=HEAD, color=AMBER, line=1.14)
D.text(s, xs + 0.42, 2.85, 5.2, 1.6,
       "An unreliable endpoint reporting nothing is not evidence of inactivity. At these six, a "
       "silence carries almost no information, and the interface has to say so at exactly those "
       "six rather than in general terms.",
       size=13.5, color=INK, line=1.28)
D.text(s, xs + 0.42, 4.45, 5.2, 1.3,
       "This is the same argument as Chapter 4's reading order, arriving from the threshold side: "
       "the number that makes a silence interpretable has to be attached to the silence, not "
       "averaged over the panel.",
       size=12.5, color=MUTED, italic=True, line=1.26)
D.source(s, "results/tables/final_thresholds.csv, reliable column")
D.notes(s, "MIN_SENS is 0.50 and the AUROC gate is 0.75. Both are stated in the code, not tuned "
           "per endpoint.")

# 9 -- what the tests pin
s = D.light()
D.head(s, "7", "What the tests pin, and what they missed", None)
D.card(s, M, 1.82, 5.85, 2.4)
D.text(s, M + 0.42, 2.06, 5.0, 0.4, "What they pin", size=17, bold=True, font=HEAD, color=DEEP)
D.text(s, M + 0.42, 2.58, 5.05, 1.5,
       "TestBackgroundPools holds five subtests: assignment is a pure function of the structure and "
       "independent of call order, every structure gets exactly one role, the bands cover exactly "
       "one hundred, and the hash band is stable and in range.\n\n"
       "TestThresholdSequenceIsAtomic holds two: the four steps stay one unit, and no member depends "
       "on a file the sequence itself rewrites.",
       size=12.5, color=INK, line=1.26)
D.card(s, M + 6.25, 1.82, W - 2 * M - 6.25, 2.4, fill=TINT)
D.text(s, M + 6.67, 2.06, 5.0, 0.4, "What they missed", size=17, bold=True, font=HEAD, color=CRIMSON)
D.text(s, M + 6.67, 2.58, 5.05, 1.5,
       "That two steps in a correctly ordered sequence can write the same field from different "
       "populations, and that the label written by the first survives the value written by the "
       "second.\n\nThe ordering was right. The semantics were not, and nothing asserted them.",
       size=12.5, color=INK, line=1.26)
D.text(s, M, 4.50, W - 2 * M, 0.5,
       "The specific test this chapter recommends",
       size=17, bold=True, font=HEAD, color=INK)
D.text(s, M, 5.05, W - 2 * M, 1.0,
       "Assert that sensitivity_basis describes the set the stored sensitivity_at_threshold was "
       "actually computed on, and that the two scripts agree on what an active is. Three endpoints, "
       "Nav1.5, SIRT1 and TAAR1, currently do not: the binder pipeline selects training actives by "
       "the label column while the specificity script selects by pChEMBL at or above 7.",
       size=13.5, color=INK, line=1.28)
D.text(s, M, 6.20, W - 2 * M, 0.5,
       "A test that pins an ordering is not a test that pins a meaning.",
       size=14, color=AMBER, italic=True)
D.notes(s, "Four of the 47 thresholds also differ between the registry and final_thresholds.csv, at "
           "D3, 5-HT1A, Nav1.5 and SIRT1. That is correct behaviour, the registry holding the later "
           "value, but nothing asserts it.")

# 10 -- closing
s = D.dark()
D.text(s, M, 1.50, 8.4, 0.5, "Chapter 5, in one statement", size=15, color=TEAL, bold=True)
D.text(s, M, 2.05, 11.2, 2.2,
       "The partition exists so that a false-positive rate can disagree with its target. The "
       "sensitivity figure had no such protection, and did not.",
       size=30, font=HEAD, color=PAPER, line=1.18)
D.text(s, M, 4.35, 11.2, 1.8,
       "Three disjoint pools stop a rate restating the quantile that produced it. Nothing played the "
       "same role for sensitivity, so a value computed on the training set could be stored under a "
       "label saying held out, and survive four documents and a submission package.",
       size=16, color=CHALK, line=1.30)
D.text(s, M, 6.30, 11.2, 0.5,
       f"Held out, the panel fires for {F['held'][0]:.3f} of the actives it should, not "
       f"{F['pub'][0]:.3f}.",
       size=14, color=AMBER, italic=True)
D.notes(s, "End here. The design principle that protected one number is the one that was missing "
           "for the other.")

D.save(OUT)
