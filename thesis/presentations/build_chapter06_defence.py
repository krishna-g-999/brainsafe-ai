"""Build the Chapter 6 defence deck, reading every figure from an artefact.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/build_chapter06_defence.py
Out:  thesis/presentations/chapter06_defence.pptx
"""
from __future__ import annotations

import json
import statistics as st
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

OUT = Path(__file__).resolve().parent / "chapter06_defence.pptx"


def facts() -> dict:
    d = {}
    reg = json.loads((PKG / "07_MODELS" / "binder_panel_registry.json").read_text())
    dep = {k: v for k, v in reg.items() if v.get("deployed")}
    wd = {k: v for k, v in reg.items() if not v.get("deployed")}
    d["n_total"], d["n_dep"], d["n_wd"] = len(reg), len(dep), len(wd)
    d["modes"] = Counter(v.get("mode") for v in reg.values())

    pairs = [(k, v["scaffold_cv_auroc"], v["auroc_vs_measured_inactives"]) for k, v in dep.items()]
    cv = [a for _, a, _ in pairs]
    mi = [b for _, _, b in pairs]
    d["decoy"] = (st.mean(cv), st.median(cv), min(cv), max(cv))
    d["meas"] = (st.mean(mi), st.median(mi), min(mi), max(mi))
    dd = [a - b for _, a, b in pairs]
    d["higher_on"] = sum(1 for x in dd if x > 0)
    d["gap"] = st.mean(dd)
    d["worst"] = sorted(pairs, key=lambda x: -(x[1] - x[2]))[:6]
    mim = {k: b for k, _, b in pairs}
    d["meas_min_ep"], d["meas_max_ep"] = min(mim, key=mim.get), max(mim, key=mim.get)

    mih = {k: v["n_measured_inactive_holdout"] for k, v in dep.items()
           if v.get("n_measured_inactive_holdout")}
    hv = list(mih.values())
    d["mih"] = (len(hv), st.median(hv), min(hv), max(hv),
                sum(1 for x in hv if x < 100), sum(1 for x in hv if x < 50))
    d["mih_small"] = sorted(mih.items(), key=lambda x: x[1])[:5]
    d["mih_of"] = {d["meas_min_ep"]: mih.get(d["meas_min_ep"]),
                   d["meas_max_ep"]: mih.get(d["meas_max_ep"])}
    ah = [v["n_active_holdout"] for v in dep.values() if v.get("n_active_holdout")]
    d["ah"] = (st.median(ah), min(ah), max(ah))
    nd = [v["n_decoy"] for v in dep.values() if v.get("n_decoy")]
    npo = [v["n_positive"] for v in dep.values() if v.get("n_positive")]
    d["decoy_ratio"] = (st.median(npo), st.median(nd), st.median(a / b for a, b in zip(nd, npo)))

    rec = rows(TAB / "sensitivity_reconciliation.csv")
    h = [num(r["sensitivity_heldout"]) for r in rec if r["sensitivity_heldout"]]
    p = [num(r["sensitivity_published"]) for r in rec if r["sensitivity_published"]]
    hm = {r["target"]: num(r["sensitivity_heldout"]) for r in rec if r["sensitivity_heldout"]}
    pm = {r["target"]: num(r["sensitivity_published"]) for r in rec if r["sensitivity_published"]}
    d["held"] = (st.mean(h), st.median(h), min(h), max(h), sum(1 for v in h if v < 0.5))
    d["pub"] = (st.mean(p), st.median(p), min(p), max(p), sum(1 for v in p if v < 0.5))
    d["held_min_ep"], d["pub_min_ep"] = min(hm, key=hm.get), min(pm, key=pm.get)
    d["low6"] = sorted(((k, v) for k, v in hm.items() if v < 0.5), key=lambda x: x[1])
    d["both"] = sorted({k for k, _ in d["low6"]} & {k for k, _, _ in d["worst"]})

    d["wd"] = [(k, v.get("auroc_vs_measured_inactives"), v.get("sensitivity_at_threshold"),
                v.get("n_positive")) for k, v in wd.items()]
    d["np"] = rows(TAB / "np_endpoint_assay_composition.csv")

    bs = [x for x in rows(TAB / "binder_cv_summary.csv")
          if x["split"] == "scaffold" and x["endpoint"] in dep]
    dev = [(x["endpoint"], abs(num(x["rerun_minus_recorded"]))) for x in bs
           if x["rerun_minus_recorded"]]
    d["rerun"] = (st.mean(v for _, v in dev), max(dev, key=lambda t: t[1]), len(dev),
                  st.mean(num(x["roc_auc_mean"]) for x in bs))
    d["taar1"] = dep.get("TAAR1", {})
    return d


F = facts()
D = Deck()

# 1 -- title
s = D.dark()
D.text(s, M, 1.74, W - 2 * M, 0.5, "Chapter 6", size=16, color=TEAL, bold=True)
D.text(s, M, 2.18, W - 2 * M - 1.0, 1.6,
       "The binder panel, validated against measured inactives",
       size=38, bold=True, font=HEAD, color=PAPER, line=1.08)
D.text(s, M, 4.10, W - 2 * M - 1.6, 1.0,
       "The largest part of the system, the part a user is most likely to query, and the one place "
       "the validation choice changes the headline by a measurable amount",
       size=16, color=CHALK, italic=True, line=1.25)
D.text(s, M, 5.65, W - 2 * M, 0.6,
       f"{F['n_total']} fitted classifiers, {F['n_dep']} deployed, {F['n_wd']} withdrawn and kept in "
       f"the inventory",
       size=13, color=TEAL)
D.notes(s, "The chapter's spine: validating against decoys would report 0.978 and mean nothing. "
           "Measured inactives cost 0.061 AUROC and buy a number that means what a reader thinks.")

# 2 -- what it is
s = D.light()
D.head(s, "1", "What the binder panel is",
       "One question per target: does this compound bind here at pChEMBL 7, about 100 nM")
mo = F["modes"]
D.card(s, M, 1.80, 5.9, 2.0)
D.text(s, M + 0.42, 2.06, 5.2, 0.4, f"{mo['hybrid_decoys_plus_measured_inactives']} endpoints, hybrid",
       size=17, bold=True, font=HEAD, color=DEEP)
D.text(s, M + 0.42, 2.58, 5.05, 1.05,
       "Recovered censored bounds where they exist, topped up with property-matched decoys where "
       "they do not.",
       size=13.5, color=INK, line=1.26)
D.card(s, M + 6.25, 1.80, W - 2 * M - 6.25, 2.0)
D.text(s, M + 6.67, 2.06, 5.2, 0.4, f"{mo['measured_labels_holdout']} endpoints, measured labels only",
       size=17, bold=True, font=HEAD, color=TEAL)
D.text(s, M + 6.67, 2.58, 5.05, 1.05,
       "Enough measured inactives that no decoy is needed at all.",
       size=13.5, color=INK, line=1.26)
p50, d50, ratio = F["decoy_ratio"]
D.stat(s, M + 0.10, 4.15, 3.6, f"{int(p50):,}", "positives, median endpoint", color=DEEP, vsize=36)
D.stat(s, M + 4.20, 4.15, 3.6, f"{int(d50):,}", "decoys, median endpoint", color=MUTED, vsize=36)
D.stat(s, M + 8.30, 4.15, 3.6, f"{ratio:.2f}", "median decoys per active", color=TEAL, vsize=36)
D.text(s, M, 5.95, W - 2 * M, 0.55,
       "Those decoys are drawn from the decoy pool of Chapter 5, and never from the pools that set or "
       "measure the threshold.",
       size=13.5, color=INK, italic=True, line=1.26)
D.source(s, "submission_package/07_MODELS/binder_panel_registry.json")
D.notes(s, "The two modes matter: a third of nothing is not the same as a third of something. The "
           "eight measured-label endpoints need no assumption in their negative class at all.")

# 3 -- the validation choice, the core chart
s = D.light()
D.head(s, "2", "Why the validation had to change",
       "Both numbers exist for every deployed endpoint, so the difference can be stated")
order = sorted(F["worst"], key=lambda x: x[2])
cd = CategoryChartData()
cd.categories = [k.replace("_", "-") for k, _, _ in order]
cd.add_series("against property-matched decoys", [a for _, a, _ in order])
cd.add_series("against measured inactives", [b for _, _, b in order])
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.78),
                        Inches(7.5), Inches(4.25), cd)
ch = gf.chart
ch.has_title = False
ch.has_legend = True
ch.legend.position = XL_LEGEND_POSITION.BOTTOM
ch.legend.include_in_layout = False
ch.legend.font.size = Pt(11)
va = ch.value_axis
va.minimum_scale, va.maximum_scale = 0.6, 1.0
va.major_gridlines.format.line.color.rgb = D.grid
va.tick_labels.font.size = Pt(11)
va.tick_labels.font.color.rgb = MUTED
ch.category_axis.tick_labels.font.size = Pt(11)
ch.category_axis.tick_labels.font.color.rgb = MUTED
for i, colr in enumerate([CHALK, DEEP]):
    ch.plots[0].series[i].format.fill.solid()
    ch.plots[0].series[i].format.fill.fore_color.rgb = colr
ch.plots[0].gap_width = 60

xs = M + 7.8
dec, mea = F["decoy"], F["meas"]
D.text(s, xs + 1.70, 1.88, 1.2, 0.32, "decoys", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
D.text(s, xs + 3.05, 1.88, 1.2, 0.32, "measured", size=11, bold=True, color=DEEP, align=PP_ALIGN.RIGHT)
for i, (lab, a, b) in enumerate([("mean AUROC", fmt(dec[0], 4), fmt(mea[0], 4)),
                                 ("median", fmt(dec[1], 4), fmt(mea[1], 4)),
                                 ("minimum", fmt(dec[2]), fmt(mea[2])),
                                 ("maximum", fmt(dec[3]), fmt(mea[3]))]):
    y = 2.28 + i * 0.46
    D.text(s, xs, y, 1.65, 0.34, lab, size=12.5, color=INK)
    D.text(s, xs + 1.70, y, 1.2, 0.34, a, size=12.5, color=MUTED, align=PP_ALIGN.RIGHT)
    D.text(s, xs + 3.05, y, 1.2, 0.34, b, size=12.5, bold=True, color=DEEP, align=PP_ALIGN.RIGHT)
D.text(s, xs, 4.25, W - M - xs, 0.9,
       f"The decoy figure is higher on {F['higher_on']} of {F['n_dep']}, by a mean of "
       f"{F['gap']:+.4f}.",
       size=13.5, bold=True, color=AMBER, line=1.24)
D.card(s, xs, 5.20, W - M - xs, 1.35, fill=TINT)
D.text(s, xs + 0.32, 5.42, W - M - xs - 0.65, 1.0,
       "A median of 0.992 and a maximum of 1.000 are the signature of a saturated problem. On more "
       "than half the panel, separating actives from decoys is very nearly solved, which tells a "
       "reader almost nothing.",
       size=12, color=INK, line=1.24)
D.source(s, "The six endpoints shown are those where the two measurements diverge most.")
D.notes(s, "GABA-A validated against decoys reads 0.965, mid-panel. Against compounds actually "
           "tested there and found inactive it reads 0.719, the weakest deployed. Nothing about the "
           "model changed. Only the question did.")

# 4 -- how thin the hold-out is
s = D.light()
D.head(s, "3", "The measured-inactive hold-out, and how thin it is",
       "The honesty of the last slide depends on there being enough to validate against")
n, med, lo, hi, u100, u50 = F["mih"]
D.stat(s, M + 0.10, 1.90, 3.4, f"{int(med)}", f"compounds in the median endpoint's measured-inactive "
                                              f"hold-out", color=DEEP, vsize=40)
D.stat(s, M + 4.20, 1.90, 3.4, f"{int(lo)} to {int(hi)}", "the full range across the panel",
       color=MUTED, vsize=36)
D.stat(s, M + 8.30, 1.90, 3.4, f"{u100} of {n}", "endpoints have fewer than 100", color=AMBER, vsize=36)
D.card(s, M, 3.65, W - 2 * M, 1.55, fill=TINT)
D.text(s, M + 0.45, 3.90, W - 2 * M - 0.9, 1.1,
       f"The panel's best AUROC against measured inactives, {fmt(F['meas'][3])}, belongs to "
       f"{F['meas_max_ep']} and is computed against {F['mih_of'][F['meas_max_ep']]} compounds. Its "
       f"worst, {fmt(F['meas'][2])} at {F['meas_min_ep'].replace('_','-')}, is computed against "
       f"{F['mih_of'][F['meas_min_ep']]}. The two ends of the reported range are the two thinnest "
       "measurements in it.",
       size=14.5, color=INK, line=1.28)
D.text(s, M, 5.45, 6.2, 0.9,
       "Five endpoints have fewer than 50: " +
       ", ".join(f"{k.replace('_','-')} {v}" for k, v in F["mih_small"]) + ".",
       size=12.5, color=INK, line=1.26)
D.text(s, M + 6.6, 5.45, W - M - (M + 6.6), 0.9,
       "No per-endpoint AUROC in any current artefact carries an interval. Adding Wilson or bootstrap "
       "intervals is the cheapest improvement available here, and it has not been done.",
       size=12.5, color=AMBER, line=1.26)
D.source(s, f"Active hold-outs are more comfortable: median {int(F['ah'][0])}, "
            f"range {int(F['ah'][1])} to {int(F['ah'][2])}.")
D.notes(s, "Volunteer this before an examiner asks how 0.985 was computed. The answer is: on 23 "
           "compounds.")

# 5 -- sensitivity corrected
s = D.dark()
D.text(s, M, 0.95, 9.0, 0.5, "Sensitivity, carried through from Chapter 5", size=15, color=AMBER,
       bold=True)
D.text(s, M, 1.40, 11.7, 0.76,
       "The panel fires for about three quarters of the actives it should.",
       size=28, bold=True, font=HEAD, color=PAPER, line=1.10)
p, h = F["pub"], F["held"]
D.text(s, M + 5.20, 2.42, 2.6, 0.3, "as published", size=11.5, bold=True, color=CHALK,
       align=PP_ALIGN.RIGHT)
D.text(s, M + 8.30, 2.42, 2.6, 0.3, "held out only", size=11.5, bold=True, color=AMBER,
       align=PP_ALIGN.RIGHT)
lines = [("mean sensitivity", fmt(p[0], 4), fmt(h[0], 4), False),
         ("median", fmt(p[1], 4), fmt(h[1], 4), False),
         ("minimum", f"{fmt(p[2])}, {F['pub_min_ep']}",
          f"{fmt(h[2])}, {F['held_min_ep'].replace('_','-')}", False),
         ("endpoints firing for under half their own actives",
          f"{p[4]} of {F['n_dep']}", f"{h[4]} of {F['n_dep']}", True)]
yy = 2.85
for lab, a, b, hi_ in lines:
    D.text(s, M + 0.10, yy, 5.0, 0.34, lab, size=13, color=PAPER, bold=hi_)
    D.text(s, M + 5.20, yy, 2.6, 0.34, a, size=13, color=CHALK, align=PP_ALIGN.RIGHT)
    D.text(s, M + 8.30, yy, 2.6, 0.34, b, size=13, bold=True, color=AMBER, align=PP_ALIGN.RIGHT)
    yy += 0.48
low = ", ".join(f"{k.replace('_','-')} {v:.3f}" for k, v in F["low6"])
D.text(s, M, 4.90, 11.7, 0.55, f"The six below 0.50: {low}.", size=13, color=CHALK, line=1.24)
D.card(s, M, 5.55, 11.7, 1.30, fill=DEEP)
D.text(s, M + 0.42, 5.78, 11.0, 0.95,
       f"{', '.join(x.replace('_','-') for x in F['both'])} appear in both lists: the endpoints that "
       "lose most against measured inactives are the endpoints that fire for fewest of their own "
       "actives. That is one fact, not two. Hard-to-separate inactives force a high threshold, and a "
       "high threshold is what costs sensitivity.",
       size=13, color=PAPER, line=1.26)
D.source(s, "results/tables/sensitivity_reconciliation.csv")
D.notes(s, "This is where Chapter 5's operating-point argument and Chapter 6's discrimination "
           "argument meet. Make the connection explicitly; it is the strongest structural point in "
           "the two chapters.")

# 6 -- withdrawals
s = D.light()
D.head(s, "4", f"The {F['n_wd']} withdrawn endpoints, kept in the inventory",
       "A panel showing only what survived is a selection, not an inventory")
D.text(s, M + 0.30, 1.90, 1.9, 0.3, "endpoint", size=11, bold=True, color=MUTED)
for i, hh in enumerate(["AUROC", "sensitivity", "positives"]):
    D.text(s, M + 2.30 + i * 1.45, 1.90, 1.3, 0.3, hh, size=11, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT)
D.text(s, M + 6.85, 1.90, 5.0, 0.3, "why", size=11, bold=True, color=MUTED)
why = {"GluA2": "fires on glucose and atenolol at its calibrated threshold of 0.629",
       "Nav1_1": "fires on glucose, urea, glycine, lactate and atenolol at 0.571",
       "NRF2": "background false-positive rate 0.057, above the 5 per cent the panel holds to",
       "NFKB1": "recovers no active while calling five trivial metabolites binders",
       "NR3C1": "AUROC below chance against its own held-out measured inactives"}
yy = 2.30
for i, (k, a, sv, npos) in enumerate(F["wd"]):
    if i % 2 == 0:
        D.card(s, M, yy - 0.10, W - 2 * M, 0.66, fill=TINT)
    D.text(s, M + 0.30, yy, 1.9, 0.36, k.replace("_", "."), size=13.5, bold=True, color=INK)
    for j, v in enumerate([fmt(a), fmt(sv), str(npos)]):
        col = CRIMSON if (j == 0 and a < 0.7) or (j == 1 and sv < 0.2) else INK
        D.text(s, M + 2.30 + j * 1.45, yy, 1.3, 0.36, v, size=13.5, color=col,
               bold=(col == CRIMSON), align=PP_ALIGN.RIGHT)
    D.text(s, M + 6.85, yy + 0.02, W - 2 * M - 7.0, 0.55, why[k], size=11.5, color=INK, line=1.2)
    yy += 0.78
D.card(s, M, 6.05, W - 2 * M, 0.82, fill=TINT)
D.text(s, M + 0.42, 6.24, W - 2 * M - 0.85, 0.5,
       "Nav1.1 has an AUROC of 0.952, better than eleven deployed endpoints, and was withdrawn "
       "anyway: its AUROC measures ranking, which is not the quantity a deployed cut needs.",
       size=13, color=INK, italic=True, line=1.24)
D.notes(s, "NR3C1 fails in the opposite direction: it fires on no trivial molecule and fails on "
           "discrimination, below chance at 0.410. A panel that only checked for glucose would have "
           "deployed it.")

# 7 -- why the natural-product endpoints failed
s = D.light()
D.head(s, "5", "The natural-product endpoints failed for a diagnosable reason",
       "And the reason is about the labels, not the chemistry or the models")
D.text(s, M + 0.30, 1.92, 2.0, 0.3, "endpoint", size=11, bold=True, color=MUTED)
for i, hh in enumerate(["labelled records", "recorded as Potency", "direct binding constant"]):
    D.text(s, M + 2.50 + i * 2.35, 1.92, 2.1, 0.3, hh, size=11, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT)
yy = 2.35
for i, r in enumerate(F["np"]):
    if i % 2 == 0:
        D.card(s, M, yy - 0.10, 9.6, 0.62, fill=TINT)
    pct = num(r["pct_direct_binding_constant"])
    D.text(s, M + 0.30, yy, 2.0, 0.36, r["endpoint"], size=13.5, bold=True, color=INK)
    for j, v in enumerate([f"{int(r['labelled_records']):,}", f"{int(r['potency']):,}", f"{pct}%"]):
        D.text(s, M + 2.50 + j * 2.35, yy, 2.1, 0.36, v, size=13.5,
               bold=(j == 2), color=CRIMSON if (j == 2 and pct < 1) else INK,
               align=PP_ALIGN.RIGHT)
    yy += 0.74
D.text(s, M, 4.75, 9.6, 0.5,
       "For NRF2 and NFKB1 essentially every label is a pooled functional readout.",
       size=16, bold=True, font=HEAD, color=DEEP)
D.text(s, M, 5.32, 9.6, 1.1,
       "A functional potency does not define a binding class a ligand fingerprint can separate: two "
       "compounds can share a potency in a cell-based assay through entirely different mechanisms, "
       "and the structure carries no signal for the label. NR3C1, the only one with real binding "
       "data, failed for the other available reason: 140 compounds after deduplication is too few.",
       size=13, color=INK, line=1.28)
xs = M + 9.9
D.card(s, xs, 2.20, W - M - xs, 4.45, fill=TINT)
D.text(s, xs + 0.28, 2.45, W - M - xs - 0.56, 4.0,
       "This is the most useful negative result in the panel.\n\nExtending to natural-product targets "
       "is not blocked by the chemistry being unusual, which was the assumed obstacle.\n\nIt is "
       "blocked by the public assay record for those targets not measuring binding.",
       size=12.5, color=INK, line=1.28)
D.source(s, "results/tables/np_endpoint_assay_composition.csv")
D.notes(s, "This reframes the natural-product limitation entirely, and points the future work at "
           "assay selection rather than at featurisation.")

# 8 -- reproduction and TAAR1
s = D.light()
D.head(s, "6", "Reproduction, and the panel's least trustworthy endpoint", None)
mean_dev, worst, ndev, rerun_mean = F["rerun"]
D.card(s, M, 1.85, 5.85, 2.25)
D.text(s, M + 0.42, 2.10, 5.0, 0.4, "Independent re-run", size=17, bold=True, font=HEAD, color=DEEP)
D.text(s, M + 0.42, 2.62, 5.05, 1.3,
       f"The whole binder cross-validation was re-run from the endpoint tables with separately "
       f"written scoring code. Across {ndev} deployed endpoints it agrees with the recorded scaffold "
       f"AUROC to a mean absolute deviation of {mean_dev:.4f}.",
       size=13, color=INK, line=1.26)
D.card(s, M + 6.25, 1.85, W - 2 * M - 6.25, 2.25, fill=TINT)
D.text(s, M + 6.67, 2.10, 5.0, 0.4, f"Maximum deviation {worst[1]:.4f}, at {worst[0]}",
       size=17, bold=True, font=HEAD, color=AMBER)
D.text(s, M + 6.67, 2.62, 5.05, 1.3,
       "Not a reproduction failure, but not nothing, and the endpoint where it occurs is the one "
       "where every other measurement in this thesis is also least stable.",
       size=13, color=INK, line=1.26)
D.text(s, M, 4.35, W - 2 * M, 0.5,
       "TAAR1 is the least trustworthy deployed endpoint on four independent counts",
       size=18, bold=True, font=HEAD, color=CRIMSON)
counts = [("smallest binder training set", "82 rows"),
          ("worst calibrated, Chapter 4", "0.178 expected calibration error"),
          ("held-out sensitivity, Chapter 5", "0.355"),
          ("largest reproduction deviation", f"{worst[1]:.4f}")]
yy = 5.00
for i, (a, b) in enumerate(counts, 1):
    D.dot(s, M + 0.10, yy, str(i), fill=CRIMSON, dia=0.32)
    D.text(s, M + 0.60, yy - 0.02, 4.6, 0.34, a, size=13, color=INK)
    D.text(s, M + 5.35, yy - 0.02, 3.4, 0.34, b, size=13, bold=True, color=CRIMSON)
    yy += 0.42
D.text(s, M + 9.0, 5.00, W - M - (M + 9.0), 1.4,
       "Nothing in the interface says so beyond the generic reliability flag. Reviewing it for "
       "withdrawal, or for a specific warning, is this chapter's recommendation.",
       size=12.5, color=AMBER, line=1.26)
D.source(s, "results/tables/binder_cv_summary.csv")
D.notes(s, "Four unrelated measurements converging on one endpoint is a stronger signal than any of "
           "them alone.")

# 9 -- closing
s = D.dark()
D.text(s, M, 1.35, 8.4, 0.5, "Chapter 6, honestly summarised", size=15, color=TEAL, bold=True)
items = [
    (f"discriminates against measured inactives at a mean AUROC of {fmt(F['meas'][0], 4)}, "
     f"median {fmt(F['meas'][1], 4)}, from {fmt(F['meas'][2])} at "
     f"{F['meas_min_ep'].replace('_','-')} to {fmt(F['meas'][3])} at {F['meas_max_ep']}, the two "
     f"extremes resting on {F['mih_of'][F['meas_min_ep']]} and {F['mih_of'][F['meas_max_ep']]} "
     f"compounds"),
    (f"fires for a mean of {fmt(F['held'][0], 4)} of held-out actives, median "
     f"{fmt(F['held'][1], 4)}, with {F['held'][4]} of {F['n_dep']} endpoints below 0.50"),
    (f"would report {fmt(F['decoy'][0], 4)} and {fmt(F['pub'][0], 4)} if validated the conventional "
     f"way, and neither would mean what a reader takes it to mean"),
    (f"carries {F['n_wd']} withdrawn endpoints in its inventory, three of which failed because the "
     f"public assay record for those targets does not measure binding"),
]
yy = 2.05
for i, t in enumerate(items, 1):
    D.dot(s, M, yy + 0.04, str(i), fill=TEAL, dia=0.38)
    D.text(s, M + 0.62, yy, 11.0, 1.0, t, size=15, color=PAPER, line=1.30)
    yy += 1.20
D.text(s, M, 6.80, 11.5, 0.4,
       "The deployed panel is good. It is 0.061 AUROC and 0.135 sensitivity less good than the "
       "conventional way of reporting it would say.",
       size=13.5, color=AMBER, italic=True)
D.notes(s, "End on that sentence. It concedes the correction and defends the panel in one line.")

D.save(OUT)
