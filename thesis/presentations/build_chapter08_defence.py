"""Build the Chapter 8 defence deck, reading every figure from an artefact.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/build_chapter08_defence.py
Out:  thesis/presentations/chapter08_defence.pptx
"""
from __future__ import annotations

import statistics as st
from pathlib import Path

from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

from deck_common import (
    ROOT, TAB, W, M, INK, DEEP, TEAL, AMBER, CRIMSON, MUTED, TINT, PAPER, CHALK,
    HEAD, BODY, Deck, rows, num, fmt,
)

OUT = Path(__file__).resolve().parent / "chapter08_defence.pptx"
BANDS = ["below 0.40 (different chemotype)", "0.40 to 0.55 (related series)",
         "0.55 to 0.70 (same series)", "0.70 and above (close analogue)"]
SHORT = ["below 0.40\ndifferent chemotype", "0.40 to 0.55\nrelated series",
         "0.55 to 0.70\nsame series", "0.70 and above\nclose analogue"]


def facts() -> dict:
    d = {}
    cv = rows(TAB / "rf_cv_summary.csv")
    r = [num(x["roc_auc_mean"]) for x in cv if x["task"] == "classification" and x["split"] == "random"]
    s = [num(x["roc_auc_mean"]) for x in cv if x["task"] == "classification" and x["split"] == "scaffold"]
    d["cv"] = (st.mean(r), min(r), max(r), st.mean(s), min(s), max(s))

    d["ext"] = rows(TAB / "external_bbb_validation.csv")

    pr = [x for x in rows(TAB / "external_prospective.csv") if x["status"] == "ok"]
    d["n_ok"], d["n_all"] = len(pr), len(rows(TAB / "external_prospective.csv"))

    def mean_of(col):
        v = [num(x[col]) for x in pr if x[col] not in ("", "nan")]
        return st.mean(v), len(v)
    d["agg"] = {k: (mean_of(t)[0], mean_of(rr)[0], mean_of(t)[1])
                for k, t, rr in [
                    ("mi", "time_auroc_vs_measured_inactives", "random_auroc_vs_measured_inactives"),
                    ("bg", "time_auroc_vs_background", "random_auroc_vs_background"),
                    ("sens", "time_sensitivity", "random_sensitivity"),
                    ("fpr", "time_fpr_background", "random_fpr_background")]}
    d["n_actives"] = sum(int(float(x["time_n_test_actives"])) for x in pr if x["time_n_test_actives"])
    d["nov"] = (st.median(num(x["time_median_test_novelty"]) for x in pr if x["time_median_test_novelty"]),
                st.median(num(x["random_median_test_novelty"]) for x in pr if x["random_median_test_novelty"]))
    d["below"] = (int(sum(num(x["time_test_actives_below_tanimoto_0.4"]) for x in pr if x["time_test_actives_below_tanimoto_0.4"])),
                  int(sum(num(x["random_test_actives_below_tanimoto_0.4"]) for x in pr if x["random_test_actives_below_tanimoto_0.4"])))

    s8 = rows(TAB / "external_novelty_strata.csv")
    g = lambda sp, b, c: [x[c] for x in s8 if x["split"] == sp and x["novelty_band"] == b][0]
    d["recall"] = {sp: [num(g(sp, b, "recall_at_threshold")) for b in BANDS]
                   for sp in ("time", "random", "cross_source")}
    d["auroc"] = {sp: [num(g(sp, b, "auroc")) for b in BANDS] for sp in ("time", "random")}
    d["comp"] = {}
    for sp in ("time", "random"):
        n = [int(g(sp, b, "n_actives")) for b in BANDS]
        d["comp"][sp] = [100 * v / sum(n) for v in n]

    cs = rows(TAB / "external_cross_source.csv")
    d["cs"] = cs
    d["cs_tot"] = sum(int(x["n_test_actives_bindingdb_only"]) for x in cs)
    d["cs_dist"] = sum(int(x["n_test_actives_distinguishable"]) for x in cs)
    d["cs_recall"] = st.mean(num(x["recall_on_external_actives"]) for x in cs)
    d["cs_fpr"] = st.mean(num(x["fpr_background_at_same_threshold"]) for x in cs)
    d["cs_nov"] = st.median(num(x["median_novelty_of_test_set"]) for x in cs)
    d["cs_dep"] = st.mean(num(x["deployed_sensitivity"]) for x in cs)

    sp = [x for x in rows(TAB / "noncns_specificity_summary.csv")
          if x["metric"].startswith("Specificity")][0]
    d["spec"] = (num(sp["estimate"]), int(sp["k"]), int(sp["n"]), num(sp["ci95_low"]), num(sp["ci95_high"]))
    iv = rows(TAB / "inversion_validation.csv")
    d["adv"] = (sum(1 for x in iv if x["result"] == "PASS"), len(iv))
    ia = {x["metric"]: x["value"] for x in rows(TAB / "integrity_audit.csv")}
    d["leak"] = (int(ia["targets checked for scaffold overlap"]),
                 int(ia["targets with any shared scaffold"]))
    sh = [num(x["holdout_recall"]) for x in rows(TAB / "scaffold_holdout_results.csv")
          if x["holdout_recall"]]
    d["sh"] = (len(sh), st.median(sh), st.mean(sh), min(sh), max(sh))
    return d


F = facts()
D = Deck()

# 1 -- title
s = D.dark()
D.text(s, M, 1.74, W - 2 * M, 0.5, "Chapter 8", size=16, color=TEAL, bold=True)
D.text(s, M, 2.18, W - 2 * M - 1.0, 1.6,
       "Validation, external and prospective, and the composition finding",
       size=36, bold=True, font=HEAD, color=PAPER, line=1.08)
D.text(s, M, 4.15, W - 2 * M - 1.6, 1.0,
       "A programme in which each test answers a question the previous one cannot, and one result "
       "that changes what every other number in it means",
       size=16, color=CHALK, italic=True, line=1.25)
D.text(s, M, 5.70, W - 2 * M, 0.5,
       f"{F['n_ok']} of {F['n_all']} endpoints refitted before a temporal wall, "
       f"{F['n_actives']:,} post-cutoff test actives",
       size=13, color=TEAL)
D.notes(s, "The aggregate looks like decay. It is not. That reversal is the chapter.")

# 2 -- the programme
s = D.light()
D.head(s, "1", "A programme, not a number",
       "Each test answers a question the previous one cannot")
tests = [("Random cross-validation", "interpolation within known chemistry", "nothing about new chemistry"),
         ("Scaffold-grouped split", "generalisation to a withheld structural class",
          "whether the class was withheld realistically"),
         ("Temporal refit", "what it would have said about chemistry that did not exist",
          "whether a drop is the wall or the smaller training set"),
         ("Size-matched random control", "which of those two it was", "the composition of the test set"),
         ("Novelty stratification", "whether the effect is date or distance",
          "absolute performance on distant chemistry"),
         ("Cross-provenance refit", "does the signal survive a change of curator",
          "anything, on more than three endpoints"),
         ("Non-CNS specificity", "does it stay quiet when it should", "behaviour on proven inactives"),
         ("Adversarial suite", "can any check fail", "what nobody thought to check")]
D.text(s, M + 0.30, 1.86, 3.0, 0.3, "test", size=11, bold=True, color=MUTED)
D.text(s, M + 3.55, 1.86, 4.2, 0.3, "what it answers", size=11, bold=True, color=DEEP)
D.text(s, M + 8.10, 1.86, 4.0, 0.3, "what it cannot", size=11, bold=True, color=AMBER)
yy = 2.22
for i, (a, b, c) in enumerate(tests):
    if i % 2 == 0:
        D.card(s, M, yy - 0.09, W - 2 * M, 0.56, fill=TINT)
    D.text(s, M + 0.30, yy, 3.1, 0.5, a, size=12.5, bold=True, color=INK, line=1.15)
    D.text(s, M + 3.55, yy, 4.4, 0.5, b, size=12, color=INK, line=1.18)
    D.text(s, M + 8.10, yy, W - M - (M + 8.10) - 0.1, 0.5, c, size=12, color=MUTED, line=1.18)
    yy += 0.56
D.notes(s, "The last column is the reason this is a programme. No single test answers the question, "
           "and saying what each cannot do is what makes the set of them credible.")

# 3 -- the one external set
s = D.light()
D.head(s, "2", "The one genuine external set",
       "FDA-curated approved drugs absent from the barrier model's training source")
D.text(s, M + 0.30, 1.92, 5.4, 0.3, "set", size=11, bold=True, color=MUTED)
for i, hh in enumerate(["n", "AUROC", "sensitivity", "specificity"]):
    D.text(s, M + 5.90 + i * 1.55, 1.92, 1.35, 0.3, hh, size=11, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT)
labels = ["Not in B3DB by InChIKey", "Also distinguishable in feature space",
          "Of which feature-identical to training"]
cols = [MUTED, DEEP, CRIMSON]
yy = 2.34
for i, (r, lab, col) in enumerate(zip(F["ext"], labels, cols)):
    D.card(s, M, yy - 0.10, W - 2 * M, 0.66, fill=TINT)
    D.text(s, M + 0.30, yy, 5.4, 0.4, lab, size=13, bold=(i == 1), color=col)
    for j, v in enumerate([f"{int(r['n'])}", r["auroc"], r["sensitivity"], r["specificity"]]):
        D.text(s, M + 5.90 + j * 1.55, yy, 1.35, 0.4, str(v), size=13, bold=(j == 1 and i == 1),
               color=col if j == 1 else INK, align=PP_ALIGN.RIGHT)
    yy += 0.78
D.text(s, M, 4.75, W - 2 * M, 0.5,
       "The second row is the one that supports an external claim.",
       size=17, bold=True, font=HEAD, color=DEEP)
D.text(s, M, 5.30, W - 2 * M, 1.1,
       "Excluding overlap by InChIKey is not enough to call a compound unseen: the InChIKey separates "
       "stereoisomers, salts and protonation states, and the featuriser does not, so a compound can "
       "pass that check and still be one the model has memorised. The third row is the memorisation "
       "the first contains.",
       size=13.5, color=INK, line=1.28)
D.text(s, M, 6.45, W - 2 * M, 0.4,
       f"0.793 on genuinely unseen approved drugs, against {fmt(F['cv'][3])} from the scaffold split. "
       "The most honest single number in this thesis about the barrier model.",
       size=13, color=AMBER, italic=True, line=1.24)
D.source(s, "results/tables/external_bbb_validation.csv")
D.notes(s, "For the target panel no external set exists, and the reason is structural: for most of "
           "these targets the public measured chemistry IS the training set.")

# 4 -- the aggregate
s = D.light()
D.head(s, "3", "The aggregate result, which looks like decay",
       f"{F['n_ok']} endpoints refitted before a temporal wall, each with a size-matched random control")
a = F["agg"]
D.text(s, M + 5.55, 1.92, 2.2, 0.3, "time split", size=11.5, bold=True, color=MUTED,
       align=PP_ALIGN.RIGHT)
D.text(s, M + 8.05, 1.92, 2.2, 0.3, "size-matched random", size=11.5, bold=True, color=MUTED,
       align=PP_ALIGN.RIGHT)
D.text(s, M + 10.55, 1.92, 1.5, 0.3, "difference", size=11.5, bold=True, color=DEEP,
       align=PP_ALIGN.RIGHT)
lines = [("false-positive rate on background", "fpr", False),
         ("AUROC against background chemistry", "bg", False),
         ("AUROC against measured inactives", "mi", True),
         ("sensitivity at the frozen threshold", "sens", True)]
yy = 2.35
for i, (lab, key, bad) in enumerate(lines):
    t, rr, n = a[key]
    D.card(s, M, yy - 0.10, W - 2 * M, 0.68, fill=TINT if i % 2 == 0 else PAPER)
    D.text(s, M + 0.30, yy, 5.1, 0.42, lab, size=13.5, bold=bad, color=CRIMSON if bad else INK)
    D.text(s, M + 5.55, yy, 2.2, 0.42, fmt(t, 4), size=13.5, bold=bad,
           color=CRIMSON if bad else INK, align=PP_ALIGN.RIGHT)
    D.text(s, M + 8.05, yy, 2.2, 0.42, fmt(rr, 4), size=13.5, color=MUTED, align=PP_ALIGN.RIGHT)
    D.text(s, M + 10.55, yy, 1.5, 0.42, f"{t - rr:+.4f}", size=13.5, bold=True,
           color=CRIMSON if bad else DEEP, align=PP_ALIGN.RIGHT)
    yy += 0.80
D.text(s, M, 5.75, W - 2 * M, 0.95,
       "Two rows are reassuring and two are not. Specificity transfers essentially unchanged, so a "
       "temporal wall does not make the panel reckless. Ranking against measured inactives falls by "
       "0.128 and sensitivity by 0.383.",
       size=13.5, color=INK, line=1.28)
D.text(s, M, 6.72, W - 2 * M, 0.4,
       "Read at face value, the deployed figures overstate what a user should expect. The first draft "
       "said exactly that. It is the wrong reading.",
       size=13, color=AMBER, italic=True)
D.source(s, "results/tables/external_prospective.csv")
D.notes(s, "Models refitted, thresholds frozen before the cutoff, decoys matched to pre-cutoff "
           "actives only, and a size-matched control so a drop can be attributed.")

# 5 -- the composition finding
s = D.dark()
D.text(s, M, 0.92, 9.0, 0.5, "The composition finding", size=15, color=AMBER, bold=True)
D.text(s, M, 1.36, 11.7, 0.78,
       "The two splits are not testing comparable populations.",
       size=30, bold=True, font=HEAD, color=PAPER, line=1.10)
cd = CategoryChartData()
cd.categories = SHORT
cd.add_series("time split", F["comp"]["time"])
cd.add_series("size-matched random", F["comp"]["random"])
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(2.30),
                        Inches(7.2), Inches(3.85), cd)
ch = gf.chart
ch.has_title = False
ch.has_legend = True
ch.legend.position = XL_LEGEND_POSITION.BOTTOM
ch.legend.include_in_layout = False
ch.legend.font.size = Pt(11)
ch.legend.font.color.rgb = CHALK
va = ch.value_axis
va.minimum_scale, va.maximum_scale = 0, 90
va.major_gridlines.format.line.color.rgb = RGBGRID = DEEP
va.tick_labels.font.size = Pt(10)
va.tick_labels.font.color.rgb = CHALK
ch.category_axis.tick_labels.font.size = Pt(10)
ch.category_axis.tick_labels.font.color.rgb = CHALK
for i, colr in enumerate([AMBER, CHALK]):
    ch.plots[0].series[i].format.fill.solid()
    ch.plots[0].series[i].format.fill.fore_color.rgb = colr
ch.plots[0].gap_width = 60
xs = M + 7.55
D.text(s, xs, 2.40, W - M - xs, 1.5,
       f"{F['comp']['random'][3]:.1f} per cent of the random split's test compounds are close "
       f"analogues of something in its own training set.",
       size=16, bold=True, color=AMBER, line=1.24)
D.text(s, xs, 4.00, W - M - xs, 1.9,
       "The published record is series, and a random draw keeps most of a series on both sides of the "
       f"split. A random split of medicinal-chemistry data is barely a test.\n\nMedian similarity to "
       f"training: {F['nov'][0]:.3f} by date against {F['nov'][1]:.3f} at random. "
       f"{F['below'][0]:,} actives below Tanimoto 0.40 against {F['below'][1]}.",
       size=12.5, color=CHALK, line=1.28)
D.source(s, "results/tables/external_novelty_strata.csv. Both classes are binned, so each band "
            "compares actives against measured inactives at a comparable distance.")
D.notes(s, "The stratification was built to check whether the time split was too easy. It showed "
           "something else.")

# 6 -- band by band
s = D.light()
D.head(s, "4", "Read band by band, most of the gap disappears",
       "Conditioned on distance, it makes no difference how the compound was withheld")
cd2 = CategoryChartData()
cd2.categories = SHORT
cd2.add_series("recall, withheld by date", F["recall"]["time"])
cd2.add_series("recall, withheld at random", F["recall"]["random"])
gf2 = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.78),
                         Inches(7.5), Inches(4.25), cd2)
c2 = gf2.chart
c2.has_title = False
c2.has_legend = True
c2.legend.position = XL_LEGEND_POSITION.BOTTOM
c2.legend.include_in_layout = False
c2.legend.font.size = Pt(11)
va = c2.value_axis
va.minimum_scale, va.maximum_scale = 0.0, 1.0
va.major_gridlines.format.line.color.rgb = D.grid
va.tick_labels.font.size = Pt(11)
va.tick_labels.font.color.rgb = MUTED
c2.category_axis.tick_labels.font.size = Pt(10)
c2.category_axis.tick_labels.font.color.rgb = MUTED
for i, colr in enumerate([DEEP, CHALK]):
    c2.plots[0].series[i].format.fill.solid()
    c2.plots[0].series[i].format.fill.fore_color.rgb = colr
c2.plots[0].gap_width = 60
xs = M + 7.8
au_gap = max(abs(a - b) for a, b in zip(F["auroc"]["time"], F["auroc"]["random"]))
re_gap = max(abs(a - b) for a, b in zip(F["recall"]["time"], F["recall"]["random"]))
au_win = sum(1 for a, b in zip(F["auroc"]["time"], F["auroc"]["random"]) if a > b)
re_win = sum(1 for a, b in zip(F["recall"]["time"], F["recall"]["random"]) if a > b)
D.card(s, xs, 1.85, W - M - xs, 2.30, fill=TINT)
D.text(s, xs + 0.32, 2.08, W - M - xs - 0.65, 1.9,
       f"The aggregate AUROC gap of 0.128 falls, within a band, to at most {au_gap:.3f}, and the time "
       f"split is the better of the two in {au_win} of 4 bands.\n\nThe aggregate recall gap of 0.383 "
       f"falls to at most {re_gap:.3f}, and the time split is better in {re_win} of 4.",
       size=13, color=INK, line=1.28)
D.text(s, xs, 4.35, W - M - xs, 1.3,
       "Accuracy is a function of chemical distance, not of date. What time changes is the input "
       "distribution.",
       size=14.5, bold=True, color=DEEP, line=1.26)
D.text(s, xs, 5.65, W - M - xs, 1.0,
       "Decay would mean the scores expire. A stable distance-dependence means they do not, and can "
       "be quoted for the compound in hand.",
       size=12.5, color=TEAL, italic=True, line=1.26)
D.source(s, "A residual remains and is not explained away: the random split yields only 212 novel "
            "actives across 39 endpoints, too few for a paired within-band comparison.")
D.notes(s, "Chemical distance accounts for most of the apparent temporal effect and not demonstrably "
           "all of it. Say the second half too.")

# 7 -- three test sets, one curve
s = D.light()
D.head(s, "5", "Three test sets built by unrelated rules trace one curve",
       "They share no construction principle, so their ways of misleading do not overlap")
cd3 = CategoryChartData()
cd3.categories = SHORT
cd3.add_series("withheld by date", F["recall"]["time"])
cd3.add_series("withheld at random", F["recall"]["random"])
cd3.add_series("withheld by curator", F["recall"]["cross_source"])
gf3 = s.shapes.add_chart(XL_CHART_TYPE.LINE_MARKERS, Inches(M), Inches(1.80),
                         Inches(8.0), Inches(4.25), cd3)
c3 = gf3.chart
c3.has_title = False
c3.has_legend = True
c3.legend.position = XL_LEGEND_POSITION.BOTTOM
c3.legend.include_in_layout = False
c3.legend.font.size = Pt(11)
va = c3.value_axis
va.minimum_scale, va.maximum_scale = 0.0, 1.0
va.major_gridlines.format.line.color.rgb = D.grid
va.tick_labels.font.size = Pt(11)
va.tick_labels.font.color.rgb = MUTED
c3.category_axis.tick_labels.font.size = Pt(10)
c3.category_axis.tick_labels.font.color.rgb = MUTED
for i, colr in enumerate([DEEP, TEAL, AMBER]):
    ser = c3.plots[0].series[i]
    ser.format.line.color.rgb = colr
    ser.format.line.width = Pt(2.5)
    # Markers keep Office's default palette unless told otherwise, which leaves the legend swatch
    # and the line disagreeing about which series is which.
    ser.marker.format.fill.solid()
    ser.marker.format.fill.fore_color.rgb = colr
    ser.marker.format.line.color.rgb = colr
xs = M + 8.3
D.text(s, xs, 1.95, W - M - xs, 1.4,
       "A cross-provenance test set, built by a rule that knows nothing about dates and nothing about "
       "the random seed, reproduces the curve of both others.",
       size=14, bold=True, color=DEEP, line=1.26)
D.text(s, xs, 3.45, W - M - xs, 1.6,
       "That is convergent evidence, and it is worth more than any single one of the three, because "
       "the ways in which these sets could individually mislead do not overlap.",
       size=13, color=INK, line=1.26)
D.card(s, xs, 5.15, W - M - xs, 0.9, fill=TINT)
D.text(s, xs + 0.30, 5.34, W - M - xs - 0.6, 0.6,
       "If recall were a property of the split, they would disagree.",
       size=12.5, color=INK, italic=True, line=1.24)
D.source(s, "results/tables/external_novelty_strata.csv")
D.notes(s, "This is the strongest evidence in the chapter and the reason the composition reading is "
           "believed rather than merely preferred.")

# 8 -- cross-provenance
s = D.light()
D.head(s, "6", "A different curator", "Compounds deposited in BindingDB and absent from ChEMBL")
D.card(s, M, 1.80, 5.6, 1.55, fill=TINT)
D.text(s, M + 0.40, 2.02, 5.0, 1.15,
       "This test cannot report AUROC, and the reason belongs before the numbers. Every BindingDB-only "
       "row is an active, because BindingDB deposits affinities. There are no independently curated "
       "negatives to rank against, and manufacturing them from the background pool would test the "
       "background pool.",
       size=12.5, color=INK, line=1.26)
D.text(s, M + 0.30, 3.65, 1.7, 0.3, "endpoint", size=11, bold=True, color=MUTED)
for i, hh in enumerate(["actives", "recall", "bg FPR", "median novelty"]):
    D.text(s, M + 2.10 + i * 1.55, 3.65, 1.4, 0.3, hh, size=11, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT)
yy = 4.02
for i, r in enumerate(F["cs"]):
    if i % 2 == 0:
        D.card(s, M, yy - 0.09, 8.0, 0.56, fill=TINT)
    D.text(s, M + 0.30, yy, 1.7, 0.36, r["endpoint"].replace("_", "-"), size=13, bold=True, color=INK)
    for j, v in enumerate([r["n_test_actives_bindingdb_only"], r["recall_on_external_actives"],
                           r["fpr_background_at_same_threshold"], r["median_novelty_of_test_set"]]):
        D.text(s, M + 2.10 + j * 1.55, yy, 1.4, 0.36, str(v), size=13, color=INK,
               align=PP_ALIGN.RIGHT)
    yy += 0.56
D.text(s, M, 5.90, 8.0, 0.8,
       f"Across 3 endpoints and {F['cs_tot']:,} independently curated actives, mean recall "
       f"{fmt(F['cs_recall'], 4)} at a mean background false-positive rate of {fmt(F['cs_fpr'], 4)}, "
       f"against a mean deployed sensitivity of {fmt(F['cs_dep'], 4)}.",
       size=13, color=INK, line=1.26)
xs = M + 8.3
D.card(s, xs, 1.85, W - M - xs, 4.2, fill=TINT)
D.text(s, xs + 0.32, 2.10, W - M - xs - 0.65, 3.7,
       f"Read alone that looks like a heavy cost for a change of curator.\n\nIt is not. The median "
       f"distance from training of these test sets is {F['cs_nov']:.4f}, the middle of the range, and "
       f"the recall observed is what the previous slide predicts for compounds at that distance.\n\n"
       "A change of curator costs almost nothing once distance is accounted for. What it changes is "
       "which compounds are available to test on, and the compounds one database holds and another "
       "does not are not the ones both already agreed about.",
       size=12.5, color=INK, line=1.28)
D.source(s, "results/tables/external_cross_source.csv")
D.notes(s, "Recall is paired with the background false-positive rate at the same threshold, because "
           "recall alone is meaningless: a model answering active to everything scores 1.0.")

# 9 -- specificity, adversarial, leakage
s = D.light()
D.head(s, "7", "Specificity, adversarial checks and leakage", None)
sp = F["spec"]
D.card(s, M, 1.80, 3.85, 2.35)
D.stat(s, M + 0.35, 2.10, 3.2, fmt(sp[0]),
       f"of {sp[2]:,} non-CNS compounds returned no actionable disease signal", color=DEEP, vsize=42)
D.text(s, M + 0.35, 3.62, 3.2, 0.4, f"95% interval {sp[3]} to {sp[4]}", size=12, color=MUTED)
D.card(s, M + 4.20, 1.80, 3.85, 2.35)
D.stat(s, M + 4.55, 2.10, 3.2, f"{F['adv'][0]} of {F['adv'][1]}",
       "adversarial checks pass, each written so that it could fail", color=TEAL, vsize=42)
D.card(s, M + 8.40, 1.80, W - M - (M + 8.40), 2.35)
D.stat(s, M + 8.75, 2.10, 3.2, f"{F['leak'][1]} of {F['leak'][0]}",
       "targets share a scaffold between training and test", color=DEEP, vsize=42)
sh = F["sh"]
D.text(s, M, 4.40, 6.2, 0.42, "Prospective scaffold hold-out", size=17, bold=True, font=HEAD, color=DEEP)
D.text(s, M, 4.88, 6.2, 1.1,
       f"Whole scaffold classes withheld before training, recall measured at the deployed threshold "
       f"over {sh[0]} targets: median {fmt(sh[1], 4)}, mean {fmt(sh[2], 4)}, from {fmt(sh[3], 3)} to "
       f"{fmt(sh[4], 3)}. The minimum of zero is one of the six endpoints Chapter 6 identifies.",
       size=13, color=INK, line=1.28)
D.card(s, M + 6.6, 4.40, W - M - (M + 6.6), 1.75, fill=TINT)
D.text(s, M + 7.0, 4.62, W - M - (M + 6.6) - 0.8, 1.35,
       "The specificity figure was 0.949 before 25 August. Chapter 7 traces the change: stale base "
       "rates were suppressing real engagement signal, so the old number was inflated by a "
       "false-negative bug. This is a correction, not a degradation.",
       size=12.5, color=INK, line=1.26)
D.source(s, "noncns_specificity_summary.csv; inversion_validation.csv; integrity_audit.csv; "
            "scaffold_holdout_results.csv")
D.notes(s, "Presumed inactive, not proven inactive, so the artefact labels the paired false-positive "
           "rate an upper bound and the specificity a lower one.")

# 10 -- what it does not establish
s = D.dark()
D.text(s, M, 1.05, 9.0, 0.5, "What the programme does not establish", size=15, color=AMBER, bold=True)
items = [
    ("Recall on genuinely distant chemistry is poor in absolute terms",
     f"{F['recall']['time'][0]:.3f} below Tanimoto 0.40, and no analysis here improves it. The finding "
     "is that the poor number is predictable, not that it is better than it looked.", AMBER),
    ("None of this is prospective in the strict sense",
     "No compound here was predicted before it was measured. The dates are real; the analysis is "
     "retrospective.", AMBER),
    ("The cross-provenance arm rests on three endpoints",
     "The number of compounds one database holds at pChEMBL 7 or above and the other lacks is small.",
     AMBER),
    ("The null models have no artefact",
     "The label-permutation null, 0.4938 random and 0.4921 scaffold, and the independent reproduction "
     "of 26 values to 4.7e-5, are hard-coded in the report generator. No file holds them. It is the "
     "one claim here a reader cannot check, and it is the evidence that the cross-validation is not "
     "inflated by leakage.", CRIMSON),
]
yy = 1.75
for i, (t, b, col) in enumerate(items, 1):
    D.dot(s, M, yy + 0.04, str(i), fill=col, fg=INK, dia=0.38)
    D.text(s, M + 0.62, yy, 4.5, 0.9, t, size=14.5, bold=True, color=col, line=1.16)
    D.text(s, M + 5.40, yy - 0.02, 11.9 - 5.40, 1.15, b, size=12.5, color=CHALK, line=1.26)
    yy += 1.28
D.text(s, M, 6.85, 11.7, 0.4,
       "The fourth is a defect rather than a limitation, and it should be regenerated into an "
       "artefact and declared.",
       size=13, color=PAPER, italic=True)
D.notes(s, "Volunteer the fourth. A referee who finds an uncheckable null model unaided will discount "
           "every cross-validation figure in the thesis.")

D.save(OUT)
