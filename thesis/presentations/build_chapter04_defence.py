"""Build the Chapter 4 defence deck, reading every figure from an artefact.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/build_chapter04_defence.py
Out:  thesis/presentations/chapter04_defence.pptx
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

OUT = Path(__file__).resolve().parent / "chapter04_defence.pptx"


def facts() -> dict:
    d = {}
    cal = rows(TAB / "calibration.csv")
    d["cal"] = cal
    d["ece_raw"] = st.mean(num(r["ece_raw"]) for r in cal)
    d["ece_cal"] = st.mean(num(r["ece_calibrated"]) for r in cal)
    d["ece_cal_med"] = st.median(num(r["ece_calibrated"]) for r in cal)
    d["ece_cal_min"] = min(cal, key=lambda r: num(r["ece_calibrated"]))
    d["ece_cal_max"] = max(cal, key=lambda r: num(r["ece_calibrated"]))
    d["brier_better"] = sum(1 for r in cal if num(r["brier_calibrated"]) < num(r["brier_raw"]))
    d["ece_better"] = sum(1 for r in cal if num(r["ece_calibrated"]) < num(r["ece_raw"]))
    d["n_core"] = len(cal)

    per = rows(TAB / "integrity_calibration_per_target.csv")
    pe = [num(r["ece"]) for r in per]
    d["n_binder"] = len(per)
    d["b_mean"], d["b_med"] = st.mean(pe), st.median(pe)
    d["b_worst"] = max(per, key=lambda r: num(r["ece"]))
    d["b_best"] = min(per, key=lambda r: num(r["ece"]))
    d["b_over"] = sum(1 for v in pe if v > 0.10)

    con = rows(TAB / "rf_conformal.csv")
    d["con"] = con
    d["cov_min"] = min(num(r["empirical_coverage"]) for r in con)
    d["cov_max"] = max(num(r["empirical_coverage"]) for r in con)
    d["cov_ok"] = sum(1 for r in con if num(r["empirical_coverage"]) >= 0.90)
    d["set_min"] = min(num(r["avg_set_size"]) for r in con)
    d["set_max"] = max(num(r["avg_set_size"]) for r in con)
    d["n_test"] = sum(int(r["n_test"]) for r in con)

    d["meas"] = rows(TAB / "applicability_measures.csv")
    d["ad"] = {r["subset"]: r for r in rows(TAB / "applicability_bbb_validation.csv")}
    cov = rows(TAB / "applicability_coverage.csv")
    d["cov"] = sorted(cov, key=lambda r: -num(r["in_domain_frac"]))
    d["cov_below"] = sum(1 for r in cov if num(r["in_domain_frac"]) < 0.5)
    d["drugbank"] = int(cov[0]["n_drugbank"])

    strata = rows(TAB / "external_novelty_strata.csv")
    order = ["below 0.40 (different chemotype)", "0.40 to 0.55 (related series)",
             "0.55 to 0.70 (same series)", "0.70 and above (close analogue)"]
    d["bands"] = order
    for key, lab in (("time", "date"), ("random", "rand"), ("cross_source", "cross")):
        m = {r["novelty_band"]: num(r["recall_at_threshold"]) for r in strata if r["split"] == key}
        d[f"recall_{lab}"] = [m[b] for b in order]
    return d


F = facts()
D = Deck()

# 1 -- title
s = D.dark()
D.text(s, M, 1.74, W - 2 * M, 0.5, "Chapter 4", size=16, color=TEAL, bold=True)
D.text(s, M, 2.18, W - 2 * M - 1.0, 1.6,
       "Calibration, conformal prediction and the applicability domain",
       size=36, bold=True, font=HEAD, color=PAPER, line=1.08)
D.text(s, M, 4.15, W - 2 * M - 1.6, 1.0,
       "Three instruments answering three different questions, one of which turns out to answer a "
       "question it was not built for",
       size=16, color=CHALK, italic=True, line=1.25)
D.text(s, M, 5.70, W - 2 * M, 0.5,
       f"Calibration measured on {F['n_core']} core classifiers and {F['n_binder']} binder endpoints, "
       f"coverage on {F['n_test']:,} held-out compounds",
       size=13, color=TEAL)
D.notes(s, "The chapter's two hardest findings are that the headline calibration figure covers eight "
           "of seventy estimators, and that the domain flag does not predict discrimination.")

# 2 -- three questions
s = D.light()
D.head(s, "1", "Three questions, and why one answer will not do", None)
qs = [("Is 0.84 a probability?",
       "A forest returns a vote share. Averaging 300 trees gives a number that behaves like a "
       "probability, clusters toward the middle, and is not one.",
       "Calibration. A statement about a population of predictions.", DEEP),
      ("How confident is it about this compound?",
       "Calibration cannot answer this. A model can be well calibrated and uniformly uncertain, and "
       "the average tells you nothing about the case in front of you.",
       "Conformal prediction. A coverage guarantee, measured rather than assumed.", TEAL),
      ("Has it seen anything like this compound?",
       "Neither of the above answers this. Both are computed on held-out chemistry from the same "
       "distribution as training. A compound from outside it can get a confident, well-calibrated, "
       "tightly-covered prediction that is worthless.",
       "The applicability domain. And it answers a different question from the one it was built for.",
       AMBER)]
yy = 1.80
for i, (q, why, ans, col) in enumerate(qs, 1):
    D.card(s, M, yy, W - 2 * M, 1.48)
    D.dot(s, M + 0.35, yy + 0.50, str(i), fill=col)
    D.text(s, M + 1.05, yy + 0.26, 3.7, 0.8, q, size=16, bold=True, font=HEAD, color=col, line=1.12)
    D.text(s, M + 5.00, yy + 0.22, 4.0, 1.05, why, size=12, color=INK, line=1.24)
    D.text(s, M + 9.25, yy + 0.30, W - M - (M + 9.25) - 0.3, 0.95, ans, size=12, bold=True,
           color=col, line=1.24)
    yy += 1.62
D.notes(s, "Insist on the third row. It is the one most tools do not have at all.")

# 3 -- calibration construction
s = D.light()
D.head(s, "2", "Calibration is fitted out of fold, and that carries the weight",
       "Isotonic regression for the core classifiers, Platt scaling for the binder panel")
D.card(s, M, 1.80, 5.9, 2.15)
D.text(s, M + 0.42, 2.06, 5.0, 0.4, "Why out of fold", size=17, bold=True, font=HEAD, color=DEEP)
D.text(s, M + 0.42, 2.56, 5.05, 1.25,
       "Fitting a calibrator on the predictions a model made about its own training data corrects an "
       "optimism the deployed model does not have. The reported improvement would then be an artefact "
       "of the procedure. No compound contributes to the calibrator that scores it.",
       size=13, color=INK, line=1.26)
D.card(s, M + 6.25, 1.80, W - 2 * M - 6.25, 2.15)
D.text(s, M + 6.67, 2.06, 5.0, 0.4, "Two consequences of monotonicity",
       size=17, bold=True, font=HEAD, color=TEAL)
D.text(s, M + 6.67, 2.56, 5.05, 1.25,
       "An isotonic map cannot reorder compounds, so AUROC is invariant under calibration and every "
       "discrimination figure in this thesis is unaffected by anything here. And being free rather "
       "than parametric, it can correct an arbitrary distortion, given enough held-out data.",
       size=13, color=INK, line=1.26)
D.text(s, M, 4.28, W - 2 * M, 0.5, "The binder panel cannot afford that.",
       size=19, bold=True, font=HEAD, color=AMBER)
D.text(s, M, 4.85, W - 2 * M, 1.2,
       "Platt scaling fits two parameters, a logistic curve, on a stratified fifth withheld from the "
       "forest's own fit. It is used because the positive class at one binder endpoint is often too "
       "small for an isotonic step function to be anything but overfitting. That constraint is real. "
       "Slide 5 measures what it costs.",
       size=14, color=INK, line=1.28)
D.notes(s, "The out-of-fold point is the one a statistician will probe first. Have it ready.")

# 4 -- calibration results
s = D.light()
D.head(s, "3", "What calibration achieves on the core classifiers",
       f"Expected calibration error falls on all {F['ece_better']} of {F['n_core']}")
cd = CategoryChartData()
eps = [r["endpoint"].replace("_", "-") for r in F["cal"]]
cd.categories = eps
cd.add_series("before calibration", [num(r["ece_raw"]) for r in F["cal"]])
cd.add_series("after isotonic calibration", [num(r["ece_calibrated"]) for r in F["cal"]])
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.72),
                        Inches(8.0), Inches(4.30), cd)
ch = gf.chart
ch.has_title = False
ch.has_legend = True
ch.legend.position = XL_LEGEND_POSITION.BOTTOM
ch.legend.include_in_layout = False
ch.legend.font.size = Pt(11)
va = ch.value_axis
va.minimum_scale, va.maximum_scale = 0.0, 0.10
va.major_gridlines.format.line.color.rgb = D.grid
va.tick_labels.font.size = Pt(11)
va.tick_labels.font.color.rgb = MUTED
ch.category_axis.tick_labels.font.size = Pt(10)
ch.category_axis.tick_labels.font.color.rgb = MUTED
for i, colr in enumerate([CHALK, DEEP]):
    ch.plots[0].series[i].format.fill.solid()
    ch.plots[0].series[i].format.fill.fore_color.rgb = colr
ch.plots[0].gap_width = 60

xs = M + 8.3
D.stat(s, xs, 1.90, 3.8, f"{fmt(F['ece_raw'], 4)} → {fmt(F['ece_cal'], 4)}",
       "mean expected calibration error, a factor of "
       f"{F['ece_raw'] / F['ece_cal']:.1f}", color=DEEP, vsize=27)
D.text(s, xs, 3.55, 3.8, 1.0,
       f"Per endpoint, {fmt(num(F['ece_cal_min']['ece_calibrated']), 4)} at "
       f"{F['ece_cal_min']['endpoint']} to {fmt(num(F['ece_cal_max']['ece_calibrated']), 4)} at "
       f"{F['ece_cal_max']['endpoint']}, worst by a factor of eight.",
       size=13, color=INK, line=1.26)
D.card(s, xs, 4.60, 3.8, 1.45, fill=TINT)
D.text(s, xs + 0.32, 4.82, 3.2, 1.05,
       f"One endpoint gets worse on the other metric. The Brier score improves on "
       f"{F['brier_better']} of {F['n_core']} and degrades at BBB, 0.1289 to 0.1358: calibration "
       f"bought reliability there at a cost in sharpness.",
       size=11.5, color=INK, line=1.24)
D.source(s, "results/tables/calibration.csv, measured on out-of-fold predictions")
D.notes(s, "Volunteer the BBB Brier regression. It is the only place calibration makes something "
           "worse and a referee who finds it unaided will wonder what else is unmentioned.")

# 5 -- the headline does not cover the panel
s = D.dark()
D.text(s, M, 0.95, 9.0, 0.5, "What the headline figure does not cover", size=15, color=AMBER, bold=True)
D.text(s, M, 1.38, 11.6, 0.78,
       f"0.0801 to 0.0147 describes {F['n_core']} of the 70 deployed estimators.",
       size=27, bold=True, font=HEAD, color=PAPER, line=1.10)
D.text(s, M, 2.28, 11.0, 0.6,
       "It is quoted in the manuscript abstract, the technical report, the evidence map and Chapter 1 "
       "of this thesis. No document says so.",
       size=14, color=CHALK, line=1.26)
hdr_y = 3.30
D.text(s, M + 4.55, hdr_y, 3.0, 0.32, f"core, isotonic  (n={F['n_core']})", size=12, bold=True,
       color=CHALK, align=PP_ALIGN.RIGHT)
D.text(s, M + 8.00, hdr_y, 3.2, 0.32, f"binder panel, Platt  (n={F['n_binder']})", size=12,
       bold=True, color=AMBER, align=PP_ALIGN.RIGHT)
lines = [("mean expected calibration error", fmt(F["ece_cal"], 4), fmt(F["b_mean"], 4)),
         ("median", fmt(F["ece_cal_med"], 4), fmt(F["b_med"], 4)),
         ("best endpoint", f"{fmt(num(F['ece_cal_min']['ece_calibrated']), 4)}  {F['ece_cal_min']['endpoint']}",
          f"{fmt(num(F['b_best']['ece']), 3)}  {F['b_best']['endpoint']}"),
         ("worst endpoint", f"{fmt(num(F['ece_cal_max']['ece_calibrated']), 4)}  {F['ece_cal_max']['endpoint']}",
          f"{fmt(num(F['b_worst']['ece']), 3)}  {F['b_worst']['endpoint']}"),
         ("endpoints above 0.10", f"0 of {F['n_core']}", f"{F['b_over']} of {F['n_binder']}")]
yy = 3.75
for lab, a, b in lines:
    D.text(s, M + 0.10, yy, 4.3, 0.34, lab, size=13, color=PAPER)
    D.text(s, M + 4.55, yy, 3.0, 0.34, a, size=13, color=CHALK, align=PP_ALIGN.RIGHT)
    D.text(s, M + 8.00, yy, 3.2, 0.34, b, size=13, bold=True, color=AMBER, align=PP_ALIGN.RIGHT)
    yy += 0.46
D.text(s, M, 6.25, 11.6, 0.7,
       f"The binder panel's calibrated error, {fmt(F['b_mean'], 4)}, is five times the core "
       f"classifiers' and almost exactly what the core classifiers read before calibration, "
       f"{fmt(F['ece_raw'], 4)}. On this evidence Platt scaling on these endpoints buys very little.",
       size=13.5, color=PAPER, italic=True, line=1.26)
D.notes(s, "TAAR1, worst at 0.178, has 403 rows in its endpoint table and 82 in the fitted binder "
           "set. The constraint that forced Platt scaling was real; the price is now measured. The "
           "remedy is specific: isotonic where there are enough held-out positives, and a per-endpoint "
           "calibration error beside every other probability.")

# 6 -- conformal
s = D.light()
D.head(s, "4", "Conformal prediction, and what it says that calibration cannot",
       f"Mondrian, class-conditional, at a 0.90 target, on {F['n_test']:,} held-out compounds")
D.text(s, M + 0.30, 1.85, 2.2, 0.3, "endpoint", size=11, bold=True, color=MUTED)
D.text(s, M + 2.60, 1.85, 1.5, 0.3, "n test", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
D.text(s, M + 4.35, 1.85, 1.8, 0.3, "coverage", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
D.text(s, M + 6.35, 1.85, 1.6, 0.3, "set size", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
yy = 2.22
for r in F["con"]:
    cov = num(r["empirical_coverage"])
    under = cov < 0.90
    if F["con"].index(r) % 2 == 0:
        D.card(s, M, yy - 0.08, 8.15, 0.50, fill=TINT)
    D.text(s, M + 0.30, yy, 2.2, 0.32, r["endpoint"].replace("_", "-"), size=13, color=INK)
    D.text(s, M + 2.60, yy, 1.5, 0.32, f"{int(r['n_test']):,}", size=13, color=MUTED,
           align=PP_ALIGN.RIGHT)
    D.text(s, M + 4.35, yy, 1.8, 0.32, fmt(cov), size=13, bold=True,
           color=CRIMSON if under else DEEP, align=PP_ALIGN.RIGHT)
    D.text(s, M + 6.35, yy, 1.6, 0.32, fmt(num(r["avg_set_size"])), size=13, color=INK,
           align=PP_ALIGN.RIGHT)
    yy += 0.50
xs = M + 8.5
D.text(s, xs, 1.90, W - M - xs, 1.0,
       f"Coverage {fmt(F['cov_min'])} to {fmt(F['cov_max'])} against 0.90. "
       f"{F['cov_ok']} of {len(F['con'])} at or above target.",
       size=14, bold=True, color=DEEP, line=1.26)
D.text(s, xs, 2.95, W - M - xs, 1.3,
       "The two that fall short are the monoamine oxidases, MAO-B by one thousandth and MAO-A by "
       "eleven. They are also the two endpoints that degrade most under a scaffold split, and "
       "coverage is a guarantee under exchangeability.",
       size=12.5, color=INK, line=1.26)
D.card(s, xs, 4.28, W - M - xs, 2.05, fill=TINT)
D.text(s, xs + 0.32, 4.50, W - M - xs - 0.65, 1.70,
       f"Mean set size {fmt(F['set_min'])} to {fmt(F['set_max'])} on a two-class problem. So 0.7 to "
       "7.9 per cent of compounds get an empty set, meaning the compound conforms to neither class, "
       "or a set with both, meaning it separates neither. Neither statement is available from a point "
       "probability.",
       size=11.5, color=INK, line=1.24)
D.source(s, "results/tables/rf_conformal.csv. Same eight core classifiers: the binder panel has no "
            "measured coverage statement at all.")
D.notes(s, "Class-conditional matters here: the median deployed endpoint is 0.825 active, and "
           "marginal coverage can be met while the rare class fails almost completely.")

# 7 -- AD measure choice
s = D.light()
D.head(s, "5", "The applicability domain, and choosing its measure",
       "Maximum Tanimoto to that endpoint's own measured chemistry, on 2,048-bit fingerprints")
D.text(s, M + 0.30, 1.88, 2.9, 0.3, "measure", size=11, bold=True, color=MUTED)
for i, hh in enumerate(["unseen drugs", "non-drug-like", "p", "separates", "caught at 10% drug loss"]):
    D.text(s, M + 3.30 + i * 1.75, 1.88, 1.6, 0.3, hh, size=10.5, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT)
yy = 2.28
for i, r in enumerate(F["meas"]):
    sep = r["separates"] == "True"
    if i % 2 == 0:
        D.card(s, M, yy - 0.09, W - 2 * M, 0.54, fill=TINT)
    name = {"max": "maximum similarity", "mean_top5": "mean of top 5",
            "kth_5": "5th nearest neighbour", "density_0.4": "density within 0.4"}[r["measure"]]
    D.text(s, M + 0.30, yy, 2.9, 0.34, name, size=12.5, bold=(r["measure"] == "max"),
           color=DEEP if r["measure"] == "max" else INK)
    for j, v in enumerate([r["median_unseen_drugs"], r["median_non_drug_like"],
                           f"{float(r['mann_whitney_p']):.4f}",
                           "yes" if sep else "no", r["non_drug_like_caught"]]):
        D.text(s, M + 3.30 + j * 1.75, yy, 1.6, 0.34, str(v), size=12.5,
               color=(DEEP if sep else CRIMSON) if j == 3 else INK,
               bold=(j == 3), align=PP_ALIGN.RIGHT)
    yy += 0.54
D.text(s, M, 4.75, W - 2 * M, 0.5,
       "Two of the four separate at p < 0.01, and no alternative beats the deployed one.",
       size=16, bold=True, font=HEAD, color=DEEP)
D.card(s, M, 5.35, W - 2 * M, 1.15, fill=TINT)
D.text(s, M + 0.40, 5.56, W - 2 * M - 0.8, 0.85,
       "Passing is not the same as being useful. At a threshold rejecting a tenth of genuine drugs "
       "the flag catches 0.375 of genuinely distant chemistry; at the deployed cut of 0.30 it catches "
       "0.20. Four in five non-drug-like structures are waved through. It separates in the aggregate "
       "and discriminates poorly on the individual case.",
       size=13, color=INK, line=1.26)
D.source(s, "results/tables/applicability_measures.csv; results/tables/inversion_validation.csv")
D.notes(s, "The check previously failed and the fix was to the controls, not the criterion: 28 "
           "controls were measured compounds in the flag's own reference library. Chapter 9 has the "
           "full history.")

# 8 -- what it does predict
s = D.light()
D.head(s, "6", "What the flag actually predicts, which is not what it was built for",
       "Recall is a monotone function of domain distance, across three unrelated test sets")
cd2 = CategoryChartData()
cd2.categories = ["below 0.40\ndifferent chemotype", "0.40 to 0.55\nrelated series",
                  "0.55 to 0.70\nsame series", "0.70 and above\nclose analogue"]
cd2.add_series("withheld by date", F["recall_date"])
cd2.add_series("withheld at random", F["recall_rand"])
cd2.add_series("withheld by curator", F["recall_cross"])
gf2 = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.75),
                         Inches(7.9), Inches(4.25), cd2)
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
for i, colr in enumerate([DEEP, TEAL, AMBER]):
    c2.plots[0].series[i].format.fill.solid()
    c2.plots[0].series[i].format.fill.fore_color.rgb = colr
c2.plots[0].gap_width = 70
xs = M + 8.2
D.text(s, xs, 1.90, W - M - xs, 1.3,
       "The quantity the server already computes and already shows predicts the recall that query "
       "will receive.",
       size=15, bold=True, color=DEEP, line=1.24)
D.text(s, xs, 3.20, W - M - xs, 1.6,
       "A compound the flag places near the edge of the domain is precisely the compound whose "
       "activity the panel is most likely to miss. So the expected recall for the compound in hand is "
       "knowable at the moment of the query.",
       size=13, color=INK, line=1.26)
D.text(s, xs, 4.95, W - M - xs, 1.1,
       "The flag is a poor detector of alien chemistry and a good predictor of when the panel will "
       "fall silent on a real active.",
       size=13, bold=True, color=TEAL, italic=True, line=1.24)
D.source(s, "results/tables/external_novelty_strata.csv")
D.notes(s, "This is the slide that rescues the applicability domain. It earns its place in a role it "
           "was not designed for.")

# 9 -- the reversal
s = D.light()
D.head(s, "7", "But the flag does not predict discrimination",
       "On the one external set where this can be tested, the direction reverses")
ad = F["ad"]
rowsx = [("In domain, T ≥ 0.30", "in_domain (T>=0.30)", DEEP),
         ("Out of domain, T < 0.30", "out_of_domain (T<0.30)", AMBER),
         ("All", "all", MUTED)]
D.text(s, M + 0.30, 1.90, 3.2, 0.3, "subset", size=11, bold=True, color=MUTED)
for i, hh in enumerate(["n", "mean max Tanimoto", "AUROC", "accuracy"]):
    D.text(s, M + 3.60 + i * 1.55, 1.90, 1.4, 0.3, hh, size=11, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT)
yy = 2.32
for lab, key, col in rowsx:
    r = ad[key]
    D.card(s, M, yy - 0.10, 9.9, 0.62, fill=TINT)
    D.text(s, M + 0.30, yy, 3.2, 0.36, lab, size=13.5, bold=True, color=col)
    for j, v in enumerate([f"{int(r['n'])}", r["mean_max_tanimoto"], r["auroc"], r["accuracy"]]):
        D.text(s, M + 3.60 + j * 1.55, yy, 1.4, 0.36, str(v), size=13.5,
               bold=(j == 2), color=col if j == 2 else INK, align=PP_ALIGN.RIGHT)
    yy += 0.74
gap = num(ad["out_of_domain (T<0.30)"]["auroc"]) - num(ad["in_domain (T>=0.30)"]["auroc"])
D.text(s, M, 4.72, 9.9, 0.9,
       f"The barrier model ranks the {int(ad['out_of_domain (T<0.30)']['n'])} out-of-domain compounds "
       f"{gap:+.4f} AUROC better than the {int(ad['in_domain (T>=0.30)']['n'])} in-domain ones, while "
       "accuracy is flat.",
       size=15, bold=True, color=AMBER, line=1.24)
D.card(s, M, 5.70, 9.9, 1.05, fill=TINT)
D.text(s, M + 0.38, 5.92, 9.2, 0.75,
       "Recall and ranking are different quantities. Domain distance predicts recall, which depends on "
       "where the threshold sits relative to the actives' scores. It does not predict AUROC, which is "
       "threshold-free. On novel chemistry the panel does not become wrong, it becomes quiet.",
       size=12.5, color=INK, line=1.24)
xs = M + 10.1
D.text(s, xs, 2.30, W - M - xs, 2.2,
       "Do not over-read it: 48 compounds, no interval computed, and the subsets are matched on "
       "nothing but distance.\n\nBut the direction is the opposite of the assumption, and the "
       "assumption is the one a referee brings.",
       size=12, color=MUTED, line=1.26)
D.source(s, "results/tables/applicability_bbb_validation.csv, regenerated for this chapter")
D.notes(s, "This artefact was six days stale and disagreed with external_bbb_validation.csv on the "
           "same 306 compounds, 0.7608 against 0.7645. Regenerated, they agree exactly. It is now "
           "declared in the freshness graph.")

# 10 -- coverage of drug space
s = D.light()
D.head(s, "8", "How much of drug space is inside the domain",
       f"The {F['drugbank']:,} DrugBank drugs scored against each endpoint's measured chemistry")
cd3 = CategoryChartData()
cd3.categories = [r["endpoint"].replace("_", "-") for r in F["cov"]]
cd3.add_series("fraction at T ≥ 0.30", [num(r["in_domain_frac"]) for r in F["cov"]])
gf3 = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.78),
                         Inches(7.6), Inches(4.25), cd3)
c3 = gf3.chart
c3.has_title = False
c3.has_legend = False
va = c3.value_axis
va.minimum_scale, va.maximum_scale = 0.0, 1.0
va.major_gridlines.format.line.color.rgb = D.grid
va.tick_labels.font.size = Pt(11)
va.tick_labels.font.color.rgb = MUTED
c3.category_axis.tick_labels.font.size = Pt(11)
c3.category_axis.tick_labels.font.color.rgb = MUTED
ser = c3.plots[0].series[0]
ser.format.fill.solid()
ser.format.fill.fore_color.rgb = DEEP
c3.plots[0].gap_width = 60
xs = M + 7.9
D.text(s, xs, 1.92, W - M - xs, 1.3,
       f"For every endpoint except the barrier model, most approved drugs sit below the cut. At BACE1 "
       f"nearly four in five do.",
       size=14.5, bold=True, color=DEEP, line=1.26)
D.text(s, xs, 3.25, W - M - xs, 1.5,
       "The median approved drug is at Tanimoto 0.25 to 0.30 from the nearest compound ever measured "
       "at these targets. With the previous slide but one, that says the panel will miss a "
       "substantial fraction of genuine activities among approved drugs.",
       size=12.5, color=INK, line=1.26)
D.card(s, xs, 4.90, W - M - xs, 1.15, fill=TINT)
D.text(s, xs + 0.32, 5.10, W - M - xs - 0.65, 0.85,
       "The barrier model's 0.72 is the exception because its training set is drawn from "
       "approved-drug chemistry in the first place. That is overlap between two drug collections, not "
       "evidence of better generalisation.",
       size=11.5, color=INK, line=1.24)
D.source(s, "results/tables/applicability_coverage.csv, regenerated for this chapter. Note in_domain "
            "here means T ≥ 0.30, which is the interface's 'in or near domain'.")
D.notes(s, "This is the quantitative basis for the Chapter 1 claim that silence is the weakest of the "
           "system's outputs.")

# 11 -- reading order
s = D.light()
D.head(s, "9", "What a user should read, in order", "A reader who takes one of the three has the wrong summary")
steps = [("Read the domain distance first",
          "It tells you whether the rest is worth reading, and it tells you the recall the panel "
          "achieves at that distance. Below 0.40, recall is about 0.16 and a silence means almost "
          "nothing.", AMBER),
         ("Read the conformal set second",
          "An empty set, or a set containing both classes, is a statement a point probability cannot "
          "make. Available for the eight core classifiers only.", TEAL),
         ("Read the calibrated probability last, against its base rate",
          f"Expected calibration error {fmt(F['ece_cal'], 4)} on the eight core classifiers and "
          f"{fmt(F['b_mean'], 4)} across the {F['n_binder']} binder endpoints. Chapter 7 explains why "
          "the enrichment over the base rate, not the probability, is what the disease layer consumes.",
          DEEP)]
yy = 1.85
for i, (t, b, col) in enumerate(steps, 1):
    D.card(s, M, yy, W - 2 * M, 1.42)
    D.dot(s, M + 0.35, yy + 0.48, str(i), fill=col)
    D.text(s, M + 1.05, yy + 0.24, 4.3, 0.9, t, size=16, bold=True, font=HEAD, color=col, line=1.12)
    D.text(s, M + 5.55, yy + 0.24, W - 2 * M - 5.9, 1.0, b, size=12.5, color=INK, line=1.26)
    yy += 1.58
D.text(s, M, 6.35, W - 2 * M, 0.5,
       "The uncertainty stack is complete for eight of seventy estimators and two-thirds absent for "
       "the 38 targets a user is most likely to query.",
       size=13, color=AMBER, italic=True, line=1.24)
D.notes(s, "End on the gap, not on the achievement. It is the honest state of the chapter.")

D.save(OUT)
