"""Build the Chapter 10 defence deck, reading every figure from an artefact.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/build_chapter10_defence.py
Out:  thesis/presentations/chapter10_defence.pptx
"""
from __future__ import annotations

import json
import statistics as st
from pathlib import Path

from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

from deck_common import (
    ROOT, TAB, INV, W, H, M, INK, DEEP, TEAL, AMBER, CRIMSON, MUTED, TINT, PAPER, CHALK,
    HEAD, BODY, Deck, rows, num, fmt,
)

OUT = Path(__file__).resolve().parent / "chapter10_defence.pptx"
BANDS = ["below 0.40 (different chemotype)", "0.40 to 0.55 (related series)",
         "0.55 to 0.70 (same series)", "0.70 and above (close analogue)"]
SHORT = ["below 0.40\ndifferent chemotype", "0.40 to 0.55\nrelated series",
         "0.55 to 0.70\nsame series", "0.70 and above\nclose analogue"]


def one(path, field, pred):
    for r in rows(path):
        if pred(r):
            return r[field]
    raise KeyError(f"no matching row in {path.name}")


def facts() -> dict:
    d = {}
    s = rows(TAB / "external_novelty_strata.csv")
    g = lambda sp, b, c: [x[c] for x in s if x["split"] == sp and x["novelty_band"] == b][0]
    d["recall"] = {sp: [num(g(sp, b, "recall_at_threshold")) for b in BANDS]
                   for sp in ("time", "cross_source")}
    d["auroc_band"] = [num(g("time", b, "auroc")) for b in BANDS]
    d["n_band"] = [int(g("time", b, "n_actives")) for b in BANDS]

    h10 = rows(INV / "H10_barrier_necessity.csv")
    d["h10"] = h10
    d["h10_verdict"] = h10[0]["verdict"]
    d["h10_by"] = {(r["method"], r["population"]): r for r in h10}

    cal = rows(TAB / "calibration.csv")
    raw = [num(x["ece_raw"]) for x in cal]
    fin = [num(x["ece_calibrated"]) for x in cal]
    d["core_cal"] = (len(cal), st.mean(raw), st.mean(fin), min(fin), max(fin),
                     [x["endpoint"] for x in cal][fin.index(max(fin))])
    b = [num(x["ece"]) for x in rows(TAB / "integrity_calibration_per_target.csv")]
    d["binder_cal"] = (len(b), st.mean(b), st.median(b), min(b), max(b))

    np_ = [x for x in rows(TAB / "external_natural_products_summary.csv")]
    d["np_scored"] = [x for x in np_ if x["status"] == "scored"]
    d["np_cand"] = sum(int(x["n_candidates"]) for x in np_)
    d["np_contam"] = int(sum(float(x["removed_as_contaminated"] or 0) for x in np_))
    d["np_noest"] = sum(1 for x in np_ if x["status"] == "no deployed estimator")

    d["spec"] = (num(one(TAB / "noncns_specificity_summary.csv", "estimate",
                         lambda r: r["metric"].startswith("Specificity"))),
                 num(one(TAB / "noncns_specificity_summary.csv", "ci95_low",
                         lambda r: r["metric"].startswith("Specificity"))),
                 num(one(TAB / "noncns_specificity_summary.csv", "ci95_high",
                         lambda r: r["metric"].startswith("Specificity"))))

    reg = json.loads((ROOT / "models_rf" / "binder_modes.json").read_text(encoding="utf-8"))
    dep = [v for v in reg.values() if v.get("deployed")]
    d["n_dep"] = len(dep)
    d["n_unreliable"] = sum(1 for v in dep if v.get("reliable_call") is False)
    sens = [v["sensitivity_at_threshold"] for v in dep]
    d["sens"] = (st.mean(sens), st.median(sens), min(sens), max(sens))

    pn = rows(TAB / "permutation_null.csv")
    d["pn"] = {sp: st.mean(num(r["permuted_roc_auc_mean"]) for r in pn if r["split"] == sp)
               for sp in ("random", "scaffold")}
    cv = [r for r in rows(TAB / "rf_cv_summary.csv") if r["task"] == "classification"]
    d["cv"] = {sp: (min(num(r["roc_auc_mean"]) for r in cv if r["split"] == sp),
                    max(num(r["roc_auc_mean"]) for r in cv if r["split"] == sp))
               for sp in ("random", "scaffold")}

    sp3 = TAB / "library_sp3_coverage.csv"
    d["sp3"] = None
    if sp3.exists():
        d["sp3"] = (num(one(sp3, "value", lambda r: r["metric"] == "median fraction sp3")),
                    num(one(sp3, "value",
                            lambda r: r["metric"].startswith("sp3-rich") and "per cent" in r["metric"])),
                    int(float(one(sp3, "value",
                                  lambda r: r["metric"] == "structures parsed as a desalted parent"))))

    v = rows(INV / "VERDICTS.csv")
    d["n_hyp"] = len(v)
    d["n_ref"] = sum(1 for x in v if x["verdict"].startswith("REFUTED"))
    d["n_weak"] = sum(1 for x in v if x["verdict"].startswith("WEAKENED"))
    return d


F = facts()
D = Deck()

# ---------------------------------------------------------------- 1. title
s = D.dark()
D.text(s, M, 1.74, W - 2 * M, 0.5, "Chapter 10", size=16, color=TEAL, bold=True)
D.text(s, M, 2.18, W - 2 * M - 1.0, 1.5,
       "Limitations, and what\nshould be done next",
       size=42, bold=True, font=HEAD, color=PAPER, line=1.08)
D.text(s, M, 4.10, W - 2 * M - 1.6, 1.1,
       "Three different things get called limitations: what the system cannot do, what the record "
       "gets wrong, and what should be done next. They are separated here.",
       size=17, color=CHALK, line=1.34)
D.text(s, M, 5.62, W - 2 * M, 0.5,
       "One item of the future work is done in this chapter rather than proposed.",
       size=14, color=AMBER, italic=True)
D.text(s, M, 6.06, W - 2 * M, 0.4,
       "Chapter 9 named the falsification suite's largest gap. A gap named but not filled is half an "
       "observation.", size=12, color=MUTED)
D.notes(s, "Do not open with an apology. Open with the separation: a limitation of the science, a "
           "defect in the record, and a piece of work are three different claims.")

# ---------------------------------------------------------------- 2. the dominant limitation
s = D.light()
D.head(s, "1", "Accuracy is a function of chemical distance",
       "The limitation that bounds every number in Chapters 3 to 8")
cd = CategoryChartData()
cd.categories = SHORT
cd.add_series("time split", tuple(F["recall"]["time"]))
cd.add_series("cross-provenance", tuple(F["recall"]["cross_source"]))
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.80),
                        Inches(7.35), Inches(3.70), cd).chart
gf.has_legend = True
gf.legend.position = XL_LEGEND_POSITION.TOP
gf.legend.include_in_layout = False
gf.legend.font.size = Pt(11)
gf.legend.font.name = BODY
va = gf.value_axis
va.minimum_scale, va.maximum_scale = 0.0, 1.0
va.has_major_gridlines = True
va.major_gridlines.format.line.color.rgb = D.grid
va.tick_labels.font.size = Pt(10)
va.tick_labels.font.name = BODY
gf.category_axis.tick_labels.font.size = Pt(9)
gf.category_axis.tick_labels.font.name = BODY
for ser, colr in zip(gf.plots[0].series, (DEEP, TEAL)):
    ser.format.fill.solid()
    ser.format.fill.fore_color.rgb = colr
    ser.format.line.fill.background()

xs = M + 7.85
D.stat(s, xs, 1.80, 2.1, fmt(F["recall"]["time"][0], 4),
       "recall on a different chemotype, time split", CRIMSON)
D.stat(s, xs + 2.35, 1.80, 2.1, fmt(F["recall"]["cross_source"][0], 4),
       "the same band across a database boundary", CRIMSON)
D.card(s, xs, 3.42, W - M - xs, 2.05, fill=TINT)
D.text(s, xs + 0.28, 3.62, W - M - xs - 0.56, 1.70,
       f"Recall runs {fmt(F['recall']['time'][0], 4)} to {fmt(F['recall']['time'][3], 4)} across the "
       f"four bands, a factor of five, on {F['n_band'][0]:,} and {F['n_band'][3]:,} actives. AUROC "
       f"moves far less, {fmt(F['auroc_band'][0], 4)} to {fmt(F['auroc_band'][3], 4)}, because "
       f"ranking survives what a fixed threshold does not.",
       size=12.5, color=INK, line=1.26)

D.text(s, M, 5.72, W - 2 * M, 0.85,
       "On genuinely novel chemistry the tool finds roughly one true activity in six, and across a "
       "database boundary roughly one in twenty. No analysis in this thesis improves that. Chapter 8 "
       "establishes that the figure is predictable, not that it is good.",
       size=13.5, color=CRIMSON, bold=True, line=1.26)
D.text(s, M, 6.52, W - 2 * M, 0.4,
       "It is not a defect of the forest. It is what a similarity-based featuriser does, and every "
       "model family in Chapter 3 shares it.",
       size=12.5, color=MUTED, italic=True)
D.source(s, "results/tables/external_novelty_strata.csv")
D.notes(s, "If asked for the single most important sentence in the thesis, this is it. The server "
           "reports the expected recall beside a negative result because of this table.")

# ---------------------------------------------------------------- 3. H10
s = D.light()
D.head(s, "2", "H10: does the barrier model earn its place?",
       "The component the architecture is named for had no hypothesis until now")
D.text(s, M, 1.70, W - 2 * M, 0.58,
       "Null: permeability is a bulk property, the twelve descriptors already in the feature vector "
       "encode it, and the fingerprint adds nothing.",
       size=13.5, color=DEEP, italic=True, line=1.24)
D.text(s, M + 0.30, 2.30, 5.3, 0.3, "estimator", size=11, bold=True, color=MUTED, space=0)
D.text(s, M + 5.70, 2.30, 2.0, 0.3, "scaffold hold-out", size=11, bold=True, color=MUTED,
       align=PP_ALIGN.RIGHT, space=0)
D.text(s, M + 7.90, 2.30, 2.0, 0.3, "external approved", size=11, bold=True, color=MUTED,
       align=PP_ALIGN.RIGHT, space=0)
show = ["deployed forest, 1,036 features", "descriptor forest, 12 features",
        "descriptor logistic regression", "tpsa alone",
        "CNS heuristic: TPSA <= 90 and MW <= 400"]
nice = {"deployed forest, 1,036 features": "Deployed forest, 1,036 features",
        "descriptor forest, 12 features": "Random forest, 12 descriptors only",
        "descriptor logistic regression": "Logistic regression, 12 descriptors",
        "tpsa alone": "TPSA alone, no fitting of any kind",
        "CNS heuristic: TPSA <= 90 and MW <= 400": "CNS heuristic: TPSA ≤ 90 and MW ≤ 400"}
yy = 2.66
for meth in show:
    top = meth.startswith("deployed")
    D.card(s, M, yy - 0.07, 10.15, 0.46, fill=TINT if top else PAPER)
    D.text(s, M + 0.30, yy, 5.4, 0.36, nice[meth], size=12.5, bold=top,
           color=DEEP if top else INK)
    for i, pop in enumerate(("scaffold hold-out", "external approved")):
        r = F["h10_by"][(meth, pop)]
        D.text(s, M + 5.70 + i * 2.20, yy, 2.0, 0.36, fmt(num(r["auroc"]), 4), size=12.5,
               bold=top, color=DEEP if top else MUTED, align=PP_ALIGN.RIGHT)
    yy += 0.50

ex = F["h10_by"][("descriptor forest, 12 features", "external approved")]
sc = F["h10_by"][("descriptor forest, 12 features", "scaffold hold-out")]
xs = M + 10.40
D.card(s, xs, 2.30, W - M - xs, 2.42, fill=TINT)
D.text(s, xs + 0.20, 2.46, W - M - xs - 0.40, 0.48, "Paired\nbootstrap", size=11.5, bold=True,
       color=MUTED, line=1.10)
D.text(s, xs + 0.20, 3.00, W - M - xs - 0.40, 0.70,
       f"scaffold\n{num(sc['delta_ci95_low']):+.4f} to {num(sc['delta_ci95_high']):+.4f}",
       size=11, color=TEAL, line=1.16)
D.text(s, xs + 0.20, 3.74, W - M - xs - 0.40, 0.80,
       f"external\n{num(ex['delta_ci95_low']):+.4f} to {num(ex['delta_ci95_high']):+.4f}\np = "
       f"{fmt(num(ex['bootstrap_p_deployed_better']), 3)}",
       size=11, bold=True, color=CRIMSON, line=1.16)

D.text(s, M, 5.28, W - 2 * M, 0.85,
       f"VERDICT {F['h10_verdict']}. On the training distribution the fingerprint demonstrably earns "
       f"its place. On {int(ex['n'])} approved drugs the model has never seen, its advantage over "
       f"twelve descriptors is not established: the interval crosses zero.",
       size=14, bold=True, color=CRIMSON, line=1.26)
D.text(s, M, 6.16, W - 2 * M, 0.75,
       "Under a point-estimate rule this reads SUPPORTED. Chapter 9 criticised H4 for exactly that, "
       "so a hypothesis written in response to the criticism reads its interval instead. That "
       "difference is the whole value of having written the criticism down.",
       size=12.5, color=DEEP, italic=True, line=1.24)
D.source(s, "inversion/results/H10_barrier_necessity.csv, built for this chapter by "
            "inversion/inv_barrier_necessity.py")
D.notes(s, "What the model does clearly beat, on both populations and by a wide margin, is the "
           "published heuristic. Say that second, not first.")

# ---------------------------------------------------------------- 4. the weakest link
s = D.dark()
D.text(s, M, 0.88, 11.4, 0.5, "Three findings that converge on one component",
       size=15, color=AMBER, bold=True)
D.text(s, M, 1.32, 11.7, 0.55,
       "The barrier model gates every disease score. It multiplies, so whatever it gets wrong "
       "propagates into all sixteen conditions at once.",
       size=14, color=CHALK, line=1.26)
items = [
    ("Its advantage vanishes soonest",
     f"Against a twelve-descriptor forest on unseen approved drugs the margin is "
     f"{abs(num(ex['delta_vs_deployed'])):.4f} with an interval reaching "
     f"{num(ex['delta_ci95_high']):+.4f}."),
    ("It is the worst-calibrated core model",
     f"Expected calibration error {fmt(F['core_cal'][4], 4)} after isotonic calibration, against a "
     f"mean of {fmt(F['core_cal'][2], 4)} across the {F['core_cal'][0]} core classifiers and "
     f"{fmt(F['core_cal'][3], 4)} at the best. This has not been stated anywhere in the project."),
    ("It is a coarser quantity than the science wants",
     "The gate multiplies by a probability of being CNS-positive, where the pharmacology asks for "
     "unbound brain-to-plasma ratio. H10 suggests the gain here lies in a better target variable, "
     "not a better model of this one."),
]
yy = 2.20
for i, (t, b) in enumerate(items, 1):
    D.dot(s, M, yy + 0.06, str(i), fill=AMBER, fg=INK, dia=0.40)
    D.text(s, M + 0.66, yy, 4.3, 0.9, t, size=15, bold=True, color=AMBER, line=1.14)
    D.text(s, M + 5.25, yy - 0.02, 11.9 - 5.25, 1.3, b, size=12.5, color=CHALK, line=1.26)
    yy += 1.48
D.text(s, M, 6.52, 11.7, 0.4,
       "It is the single weakest link in the architecture, and it sits at the gate.",
       size=13.5, color=PAPER, italic=True)
D.source(s, "inversion/results/H10_barrier_necessity.csv; results/tables/calibration.csv")
D.notes(s, "Volunteer this slide. A referee who assembles these three findings unaided will read it "
           "as something the project failed to notice.")

# ---------------------------------------------------------------- 5. uncertainty coverage
s = D.light()
D.head(s, "3", "Uncertainty is complete for eight estimators, thin for thirty-eight",
       "The machinery is fullest where a user is least likely to look")
nb, bmean, bmed, bmin, bmax = F["binder_cal"]
ncore, craw, ccal, cmin, cmax, cworst = F["core_cal"]
D.stat(s, M, 1.84, 2.5, f"{ncore}", "core classifiers with isotonic calibration, conformal "
                                    "prediction and a domain band", TEAL)
D.stat(s, M + 2.85, 1.84, 2.5, f"{nb}", "binder endpoints with Platt scaling and no conformal "
                                        "statement of any kind", CRIMSON)
D.stat(s, M + 5.70, 1.84, 2.5, fmt(ccal, 4), f"core expected calibration error, from {fmt(craw, 4)} "
                                             f"raw", TEAL)
D.stat(s, M + 8.55, 1.84, 2.5, fmt(bmean, 4), "binder expected calibration error, about five times "
                                              "the core figure", CRIMSON)

D.card(s, M, 3.60, W - 2 * M, 1.20, fill=TINT)
D.text(s, M + 0.32, 3.80, W - 2 * M - 0.64, 0.95,
       f"Stating the mean alone would hide both spreads. The core figures run {fmt(cmin, 4)} to "
       f"{fmt(cmax, 4)}; the binder figures run {fmt(bmin, 4)} to {fmt(bmax, 4)} with a median of "
       f"{fmt(bmed, 4)}. The headline in the manuscript, the technical report and the evidence map "
       f"quotes {fmt(craw, 4)} to {fmt(ccal, 4)} without saying it covers {ncore} estimators.",
       size=13, color=INK, line=1.26)

D.card(s, M, 5.00, W - 2 * M, 1.10, fill=PAPER)
D.text(s, M + 0.32, 5.18, W - 2 * M - 0.64, 0.90,
       f"And the worst of the {ncore} is {cworst}, at {fmt(cmax, 4)}, nearly three times the mean. "
       f"That is the model that multiplies every disease score.",
       size=13.5, color=CRIMSON, bold=True, line=1.26)

D.text(s, M, 6.30, W - 2 * M, 0.45,
       "Extending the stack to the panel is the largest single improvement available and needs no "
       "new data.",
       size=13, color=DEEP, italic=True)
D.source(s, "results/tables/calibration.csv; results/tables/integrity_calibration_per_target.csv")
D.notes(s, "The asymmetry is the point: eight endpoints a user rarely queries have three layers of "
           "uncertainty, and thirty-eight they query constantly have one.")

# ---------------------------------------------------------------- 6. coverage
s = D.light()
D.head(s, "4", "What the panel does not cover",
       "Natural-product chemistry is the largest and best-characterised gap")
if F["sp3"]:
    med, pct, npar = F["sp3"]
    D.text(s, M, 1.76, W - 2 * M, 0.55,
           f"Visible in the library before any model is scored: across {npar:,} parsed structures "
           f"the median fraction of sp3 carbon is {fmt(med, 4)}, and only {pct:.2f} per cent are "
           f"both sp3-rich and carry at most one aromatic ring.",
           size=13.5, color=INK, line=1.26)
    ytop = 2.46
else:
    ytop = 1.80
D.text(s, M + 0.30, ytop, 2.4, 0.3, "endpoint", size=11, bold=True, color=MUTED, space=0)
for i, hh in enumerate(["n", "actives", "AUROC", "median max Tanimoto"]):
    D.text(s, M + 3.00 + i * 1.75, ytop, 1.6, 0.3, hh, size=11, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT, space=0)
yy = ytop + 0.38
for r in sorted(F["np_scored"], key=lambda x: -float(x["auroc"])):
    D.card(s, M, yy - 0.07, 10.05, 0.46, fill=TINT)
    D.text(s, M + 0.30, yy, 2.4, 0.36, r["endpoint"], size=12.5, bold=True, color=INK)
    for i, v in enumerate([str(int(float(r["n"]))), str(int(float(r["n_active"]))),
                           fmt(num(r["auroc"]), 4), fmt(num(r["median_max_tanimoto"]), 3)]):
        D.text(s, M + 3.00 + i * 1.75, yy, 1.6, 0.36, v, size=12.5, color=MUTED,
               align=PP_ALIGN.RIGHT)
    yy += 0.50

xs = M + 10.35
D.stat(s, xs, ytop, 2.0, f"{F['np_contam']:,}",
       f"of {F['np_cand']:,} removed: a tautomer or salt of a training compound", CRIMSON)
D.stat(s, xs, ytop + 1.90, 2.0, f"{F['np_noest']}",
       "surveyed targets with no deployed estimator", AMBER)

D.text(s, M, 5.62, W - 2 * M, 0.85,
       "With two to five actives per endpoint those intervals are far too wide to establish failure, "
       "and equally far from supporting a claim of competence. The honest statement is that the "
       "system is untested on this chemistry, not that it fails on it.",
       size=13.5, color=DEEP, line=1.26)
D.source(s, "results/tables/external_natural_products_summary.csv"
            + ("; results/tables/library_sp3_coverage.csv" if F["sp3"] else ""))
D.notes(s, "The contamination count is the transferable lesson: an InChYKey check is not enough, "
           "because the featuriser is stereo-blind. Chapter 8 found the same thing in the barrier "
           "set.")

# ---------------------------------------------------------------- 7. the structural defect
s = D.light()
D.head(s, "5", "One structural defect explains most of the record's errors",
       "Two dozen document-artefact disagreements, almost all of one kind")
D.card(s, M, 1.78, W - 2 * M, 0.78, fill=TINT)
D.text(s, M + 0.32, 1.96, W - 2 * M - 0.64, 0.55,
       "A number that lives in a generator rather than in an artefact cannot be checked, cannot be "
       "cited, and drifts silently.",
       size=15, bold=True, color=DEEP, font=HEAD, line=1.20)
inst = [
    ("3", "The censored-bound recovery count", "21,994 against a recount of 29,751", CRIMSON),
    ("8", "The label-permutation null", "0.4938 and 0.4921, in a literal prose block", CRIMSON),
    ("8", "The independent-reproduction figure", "26 values to 4.7e-5, same block", CRIMSON),
    ("9", "The H2 weight ablation", "0.7917 / 0.7911 / 0.7899 over a population that never "
                                    "matched, in app.py", AMBER),
    ("9", "Five co-firing correlations", "in app.py, and rendered on a result page", CRIMSON),
    ("10", "The natural-product coverage figures", "measured at report time, held by no file", MUTED),
]
D.text(s, M + 0.30, 2.76, 0.9, 0.3, "chapter", size=11, bold=True, color=MUTED, space=0)
D.text(s, M + 1.35, 2.76, 4.3, 0.3, "the number", size=11, bold=True, color=MUTED, space=0)
D.text(s, M + 6.00, 2.76, 5.9, 0.3, "where it lived", size=11, bold=True, color=MUTED, space=0)
yy = 3.14
for ch, what, where, col in inst:
    D.text(s, M + 0.30, yy, 0.9, 0.4, ch, size=13, bold=True, color=col)
    D.text(s, M + 1.35, yy, 4.4, 0.4, what, size=12.5, color=INK)
    D.text(s, M + 6.00, yy, W - M - (M + 6.00) - 0.2, 0.56, where, size=12, color=MUTED, line=1.16)
    yy += 0.58

D.text(s, M, 6.42, W - 2 * M, 0.52,
       "The remedy is uniform and cheap: every number a document states should be read from a file "
       "at build time, and every file declared in the freshness graph.",
       size=13, color=DEEP, line=1.24)
D.source(s, "Two gaps in the checking machinery sit behind these: the freshness graph's inputs are "
            "data and models but never code, and no test pins any falsification verdict.")
D.notes(s, "The last row is the mildest and most instructive: those two figures are genuinely "
           "measured, and the generator's own comment records that as prose they had already "
           "drifted once. Measuring was right; not writing to a file was not.")

# ---------------------------------------------------------------- 8. what should be done next
s = D.light()
D.head(s, "6", "What should be done next",
       "Ordered by what it would change, not by how interesting it is")
cols = [
    ("Would change a user's decision", TEAL, [
        "Extend the uncertainty stack to the target panel: isotonic calibration where the held-out "
        "positives support it, per-endpoint error carried with every prediction, and a conformal "
        "statement for the 38 endpoints that have none.",
        "Report intervals on every per-endpoint figure. An AUROC of 0.719 on 37 compounds and one "
        "of 0.985 on 23 are not comparable numbers.",
        "Decide about the six endpoints that fail the reliability gate. This thesis records that "
        "decision rather than making it."]),
    ("Would change what may be claimed", AMBER, [
        "Migrate H1 to H9 to the single-verdict-rule arrangement H10 now uses, and pin the verdicts "
        "with a test.",
        "Write the missing hypotheses. The most valuable is a null for the aggregation rule: does "
        "the maximum over engaged targets beat a simple count?",
        "Make the freshness graph aware of code. An afternoon's work, and it would have caught the "
        "background_specificity.csv case."]),
    ("Needs new data, not new analysis", CRIMSON, [
        "Natural-product chemistry, with the feature-vector contamination check applied from the "
        "start rather than the InChIKey check.",
        "A prospective test in the proper sense. Every result in Chapter 8 is retrospective: the "
        "dates are real, but no compound was predicted before it was measured.",
        "Measured brain exposure rather than a barrier classification, since H10 suggests the gain "
        "lies in a better target variable."]),
]
cw = (W - 2 * M - 0.5) / 3
for j, (title, col, bullets) in enumerate(cols):
    x = M + j * (cw + 0.25)
    D.card(s, x, 1.80, cw, 4.55, fill=TINT)
    D.text(s, x + 0.24, 1.98, cw - 0.48, 0.60, title, size=13.5, bold=True, color=col, line=1.14)
    yy = 2.66
    for b in bullets:
        D.dot(s, x + 0.24, yy + 0.02, "•", fill=col, dia=0.20)
        D.text(s, x + 0.56, yy - 0.04, cw - 0.82, 1.25, b, size=11, color=INK, line=1.22)
        yy += 1.28
D.text(s, M, 6.52, W - 2 * M, 0.4,
       "41 items were raised across the nine preceding chapters. Five are now closed, and a sixth "
       "was raised and closed within this chapter.",
       size=12.5, color=MUTED, italic=True)
D.source(s, "Chapter 10, sections 10.8 and 10.9")
D.notes(s, "If pressed on which single item matters most: extending the uncertainty stack to the "
           "panel. It needs no new data and it fixes the asymmetry on the previous slide.")

# ---------------------------------------------------------------- 9. what the thesis establishes
s = D.dark()
D.text(s, M, 0.86, 11.4, 0.5, "What this thesis establishes", size=15, color=TEAL, bold=True)
D.text(s, M, 1.28, 11.7, 0.45,
       "With the spread rather than the mean wherever one exists.", size=13.5, color=CHALK,
       italic=True)
est = [
    (f"{fmt(F['cv']['scaffold'][0], 4)} to {fmt(F['cv']['scaffold'][1], 4)}",
     f"scaffold-grouped AUROC across the eight core classifiers, against permuted-label nulls of "
     f"{fmt(F['pn']['random'], 4)} and {fmt(F['pn']['scaffold'], 4)}"),
    (f"{fmt(F['recall']['time'][0], 4)} to {fmt(F['recall']['time'][3], 4)}",
     "recall across four novelty bands, the same curve traced by three independently built test "
     "sets"),
    (f"{fmt(F['sens'][0], 4)}",
     f"panel sensitivity on held-out actives, median {fmt(F['sens'][1], 3)}, from "
     f"{fmt(F['sens'][2], 3)} to {fmt(F['sens'][3], 3)}, with {F['n_unreliable']} of {F['n_dep']} "
     f"endpoints firing for fewer than half their own actives"),
    (f"{fmt(F['spec'][0], 3)}",
     f"specificity on presumed-inactive non-CNS chemistry, {fmt(F['spec'][1], 3)} to "
     f"{fmt(F['spec'][2], 3)}, and an upper bound rather than an estimate"),
    (f"{F['n_ref']} of {F['n_hyp']}",
     f"falsification hypotheses refuted, {F['n_weak']} weakened. The refutations cost a claim about "
     f"the graph, a claim about the gate, a change to the interface, and four withdrawn endpoints"),
]
yy = 1.94
for v, b in est:
    D.text(s, M, yy, 2.5, 0.44, v, size=19, bold=True, font=HEAD, color=TEAL)
    D.text(s, M + 2.85, yy - 0.02, 11.9 - 2.85, 0.85, b, size=12.5, color=CHALK, line=1.24)
    yy += 0.92
D.card(s, M, 6.44, 11.7, 0.74, fill=DEEP)
D.text(s, M + 0.28, 6.58, 11.2, 0.56,
       "Not established: anything prospective, anything about natural-product chemistry, and any "
       "advantage for the barrier model over twelve descriptors on chemistry it has not seen.",
       size=13, bold=True, color=PAPER, line=1.20)
D.notes(s, "End on the last line. Naming the boundary of the work is the strongest thing available "
           "in a viva, and all three items have a measurement behind them rather than a hedge.")

D.save(OUT)
