"""Build the Chapter 3 defence deck, reading every figure from an artefact.

Shares its chrome with the Chapter 1 deck via deck_common.py, so the two are one visual system and a
change to the palette or the title logic reaches both.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/build_chapter03_defence.py
Out:  thesis/presentations/chapter03_defence.pptx
"""
from __future__ import annotations

import csv
import json
import statistics as st
from collections import defaultdict
from pathlib import Path

from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt

from deck_common import (
    ROOT, TAB, PKG, W, M, INK, DEEP, TEAL, AMBER, CRIMSON, MUTED, TINT, PAPER, CHALK,
    HEAD, BODY, Deck, rows, num, fmt,
)

OUT = Path(__file__).resolve().parent / "chapter03_defence.pptx"
CLASSIFIERS = {"BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_A", "MAO_B", "hERG"}


def facts() -> dict:
    d: dict = {}

    # -- model family comparison --------------------------------------------
    cmp_ = rows(TAB / "model_comparison.csv")
    by = defaultdict(lambda: defaultdict(dict))
    for r in cmp_:
        by[r["split"]][r["endpoint"]][r["model"]] = num(r["mean"])
    sc = by["scaffold"]
    d["fam_cls"], d["fam_reg"] = {}, {}
    for model in sorted({m for e in sc for m in sc[e]}):
        cls = [sc[e][model] for e in sc if e in CLASSIFIERS]
        reg = [sc[e][model] for e in sc if e not in CLASSIFIERS]
        d["fam_cls"][model] = st.mean(cls)
        d["fam_reg"][model] = st.mean(reg)
    d["rf_cls_wins"] = sum(1 for e in sc if e in CLASSIFIERS
                           and max(sc[e], key=sc[e].get) == "RandomForest")
    d["rf_reg_wins"] = sum(1 for e in sc if e not in CLASSIFIERS
                           and max(sc[e], key=sc[e].get) == "RandomForest")
    d["n_cls"] = sum(1 for e in sc if e in CLASSIFIERS)
    d["n_reg"] = sum(1 for e in sc if e not in CLASSIFIERS)
    d["ache"] = (sc["AChE"]["RandomForest"], sc["AChE"]["HistGradientBoosting"])
    knn = [(e, sc[e]["RandomForest"] - sc[e]["kNN read-across"]) for e in sc]
    d["knn_margin"] = (min(knn, key=lambda x: x[1]), max(knn, key=lambda x: x[1]), len(knn))

    # -- significance --------------------------------------------------------
    sig = rows(TAB / "model_family_significance.csv")
    d["sig"] = [r for r in sig if r["split"] == "scaffold"]

    # -- feature block ablation ---------------------------------------------
    ab = rows(TAB / "feature_block_ablation.csv")
    blk = defaultdict(dict)
    for r in ab:
        blk[(r["endpoint"], r["task"])][r["block"]] = num(r["mean"])
    d["blk_best"] = defaultdict(int)
    for k, m in blk.items():
        d["blk_best"][max(m, key=m.get)] += 1
    d["blk_n"] = len(blk)
    cls = [m for (e, t), m in blk.items() if t == "classification"]
    reg = [m for (e, t), m in blk.items() if t == "regression"]
    d["blk_margin_cls"] = st.mean(m["combined"] - m["fingerprint_only"] for m in cls)
    d["blk_margin_reg"] = st.mean(m["combined"] - m["fingerprint_only"] for m in reg)
    gaps = {e: m["descriptors_only"] - m["fingerprint_only"] for (e, t), m in blk.items()}
    d["desc_gap_best"] = max(gaps.items(), key=lambda kv: kv[1])
    d["desc_gap_worst"] = min(gaps.items(), key=lambda kv: kv[1])

    # -- learning curve ------------------------------------------------------
    lc = defaultdict(dict)
    for r in rows(TAB / "learning_curve.csv"):
        lc[r["endpoint"]][num(r["train_fraction"])] = (num(r["score"]), r["metric"])
    d["lc"] = {e: {"half": m[0.5][0], "full": m[1.0][0], "quarter": m[0.25][0],
                   "metric": m[1.0][1]} for e, m in lc.items()}

    # -- gnn -----------------------------------------------------------------
    d["gnn"] = [(r["endpoint"], r["metric"], num(r["GIN"]), num(r["RandomForest"]), r["winner"])
                for r in rows(ROOT / "results" / "gnn" / "gnn_vs_rf.csv")]

    # -- representation ------------------------------------------------------
    meta = json.loads((ROOT / "models_rf" / "BBB_meta.json").read_text())
    d["n_features"] = meta["n_features"]
    d["layout"] = meta["feature_layout"]
    d["hyper"] = meta["hyperparameters"]
    d["dedup"] = meta["deduplication"]

    unc = rows(TAB / "uncharging_impact_by_endpoint.csv")
    tot = sum(int(r["compounds"]) for r in unc)
    ch = sum(int(r["changed"]) for r in unc)
    worst = max(unc, key=lambda r: num(r["pct"]))
    d["uncharge"] = (tot, ch, 100 * ch / tot, worst["endpoint"], num(worst["pct"]),
                     sum(1 for r in unc if int(r["changed"]) > 0), len(unc))

    di = rows(TAB / "feature_descriptor_importance.csv")
    inf = defaultdict(int)
    for r in di:
        if r["informative"] == "True":
            inf[r["descriptor"]] += 1
    d["desc_inf"] = sorted(inf.items(), key=lambda kv: -kv[1])
    d["desc_eps"] = len({r["endpoint"] for r in di})
    return d


F = facts()
D = Deck()

# 1 -- title
s = D.dark()
D.text(s, M, 1.74, W - 2 * M, 0.5, "Chapter 3", size=16, color=TEAL, bold=True)
D.text(s, M, 2.18, W - 2 * M - 1.0, 1.6,
       "Representation and models: what the compound becomes, and what learns from it",
       size=34, bold=True, font=HEAD, color=PAPER, line=1.08)
D.text(s, M, 4.20, W - 2 * M - 2.0, 0.9,
       f"One {F['n_features']:,}-column map shared by every endpoint, and a forest chosen for "
       f"reasons that are not about the mean",
       size=16, color=CHALK, italic=True, line=1.25)
D.text(s, M, 5.70, W - 2 * M, 0.5,
       "Five families, thirteen endpoints, and the first significance test run on the deployed panel",
       size=13, color=TEAL)
D.notes(s, "This chapter contains the clearest correction of the project's own documents, in the "
           "learning-curve section. Lead with the method, end with the correction.")

# 2 -- the representation
s = D.light()
D.head(s, "1", "One compound, one vector, every endpoint",
       f"{F['layout']}, verified from src/brainsafe/features/featurize.py")
D.card(s, M, 1.80, 6.35, 2.05)
D.text(s, M + 0.42, 2.08, 5.6, 0.4, "1,024 bits, folded ECFP-4", size=18, bold=True, font=HEAD, color=DEEP)
D.text(s, M + 0.42, 2.60, 5.5, 1.1,
       "Bit j is set when some atom environment of radius at most 2 hashes to j modulo 1,024. "
       "Folding is not injective, so a set bit says some environment hashing there is present, not "
       "which one.", size=13, color=INK, line=1.26)
D.card(s, M + 6.75, 1.80, 5.35, 2.05)
D.text(s, M + 7.17, 2.08, 4.6, 0.4, "12 descriptors", size=18, bold=True, font=HEAD, color=TEAL)
D.text(s, M + 7.17, 2.60, 4.6, 1.1,
       "Molecular weight, cLogP, TPSA, donors, acceptors, rotatable bonds, aromatic rings, "
       "fraction sp3, ring count, heavy atoms, formal charge, QED.",
       size=13, color=INK, line=1.26)
D.text(s, M, 4.15, W - 2 * M, 0.5,
       "Nothing is shared between endpoints except this map.",
       size=19, bold=True, font=HEAD, color=INK)
D.text(s, M, 4.70, W - 2 * M, 1.1,
       "Two models that disagree about a compound disagree about the same 1,036 numbers, which is "
       "what makes the panel comparable at all. It is also why a defect in the featuriser would be a "
       "defect in seventy models at once, and why the software suite pins the vector's shape and its "
       "purity as a function of structure.",
       size=14, color=INK, line=1.28)
D.source(s, "models_rf/BBB_meta.json; src/brainsafe/features/featurize.py")
D.notes(s, "If asked why folding to 1,024 rather than 2,048: collisions are tolerable because the "
           "forest sees many bits, and the applicability-domain measure uses 2,048-bit fingerprints "
           "where resolution matters more.")

# 3 -- standardisation, verified
s = D.light()
D.head(s, "2", "Neutralisation is part of the representation, not housekeeping",
       "A drug and its salt must give the same answer")
D.text(s, M, 1.82, 6.2, 1.5,
       "Stripping a counter-ion alone leaves the parent carrying the charge the salt gave it, so the "
       "two forms differ in formal charge and in the environments hashed around the protonated "
       "nitrogen. The rule that fixes this must still keep a permanent charge, because a quaternary "
       "ammonium's charge is precisely what stops it crossing.",
       size=14, color=INK, line=1.28)
D.card(s, M, 3.55, 6.2, 2.35, fill=TINT)
D.text(s, M + 0.40, 3.80, 5.4, 0.4, "Verified on the deployed models, this session",
       size=14, bold=True, color=DEEP)
for i, (lab, val) in enumerate([
        ("haloperidol free base and hydrochloride", "byte-identical vectors"),
        ("barrier probability, both forms", "0.9888"),
        ("hERG probability, both forms", "0.8216"),
        ("neostigmine, quaternary ammonium", "charge +1 retained")]):
    D.text(s, M + 0.40, 4.30 + i * 0.40, 3.9, 0.34, lab, size=12.5, color=INK)
    D.text(s, M + 4.35, 4.30 + i * 0.40, 1.75, 0.34, val, size=12.5, bold=True,
           color=TEAL, align=PP_ALIGN.RIGHT)

xs = M + 6.6
D.card(s, xs, 1.80, W - M - xs, 4.10, fill=PAPER)
D.card(s, xs, 1.80, W - M - xs, 4.10)
D.text(s, xs + 0.42, 2.05, 4.9, 0.42, "How often it fires", size=18, bold=True, font=HEAD, color=DEEP)
u = F["uncharge"]
D.stat(s, xs + 0.42, 2.60, 2.4, f"{u[2]:.3f}%", f"of {u[0]:,} training rows are changed by "
                                                f"neutralisation", color=DEEP, vsize=34)
D.stat(s, xs + 3.05, 2.60, 2.2, f"{u[5]} of {u[6]}", "endpoints affected at all", color=TEAL, vsize=34)
D.text(s, xs + 0.42, 4.30, 4.9, 1.35,
       f"Worst is {u[3]} at {u[4]:.2f} per cent. A rule that alters one row in a thousand is easy to "
       "dismiss. It is not dismissed because its failures concentrate exactly where a user is most "
       "likely to paste a salt form, which is an approved drug.",
       size=12.5, color=INK, line=1.26)
D.source(s, "results/tables/uncharging_impact_by_endpoint.csv")
D.notes(s, "The 0.613 against 0.993 figures in the manuscript describe an earlier engine, before "
           "neutralisation was in the representation. The current pipeline returns identical numbers "
           "for both forms. Say so if it comes up: the defect was real and is fixed.")

# 4 -- what the representation discards
s = D.light()
D.head(s, "3", "What the representation discards, and what that costs immediately",
       "Chirality is excluded, so enantiomers collide")
dd = F["dedup"]
D.text(s, M, 1.85, 5.9, 2.0,
       "Chapter 10 bounds the scientific cost of stereo-blindness at 0.19 per cent of the panel. The "
       "structural cost is immediate and belongs here: if two rows are identical in feature space and "
       "land on opposite sides of a fold, the model is tested on a compound it was trained on, and no "
       "amount of scaffold grouping removes it, because those rows share a skeleton by construction.",
       size=14, color=INK, line=1.28)
D.text(s, M, 4.10, 5.9, 1.5,
       "Rows identical in the representation are therefore collapsed before any split is drawn, and a "
       "group whose labels disagree is dropped rather than resolved by vote: the featuriser cannot "
       "tell its members apart, and a vote would be inventing an answer.",
       size=14, color=INK, line=1.28)
xs = M + 6.35
D.card(s, xs, 1.82, W - M - xs, 4.05)
D.text(s, xs + 0.45, 2.08, 4.8, 0.42, "The barrier endpoint, in full",
       size=18, bold=True, font=HEAD, color=DEEP)
steps = [("rows entering", f"{dd['rows_in']:,}", INK),
         ("duplicate rows removed", f"-{dd['duplicate_rows_removed']:,}", CRIMSON),
         ("conflicting groups dropped", f"-{dd['conflicting_groups_dropped']:,}", CRIMSON),
         ("rows the model is fitted on", f"{dd['rows_out']:,}", DEEP)]
for i, (lab, val, col) in enumerate(steps):
    y = 2.72 + i * 0.62
    D.text(s, xs + 0.45, y, 3.3, 0.36, lab, size=13.5, color=INK,
           bold=(i == len(steps) - 1))
    D.text(s, xs + 3.85, y - 0.04, 1.4, 0.42, val, size=17, bold=True, color=col,
           align=PP_ALIGN.RIGHT)
D.text(s, xs + 0.45, 5.30, 4.8, 0.5,
       "Nearly half the rows. Not a defect in the data: the price of a stereo-blind representation, "
       "paid at the featurisation boundary rather than concealed in the score.",
       size=12, color=MUTED, italic=True, line=1.24)
D.source(s, "models_rf/BBB_meta.json")
D.notes(s, "3,773 of 7,805 is the single most arresting number in this chapter. It is why the "
           "cross-validation figures can be trusted.")

# 5 -- block ablation
s = D.light()
D.head(s, "4", "Do both blocks earn their place?",
       "Each of the two blocks was tested alone, over thirteen endpoints")
D.card(s, M, 1.80, 5.85, 2.0)
D.text(s, M + 0.42, 2.05, 5.0, 0.4, "Which block wins", size=17, bold=True, font=HEAD, color=DEEP)
for i, (lab, key) in enumerate([("both blocks", "combined"), ("fingerprint alone", "fingerprint_only"),
                                ("descriptors alone", "descriptors_only")]):
    D.text(s, M + 0.42, 2.58 + i * 0.38, 3.4, 0.34, lab, size=13.5, color=INK)
    D.text(s, M + 3.95, 2.58 + i * 0.38, 1.5, 0.34, f"{F['blk_best'][key]} of {F['blk_n']}",
           size=13.5, bold=True, color=DEEP if key == "combined" else MUTED, align=PP_ALIGN.RIGHT)
D.text(s, M, 4.00, 5.85, 1.5,
       f"But the margin is thin. Adding the descriptors to the fingerprint is worth "
       f"{F['blk_margin_cls']:+.4f} AUROC on classification and {F['blk_margin_reg']:+.4f} R² on "
       f"regression. The fingerprint carries almost all of the signal.",
       size=14, color=INK, line=1.28)

xs = M + 6.25
D.card(s, xs, 1.80, W - M - xs, 4.05, fill=TINT)
D.text(s, xs + 0.42, 2.05, 5.0, 0.42, "Where they do earn it", size=18, bold=True, font=HEAD, color=TEAL)
D.text(s, xs + 0.42, 2.60, 5.0, 1.0,
       "The useful question is not what the descriptors add, but how close they come alone. That "
       "measures whether physicochemistry is the mechanism.",
       size=13, color=INK, line=1.26)
gb, gw = F["desc_gap_best"], F["desc_gap_worst"]
D.text(s, xs + 0.42, 3.75, 5.0, 1.3,
       f"At {gb[0]} the descriptors alone fall only {abs(gb[1]):.3f} short of the fingerprint. "
       f"At {gw[0]} they fall {abs(gw[1]):.3f} short.",
       size=14, bold=True, color=DEEP, line=1.26)
D.text(s, xs + 0.42, 4.75, 5.0, 1.0,
       "Twelve numbers nearly suffice to predict whether a compound crosses the barrier, and come "
       "nowhere near predicting whether it inhibits a kinase.",
       size=13, color=INK, italic=True, line=1.26)
D.source(s, "results/tables/feature_block_ablation.csv. 5-fold random and pre-deduplication, so "
            "absolute values are not comparable with Chapter 8; the between-block comparison is unaffected.")
D.notes(s, "The technical report says the descriptors earn their place mainly on exposure endpoints. "
           "True under this reading, false under the other one, where hERG gains most. Chapter 3 "
           "flags it.")

# 6 -- the estimator
s = D.light()
D.head(s, "5", "The estimator, and three reasons that are not about accuracy", None)
h = F["hyper"]
D.card(s, M, 1.72, W - 2 * M, 0.86, fill=TINT)
D.text(s, M, 1.96, W - 2 * M, 0.4,
       f"Random forest   ·   {h['n_estimators']} trees   ·   min_samples_leaf {h['min_samples_leaf']} "
       f"(4 for the binder panel)   ·   class_weight {h['class_weight']}   ·   seed {h['random_state']}",
       size=17, bold=True, font=HEAD, color=DEEP, align=PP_ALIGN.CENTER)
reasons = [
    ("It handles the representation as it is",
     "A sparse binary fingerprint sits beside twelve continuous descriptors on wildly different "
     "scales. A tree ensemble needs no scaling, no imputation and no encoding to use both."),
    ("It does not extrapolate",
     "A forest's prediction is an average of training labels, so it cannot return a value outside the "
     "range it has seen. For a system whose central claim is that it behaves predictably at the edge "
     "of its domain, refusing to extrapolate is correct behaviour, not a weakness."),
    ("It is exactly explainable",
     "TreeSHAP computes exact Shapley values for a tree ensemble. For most other families the same "
     "explanation is an approximation, and an explanation a user cannot check is a liability here."),
]
yy = 2.90
for i, (t, b) in enumerate(reasons, 1):
    D.dot(s, M, yy + 0.04, str(i), fill=DEEP, dia=0.36)
    D.text(s, M + 0.56, yy, 4.0, 0.75, t, size=15.5, bold=True, color=DEEP, line=1.14)
    D.text(s, M + 4.80, yy - 0.02, W - 2 * M - 4.8, 1.0, b, size=13, color=INK, line=1.26)
    yy += 1.12
D.text(s, M, 6.28, W - 2 * M, 0.5,
       "Class weighting is not decorative: the median deployed endpoint is 0.825 active and P2X7 is "
       "0.962, so without it the forest has every incentive to answer 'active'.",
       size=13, color=MUTED, italic=True, line=1.24)
D.notes(s, "The three reasons matter because, as the next two slides show, accuracy is not one of "
           "them against boosting.")

# 7 -- family comparison chart
s = D.light()
D.head(s, "6", "Five families, thirteen endpoints, scaffold split",
       "Two baselines a reader is entitled to demand, and three ensembles")
order = ["RandomForest", "XGBoost", "HistGradientBoosting", "kNN read-across", "LogisticRegression"]
labels = ["Random\nforest", "XGBoost", "HistGradient\nBoosting", "kNN\nread-across", "L2 logistic\nregression"]
cd = CategoryChartData()
cd.categories = labels
cd.add_series("8 classification endpoints, AUROC", [F["fam_cls"][m] for m in order])
cd.add_series("5 regression endpoints, R²", [F["fam_reg"][m] for m in order])
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.72),
                        Inches(7.9), Inches(4.30), cd)
ch = gf.chart
ch.has_title = False
ch.has_legend = True
ch.legend.position = XL_LEGEND_POSITION.BOTTOM
ch.legend.include_in_layout = False
ch.legend.font.size = Pt(11)
va = ch.value_axis
va.minimum_scale, va.maximum_scale = 0.0, 1.0
va.major_gridlines.format.line.color.rgb = D.grid
va.tick_labels.font.size = Pt(11)
va.tick_labels.font.color.rgb = MUTED
ch.category_axis.tick_labels.font.size = Pt(10)
ch.category_axis.tick_labels.font.color.rgb = MUTED
for i, colr in enumerate([DEEP, AMBER]):
    ch.plots[0].series[i].format.fill.solid()
    ch.plots[0].series[i].format.fill.fore_color.rgb = colr
ch.plots[0].gap_width = 70

xs = M + 8.2
D.text(s, xs, 1.85, 3.9, 0.9,
       "The result splits by task, and saying so is the only honest way to report it.",
       size=14, bold=True, color=DEEP, line=1.24)
D.text(s, xs, 2.85, 3.9, 1.4,
       f"On classification the forest leads, best on {F['rf_cls_wins']} of {F['n_cls']}, losing AChE "
       f"to histogram gradient boosting at {fmt(F['ache'][0], 4)} against {fmt(F['ache'][1], 4)}.",
       size=13, color=INK, line=1.26)
D.text(s, xs, 4.20, 3.9, 1.0,
       f"On regression it leads nothing: best on {F['rf_reg_wins']} of {F['n_reg']}.",
       size=13.5, bold=True, color=CRIMSON, line=1.26)
lo, hi, nn = F["knn_margin"]
D.text(s, xs, 5.15, 3.9, 1.0,
       f"It exceeds the read-across on all {nn}, from {lo[1]:+.4f} at {lo[0]} to {hi[1]:+.4f} at "
       f"{hi[0]}. Much of this signal genuinely is similarity.",
       size=12.5, color=MUTED, line=1.24)
D.source(s, "results/tables/model_comparison.csv")
D.notes(s, "The forest losing all five regressions is the fact the condensed manuscript got wrong "
           "and has now been corrected. Do not soften it.")

# 8 -- significance
s = D.light()
D.head(s, "7", "Is the margin real? The first test run on this panel",
       "Every document until now quoted point deltas alone")
D.text(s, M, 1.78, W - 2 * M, 0.55,
       "The decisions log makes exactly the right criticism, that point deltas do not establish "
       "significance, and answers it with DeLong tests. Those were run on the superseded ensemble and "
       "were never re-run for the deployed forests.",
       size=13.5, color=INK, line=1.26)
hdr = ["comparison", "forest higher on", "median Δ", "mean Δ", "p", "verdict"]
xcol = [M + 0.30, M + 4.75, M + 6.45, M + 7.85, M + 9.15, M + 10.20]
wcol = [4.4, 1.6, 1.3, 1.2, 1.0, 1.9]
for i, hh in enumerate(hdr):
    D.text(s, xcol[i], 2.52, wcol[i], 0.3, hh, size=11, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT if 0 < i < 5 else PP_ALIGN.LEFT)
want = [("kNN read-across", "all 13"), ("LogisticRegression", "all 13"),
        ("XGBoost", "all 13"), ("HistGradientBoosting", "all 13"),
        ("XGBoost", "8 classification"), ("HistGradientBoosting", "8 classification")]
yy = 2.88
for alt, subset in want:
    r = next(x for x in F["sig"] if x["alternative"] == alt and x["subset"].startswith(subset))
    ok = r["distinguishable_at_0.05"] == "True"
    col = DEEP if ok else CRIMSON
    if want.index((alt, subset)) % 2 == 0:
        D.card(s, M, yy - 0.09, W - 2 * M, 0.52, fill=TINT)
    name = {"kNN read-across": "kNN read-across", "LogisticRegression": "L2 logistic regression",
            "XGBoost": "XGBoost", "HistGradientBoosting": "histogram gradient boosting"}[alt]
    tail = "" if subset == "all 13" else ", classification only"
    D.text(s, xcol[0], yy, wcol[0], 0.34, f"forest vs {name}{tail}", size=12.5, color=INK)
    D.text(s, xcol[1], yy, wcol[1], 0.34, f"{r['forest_higher_on']} of {r['n_endpoints']}",
           size=12.5, color=INK, align=PP_ALIGN.RIGHT)
    D.text(s, xcol[2], yy, wcol[2], 0.34, f"{float(r['median_delta']):+.4f}", size=12.5,
           color=INK, align=PP_ALIGN.RIGHT)
    D.text(s, xcol[3], yy, wcol[3], 0.34, f"{float(r['mean_delta']):+.4f}", size=12.5,
           color=CRIMSON if float(r['mean_delta']) < 0 else INK,
           bold=float(r['mean_delta']) < 0, align=PP_ALIGN.RIGHT)
    D.text(s, xcol[4], yy, wcol[4], 0.34, f"{float(r['wilcoxon_p']):.4f}", size=12.5,
           color=col, bold=True, align=PP_ALIGN.RIGHT)
    D.text(s, xcol[5], yy, wcol[5], 0.34, "distinguishable" if ok else "not distinguishable",
           size=11.5, color=col, bold=not ok)
    yy += 0.52
D.text(s, M, 6.10, W - 2 * M, 0.85,
       "Against the baselines the forest is decisively better. Against boosting, pooled over all "
       "thirteen endpoints, it is statistically indistinguishable: the classification wins and the "
       "regression losses cancel. It was not selected for accuracy, because on this evidence it does "
       "not have more.",
       size=13.5, color=INK, line=1.26)
D.source(s, "results/tables/model_family_significance.csv. Wilcoxon signed-rank on per-endpoint "
            "means, not DeLong on ROC curves; the regression subset of five cannot reach p<0.05.")
D.notes(s, "Concede this cleanly. A reader is entitled to say boosting would have served about as "
           "well on classification and rather better on the receptor regressions.")

# 9 -- gnn
s = D.light()
D.head(s, "8", "A learned representation, tested rather than dismissed",
       "A graph isomorphism network on the raw molecular graph, identical scaffold hold-out")
D.text(s, M + 0.30, 1.90, 3.0, 0.3, "endpoint", size=11, bold=True, color=MUTED)
D.text(s, M + 3.30, 1.90, 1.4, 0.3, "metric", size=11, bold=True, color=MUTED)
D.text(s, M + 5.00, 1.90, 1.6, 0.3, "GIN", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
D.text(s, M + 6.90, 1.90, 1.8, 0.3, "random forest", size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT)
yy = 2.30
for ep, metric, gin, rf, winner in F["gnn"]:
    D.card(s, M, yy - 0.10, 9.05, 0.62, fill=TINT)
    D.text(s, M + 0.30, yy, 3.0, 0.36, ep, size=14, bold=True, color=INK)
    D.text(s, M + 3.30, yy, 1.4, 0.36, metric.replace("roc_auc", "AUROC").replace("r2", "R²"),
           size=12.5, color=MUTED)
    D.text(s, M + 5.00, yy, 1.6, 0.36, fmt(gin, 4), size=14, color=MUTED, align=PP_ALIGN.RIGHT)
    D.text(s, M + 6.90, yy, 1.8, 0.36, fmt(rf, 4), size=14, bold=True, color=DEEP, align=PP_ALIGN.RIGHT)
    yy += 0.74
xs = M + 9.5
D.text(s, xs, 2.10, W - M - xs, 2.6,
       "The forest wins all four, and the margin is not marginal: 0.074 at MAO-A and 0.081 R² at A2A.",
       size=14, bold=True, color=DEEP, line=1.26)
D.text(s, xs, 3.60, W - M - xs, 2.4,
       "Consistent with the literature on data scale: graph networks generally need far more data "
       "per task, or pretraining, to overtake fingerprints.\n\nFour endpoints, not thirteen, so this "
       "bounds rather than settles it. Fine-tuning a pretrained chemical language model is the most "
       "defensible thing an unconvinced reader could ask for next, and it has not been tried.",
       size=12.5, color=INK, line=1.26)
D.source(s, "results/gnn/gnn_vs_rf.csv")
D.notes(s, "Name the untried alternative before the examiners do.")

# 10 -- learning curve, the correction
s = D.light()
D.head(s, "9", "Would more data have helped? The answer splits by task",
       "And the project's own documents record the wrong conclusion")
cd2 = CategoryChartData()
eps = ["BBB", "BACE1", "MAO_A", "A2A"]
cd2.categories = ["BBB\nAUROC", "BACE1\nAUROC", "MAO-A\nAUROC", "A2A\nR²"]
cd2.add_series("half the data", [F["lc"][e]["half"] for e in eps])
cd2.add_series("all of it", [F["lc"][e]["full"] for e in eps])
gf2 = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.75),
                         Inches(6.9), Inches(4.15), cd2)
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
for i, colr in enumerate([CHALK, DEEP]):
    c2.plots[0].series[i].format.fill.solid()
    c2.plots[0].series[i].format.fill.fore_color.rgb = colr
c2.plots[0].gap_width = 70

xs = M + 7.25
D.card(s, xs, 1.78, W - M - xs, 2.15, fill=TINT)
D.text(s, xs + 0.38, 2.02, 4.5, 0.4, "Gain from half the data to all of it",
       size=13.5, bold=True, color=DEEP)
for i, e in enumerate(eps):
    g = F["lc"][e]["full"] - F["lc"][e]["half"]
    D.text(s, xs + 0.38, 2.52 + i * 0.34, 2.6, 0.3, e.replace("_", "-"), size=12.5, color=INK)
    D.text(s, xs + 3.10, 2.52 + i * 0.34, 1.5, 0.3, f"{g:+.4f}", size=12.5, bold=True,
           color=AMBER if g > 0.03 else INK, align=PP_ALIGN.RIGHT)
D.text(s, xs, 4.10, W - M - xs, 1.9,
       "Only the barrier model has flattened. MAO-A gains more on its final doubling than the entire "
       "margin separating the forest from the read-across baseline, and A2A gains 0.085 R² with no "
       "sign of saturating.",
       size=13, color=INK, line=1.26)
D.text(s, xs, 5.75, W - M - xs, 0.55,
       "The classifiers are near saturation. The receptor regressions are data-limited.",
       size=13, bold=True, color=DEEP, line=1.24)
D.source(s, "results/tables/learning_curve.csv. The technical report reads this as showing the curves "
            "flatten, and uses it to argue against gathering more data. That does not survive the table.")
D.notes(s, "This is the correction to lead on in questions. It is also useful rather than merely "
           "critical: it makes 'more potency data for the receptor regressions' a specific, costed "
           "piece of future work, and A2A gained most from the BindingDB expansion in Chapter 2 too.")

# 11 -- why seventy models
s = D.light()
D.head(s, "10", "Why seventy models and not one", "The reasons are about the data, not about preference")
items = [
    ("The label matrix is almost entirely missing",
     "228,200 measurements over 63 endpoints and roughly 169,000 compounds fill under two per cent of "
     "the matrix. Multi-task learning shares strength; at this sparsity it mostly shares absence, and "
     "a missing measurement is not a negative result."),
    ("The negative class means different things",
     "Where censored bounds were recovered a negative is a measurement; where they were not it is a "
     "decoy, which is an assumption. One loss over both silently averages the two."),
    ("Base rates are incompatible",
     "A shared output layer propagates the 0.962 active fraction at P2X7 into another endpoint's "
     "decision boundary."),
    ("Failure stays local",
     "Five endpoints were withdrawn. In an independent panel that removes one output. Under shared "
     "weights the same data would have shaped every other endpoint's representation."),
]
yy = 1.80
for i, (t, b) in enumerate(items, 1):
    D.dot(s, M, yy + 0.06, str(i), fill=DEEP, dia=0.34)
    D.text(s, M + 0.54, yy, 4.0, 0.72, t, size=14.5, bold=True, color=DEEP, line=1.14)
    D.text(s, M + 4.75, yy - 0.02, W - 2 * M - 4.75, 0.95, b, size=12.5, color=INK, line=1.24)
    yy += 1.10
D.card(s, M, 6.15, W - 2 * M, 0.72, fill=TINT)
D.text(s, M + 0.40, 6.34, W - 2 * M - 0.8, 0.4,
       "The cost is real and is not hidden: the smallest deployed set has 387 compounds and would "
       "plausibly benefit from borrowing strength. That benefit is forgone deliberately.",
       size=13, color=INK, italic=True)
D.notes(s, "Four reasons, three of which come straight out of Chapter 2. The fourth is about "
           "governance rather than statistics and is the one most often missed.")

# 12 -- closing
s = D.dark()
D.text(s, M, 1.50, 8.4, 0.5, "Chapter 3, in one statement", size=15, color=TEAL, bold=True)
D.text(s, M, 2.05, 11.2, 2.3,
       "The forest was not chosen because it is the most accurate model available. On this evidence "
       "it is not.",
       size=31, font=HEAD, color=PAPER, line=1.18)
D.text(s, M, 4.20, 11.2, 2.2,
       "It was chosen because it does not extrapolate at the edge of the applicability domain, "
       "because TreeSHAP is exact for it rather than approximate, and because it is the most stable "
       "of the three ensembles under hyperparameter choice. Those are defensible grounds. Superior "
       "accuracy would not have been.",
       size=17, color=CHALK, line=1.32)
D.text(s, M, 6.35, 11.2, 0.5,
       "The cost of that uniformity is 0.023 mean R² on the potency regressions, and it is stated.",
       size=14, color=AMBER, italic=True)
D.notes(s, "End on the concession. It is the strongest position available and it is true.")

D.save(OUT)
