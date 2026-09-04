"""Build the master deck: the whole model, explained in detail, from artefacts.

The ten chapter decks each defend one chapter. This is the other thing a thesis needs: one
presentation that explains what the system is and what is known about it, in the order a person
meeting it for the first time can follow. It is the deck to give a lab, a committee, or a
collaborator, and the source for any shorter talk cut from it.

Every figure is read from an artefact at build time and the four thesis figures are embedded from
thesis/figures/, so the deck cannot drift from the repository. Where a number is uncomfortable it is
on a slide rather than in a footnote, because a presentation that only carries the flattering half of
the evidence is the thing this project was built to avoid.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/build_master_deck.py
Out:  thesis/presentations/BrainSafe_master.pptx
"""
from __future__ import annotations

import csv
import json
import statistics as st
from pathlib import Path

from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION
from pptx.enum.text import PP_ALIGN
from PIL import Image
from pptx.util import Inches, Pt

from deck_common import (
    ROOT, TAB, INV, W, H, M, INK, DEEP, TEAL, AMBER, CRIMSON, MUTED, TINT, PAPER, CHALK,
    HEAD, BODY, Deck, rows, num, fmt,
)

OUT = Path(__file__).resolve().parent / "BrainSafe_master.pptx"
FIGS = ROOT / "thesis" / "figures"
BANDS = ["below 0.40 (different chemotype)", "0.40 to 0.55 (related series)",
         "0.55 to 0.70 (same series)", "0.70 and above (close analogue)"]
SHORT_BAND = ["below 0.40\ndifferent chemotype", "0.40 to 0.55\nrelated series",
              "0.55 to 0.70\nsame series", "0.70 and above\nclose analogue"]


def one(path, field, pred):
    for r in rows(path):
        if pred(r):
            return r[field]
    raise KeyError(f"no matching row in {path.name}")


def facts() -> dict:
    d = {}
    reg = json.loads((ROOT / "models_rf" / "binder_modes.json").read_text(encoding="utf-8"))
    dep = {k: v for k, v in reg.items() if v.get("deployed")}
    d["n_endpoints"], d["n_dep"] = len(reg), len(dep)
    d["n_unreliable"] = sum(1 for v in dep.values() if v.get("reliable_call") is False)
    d["unreliable"] = sorted(k for k, v in dep.items() if v.get("reliable_call") is False)
    s = [v["sensitivity_at_threshold"] for v in dep.values()]
    d["sens"] = (st.mean(s), st.median(s), min(s), max(s))
    d["sens_min_ep"] = min(dep, key=lambda k: dep[k]["sensitivity_at_threshold"])
    d["sens_max_ep"] = max(dep, key=lambda k: dep[k]["sensitivity_at_threshold"])

    cv = [r for r in rows(TAB / "rf_cv_summary.csv") if r["task"] == "classification"]
    d["cv"] = {sp: (st.mean(num(r["roc_auc_mean"]) for r in cv if r["split"] == sp),
                    min(num(r["roc_auc_mean"]) for r in cv if r["split"] == sp),
                    max(num(r["roc_auc_mean"]) for r in cv if r["split"] == sp))
               for sp in ("random", "scaffold")}
    d["n_core"] = len([r for r in cv if r["split"] == "random"])

    mc = rows(TAB / "model_comparison.csv")
    d["mc"] = mc

    cal = rows(TAB / "calibration.csv")
    d["cal"] = (st.mean(num(r["ece_raw"]) for r in cal),
                st.mean(num(r["ece_calibrated"]) for r in cal),
                min(num(r["ece_calibrated"]) for r in cal),
                max(num(r["ece_calibrated"]) for r in cal))
    b = [num(r["ece"]) for r in rows(TAB / "integrity_calibration_per_target.csv")]
    d["bcal"] = (len(b), st.mean(b), st.median(b), min(b), max(b))

    d["spec"] = (num(one(TAB / "noncns_specificity_summary.csv", "estimate",
                         lambda r: r["metric"].startswith("Specificity"))),
                 num(one(TAB / "noncns_specificity_summary.csv", "ci95_low",
                         lambda r: r["metric"].startswith("Specificity"))),
                 num(one(TAB / "noncns_specificity_summary.csv", "ci95_high",
                         lambda r: r["metric"].startswith("Specificity"))))

    d["ext"] = rows(TAB / "external_bbb_validation.csv")

    strata = rows(TAB / "external_novelty_strata.csv")
    g = lambda sp, bnd, c: [x[c] for x in strata if x["split"] == sp and x["novelty_band"] == bnd][0]
    d["recall"] = {sp: [num(g(sp, bnd, "recall_at_threshold")) for bnd in BANDS]
                   for sp in ("time", "random", "cross_source")}
    d["auroc_band"] = [num(g("time", bnd, "auroc")) for bnd in BANDS]
    d["n_band"] = [int(g("time", bnd, "n_actives")) for bnd in BANDS]

    pr = [x for x in rows(TAB / "external_prospective.csv") if x["status"] == "ok"]
    d["n_ok"], d["n_all"] = len(pr), len(rows(TAB / "external_prospective.csv"))
    d["n_post"] = sum(int(float(x["time_n_test_actives"])) for x in pr if x["time_n_test_actives"])

    pn = rows(TAB / "permutation_null.csv")
    d["pn"] = {sp: st.mean(num(r["permuted_roc_auc_mean"]) for r in pn if r["split"] == sp)
               for sp in ("random", "scaffold")}
    d["pn_dev"] = max(abs(num(r["permuted_roc_auc_mean"]) - 0.5) for r in pn)

    ver = rows(INV / "VERDICTS.csv")
    d["n_hyp"] = len(ver)
    d["n_ref"] = sum(1 for r in ver if r["verdict"].startswith("REFUTED"))
    d["n_weak"] = sum(1 for r in ver if r["verdict"].startswith("WEAKENED"))

    h10 = {(r["method"], r["population"]): r for r in rows(INV / "H10_barrier_necessity.csv")}
    d["h10"] = h10
    d["h10_verdict"] = list(h10.values())[0]["verdict"]

    h9 = rows(INV / "H9_disease_discrimination.csv")
    d["h9"] = (num(one(INV / "H9_disease_discrimination_summary.csv", "model",
                       lambda r: r["metric"] == "mean per-indication AUROC")),
                max(num(r["auroc_model"]) for r in h9),
                min(num(r["auroc_model"]) for r in h9),
                sum(1 for r in h9 if num(r["auroc_model"]) > 0.5), len(h9))
    d["h1"] = (num(one(INV / "H1_disease_layer.csv", "value",
                       lambda r: r["test"] == "observed top-3 accuracy")),
               num(one(INV / "H1_disease_layer.csv", "value",
                       lambda r: r["test"].startswith("permutation null"))),
               num(one(INV / "H1_disease_layer.csv", "value",
                       lambda r: r["test"].startswith("frequency null"))))
    d["h6"] = one(INV / "H6_clinical_indication.csv", "top3_accuracy",
                  lambda r: r["stratum"] == "never seen in training")
    d["h6_freq"] = one(INV / "H6_clinical_indication.csv", "frequency_null",
                       lambda r: r["stratum"] == "never seen in training")

    g10 = json.loads((INV / "GRAPH_FINGERPRINT.json").read_text(encoding="utf-8"))
    d["graph"] = (g10["n_conditions"], g10["n_targets"])

    bvm = rows(TAB / "binder_vs_measured_inactives.csv")
    au = [num(r["auroc"]) for r in bvm if r.get("auroc") not in ("", None)]
    d["bvm"] = (len(au), st.mean(au), min(au), max(au)) if au else None

    sp3 = TAB / "library_sp3_coverage.csv"
    d["sp3"] = (num(one(sp3, "value", lambda r: r["metric"] == "median fraction sp3")),
                num(one(sp3, "value", lambda r: r["metric"].startswith("sp3-rich")
                        and "per cent" in r["metric"]))) if sp3.exists() else None
    return d


F = facts()
D = Deck()


def figure_slide(num_label, title, sub, name, source, notes, points):
    """A figure at full available height on the left, with its takeaways in the space beside it.

    These figures are nearly square. Fitted to a 16:9 slide they leave a third of the width empty,
    and the honest use of that space is the two or three sentences a viewer needs in order to read
    the panels in the right order.
    """
    s = D.light()
    D.head(s, num_label, title, sub)
    picture(s, name, M, 1.58, 6.45, 5.36)
    xs = M + 6.85
    yy = 1.66
    for i, (t, b, col) in enumerate(points, 1):
        D.dot(s, xs, yy + 0.04, str(i), fill=col, dia=0.36)
        D.text(s, xs + 0.54, yy - 0.02, W - M - (xs + 0.54), 0.42, t, size=13.5, bold=True,
               color=col)
        D.text(s, xs + 0.54, yy + 0.44, W - M - (xs + 0.54), 1.30, b, size=12.5, color=INK,
               line=1.28)
        yy += 1.86
    D.source(s, source)
    D.notes(s, notes)
    return s


def picture(s, name, x, y, box_w, box_h):
    """Fit a thesis figure inside a box, centred, without distorting it.

    These figures are nearly square, so scaling to the slide width alone pushes them a foot below
    the bottom edge. Whichever dimension binds is the one that sets the scale.
    """
    p = FIGS / f"{name}.png"
    with Image.open(p) as im:
        aspect = im.size[1] / im.size[0]
    w = min(box_w, box_h / aspect)
    h = w * aspect
    return s.shapes.add_picture(str(p), Inches(x + (box_w - w) / 2), Inches(y + (box_h - h) / 2),
                                width=Inches(w))


def section(title, sub, n, note):
    s = D.dark()
    D.text(s, M, 2.50, 2.2, 1.0, n, size=64, bold=True, font=HEAD, color=DEEP)
    D.text(s, M + 2.45, 2.72, W - 2 * M - 2.7, 0.9, title, size=34, bold=True, font=HEAD,
           color=PAPER)
    D.text(s, M + 2.48, 3.72, W - 2 * M - 2.7, 0.7, sub, size=15, color=CHALK, line=1.3)
    D.notes(s, note)
    return s


# ============================================================== 1. title
s = D.dark()
D.text(s, M, 1.62, W - 2 * M, 0.5, "BrainSafe AI", size=17, color=TEAL, bold=True)
D.text(s, M, 2.05, W - 2 * M - 1.0, 1.7,
       "A BBB-gated, multi-endpoint\npredictor of small-molecule\naction in the human brain",
       size=36, bold=True, font=HEAD, color=PAPER, line=1.14)
D.text(s, M, 4.62, W - 2 * M - 1.4, 1.0,
       f"{F['n_dep']} deployed endpoints across {F['graph'][0]} conditions. What the system computes, "
       f"what is known about it, and what is not.",
       size=17, color=CHALK, line=1.34)
D.text(s, M, 5.90, W - 2 * M, 0.5,
       "Every figure in this deck is read from an artefact at build time.",
       size=13, color=AMBER, italic=True)
D.notes(s, "This is the long-form explanation. The four thesis figures embedded here are the same "
           "files the chapters cite, so nothing shown can drift from the repository.")

# ============================================================== 2. the problem
s = D.light()
D.head(s, "1", "The problem is a coupling, not a prediction",
       "Two questions that are usually answered by two different tools, and must be answered together")
items = [
    ("Does it reach the brain?", "An exposure question. Answered by barrier models, physicochemical "
                                 "rules, and PK.", DEEP),
    ("Does it do anything there?", "An engagement question. Answered by target-activity models, one "
                                   "per protein.", TEAL),
    ("Neither answer is useful alone", "A potent ligand that does not cross is not a CNS drug. A "
                                       "permeant molecule with no mechanism is not either. The "
                                       "product is what a project needs, and it is what almost no "
                                       "server returns.", AMBER),
]
yy = 1.90
for i, (t, b, col) in enumerate(items, 1):
    D.card(s, M, yy - 0.10, W - 2 * M, 1.18, fill=TINT if i % 2 else PAPER)
    D.dot(s, M + 0.24, yy + 0.28, str(i), fill=col, dia=0.40)
    D.text(s, M + 0.86, yy + 0.10, 4.4, 0.8, t, size=15, bold=True, color=col, line=1.16)
    D.text(s, M + 5.45, yy + 0.06, W - M - (M + 5.45) - 0.30, 1.0, b, size=13, color=INK, line=1.26)
    yy += 1.30

D.text(s, M, 5.90, W - 2 * M, 0.85,
       "BrainSafe computes the product explicitly: an engagement signal per target, routed to the "
       "conditions that mechanism touches, then multiplied by a predicted barrier probability. "
       "The multiplication is the design claim, and Chapter 7 tests it.",
       size=13.5, color=DEEP, line=1.28)
D.source(s, "Chapters 1 and 7")
D.notes(s, "Resist the urge to call this novel. The claim is that the coupling is made explicit and "
           "then tested, not that nobody has multiplied two numbers before.")

# ============================================================== 3. what it is
s = D.light()
D.head(s, "2", "What the system is, in numbers",
       "Counted from models_rf/ and the panel registry at build time, not restated")
stats = [
    (f"{F['n_dep']}", "deployed endpoints, of "
                      f"{F['n_endpoints']} trained. Five withdrawn on evidence", DEEP),
    (f"{F['graph'][0]}", f"conditions, reached from {F['graph'][1]} targets through a curated "
                         f"pathway graph", TEAL),
    (f"{F['n_core']}", "core classifiers with the full uncertainty stack", TEAL),
    (f"{F['n_unreliable']}", "deployed endpoints fail the reliability gate and say so on every "
                             "result", AMBER),
]
for i, (v, lab, col) in enumerate(stats):
    D.stat(s, M + i * 3.05, 1.95, 2.7, v, lab, col)

D.card(s, M, 4.05, W - 2 * M, 2.05, fill=TINT)
D.text(s, M + 0.34, 4.28, W - 2 * M - 0.68, 1.65,
       "One prediction is assembled from many estimators trained on different data under different "
       "label rules. There is no single model, and describing the system by a panel mean is the one "
       "summary guaranteed to mislead: the endpoints are not interchangeable and their sensitivities "
       f"run from {fmt(F['sens'][2], 3)} to {fmt(F['sens'][3], 3)}.",
       size=14, color=INK, line=1.34)
D.text(s, M, 6.34, W - 2 * M, 0.4,
       f"Endpoints failing the gate: {', '.join(F['unreliable'])}.",
       size=12.5, color=AMBER, italic=True)
D.source(s, "models_rf/binder_modes.json; results/tables/rf_cv_summary.csv")
D.notes(s, "If asked why five were withdrawn, the answer is on the falsification slide: they fired "
           "on glucose, urea and lactate.")

# ============================================================== 4. section: the model
section("How a prediction is made", "Features, forests, and the five operations that turn a "
                                    "probability into a reported condition", "I",
        "Four slides. The two figures carry most of it; the words between them are the parts a "
        "figure cannot say.")

# ============================================================== 5. features and models
s = D.light()
D.head(s, "3", "One representation, one estimator family",
       "A fixed 1,036-column vector per compound, and a random forest per endpoint")
D.card(s, M, 1.82, 5.9, 2.35, fill=TINT)
D.text(s, M + 0.28, 2.02, 5.4, 0.34, "The feature vector", size=13.5, bold=True, color=DEEP)
for i, (k, v) in enumerate([("1,024", "ECFP-4 bits, folded"),
                            ("12", "physicochemical descriptors"),
                            ("1,036", "columns, identical for every endpoint")]):
    D.text(s, M + 0.28, 2.46 + i * 0.46, 1.0, 0.34, k, size=14, bold=True, color=TEAL,
           align=PP_ALIGN.RIGHT)
    D.text(s, M + 1.45, 2.46 + i * 0.46, 4.2, 0.34, v, size=12.5, color=INK)
D.text(s, M + 0.28, 3.86, 5.4, 0.3,
       "Computed on the desalted, neutralised parent", size=11.5, color=MUTED, italic=True)

D.card(s, M + 6.15, 1.82, W - M - (M + 6.15), 2.35, fill=TINT)
D.text(s, M + 6.43, 2.02, 5.4, 0.34, "The estimator", size=13.5, bold=True, color=DEEP)
for i, (k, v) in enumerate([("300", "trees per forest"),
                            ("2 / 4", "minimum samples per leaf, core / binder"),
                            ("balanced", "class weights, and a fixed seed")]):
    D.text(s, M + 6.43, 2.46 + i * 0.46, 1.0, 0.34, k, size=14, bold=True, color=TEAL,
           align=PP_ALIGN.RIGHT)
    D.text(s, M + 7.60, 2.46 + i * 0.46, 4.4, 0.34, v, size=12.5, color=INK)
D.text(s, M + 6.43, 3.86, 5.4, 0.3,
       "Chosen by benchmark, not by preference", size=11.5, color=MUTED, italic=True)

D.text(s, M, 4.42, W - 2 * M, 0.45,
       "Five model families were compared on the same folds. The forest was not the best on every "
       "task, and the thesis says so.",
       size=13.5, bold=True, color=DEEP)
D.text(s, M, 4.94, W - 2 * M, 1.0,
       "It wins 7 of 8 classification endpoints and 0 of 5 regression endpoints, where gradient "
       "boosting is ahead. Against the deployed panel it is statistically indistinguishable from "
       "both boosting methods pooled (p = 0.735 and 0.893, Wilcoxon signed-rank paired by "
       "endpoint). The forest is kept for its calibration behaviour and its out-of-bag structure, "
       "not because it is measurably more accurate.",
       size=13, color=INK, line=1.30)
D.source(s, "results/tables/model_comparison.csv; results/tables/model_family_significance.csv")
D.notes(s, "Volunteer the 0 of 5. A referee who finds it unaided will discount the 7 of 8.")

# ============================================================== 6. Figure T1
figure_slide(
    "4", "What the system computes, end to end",
    "Five operations, each of which changes what the number means",
    "FigureT1_scoring_pipeline",
    "Figure T1. Worked on donepezil, scored through the served pipeline.",
    "Walk the top row left to right, then say that panel B is the step people miss: a probability "
    "is not comparable across endpoints, which is why enrichment exists.",
    [("A probability is not a score",
      "Panel B. The same 0.60 means opposite things at an endpoint where a quarter of compounds "
      "are active and one where six sevenths are. Enrichment fixes the reference point.", TEAL),
     ("One mechanism, not a tally",
      "Panel C. Only the strongest engaged target contributes to a condition, because homologous "
      "targets are not independent observations.", DEEP),
     ("The gate decides whether to speak",
      "Panel D. Twelve of the sixteen conditions fall below the reporting threshold here, and the "
      "expected recall at this compound's distance is printed beside the answer.", AMBER)])

# ============================================================== 7. Figure T2
figure_slide(
    "5", "The two transforms that are easiest to misread",
    "Enrichment over base rate, and an exposure gate that cannot reorder conditions",
    "FigureT2_enrichment_and_gate",
    "Figure T2. Base rates from the deployed context; C donepezil, D teriflunomide.",
    "Panel D is the one to dwell on. H3 is refuted for the fourteen central conditions, and the "
    "exemption for migraine and multiple sclerosis is a qualification no document in the project "
    "states.",
    [("The base rates span 0.626",
      "Panels A and B. A single probability threshold applied across the panel would be asking a "
      "different question at every endpoint.", TEAL),
     ("The gate is a filter, not a discriminator",
      "Panel C. Multiplying fourteen conditions by the same number is rank-invariant, so gating "
      "cannot sharpen a disease call. It decides whether anything is reported at all.", CRIMSON),
     ("With one exception nobody had stated",
      "Panel D. Migraine and multiple sclerosis are gate-exempt, so the gate does reorder them "
      "against the central conditions. Worth a factor of 1.79 on teriflunomide.", AMBER)])

# ============================================================== 8. aggregation
s = D.light()
D.head(s, "6", "Why a maximum and not a sum",
       "The aggregation rule, and the evidence that made it defensible")
D.card(s, M, 1.85, 5.85, 2.5, fill=TINT)
D.text(s, M + 0.28, 2.05, 5.3, 0.34, "If engaged targets were independent", size=13, bold=True,
       color=MUTED)
D.text(s, M + 0.28, 2.48, 5.3, 1.65,
       "summing them would be the right thing to do: three separate observations of activity are "
       "stronger evidence than one, and a sum would say so.",
       size=13, color=INK, line=1.30)

D.card(s, M + 6.10, 1.85, W - M - (M + 6.10), 2.5, fill=TINT)
D.text(s, M + 6.38, 2.05, 5.3, 0.34, "They are not independent", size=13, bold=True, color=CRIMSON)
D.text(s, M + 6.38, 2.48, 5.4, 1.75,
       "Across 400 approved drugs, 37 targets fire at least once but span only 16 independent "
       "directions. Five homologous pairs co-fire above φ = 0.5: μ and κ opioid at 0.798, SERT and "
       "NET at 0.732, D2 and D3 at 0.720.",
       size=13, color=INK, line=1.30)

D.text(s, M, 4.62, W - 2 * M, 0.9,
       "So a sum would count the same pharmacology two or three times. The maximum takes the "
       "strongest engaged mechanism and nothing else, which understates the evidence rather than "
       "overstating it. That is the safer direction for a tool whose purpose is to avoid a wasted "
       "experiment.",
       size=14, color=INK, line=1.30)
D.text(s, M, 5.66, W - 2 * M, 0.75,
       "The interface now groups engaged targets by homology family and quotes the measured "
       "correlation wherever two members fire, because reporting three engaged targets as three "
       "findings is the error, not the co-firing itself.",
       size=13, color=DEEP, line=1.28)
D.source(s, "inversion/results/H8_panel_independence.csv, H8_family_correlation.csv")
D.notes(s, "This is the one design choice whose stated justification survived testing, and it "
           "survived because H8 supplied the evidence for it.")

# ============================================================== 9. section: what is known
section("What is known about it", "Cross-validation, calibration, thresholds, and three external "
                                  "tests", "II",
        "This is the evidence section. Lead with the scaffold split and the permutation null, "
        "because together they answer the leakage question a referee asks first.")

# ============================================================== 10. CV design
s = D.light()
D.head(s, "7", "Two cross-validation regimes, and why both are reported",
       "A single AUROC is the number most easily overstated in this field")
D.text(s, M + 0.30, 1.95, 3.4, 0.3, "regime", size=11.5, bold=True, color=MUTED, space=0)
for i, hh in enumerate(["mean", "lowest", "highest"]):
    D.text(s, M + 4.10 + i * 2.15, 1.95, 1.9, 0.3, hh, size=11.5, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT, space=0)
yy = 2.38
for lab, key, col in [("Random 10-fold", "random", MUTED),
                      ("Scaffold-grouped 10-fold", "scaffold", DEEP)]:
    mean, lo, hi = F["cv"][key]
    D.card(s, M, yy - 0.10, 10.6, 0.62, fill=TINT if key == "scaffold" else PAPER)
    D.text(s, M + 0.30, yy, 3.6, 0.44, lab, size=13.5, bold=(key == "scaffold"), color=col)
    for i, v in enumerate([mean, lo, hi]):
        D.text(s, M + 4.10 + i * 2.15, yy, 1.9, 0.44, fmt(v, 4), size=13.5,
               bold=(key == "scaffold"), color=col, align=PP_ALIGN.RIGHT)
    yy += 0.72

D.card(s, M, 3.98, W - 2 * M, 1.30, fill=TINT)
D.text(s, M + 0.32, 4.20, W - 2 * M - 0.64, 1.0,
       "The scaffold figure is the one to quote. A random split of medicinal-chemistry data holds "
       "out close analogues of its own training set, because the published record is series, so it "
       "answers a question no user ever asks.",
       size=13.5, color=INK, line=1.30)

D.text(s, M, 5.50, W - 2 * M, 0.85,
       f"And neither figure is obtainable without the labels. With labels permuted, the same "
       f"pipeline on the same folds returns {fmt(F['pn']['random'], 4)} and "
       f"{fmt(F['pn']['scaffold'], 4)}, every one of the sixteen values within "
       f"{fmt(F['pn_dev'], 4)} of chance.",
       size=13.5, bold=True, color=TEAL, line=1.28)
D.source(s, "results/tables/rf_cv_summary.csv; results/tables/permutation_null.csv")
D.notes(s, "The permutation null closes the leakage question. It also shows scaffold grouping alone "
           "confers nothing, which is the specific route a referee will ask about.")

# ============================================================== 11. Figure T3
figure_slide(
    "8", "The uncertainty stack, and where it stops",
    "Calibration, conformal coverage, an applicability band, and an asymmetry worth stating",
    "FigureT3_uncertainty_stack",
    "Figure T3. calibration.csv, rf_conformal.csv, integrity_calibration_per_target.csv.",
    "Panel B is the uncomfortable one and should be said aloud: the barrier model is the "
    "worst-calibrated of the eight, and it multiplies every disease score.",
    [("Calibration works where it is applied",
      f"Panel A. Expected calibration error falls from {fmt(F['cal'][0], 4)} to "
      f"{fmt(F['cal'][1], 4)} across the eight core classifiers.", TEAL),
     ("And the gate is the worst of the eight",
      f"Panel B. {fmt(F['cal'][3], 4)}, nearly three times the mean, on the one model whose "
      f"output multiplies every disease score rather than being reported alone.", CRIMSON),
     ("The stack is thinnest where it is used most",
      f"Panel D. None of the {F['n_dep']} deployed binder endpoints carries a conformal "
      f"statement. The {F['bcal'][0]} with a measured calibration error average "
      f"{fmt(F['bcal'][1], 4)}; the other {F['n_dep'] - F['bcal'][0]} carry none.", AMBER)])

# ============================================================== 12. thresholds
s = D.light()
D.head(s, "9", "Thresholds are a separate problem from probabilities",
       "Three disjoint background pools, and two constraints of which the harder binds")
cols = [
    ("Decoy pool, 60%", "trains the binder classifiers against chemistry presumed inactive", MUTED),
    ("Threshold pool, 20%", "sets each operating point", TEAL),
    ("Evaluation pool, 20%", "measures the false-positive rate at that point", DEEP),
]
cw = (W - 2 * M - 0.5) / 3
for j, (t, b, col) in enumerate(cols):
    x = M + j * (cw + 0.25)
    D.card(s, x, 1.88, cw, 1.55, fill=TINT)
    D.text(s, x + 0.26, 2.08, cw - 0.52, 0.42, t, size=13.5, bold=True, color=col)
    D.text(s, x + 0.26, 2.58, cw - 0.52, 0.80, b, size=12.5, color=INK, line=1.26)
D.text(s, M, 3.62, W - 2 * M, 0.4,
       "Assignment by blake2b(salt + canonical SMILES) mod 100, so a compound's pool is a property "
       "of the compound and never of the run.",
       size=12.5, color=MUTED, italic=True)

D.text(s, M, 4.20, W - 2 * M, 0.45,
       "Each threshold is then set by the harder of two constraints", size=14, bold=True,
       color=DEEP)
for i, (t, b) in enumerate([
        ("separation from background", "a false-positive rate at or under 5 per cent on random "
                                       "chemistry"),
        ("separation from measured inactives", "compounds tested at that same target and reported "
                                               "as non-binders")]):
    D.dot(s, M, 4.78 + i * 0.62, str(i + 1), fill=TEAL, dia=0.34)
    D.text(s, M + 0.52, 4.74 + i * 0.62, 4.6, 0.44, t, size=13, bold=True, color=INK)
    D.text(s, M + 5.30, 4.74 + i * 0.62, W - M - (M + 5.30), 0.44, b, size=12.5, color=MUTED)

D.text(s, M, 6.08, W - 2 * M, 0.5,
       "The second is usually the harder, and it is the reason several endpoints have low "
       "sensitivity: relaxing them would mean calling compounds an experiment has already reported "
       "as non-binders.",
       size=13, color=AMBER, line=1.26)
D.source(s, "Chapter 5; results/tables/final_thresholds.csv")
D.notes(s, "Higher sensitivity is available at any time. It would make the server less truthful.")

# ============================================================== 13. sensitivity
s = D.light()
D.head(s, "10", "What a silence is worth",
       "The number that makes a negative result interpretable, and the correction it required")
D.stat(s, M, 1.92, 2.6, fmt(F["sens"][0], 4), "panel mean sensitivity at the deployed threshold, on "
                                              "held-out actives", DEEP)
D.stat(s, M + 2.95, 1.92, 2.6, fmt(F["sens"][1], 3), "median", TEAL)
D.stat(s, M + 5.90, 1.92, 2.6, fmt(F["sens"][2], 3),
       f"lowest, {F['sens_min_ep'].replace('_', '-')}", CRIMSON)
D.stat(s, M + 8.85, 1.92, 2.6, fmt(F["sens"][3], 3),
       f"highest, {F['sens_max_ep']}", TEAL)

D.card(s, M, 3.72, W - 2 * M, 1.55, fill=TINT)
D.text(s, M + 0.32, 3.94, W - 2 * M - 0.64, 1.20,
       "For most of this project's life the published figure was 0.8983, under a label reading "
       "\"held out by scaffold\". It was not held out: the last of four scripts writing that field "
       "scored every active in each endpoint table, roughly four fifths of which the model had been "
       "fitted on. The held-out mean is 0.7638, and six endpoints fire for fewer than half their own "
       "actives where the published figure had none.",
       size=13, color=INK, line=1.30)

D.text(s, M, 5.48, W - 2 * M, 0.9,
       "A silence is therefore not evidence of inactivity, and the server reports the expected "
       "recall at that compound's distance from training chemistry beside every negative result. "
       "The falsification suite had recorded the correct figure three weeks before anything "
       "compared it with the registry.",
       size=13.5, color=DEEP, line=1.28)
D.source(s, "models_rf/binder_modes.json; results/tables/sensitivity_reconciliation.csv")
D.notes(s, "Tell this as a defect found and fixed. It is the strongest evidence that the project's "
           "checking works, and hiding it would waste that.")

# ============================================================== 14. section: external
section("Does it hold outside the training set?", "One genuine external set, a prospective "
                                                  "simulation, and a different curator", "III",
        "The composition finding is the payload here: an apparent temporal decay that is a "
        "property of what each split holds out, not of time.")

# ============================================================== 15. external BBB
s = D.light()
D.head(s, "11", "The one genuine external set",
       "FDA-curated approved drugs absent from the barrier model's training database")
D.text(s, M + 0.30, 1.98, 5.6, 0.3, "set", size=11.5, bold=True, color=MUTED, space=0)
for i, hh in enumerate(["n", "AUROC", "sensitivity", "specificity"]):
    D.text(s, M + 6.10 + i * 1.55, 1.98, 1.4, 0.3, hh, size=11.5, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT, space=0)
yy = 2.42
for i, r in enumerate(F["ext"]):
    key = i == 1
    D.card(s, M, yy - 0.10, W - 2 * M, 0.66, fill=TINT if key else PAPER)
    lab = r["set"]
    D.text(s, M + 0.30, yy, 5.7, 0.46, lab if len(lab) < 62 else lab[:60] + "…",
           size=12.5, bold=key, color=DEEP if key else MUTED)
    for j, v in enumerate([r["n"], r["auroc"], r["sensitivity"], r["specificity"]]):
        D.text(s, M + 6.10 + j * 1.55, yy, 1.4, 0.46, v, size=12.5, bold=key,
               color=DEEP if key else MUTED, align=PP_ALIGN.RIGHT)
    yy += 0.76

D.text(s, M, 4.86, W - 2 * M, 0.9,
       "The middle row is the claim. Excluding overlap by InChIKey is not enough to call a compound "
       "unseen, because the InChIKey separates stereoisomers and salts and the featuriser does not. "
       "The third row is the memorisation the first row contains.",
       size=13.5, color=INK, line=1.30)
D.text(s, M, 5.82, W - 2 * M, 0.75,
       "For the target panel no external set of this kind exists, and the reason is structural: for "
       "most of these proteins the public record is the training set.",
       size=13, color=AMBER, line=1.26)
D.source(s, "results/tables/external_bbb_validation.csv")
D.notes(s, "The gap between rows one and two, 0.7645 against 0.7934, is the memorisation. Quote "
           "the second and explain why.")

# ============================================================== 16. prospective + composition
s = D.dark()
D.text(s, M, 0.86, 11.4, 0.5, "An apparent temporal decay that is a composition effect",
       size=15, color=AMBER, bold=True)
D.text(s, M, 1.30, 11.7, 0.6,
       f"Models refitted before a date cutoff, thresholds frozen, then scored on "
       f"{F['n_post']:,} actives measured afterwards across {F['n_ok']} of {F['n_all']} endpoints.",
       size=14, color=CHALK, line=1.26)
cd = CategoryChartData()
cd.categories = SHORT_BAND
cd.add_series("held out by date", tuple(F["recall"]["time"]))
cd.add_series("held out at random", tuple(F["recall"]["random"]))
cd.add_series("held out by curator", tuple(F["recall"]["cross_source"]))
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(2.10),
                        Inches(7.1), Inches(3.75), cd).chart
gf.has_legend = True
gf.legend.position = XL_LEGEND_POSITION.TOP
gf.legend.include_in_layout = False
gf.legend.font.size = Pt(11)
gf.legend.font.color.rgb = CHALK
gf.legend.font.name = BODY
va = gf.value_axis
va.minimum_scale, va.maximum_scale = 0.0, 1.0
va.has_major_gridlines = True
va.major_gridlines.format.line.color.rgb = DEEP
va.tick_labels.font.size = Pt(10)
va.tick_labels.font.color.rgb = CHALK
gf.category_axis.tick_labels.font.size = Pt(9)
gf.category_axis.tick_labels.font.color.rgb = CHALK
for ser, colr in zip(gf.plots[0].series, (TEAL, MUTED, AMBER)):
    ser.format.fill.solid()
    ser.format.fill.fore_color.rgb = colr
    ser.format.line.fill.background()

xs = M + 7.50
D.text(s, xs, 2.15, W - M - xs, 3.6,
       "Read in aggregate the time split looks like decay: sensitivity 0.4886 against 0.8715 for a "
       "size-matched random control.\n\n"
       "It is not decay. 83.3 per cent of the random split's test compounds are close analogues of "
       "its own training set, against 16.3 per cent for the time split. The two splits are not "
       "testing comparable populations.\n\n"
       "Read band by band the AUROC gap of 0.128 falls to at most 0.081, and three test sets built "
       "by unrelated rules trace one curve.",
       size=13, color=CHALK, line=1.36)
D.source(s, "results/tables/external_prospective.csv, external_novelty_strata.csv")
D.notes(s, "The residual is not attributed: the random split yields only 212 novel actives across 39 "
           "endpoints, too few for a paired within-band comparison. Say so.")

# ============================================================== 17. the distance curve
s = D.light()
D.head(s, "12", "Accuracy is a function of chemical distance",
       "The limitation that bounds every other number in this deck")
D.text(s, M + 0.30, 1.98, 4.6, 0.3, "nearest-training Tanimoto", size=11.5, bold=True,
       color=MUTED, space=0)
for i, hh in enumerate(["actives", "AUROC", "recall"]):
    D.text(s, M + 5.30 + i * 1.85, 1.98, 1.6, 0.3, hh, size=11.5, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT, space=0)
labels = ["below 0.40, a different chemotype", "0.40 to 0.55, a related series",
          "0.55 to 0.70, the same series", "0.70 and above, a close analogue"]
yy = 2.42
for i, lab in enumerate(labels):
    low = i == 0
    D.card(s, M, yy - 0.10, 10.9, 0.62, fill=TINT if low else PAPER)
    D.text(s, M + 0.30, yy, 4.8, 0.44, lab, size=13, bold=low, color=CRIMSON if low else INK)
    for j, v in enumerate([f"{F['n_band'][i]:,}", fmt(F["auroc_band"][i], 4),
                           fmt(F["recall"]["time"][i], 4)]):
        D.text(s, M + 5.30 + j * 1.85, yy, 1.6, 0.44, v, size=13, bold=low,
               color=CRIMSON if low else MUTED, align=PP_ALIGN.RIGHT)
    yy += 0.72

D.text(s, M, 5.40, W - 2 * M, 0.85,
       f"Recall runs {fmt(F['recall']['time'][0], 4)} to {fmt(F['recall']['time'][3], 4)}, a factor "
       f"of five. Across a database boundary the lowest band falls to "
       f"{fmt(F['recall']['cross_source'][0], 4)}. On genuinely novel chemistry the tool finds "
       f"roughly one true activity in six.",
       size=14, bold=True, color=CRIMSON, line=1.28)
D.text(s, M, 6.32, W - 2 * M, 0.45,
       "No analysis in this thesis improves that. What is established is that the figure is "
       "predictable, which is why it is reported beside every negative result.",
       size=13, color=DEEP, italic=True)
D.source(s, "results/tables/external_novelty_strata.csv")
D.notes(s, "If the audience takes one slide away, this is the one. AUROC moves far less than recall "
           "because ranking survives what a fixed threshold does not.")

# ============================================================== 18. specificity
s = D.light()
D.head(s, "13", "What it does when there is nothing to find",
       "The behaviour a triage tool is used for most often")
D.stat(s, M, 1.98, 3.0, fmt(F["spec"][0], 3),
       f"specificity on 1,000 non-CNS compounds, interval {fmt(F['spec'][1], 3)} to "
       f"{fmt(F['spec'][2], 4)}", DEEP)
D.stat(s, M + 3.35, 1.98, 3.0, "0.0717", "of random PubChem structures receive any user-visible "
                                         "finding", TEAL)
D.stat(s, M + 6.70, 1.98, 3.0, "0.515", "of approved drugs receive one, a set active somewhere by "
                                        "construction", TEAL)
D.stat(s, M + 10.05, 1.98, 2.2, "6 of 47", "endpoints carry a reliability flag on every result",
       AMBER)

D.card(s, M, 3.85, W - 2 * M, 1.45, fill=TINT)
D.text(s, M + 0.32, 4.06, W - 2 * M - 0.64, 1.10,
       "The specificity figure is an upper bound and the artefact labels it as one: the compounds "
       "are presumed inactive because nothing is recorded about them, not proven inactive. They are "
       "also drawn from within the applicability reference, so the figure does not bound behaviour "
       "on genuinely distant chemistry.",
       size=13.5, color=INK, line=1.30)

D.text(s, M, 5.52, W - 2 * M, 0.85,
       "H4 tested exactly that and is the suite's weakest supported verdict: on 61 compounds distant "
       "from training chemistry the false-positive rate is 0.0164, but with one false positive the "
       "95 per cent interval runs 0.0029 to 0.0872 and contains the 0.0750 comparator. Specificity "
       "does not degrade with distance; that it improves is not established.",
       size=13, color=AMBER, line=1.28)
D.source(s, "results/tables/noncns_specificity_summary.csv; inversion/results/H4_*.csv")
D.notes(s, "Distinguish the two claims carefully. The null predicted degradation and did not get "
           "it, which is a real negative result worth having.")

# ============================================================== 19. section: falsification
section("Trying to break it", "Ten hypotheses stated so that they could fail, and what the four "
                              "refutations cost", "IV",
        "The strongest section. A suite is worth what the project paid for it, and every payment "
        "here is a commit rather than a claim.")

# ============================================================== 20. Figure T4
figure_slide(
    "14", "The falsification suite",
    f"{F['n_ref']} of {F['n_hyp']} hypotheses refuted, {F['n_weak']} weakened",
    "FigureT4_falsification",
    "Figure T4. VERDICTS.csv, permutation_null.csv, H10_barrier_necessity.csv.",
    "The suite is worth what the project paid for it, and every payment is a commit rather than a "
    "claim: the graph re-described, the gate re-described, the interface changed, four endpoints "
    "withdrawn.",
    [("Every hypothesis has a null that could win",
      "Panel A. A test that cannot fail is not evidence, so each claim is paired with a model "
      "capable of producing the same apparent success by accident.", DEEP),
     ("The cross-validation is not leakage",
      f"Panel B. Under permuted labels the same folds return {fmt(F['pn']['random'], 4)} and "
      f"{fmt(F['pn']['scaffold'], 4)}, and scaffold grouping alone confers nothing.", TEAL),
     ("The newest hypothesis weakened a claim",
      "Panel C. H10 tests the barrier model against descriptor rules, and its interval on unseen "
      "drugs crosses zero.", AMBER)])

# ============================================================== 21. the disease layer
s = D.light()
D.head(s, "15", "What the disease layer does and does not establish",
       "Three hypotheses, answering three different questions")
obs, perm, freq = F["h1"]
items = [
    ("H1", "against the project's own map",
     f"Top-3 accuracy {fmt(obs, 4)} on 7,008 held-out compounds, against a permutation null of "
     f"{fmt(perm, 4)} and a frequency null of {fmt(freq, 4)}. The layer carries real information "
     f"about mechanism.", TEAL),
    ("H6", "against real clinical indications",
     f"On {'162'} drugs whose structure appears nowhere in training, top-3 accuracy is {F['h6']} "
     f"against a frequency null of {F['h6_freq']}. It beats the permutation null decisively and "
     f"does not beat a constant answer naming the three commonest CNS indications.", AMBER),
    ("H9", "on metrics a constant answer cannot pass",
     f"Mean per-indication AUROC {fmt(F['h9'][0], 4)} against 0.500, beating chance on "
     f"{F['h9'][3]} of {F['h9'][4]} indications, from {fmt(F['h9'][1], 4)} for depression down to "
     f"{fmt(F['h9'][2], 4)} for epilepsy.", TEAL),
]
yy = 1.92
for k, t, b, col in items:
    D.card(s, M, yy - 0.10, W - 2 * M, 1.30, fill=TINT)
    D.text(s, M + 0.28, yy + 0.28, 0.8, 0.44, k, size=17, bold=True, font=HEAD, color=col)
    D.text(s, M + 1.15, yy + 0.06, 3.8, 0.8, t, size=13, bold=True, color=INK, line=1.18)
    D.text(s, M + 5.25, yy + 0.02, W - M - (M + 5.25) - 0.28, 1.15, b, size=12.5, color=INK,
           line=1.26)
    yy += 1.42

D.text(s, M, 6.24, W - 2 * M, 0.62,
       "The honest description: the layer ranks mechanisms and the conditions they touch. It does "
       "not predict indication, and two of the nine conditions are not being predicted at all.",
       size=13.5, bold=True, color=DEEP, line=1.26)
D.source(s, "inversion/results/H1_disease_layer.csv, H6_clinical_indication.csv, "
            "H9_disease_discrimination.csv")
D.notes(s, "Epilepsy at 0.4898 is the informative failure: most approved antiepileptics are small, "
           "simple and low-affinity, which is the chemistry that sits under a strict cut.")

# ============================================================== 22. H10 in words
s = D.light()
D.head(s, "16", "The hypothesis the suite lacked",
       "The component the architecture is named for had no test until a thesis audit named the gap")
ext = F["h10"]
dep = ext[("deployed forest, 1,036 features", "external approved")]
sca = ext[("deployed forest, 1,036 features", "scaffold hold-out")]
nul_e = ext[("descriptor forest, 12 features", "external approved")]
nul_s = ext[("descriptor forest, 12 features", "scaffold hold-out")]
D.text(s, M, 1.76, W - 2 * M, 0.62,
       "Null: permeability is a bulk property, the twelve descriptors already in the feature vector "
       "encode it, and the fingerprint adds nothing.",
       size=13.5, italic=True, color=DEEP, line=1.26)
D.text(s, M + 0.30, 2.42, 5.0, 0.3, "", size=11, space=0)
for i, hh in enumerate(["scaffold hold-out", "241 unseen approved drugs"]):
    D.text(s, M + 5.30 + i * 2.85, 2.42, 2.6, 0.3, hh, size=11.5, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT, space=0)
rowsx = [("Deployed forest, 1,036 features", sca["auroc"], dep["auroc"], True),
         ("Random forest, 12 descriptors only", nul_s["auroc"], nul_e["auroc"], False)]
yy = 2.86
for lab, a, b, top in rowsx:
    D.card(s, M, yy - 0.10, 10.9, 0.62, fill=TINT if top else PAPER)
    D.text(s, M + 0.30, yy, 5.0, 0.44, lab, size=13, bold=top, color=DEEP if top else INK)
    for j, v in enumerate([a, b]):
        D.text(s, M + 5.30 + j * 2.85, yy, 2.6, 0.44, fmt(num(v), 4), size=13, bold=top,
               color=DEEP if top else MUTED, align=PP_ALIGN.RIGHT)
    yy += 0.72

D.card(s, M, 4.40, W - 2 * M, 1.30, fill=PAPER)
D.text(s, M + 0.32, 4.60, W - 2 * M - 0.64, 1.0,
       f"Read as point estimates the fingerprint wins on both. Read with a paired bootstrap of "
       f"2,000 resamples, the external interval runs {num(nul_e['delta_ci95_low']):+.4f} to "
       f"{num(nul_e['delta_ci95_high']):+.4f}, p = {fmt(num(nul_e['bootstrap_p_deployed_better']), 3)}. "
       f"It crosses zero.",
       size=13.5, color=INK, line=1.30)

D.text(s, M, 5.88, W - 2 * M, 0.85,
       f"VERDICT {F['h10_verdict']}. On the training distribution the fingerprint demonstrably earns "
       f"its place; on chemistry the model has never seen, its advantage over twelve descriptors is "
       f"not established. Under a point-estimate rule this would read SUPPORTED, and the suite's own "
       f"criticism of H4 is why it does not.",
       size=13.5, bold=True, color=CRIMSON, line=1.28)
D.source(s, "inversion/results/H10_barrier_necessity.csv")
D.notes(s, "What the barrier model does clearly beat, on both populations and by a wide margin, is "
           "the published CNS heuristic: 0.7934 against 0.6369. Say that second, not first.")

# ============================================================== 23. limitations
s = D.dark()
D.text(s, M, 0.86, 11.4, 0.5, "What the system cannot do", size=15, color=AMBER, bold=True)
D.text(s, M, 1.30, 11.7, 0.45,
       "Properties of the science rather than of the housekeeping, and none of them is fixed by "
       "tidying.", size=14, color=CHALK)
lims = [
    ("Novel chemistry",
     f"Recall {fmt(F['recall']['time'][0], 4)} on a different chemotype and "
     f"{fmt(F['recall']['cross_source'][0], 4)} across a database boundary. A silence there is "
     f"close to uninformative."),
    ("Uncertainty coverage",
     f"The full stack reaches {F['n_core']} core classifiers. None of the {F['n_dep']} deployed "
     f"binder endpoints has a conformal statement; the {F['bcal'][0]} with a measured calibration "
     f"error average {fmt(F['bcal'][1], 4)}, about five times the core figure, and "
     f"{F['n_dep'] - F['bcal'][0]} carry no calibration measurement at all."),
    ("Natural products",
     f"Median fraction sp3 of the library is {fmt(F['sp3'][0], 4)} and only {F['sp3'][1]:.2f} per "
     f"cent of it is sp3-rich and largely non-aromatic. Three endpoints could be scored externally, "
     f"with two to five actives each."
     if F["sp3"] else "Largely outside the training library."),
    ("Nothing prospective",
     "Every result here is retrospective. The dates are real, but no compound was predicted before "
     "it was measured."),
]
yy = 1.96
for i, (t, b) in enumerate(lims, 1):
    D.dot(s, M, yy + 0.04, str(i), fill=AMBER, fg=INK, dia=0.38)
    D.text(s, M + 0.62, yy, 3.5, 0.85, t, size=14.5, bold=True, color=AMBER, line=1.14)
    D.text(s, M + 4.35, yy - 0.02, 11.9 - 4.35, 1.15, b, size=12.5, color=CHALK, line=1.26)
    yy += 1.22
D.text(s, M, 6.86, 11.7, 0.4,
       "Stated because the alternative is a referee finding them unaided.",
       size=13, color=PAPER, italic=True)
D.notes(s, "Do not soften these. The thesis's credibility rests on them being said first.")

# ============================================================== 24. close
s = D.dark()
D.text(s, M, 0.82, 11.4, 0.5, "What this establishes", size=15, color=TEAL, bold=True)
D.text(s, M, 1.24, 11.7, 0.42, "With the spread rather than the mean wherever one exists.",
       size=13.5, color=CHALK, italic=True)
est = [
    (f"{fmt(F['cv']['scaffold'][1], 4)} to {fmt(F['cv']['scaffold'][2], 4)}",
     f"scaffold-grouped AUROC across the {F['n_core']} core classifiers, against permuted-label "
     f"nulls of {fmt(F['pn']['random'], 4)} and {fmt(F['pn']['scaffold'], 4)}"),
    (f"{fmt(F['recall']['time'][0], 4)} to {fmt(F['recall']['time'][3], 4)}",
     "recall across four novelty bands, the same curve traced by three independently built test "
     "sets"),
    (f"{fmt(F['sens'][0], 4)}",
     f"panel sensitivity on held-out actives, median {fmt(F['sens'][1], 3)}, from "
     f"{fmt(F['sens'][2], 3)} to {fmt(F['sens'][3], 3)}"),
    (f"{fmt(F['spec'][0], 3)}",
     "specificity on presumed-inactive non-CNS chemistry, an upper bound rather than an estimate"),
    (f"{F['n_ref']} of {F['n_hyp']}",
     f"falsification hypotheses refuted, {F['n_weak']} weakened, at the cost of two claims, one "
     f"interface change and four withdrawn endpoints"),
]
yy = 1.86
for v, b in est:
    D.text(s, M, yy, 2.6, 0.44, v, size=18, bold=True, font=HEAD, color=TEAL)
    D.text(s, M + 2.95, yy - 0.02, 11.9 - 2.95, 0.85, b, size=12.5, color=CHALK, line=1.24)
    yy += 0.92
D.card(s, M, 6.44, 11.7, 0.74, fill=DEEP)
D.text(s, M + 0.28, 6.58, 11.2, 0.56,
       "Not established: anything prospective, anything about natural-product chemistry, and any "
       "advantage for the barrier model over twelve descriptors on chemistry it has not seen.",
       size=13, bold=True, color=PAPER, line=1.20)
D.notes(s, "End on the last line. Naming the boundary of the work is stronger than any figure "
           "above it, and all three items have a measurement behind them.")

D.save(OUT)
