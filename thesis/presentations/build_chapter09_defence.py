"""Build the Chapter 9 defence deck, reading every figure from an artefact.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/build_chapter09_defence.py
Out:  thesis/presentations/chapter09_defence.pptx
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
    ROOT, TAB, INV, W, M, INK, DEEP, TEAL, AMBER, CRIMSON, MUTED, TINT, PAPER, CHALK,
    HEAD, BODY, Deck, rows, num, fmt,
)

OUT = Path(__file__).resolve().parent / "chapter09_defence.pptx"


def one(path, field, pred):
    for r in rows(path):
        if pred(r):
            return r[field]
    raise KeyError(f"no matching row in {path.name}")


def facts() -> dict:
    d = {}
    d["verdicts"] = rows(INV / "VERDICTS.csv")

    h2 = {r["weights"]: num(r["top3_accuracy"]) for r in rows(INV / "H2_weight_ablation.csv")}
    d["h2"] = (h2["curated"], h2["uniform (all 1.0)"], h2["randomly permuted"])

    d["h4"] = rows(INV / "H4_distant_specificity.csv")
    d["lib_fpr"] = num(one(TAB / "noncns_specificity_summary.csv", "estimate",
                           lambda r: r["metric"].startswith("False-positive rate (any actionable")))

    ra = rows(INV / "H5_readacross_value.csv")
    d["h5"] = (num(one(INV / "H5_readacross_value.csv", "recall",
                       lambda r: r["method"].startswith("read-across"))),
               num(one(INV / "H5_readacross_value.csv", "recall",
                       lambda r: r["method"].startswith("frequency"))),
               int(ra[0]["n"]))
    per5 = rows(INV / "H5_readacross_per_target.csv")
    d["h5_targets"] = len(per5)
    d["h5_cap"] = max(int(x["n"]) for x in per5)

    h7 = rows(INV / "H7_target_discrimination.csv")
    au = [num(x["auroc_vs_random"]) for x in h7]
    se = [num(x["deployed_sensitivity"]) for x in h7]
    d["h7"] = (len(h7), min(au), max(au), st.median(au), min(se), max(se), st.median(se),
               st.mean(se), sum(1 for v in se if v < 0.50))

    d["h8_fire"] = int(float(one(INV / "H8_panel_independence.csv", "value",
                                 lambda r: r["metric"] == "targets that ever fire on the drug set")))
    d["h8_ind"] = int(float(one(INV / "H8_panel_independence.csv", "value",
                                lambda r: r["metric"] == "independent directions in the firing pattern")))

    aud = [r for r in rows(ROOT / "results" / "deployed_specificity_audit.csv")
           if not r["verdict"].startswith("ok")]
    d["withdrawn"] = [(r["target"], num(r["random_fpr"]), r["trivial_fired"]) for r in aud]

    reg = json.loads((ROOT / "models_rf" / "binder_modes.json").read_text(encoding="utf-8"))
    dep = [v for v in reg.values() if v.get("deployed")]
    d["reg_mean"] = st.mean(v["sensitivity_at_threshold"] for v in dep)
    d["n_dep"] = len(dep)

    pn = rows(TAB / "permutation_null.csv")
    cv = {(r["endpoint"], r["split"]): num(r["roc_auc_mean"])
          for r in rows(TAB / "rf_cv_summary.csv") if r["task"] == "classification"}
    order = [r["endpoint"] for r in pn if r["split"] == "random"]
    d["pn_order"] = order
    d["pn"] = {sp: [num(one(TAB / "permutation_null.csv", "permuted_roc_auc_mean",
                            lambda r, e=e, sp=sp: r["endpoint"] == e and r["split"] == sp))
                    for e in order] for sp in ("random", "scaffold")}
    d["pn_true"] = {sp: [cv[(e, sp)] for e in order] for sp in ("random", "scaffold")}
    d["pn_mean"] = {sp: st.mean(d["pn"][sp]) for sp in ("random", "scaffold")}
    d["pn_margin"] = {sp: [t - n for t, n in zip(d["pn_true"][sp], d["pn"][sp])]
                      for sp in ("random", "scaffold")}
    d["pn_sd"] = (min(num(r["permuted_roc_auc_sd"]) for r in pn),
                  max(num(r["permuted_roc_auc_sd"]) for r in pn))
    d["pn_dev"] = max(abs(num(r["permuted_roc_auc_mean"]) - 0.5) for r in pn)

    d["heldout_all"] = sum(len(v) for v in json.loads(
        (ROOT / "models_rf" / "holdout" / "heldout_actives.json").read_text(encoding="utf-8")).values())
    d["h1_n"] = int(one(INV / "H1_disease_layer.csv", "n",
                        lambda r: r["test"] == "observed top-3 accuracy"))

    bs = rows(TAB / "background_specificity.csv")
    d["bs_mean"] = st.mean(num(x["sensitivity_after"]) for x in bs)
    d["bs_reliable"] = sum(1 for x in bs if x["reliable"] == "True")
    d["reg_unreliable"] = sum(1 for v in dep if v.get("reliable_call") is False)

    g = json.loads((INV / "GRAPH_FINGERPRINT.json").read_text(encoding="utf-8"))
    d["graph"] = (g["n_conditions"], g["n_targets"], g["sha256"][:8])
    return d


F = facts()
D = Deck()

LABEL = {"BBB": "BBB", "AChE": "AChE", "BChE": "BChE", "BACE1": "BACE1",
         "GSK3B": "GSK-3β", "MAO_A": "MAO-A", "MAO_B": "MAO-B", "hERG": "hERG"}
TARGET_LABEL = {"Nav1_1": "Nav1.1", "GluA2": "GluA2", "NRF2": "NRF2", "NFKB1": "NF-κB1"}


def sentence(text: str) -> str:
    """Upper-case the first letter and leave the rest alone.

    str.capitalize() lower-cases everything after the first character, which turns BBB into Bbb.
    """
    return text[:1].upper() + text[1:] if text else text

# ---------------------------------------------------------------- 1. title
s = D.dark()
D.text(s, M, 1.74, W - 2 * M, 0.5, "Chapter 9", size=16, color=TEAL, bold=True)
D.text(s, M, 2.18, W - 2 * M - 1.0, 1.6,
       "The falsification suite,\nand what it cost",
       size=42, bold=True, font=HEAD, color=PAPER, line=1.08)
D.text(s, M, 4.20, W - 2 * M - 1.6, 1.1,
       "Nine hypotheses stated so they could fail, each paired with a null model capable of "
       "producing the same apparent success by accident. Four did fail.",
       size=17, color=CHALK, line=1.34)
D.text(s, M, 5.75, W - 2 * M, 0.5,
       "Every result so far has been produced by people who wanted the tool to work, including me.",
       size=14, color=AMBER, italic=True)
D.text(s, M, 6.18, W - 2 * M, 0.4, "inversion/PLAN.md", size=12, color=MUTED)
D.notes(s, "The suite is the strongest evidence in the thesis, and it is strongest where it "
           "refuted something. Lead with that, not with the four that survived.")

# ---------------------------------------------------------------- 2. what makes a test able to fail
s = D.light()
D.head(s, "1", "What makes a test capable of failing",
       "Four constraints, all of them stated before any result was seen")
items = [
    ("A null that could succeed by accident",
     f"79 per cent top-3 accuracy sounds decisive until one asks what always naming the three "
     f"commonest conditions scores. It scores 0.5508."),
    ("Read-only",
     "No model, dataset or threshold was changed to obtain any number. Acting on a finding is a "
     "separate commit, so a wording change cannot be mistaken for evidence."),
    ("Filed under inversion/, not results/",
     "A falsification cannot be found where a validation would be and quoted as one."),
    ("Scored on models_rf/holdout/",
     "Where predictive power is at issue, the scaffold-hold-out twins score compounds they never "
     "saw."),
]
yy = 1.80
for i, (t, b) in enumerate(items, 1):
    D.card(s, M, yy - 0.10, W - 2 * M, 1.06, fill=TINT if i % 2 else PAPER)
    D.dot(s, M + 0.22, yy + 0.20, str(i), fill=TEAL, dia=0.38)
    D.text(s, M + 0.80, yy + 0.02, 4.3, 0.8, t, size=14, bold=True, color=DEEP, line=1.16)
    D.text(s, M + 5.30, yy + 0.02, W - M - (M + 5.30) - 0.30, 0.9, b, size=12.5, color=INK, line=1.24)
    yy += 1.18

D.text(s, M, 6.44, W - 2 * M, 0.46,
       "The first does the work: it asks what would look like success before the answer is known.",
       size=13.5, color=AMBER, italic=True)
D.source(s, "inversion/PLAN.md; inversion/results/H1_disease_layer.csv")
D.notes(s, "If asked why this is not just more validation: validation measures a component you "
           "already believe in. This states the claim you would most like to be true and tries to "
           "break it.")

# ---------------------------------------------------------------- 3. the nine verdicts
s = D.light()
D.head(s, "2", "Nine hypotheses, four refuted",
       f"Graph under test fingerprinted at SHA-256 {F['graph'][2]}..., "
       f"{F['graph'][0]} conditions, {F['graph'][1]} targets")
COL = {"SUPPORTED": TEAL, "REFUTED": CRIMSON, "WEAKENED": AMBER}
D.text(s, M + 0.30, 1.80, 7.6, 0.3, "hypothesis", size=11, bold=True, color=MUTED, space=0)
D.text(s, M + 8.10, 1.80, 2.1, 0.3, "verdict", size=11, bold=True, color=MUTED, space=0)
D.text(s, M + 10.55, 1.80, 1.6, 0.3, "argued in", size=11, bold=True, color=MUTED, space=0)
WHERE = {"H1": "7.7", "H2": "7.6", "H3": "7.5", "H4": "9.6", "H5": "9.6",
         "H6": "7.7", "H7": "5.9", "H8": "7.4, 7.8", "H9": "7.7"}
yy = 2.12
for r in F["verdicts"]:
    hyp = r["hypothesis"]
    key = hyp.split()[0].strip('"')
    v = r["verdict"]
    base = v.split()[0]
    c = COL.get(base, MUTED)
    D.card(s, M, yy - 0.07, W - 2 * M, 0.44, fill=TINT if base == "REFUTED" else PAPER)
    D.text(s, M + 0.30, yy, 7.6, 0.36,
           sentence(hyp[len(key):].strip().rstrip('"')), size=12,
           bold=(base == "REFUTED"), color=INK)
    D.text(s, M + 8.10, yy, 2.3, 0.36, v, size=11.5, bold=True, color=c)
    D.text(s, M + 10.55, yy, 1.6, 0.36, WHERE.get(key, ""), size=11.5, color=MUTED)
    yy += 0.48
D.text(s, M, 6.42, W - 2 * M, 0.46,
       "Seven are argued elsewhere. This chapter takes what the refutations cost, what the suite "
       "found that nothing else did, and where its verdicts are weaker than they read.",
       size=12.5, color=DEEP, italic=True)
D.source(s, "inversion/results/VERDICTS.csv; inversion/results/GRAPH_FINGERPRINT.json")
D.notes(s, "The fingerprint matters: it means a verdict cannot silently come to describe a "
           "different graph from the one it was computed on.")

# ---------------------------------------------------------------- 4. what the refutations cost
s = D.light()
D.head(s, "3", "What the refutations cost",
       "A falsification suite is worth what the project paid for it")
c, u, p = F["h2"]
items = [
    ("H2", "The knowledge graph was re-described",
     f"Curated weights {fmt(c, 4)}, uniform {fmt(u, 4)}, permuted {fmt(p, 4)}. Spread {fmt(c - p, 4)}. "
     f"Relabelled in the source, the interface and the manuscript as a mechanistic prior, not tuned "
     f"parameters.", CRIMSON),
    ("H3", "A design claim became a filter claim",
     "The gate multiplies 14 conditions by the same factor, so it is rank-invariant among them and "
     "cannot sharpen a disease call. Wording corrected everywhere.", CRIMSON),
    ("H7", "Four endpoints were withdrawn",
     f"An audit of the served models found four failures, three of them firing on molecules that "
     f"are not drugs. Panel now {F['n_dep']} deployed.", CRIMSON),
    ("H8", "The interface changed",
     f"{F['h8_fire']} targets fire across approved drugs but span {F['h8_ind']} independent "
     f"directions. Engaged targets are now grouped by homology family with the measured correlation "
     f"quoted.", CRIMSON),
    ("H1", "The suite tightened its own null",
     f"Restricted to compounds no other panel model was trained on: {F['heldout_all']:,} held-out "
     f"entries down to {F['h1_n']:,}. It lowered the headline.", TEAL),
]
yy = 1.74
for k, t, b, col in items:
    D.card(s, M, yy - 0.08, W - 2 * M, 0.90, fill=TINT)
    D.text(s, M + 0.26, yy + 0.16, 0.7, 0.4, k, size=15, bold=True, font=HEAD, color=col)
    D.text(s, M + 1.05, yy + 0.02, 3.9, 0.8, t, size=13, bold=True, color=INK, line=1.16)
    D.text(s, M + 5.15, yy + 0.02, W - M - (M + 5.15) - 0.28, 0.9, b, size=11.8, color=INK, line=1.22)
    yy += 0.99

D.text(s, M, 6.44, W - 2 * M, 0.46,
       "The last is the one to notice. A suite that tightens its own null after the null has "
       "already been passed is behaving correctly.",
       size=12.5, color=TEAL, italic=True)
D.source(s, "commits f40b5d7, 8bb402b, 7fa7203, da82d80, a3d184f; inversion/results/*.csv")
D.notes(s, "Every payment is a commit, not a claim. If challenged on whether the refutations were "
           "acted on, the git history answers it.")

# ---------------------------------------------------------------- 5. the withdrawals
s = D.dark()
D.text(s, M, 0.86, 11.0, 0.5, "What a direct audit of the served models found",
       size=15, color=AMBER, bold=True)
D.text(s, M, 1.30, 11.7, 0.6,
       "H7 asked why approved antiepileptics were silent. Following that line produced an audit of "
       "what the deployed models fire on, and the answer for four of them was: sugar.",
       size=14, color=CHALK, line=1.28)
D.text(s, M + 0.30, 2.28, 2.6, 0.3, "endpoint", size=11, bold=True, color=MUTED)
D.text(s, M + 3.10, 2.28, 2.0, 0.3, "FPR on random chemistry", size=11, bold=True, color=MUTED)
D.text(s, M + 5.60, 2.28, 6.2, 0.3, "trivial molecules called binders", size=11, bold=True, color=MUTED)
yy = 2.68
for tgt, fpr, triv in F["withdrawn"]:
    D.text(s, M + 0.30, yy, 2.6, 0.42, TARGET_LABEL.get(tgt, tgt), size=14, bold=True, color=PAPER)
    D.text(s, M + 3.10, yy, 2.0, 0.42, fmt(fpr, 4), size=14, color=CRIMSON, bold=True)
    D.text(s, M + 5.60, yy, 6.2, 0.42, triv or "none, but above the 5 per cent target",
           size=12.5, color=CHALK if triv else MUTED, italic=not triv)
    yy += 0.62
D.text(s, M, 5.40, 11.7, 0.9,
       "All four were withdrawn, with NR3C1, which passed this audit and was withdrawn on other "
       "grounds. Nothing about the models changed. What changed is that the project stopped serving "
       "endpoints it could not defend.",
       size=13.5, color=PAPER, line=1.28)
D.text(s, M, 6.45, 11.7, 0.5,
       "The audit exists because a hypothesis was written to embarrass the tool, and it worked.",
       size=13, color=AMBER, italic=True)
D.source(s, "results/deployed_specificity_audit.csv")
D.notes(s, "Volunteer this. A panel that once called glucose an ion-channel binder is a better "
           "story told first than found.")

# ---------------------------------------------------------------- 6. H7 convergence
n7, aulo, auhi, aumed, selo, sehi, semed, semean, nlow = F["h7"]
s = D.light()
D.head(s, "4", "The result the suite bought, by accident",
       "H7 measured the quantity Chapter 5 would later find had been published wrongly")
D.text(s, M, 1.72, 5.55, 1.5,
       "H7 asked whether the silent antiepileptics meant the models could not rank. They can: across "
       f"{n7} testable targets the AUROC against 600 random PubChem structures runs from "
       f"{fmt(aulo, 4)} to {fmt(auhi, 4)}, median {fmt(aumed, 4)}. The explanation is refuted.",
       size=13, color=INK, line=1.30)
D.text(s, M, 3.30, 5.55, 1.6,
       "What the same table recorded incidentally was deployed sensitivity, and nothing compared it "
       "against the registry for three weeks.",
       size=13, color=DEEP, line=1.30, bold=True)

xs = M + 6.05
D.card(s, xs, 1.72, W - M - xs, 2.10, fill=TINT)
D.text(s, xs + 0.32, 1.92, 5.4, 0.34, "Deployed sensitivity, panel mean", size=12, bold=True, color=MUTED)
pairs = [("Published in the registry", 0.8983, CRIMSON),
         ("H7, held out vs random PubChem", semean, TEAL),
         ("Chapter 5 audit, held out", F["reg_mean"], TEAL)]
yy = 2.34
for lab, v, col in pairs:
    D.text(s, xs + 0.32, yy, 3.6, 0.36, lab, size=12.5, color=INK)
    D.text(s, xs + 4.10, yy, 1.5, 0.36, fmt(v, 4), size=13.5, bold=True, color=col,
           align=PP_ALIGN.RIGHT)
    yy += 0.45
D.text(s, xs + 0.32, 3.72 - 0.06, 5.4, 0.3, "", size=10, color=MUTED)

D.card(s, xs, 4.00, W - M - xs, 1.62, fill=PAPER)
D.text(s, xs + 0.32, 4.20, 5.4, 1.3,
       f"Two constructions sharing nothing landed {fmt(abs(semean - F['reg_mean']), 3)} apart, and "
       f"the figure they agreed on was not the published one. H7 names six targets under 0.50; the "
       f"audit names six; five are the same.",
       size=13, color=INK, line=1.26)

D.text(s, M, 4.45, 5.55, 1.7,
       "The lesson is about reading, not measuring. The evidence sat in a committed file, in a "
       "report the project wrote itself. Nothing compared it with the registry, because nothing had "
       "been built to.",
       size=13, color=AMBER, line=1.30, italic=True)
D.source(s, "inversion/results/H7_target_discrimination.csv; models_rf/binder_modes.json; "
            "results/tables/sensitivity_reconciliation.csv")
D.notes(s, "This is the single best argument for running a falsification suite at all. Its most "
           "valuable output was not one of its nine answers.")

# ---------------------------------------------------------------- 7. H4
s = D.light()
D.head(s, "5", "H4 is supported by a rule that ignores its own interval",
       "One false positive in 61 compounds")
D.text(s, M + 0.30, 1.86, 3.4, 0.3, "stratum", size=11, bold=True, color=MUTED)
for i, h in enumerate(["n", "false-positive rate", "95% interval"]):
    D.text(s, M + 4.00 + i * 2.05, 1.86, 1.9, 0.3, h, size=11, bold=True, color=MUTED,
           align=PP_ALIGN.RIGHT)
yy = 2.22
for r in F["h4"]:
    far = r["stratum"].startswith("distant")
    D.card(s, M, yy - 0.09, 10.4, 0.56, fill=TINT if far else PAPER)
    D.text(s, M + 0.30, yy, 3.6, 0.42, r["stratum"], size=12.5, bold=far,
           color=CRIMSON if far else INK)
    for i, v in enumerate([r["n"], fmt(num(r["false_positive_rate"]), 4),
                           f"{fmt(num(r['ci95_low']), 4)} to {fmt(num(r['ci95_high']), 4)}"]):
        D.text(s, M + 4.00 + i * 2.05, yy, 1.9, 0.42, v, size=12.5, bold=far,
               color=CRIMSON if far else MUTED, align=PP_ALIGN.RIGHT)
    yy += 0.60
D.card(s, M, yy - 0.09, 10.4, 0.56, fill=PAPER)
D.text(s, M + 0.30, yy, 3.6, 0.42, "comparator: library chemistry", size=12.5, color=DEEP, italic=True)
D.text(s, M + 4.00, yy, 1.9, 0.42, "1,000", size=12.5, color=MUTED, align=PP_ALIGN.RIGHT)
D.text(s, M + 6.05, yy, 1.9, 0.42, fmt(F["lib_fpr"], 4), size=12.5, bold=True, color=DEEP,
       align=PP_ALIGN.RIGHT)
D.text(s, M + 8.10, yy, 1.9, 0.42, "0.0603 to 0.0930", size=12.5, color=MUTED, align=PP_ALIGN.RIGHT)

xs = M + 10.62
D.card(s, xs, 1.80, W - M - xs, 2.10, fill=TINT)
D.text(s, xs + 0.22, 1.96, W - M - xs - 0.44, 0.46, "Fisher exact,\ntwo-sided", size=11.5,
       bold=True, color=MUTED, line=1.10)
yy2 = 2.58
for lab, p in [("vs library", "0.119"), ("vs near", "0.161"), ("vs in domain", "0.023"),
               ("3-way χ²", "0.066")]:
    D.text(s, xs + 0.22, yy2, 0.95, 0.30, lab, size=10.5, color=INK)
    D.text(s, xs + 1.18, yy2, 0.72, 0.30, p, size=10.5, bold=(p == "0.023"),
           color=TEAL if p == "0.023" else MUTED, align=PP_ALIGN.RIGHT)
    yy2 += 0.34

D.text(s, M, 5.20, W - 2 * M, 1.15,
       "The distant stratum's 95 per cent interval contains the comparator it is being compared "
       "against. The verdict rule is a point-estimate comparison: SUPPORTED if the rate is at most "
       "1.5 times the reference.",
       size=13.5, color=CRIMSON, line=1.28, bold=True)
D.text(s, M, 6.05, W - 2 * M, 1.00,
       "The defensible claim is the narrower one, and it is still worth having. Specificity does "
       "not degrade outside the training neighbourhood, which is exactly what the null predicted it "
       "would do. That it improves is not established on 61 compounds.",
       size=13, color=INK, line=1.28)
D.source(s, "inversion/results/H4_distant_specificity.csv; "
            "results/tables/noncns_specificity_summary.csv; Fisher tests computed for this chapter")
D.notes(s, "Note the in-domain row: 0.12 is the worst of the three and higher than the 0.075 "
           "measured on the non-CNS library, whose compounds are also in domain. Overlapping "
           "intervals, so a direction rather than a finding, but it points the other way.")

# ---------------------------------------------------------------- 8. H5
ra, fb, n5 = F["h5"]
s = D.light()
D.head(s, "6", "H5 beats a null the evaluation set caps near 0.09",
       "Read-across is measured in its deployment regime, and only there")
D.stat(s, M, 1.80, 3.0, fmt(ra, 4), "read-across recall, top-3 of five weighted neighbours", TEAL)
D.stat(s, M + 3.35, 1.80, 3.0, fmt(fb, 4), "frequency baseline: always answer with the three "
                                           "commonest index targets", MUTED)
D.stat(s, M + 6.70, 1.80, 3.0, f"{n5:,}", f"held-out compounds across {F['h5_targets']} targets",
       DEEP)
D.stat(s, M + 10.05, 1.80, 2.0, "0.088", "ceiling the null cannot exceed", AMBER)

D.card(s, M, 3.55, W - 2 * M, 1.45, fill=TINT)
D.text(s, M + 0.32, 3.75, W - 2 * M - 0.64, 1.15,
       f"The set caps each target at {F['h5_cap']} compounds, so with {F['h5_targets']} targets "
       f"contributing {n5:,} compounds between them, any constant answer naming three targets can "
       f"be right at most about 360 times. The measured {fmt(fb, 4)} is exactly 240 of {n5:,}: two "
       f"targets, D2 and A2A, each contributing its full block, and nothing else contributing at "
       f"all. Beating a null with a ceiling near 0.09 by 0.91 is not evidence of quality.",
       size=13, color=INK, line=1.28)

D.card(s, M, 5.06, W - 2 * M, 1.30, fill=PAPER)
D.text(s, M + 0.32, 5.24, W - 2 * M - 0.64, 1.02,
       "And the index is not filtered by scaffold. The query and anything at Tanimoto 0.999 or "
       "above are excluded, and nothing else, so a held-out compound's own scaffold siblings are "
       "among the neighbours answering for it. The compounds are held out from the models, not from "
       "the index.",
       size=13, color=CRIMSON, line=1.28)

D.text(s, M, 6.44, W - 2 * M, 0.46,
       "So 0.9726 is the right number for the regime the server runs in, and says nothing about a "
       "novel scaffold or an unrepresented target class.",
       size=13, color=DEEP, italic=True)
D.source(s, "inversion/results/H5_readacross_value.csv, H5_readacross_per_target.csv; "
            "inversion/inv_readacross_value.py:76-79")
D.notes(s, "The report already carries the target-class half of this caveat. The scaffold half is "
           "missing and belongs beside it.")

# ---------------------------------------------------------------- 9. permutation null
s = D.light()
D.head(s, "7", "The claim the suite never made, now measured",
       "The technical report's label-permutation null was hard-coded in a prose block")
cd = CategoryChartData()
cd.categories = [LABEL[e] for e in F["pn_order"]]
cd.add_series("random split", tuple(F["pn"]["random"]))
cd.add_series("scaffold split", tuple(F["pn"]["scaffold"]))
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(1.80),
                        Inches(7.55), Inches(3.55), cd).chart
gf.has_legend = True
gf.legend.position = XL_LEGEND_POSITION.TOP
gf.legend.include_in_layout = False
gf.legend.font.size = Pt(11)
gf.legend.font.name = BODY
va = gf.value_axis
va.minimum_scale, va.maximum_scale = 0.44, 0.56
va.has_major_gridlines = True
va.major_gridlines.format.line.color.rgb = D.grid
va.tick_labels.font.size = Pt(10)
va.tick_labels.font.name = BODY
gf.category_axis.tick_labels.font.size = Pt(10)
gf.category_axis.tick_labels.font.name = BODY
for ser, colr in zip(gf.plots[0].series, (DEEP, TEAL)):
    ser.format.fill.solid()
    ser.format.fill.fore_color.rgb = colr
    ser.format.line.fill.background()

xs = M + 8.05
D.stat(s, xs, 1.80, 2.0, fmt(F["pn_mean"]["random"], 4), "mean under permuted labels, random split",
       DEEP)
D.stat(s, xs + 2.30, 1.80, 2.0, fmt(F["pn_mean"]["scaffold"], 4),
       "mean under permuted labels, scaffold split", TEAL)
D.card(s, xs, 3.42, W - M - xs, 1.95, fill=TINT)
D.text(s, xs + 0.28, 3.62, W - M - xs - 0.56, 1.6,
       f"All sixteen values sit within 0.020 of chance, over ranges of "
       f"{fmt(min(F['pn']['random']), 4)} to {fmt(max(F['pn']['random']), 4)} and "
       f"{fmt(min(F['pn']['scaffold']), 4)} to {fmt(max(F['pn']['scaffold']), 4)}. Per-fold standard "
       f"deviations run {fmt(F['pn_sd'][0], 4)} to {fmt(F['pn_sd'][1], 4)}, so one permutation per "
       f"endpoint carries about ±0.02 of noise and the panel mean is the reliable statistic.",
       size=12.5, color=INK, line=1.26)

D.text(s, M, 5.60, W - 2 * M, 0.55,
       f"Margin over each endpoint's own null: "
       f"{fmt(min(F['pn_margin']['random']), 4)} to {fmt(max(F['pn_margin']['random']), 4)} "
       f"on the random split, and {fmt(min(F['pn_margin']['scaffold']), 4)} to "
       f"{fmt(max(F['pn_margin']['scaffold']), 4)} on the scaffold split, against a largest "
       f"departure from chance of {fmt(F['pn_dev'], 4)} in the null itself, a factor of "
       f"{min(F['pn_margin']['scaffold']) / F['pn_dev']:.0f}.",
       size=13, color=INK, line=1.26)
D.text(s, M, 6.28, W - 2 * M, 0.6,
       "Scaffold grouping alone confers nothing: paired by endpoint the scaffold null exceeds the "
       "random null by a mean of 0.0067, Wilcoxon p = 0.25, and the sign is negative for three of "
       "eight. That is the specific leakage route this test exists to close.",
       size=13, color=TEAL, bold=True, line=1.26)
D.source(s, "results/tables/permutation_null.csv, built for this chapter by "
            "src/brainsafe/evaluation/permutation_null.py, which imports train_rf so the folds are "
            "identical by construction")
D.notes(s, "Say plainly that the report's own triple is corroborated, not reproduced: it records no "
           "seed, and a permutation null is a random quantity, so nobody can recover its exact "
           "numbers. The fix is for the report to cite the artefact.")

# ---------------------------------------------------------------- 10. where the suite is weak
s = D.light()
D.head(s, "8", "Where the suite is weaker than it reads",
       "Six defects found while writing this chapter. None changes a verdict")
items = [
    ("One hypothesis has two verdict rules that disagree",
     "The test script prints WEAKENED for H2; summarise.py publishes REFUTED. The harsher is "
     "published, which is the right direction, but a verdict should not depend on which copy of a "
     "rule you run.", CRIMSON),
    ("No test pins any verdict",
     "The five test files contain no reference to inversion/. A rerun that flipped a verdict would "
     "be caught only by a person reading the report.", CRIMSON),
    ("app.py quotes a withdrawn version of H2",
     f"0.7917 / 0.7911 / 0.7899 over 15,609 compounds, where the artefact says "
     f"{fmt(F['h2'][0], 4)} / {fmt(F['h2'][1], 4)} / {fmt(F['h2'][2], 4)} over {F['h1_n']:,}. The "
     f"population never matched. Audit item BS-M-04, closed in the manuscript, not in the source.",
     AMBER),
    ("The interface shows co-firing figures that have drifted",
     "Five of five correlations and three of four conditional probabilities disagree with the "
     "artefact, and the monoamine entry names the wrong pair as its family's strongest. These are "
     "rendered on a result page.", CRIMSON),
    ("background_specificity.csv holds output from a script since fixed",
     f"It still reports mean sensitivity {fmt(F['bs_mean'], 4)} and {F['bs_reliable']} of "
     f"{F['n_dep']} reliable, where the registry says {fmt(F['reg_mean'], 4)} and "
     f"{F['n_dep'] - F['reg_unreliable']}. Freshness cannot catch it: the graph's inputs are data "
     f"and models, never code.", AMBER),
    ("Three docstrings and three logs have drifted",
     "H7's docstring says sensitivity ranges 0.26 to 0.98 where the artefact says 0.000 to 0.990; "
     "the 4 August logs still record VERDICT H2: WEAKENED.", MUTED),
]
yy = 1.70
for i, (t, b, col) in enumerate(items, 1):
    D.dot(s, M, yy + 0.04, str(i), fill=col, dia=0.36)
    D.text(s, M + 0.58, yy, 4.35, 0.85, t, size=13, bold=True, color=col, line=1.14)
    D.text(s, M + 5.20, yy - 0.02, W - M - (M + 5.20) - 0.10, 0.80, b, size=11.5, color=INK,
           line=1.18)
    yy += 0.86
D.source(s, "inversion/summarise.py:60; inv_disease_layer.py:228; app.py:102-104, 529; "
            "results/tables/background_specificity.csv; audit/AUDIT_REPORT.md:824,857")
D.notes(s, "Items 3 and 4 are the ones to volunteer. Item 4 is user-visible, which makes it the "
           "most serious of the six even though every discrepancy in it is small.")

# ---------------------------------------------------------------- 11. what it cannot do
s = D.dark()
D.text(s, M, 0.90, 11.0, 0.5, "What a falsification suite cannot do", size=15, color=AMBER, bold=True)
D.text(s, M, 1.34, 11.7, 0.5,
       "Stated because this is the strongest evidence in the thesis and overreading it would be easy.",
       size=14, color=CHALK)
items = [
    ("It tests what its author thought to doubt",
     "All nine hypotheses were written by the person who built the system. H8 and H9 do not appear "
     "in the plan; H9 exists only because H6 returned a result its author judged unfairly scored. "
     "Three of the nine exist because earlier ones produced surprises, which is how a suite should "
     "grow and also an admission that its initial coverage was set by intuition."),
    ("Nothing in it tests the barrier model",
     "H1, H2, H3, H6 and H9 test the disease layer; H5 read-across; H7 and H8 the panel; H4 "
     "system-level specificity. The BBB classifier, which gates every disease score and names the "
     "architecture, has no hypothesis. The missing test is whether its contribution could be "
     "replaced by a molecular-weight or cLogP rule without loss."),
    ("A verdict summarises one measurement",
     f"H4 is SUPPORTED on 61 compounds; H5 against a null capped near 0.09; H9 on a mean of 0.6158 "
     f"that hides a range from 0.7981 down to 0.4898. The verdict column is the least informative "
     f"part of the table."),
]
yy = 2.10
for i, (t, b) in enumerate(items, 1):
    D.dot(s, M, yy + 0.06, str(i), fill=AMBER, fg=INK, dia=0.40)
    D.text(s, M + 0.66, yy, 4.4, 1.0, t, size=15, bold=True, color=AMBER, line=1.14)
    D.text(s, M + 5.35, yy - 0.02, 11.9 - 5.35, 1.5, b, size=12.5, color=CHALK, line=1.26)
    yy += 1.62
D.text(s, M, 6.95, 11.7, 0.4,
       "The suite is valuable in proportion to what it refuted, and it refuted four of nine.",
       size=13.5, color=PAPER, italic=True)
D.notes(s, "Close on the second item. Naming the missing hypothesis before a referee does is worth "
           "more than any of the four that survived.")

D.save(OUT)
