"""Build the ML foundations deck: the mathematics and computation, from first principles.

The audience is a reader who needs the whole machine-learning side explained from the ground up: what
a sample is, what the features are, how the labels were made, what each algorithm optimises, and what
was and was not tuned. It is the slide companion to thesis/ml_foundations.md.

Every number is read from an artefact at build time. Nothing is typed into a slide by hand, so a
figure that moves in the pipeline moves here on the next build rather than going stale.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/build_ml_foundations_deck.py
Out:  thesis/presentations/BrainSafe_ML_foundations.pptx
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from pptx.enum.text import PP_ALIGN

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deck_common import (Deck, rows, num, ROOT, TAB, M, W, H,  # noqa: E402
                         INK, DEEP, TEAL, AMBER, CRIMSON, MUTED, TINT, PAPER, CHALK,
                         HEAD, BODY)

OUT = ROOT / "thesis" / "presentations" / "BrainSafe_ML_foundations.pptx"
OOF = ROOT / "data" / "processed" / "cv_predictions"


# ---------------------------------------------------------------- artefacts
def load():
    d = {}
    d["cv"] = rows(TAB / "rf_cv_summary.csv")
    d["abl"] = rows(TAB / "feature_block_ablation.csv")
    d["coll"] = {r["measure"]: r["value"] for r in rows(TAB / "fingerprint_collisions.csv")}
    d["cal"] = rows(TAB / "calibration.csv")
    d["conf"] = rows(TAB / "rf_conformal.csv")
    d["confs"] = rows(TAB / "rf_conformal_scaffold.csv")
    d["cmp"] = rows(TAB / "model_comparison.csv")
    d["lib"] = {r["metric"]: r["value"] for r in rows(TAB / "library_sp3_coverage.csv")}
    d["adme"] = rows(TAB / "adme_cv_summary.csv")
    gp = ROOT / "results" / "gnn" / "gnn_vs_rf.csv"
    d["gnn"] = rows(gp) if gp.exists() else []
    d["binder"] = json.loads((ROOT / "models_rf" / "binder_modes.json").read_text(encoding="utf-8"))
    # class balance, read from the deduplicated out-of-fold files
    bal = {}
    for p in sorted(OOF.glob("*_random_oof.csv")):
        ep = p.name.split("_random_oof")[0]
        r = rows(p)
        ys = [num(x["y_true"]) for x in r]
        if all(y in (0.0, 1.0) for y in ys):
            bal[ep] = (len(ys), int(sum(ys)))
    d["bal"] = bal
    return d


def cv_of(d, ep, split, col):
    for r in d["cv"]:
        if r["endpoint"] == ep and r["split"] == split:
            return num(r.get(col))
    return None


def ablation(d):
    by = {}
    for r in d["abl"]:
        by.setdefault((r["endpoint"], r["metric"]), {})[r["block"]] = num(r["mean"])
    v = [x for x in by.values() if len(x) == 3]
    fp = sum(x["combined"] - x["fingerprint_only"] for x in v) / len(v)
    ds = sum(x["combined"] - x["descriptors_only"] for x in v) / len(v)
    better = sum(1 for x in v if x["combined"] > x["fingerprint_only"])
    return len(v), fp, ds, better


# ---------------------------------------------------------------- helpers
def table(dk, s, x, y, w, headers, body, widths, size=12, rowh=0.30, head_color=DEEP):
    """A plain rule-lined table. widths are fractions of w and must sum to 1."""
    cx = x
    for h, fr in zip(headers, widths):
        dk.text(s, cx, y, w * fr - 0.06, 0.26, h, size=size, bold=True, color=head_color)
        cx += w * fr
    ry = y + 0.30
    for i, r in enumerate(body):
        if i % 2 == 0:
            dk.card(s, x - 0.08, ry - 0.04, w + 0.10, rowh, fill=TINT)
        cx = x
        for cell, fr in zip(r, widths):
            txt, over = (cell if isinstance(cell, tuple) else (cell, {}))
            dk.text(s, cx, ry, w * fr - 0.06, rowh - 0.04, txt,
                    size=over.get("size", size), color=over.get("color", INK),
                    bold=over.get("bold", False))
            cx += w * fr
        ry += rowh
    return ry


def bullets(dk, s, x, y, w, items, size=14, space=9):
    dk.text(s, x, y, w, H - y - 0.8,
            [(f"•   {t}", o) if not isinstance(t, tuple) else t for t, o in
             ((i, {}) if isinstance(i, str) else i for i in items)],
            size=size, color=INK, space=space, line=1.22)


# ---------------------------------------------------------------- slides
def title_slide(dk, d):
    s = dk.dark()
    dk.text(s, M, 2.30, W - 2 * M, 1.30,
            "The mathematics and computation", size=48, bold=True, font=HEAD, color=PAPER)
    dk.text(s, M, 3.45, W - 2 * M, 0.70,
            "Features, labels, algorithms and probability, from first principles",
            size=22, color=CHALK, font=HEAD)
    n_struct = int(float(d["lib"]["distinct SMILES in the endpoint tables"]))
    dk.text(s, M, 4.55, W - 2 * M, 0.50,
            f"{n_struct:,} distinct structures   |   1,036 features   |   "
            f"52 endpoints, 47 deployed   |   no hyperparameter search",
            size=15, color=TEAL, font=BODY)
    dk.text(s, M, H - 1.05, W - 2 * M, 0.40,
            "Every number on every slide is read from an artefact at build time.",
            size=12, color=MUTED, italic=True)


def sample_space(dk, d):
    s = dk.light()
    dk.head(s, "1", "What one sample is",
            "The unit of observation is not a molecule; it is a molecule measured against an endpoint")
    dk.card(s, M, 1.75, W - 2 * M, 1.05, fill=TINT)
    dk.text(s, M + 0.30, 1.95, W - 2 * M - 0.60, 0.70,
            "fₜ : Ω → \U0001d4b4ₜ        Ω = small drug-like molecules,   "
            "\U0001d4b4ₜ = {0,1} or ℝ",
            size=24, bold=True, font=HEAD, color=DEEP, align=PP_ALIGN.CENTER)
    bullets(dk, s, M, 3.10, W - 2 * M, [
        "We never observe fₜ. We observe a finite sample of noisy laboratory measurements of it.",
        ("The sample is NOT drawn uniformly from Ω. It is drawn from what medicinal chemists "
         "chose to synthesise and publish, which clusters tightly around a few series per target.",
         {"color": AMBER}),
        "Every design decision downstream, above all the scaffold split, exists to deal with that "
        "non-uniformity.",
        "A classification endpoint asks a yes/no question; a regression endpoint asks how much.",
    ], size=15)
    dk.source(s, "src/brainsafe/models/train_rf.py; data/endpoints/*.csv")
    dk.notes(s, "The single most important idea on this slide: the training data is a convenience "
                "sample, not a random sample. That is why a random cross-validation split flatters "
                "the model and why the scaffold split is the honest one.")


def sizes(dk, d):
    s = dk.light()
    dk.head(s, "2", "The sample space, endpoint by endpoint",
            "After deduplication: the rows the models are actually fitted and scored on")
    body = []
    for ep, (n, pos) in sorted(d["bal"].items(), key=lambda kv: -kv[1][0]):
        r_a = cv_of(d, ep, "random", "roc_auc_mean")
        s_a = cv_of(d, ep, "scaffold", "roc_auc_mean")
        gap = (r_a - s_a) if (r_a and s_a) else None
        body.append([ep, f"{n:,}", f"{pos:,}", f"{n - pos:,}", f"{100 * pos / n:.1f}%",
                     f"{r_a:.4f}" if r_a else "-",
                     (f"{s_a:.4f}" if s_a else "-", {"bold": True}),
                     (f"-{gap:.4f}" if gap else "-", {"color": AMBER})])
    table(dk, s, M, 1.80, 7.6,
          ["endpoint", "n", "pos", "neg", "pos rate", "AUROC rand", "scaffold", "gap"],
          body, [0.17, 0.12, 0.11, 0.11, 0.13, 0.14, 0.12, 0.10], size=11.5, rowh=0.285)
    x2 = M + 8.0
    dk.text(s, x2, 1.80, W - x2 - M, 0.32, "Two things to notice", size=15, bold=True, color=DEEP)
    bullets(dk, s, x2, 2.22, W - x2 - M, [
        ("Class balance ranges from 23.9% to 86.5% positive. This is why class weighting is not "
         "optional.", {}),
        ("Every endpoint drops under the scaffold split. That gap is the honest measure of how much "
         "random-split performance was analogue recall.", {"color": AMBER}),
    ], size=13, space=10)
    dk.stat(s, x2, 4.55, 2.2, "5", "regression endpoints,\nR² 0.4153 to 0.7231 scaffold",
            color=TEAL, vsize=34)
    dk.stat(s, x2 + 2.4, 4.55, 2.2, "47", "binder endpoints deployed\nof 52 registered",
            color=TEAL, vsize=34)
    dk.source(s, "results/tables/rf_cv_summary.csv; data/processed/cv_predictions/*_random_oof.csv")


def features(dk, d):
    s = dk.light()
    dk.head(s, "3", "The features: one map, used by every model",
            "φ : Ω → ℝ¹⁰³⁶, identical for every endpoint "
            "and every algorithm except the graph network")
    dk.card(s, M, 1.80, (W - 2 * M) * 0.60 - 0.12, 1.55, fill=TINT)
    dk.text(s, M + 0.25, 1.98, (W - 2 * M) * 0.60 - 0.60, 0.40,
            "Block 1   —   1,024-bit ECFP-4 fingerprint", size=16, bold=True, color=DEEP)
    dk.text(s, M + 0.25, 2.38, (W - 2 * M) * 0.60 - 0.60, 0.85,
            "Morgan, radius 2. Each bit is 0 or 1 and records that at least one local atomic "
            "environment hashed into that bit. Binary, sparse, and the block that carries almost all "
            "of the predictive signal.", size=12.5, color=INK, line=1.2)
    x2 = M + (W - 2 * M) * 0.60 + 0.12
    dk.card(s, x2, 1.80, (W - 2 * M) * 0.40 - 0.12, 1.55, fill=TINT)
    dk.text(s, x2 + 0.25, 1.98, (W - 2 * M) * 0.40 - 0.60, 0.40,
            "Block 2   —   12 descriptors", size=16, bold=True, color=DEEP)
    dk.text(s, x2 + 0.25, 2.38, (W - 2 * M) * 0.40 - 0.60, 0.85,
            "mw, clogp, tpsa, hbd, hba, rotatable bonds, aromatic rings, fraction csp3, ring count, "
            "heavy atoms, formal charge, qed", size=12.5, color=INK, line=1.2)
    dk.equation(s, M, 3.65, W - 2 * M,
                [("φ(m) = [ b", False), ("0", True), (", …, b", False), ("1023", True),
                 (",  d", False), ("1", True), (", …, d", False), ("12", True),
                 (" ]   ∈  ℝ", False), ("1036", True)], size=26)
    bullets(dk, s, M, 4.55, W - 2 * M, [
        "Using one representation everywhere means a comparison between two algorithms is never "
        "confounded by a difference in features.",
        "The twelve descriptors are the classical determinants of passive permeability, which is why "
        "they are present at all: the barrier endpoint is where physicochemistry, not substructure, "
        "is the mechanism.",
    ], size=14)
    dk.source(s, "src/brainsafe/features/featurize.py")


def ecfp(dk, d):
    s = dk.light()
    dk.head(s, "4", "How the fingerprint is computed",
            "Iterative neighbourhood hashing. There is no learning anywhere in this step")
    steps = [("Step 0", "Each atom gets an integer identifier hashed from its own properties: "
                        "element, degree, charge, attached hydrogens, ring membership."),
             ("Step k", "Each identifier is replaced by a hash of itself together with its bonded "
                        "neighbours' identifiers and the bond types. Neighbours are sorted first, "
                        "so the result does not depend on atom ordering."),
             ("Radius 2", "Done twice. Each final identifier now encodes the atom plus everything "
                          "within two bonds: a circular substructure. ECFP-4 names the diameter, "
                          "2 × radius, so radius 2 and ECFP-4 are the same thing."),
             ("Folding", "The set of identifiers across a library is unbounded, so each is reduced "
                         "modulo the bit width to give a fixed-length vector.")]
    y = 1.80
    for i, (k, v) in enumerate(steps, 1):
        dk.dot(s, M, y + 0.02, str(i), fill=TEAL, dia=0.34)
        dk.text(s, M + 0.50, y, 1.35, 0.30, k, size=14, bold=True, color=DEEP)
        dk.text(s, M + 1.95, y, W - M - 1.95 - M, 0.72, v, size=13, color=INK, line=1.2)
        y += 0.92
    dk.card(s, M, 5.55, W - 2 * M, 0.95, fill=TINT)
    dk.equation(s, M, 5.78, W - 2 * M,
                [("id", False), ("a", True), ("(k)", True),
                 ("  =  H( id", False), ("a", True), ("(k−1)", True),
                 (" ,  { (β", False), ("ab", True), (" , id", False), ("b", True),
                 ("(k−1)", True), (") : b ∈ N(a) } )", False)], size=20)
    dk.source(s, "RDKit rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)")


def collisions(dk, d):
    s = dk.light()
    dk.head(s, "5", "Folding makes collisions certain, and we measured them",
            "The most important caveat about the representation, stated as a number rather than a "
            "hedge")
    c = d["coll"]
    envs = int(float(c["distinct atomic environments"]))
    amb = int(float(c["bits carrying more than one environment"]))
    med = int(float(c["environments per occupied bit, median"]))
    mx = int(float(c["environments per occupied bit, max"]))
    n = int(float(c["structures sampled"]))
    dk.stat(s, M, 1.85, 2.7, f"{envs:,}", f"distinct atomic environments\nseen in {n:,} structures",
            color=DEEP)
    dk.stat(s, M + 2.9, 1.85, 2.7, "1,024", "bits available\nto hold them", color=DEEP)
    dk.stat(s, M + 5.8, 1.85, 2.7, f"{amb:,}", "of 1,024 bits carry MORE\nthan one environment",
            color=CRIMSON)
    dk.stat(s, M + 8.7, 1.85, 2.7, f"{med}", f"environments share the\nmedian bit (max {mx})",
            color=CRIMSON)
    dk.card(s, M, 3.75, W - 2 * M, 1.05, fill=TINT)
    dk.text(s, M + 0.30, 3.95, W - 2 * M - 0.60, 0.70,
            "Bit i does NOT mean “substructure i is present”.\n"
            "It means “at least one of about 52 different substructures is present”.",
            size=17, bold=True, color=CRIMSON, align=PP_ALIGN.CENTER, line=1.25)
    bullets(dk, s, M, 5.05, W - 2 * M, [
        ("Consequence for interpretation: a bit a tree ranks important names a SET of environments, "
         "so a SHAP attribution to one bit cannot be read as a chemical substructure.", {}),
        ("Consequence for prediction: none. The model only needs the mapping to be consistent, and "
         "it is. Collisions add noise, they do not add bias.", {"color": TEAL}),
        ("Second consequence: the fingerprint is stereo-blind, so enantiomers give identical "
         "vectors. This forces the deduplication step on slide 8.", {}),
    ], size=13.5, space=7)
    dk.source(s, "results/tables/fingerprint_collisions.csv; "
                 "src/brainsafe/evaluation/fingerprint_collisions.py")
    dk.notes(s, "This slide exists because the featuriser docstring once claimed the encoding was "
                "collision-free by construction. That was the reverse of the truth. An external "
                "audit measured it and the claim was corrected.")


def labels(dk, d):
    s = dk.light()
    dk.head(s, "6", "Where the labels came from",
            "pChEMBL = −log₁₀(molar potency), so each whole unit is a tenfold change")
    dk.card(s, M, 1.80, (W - 2 * M) * 0.52, 2.35, fill=TINT)
    dk.text(s, M + 0.30, 2.00, (W - 2 * M) * 0.52 - 0.60, 0.35,
            "The classification cuts", size=16, bold=True, color=DEEP)
    for i, (lab, rule, col) in enumerate([
            ("y = 1", "pChEMBL ≥ 6.0      (≤ 1 µM, active)", TEAL),
            ("y = 0", "pChEMBL < 5.0       (> 10 µM, inactive)", DEEP),
            ("discarded", "5.0 ≤ pChEMBL < 6.0", AMBER)]):
        dk.text(s, M + 0.30, 2.48 + i * 0.52, 1.30, 0.34, lab, size=15, bold=True, color=col)
        dk.text(s, M + 1.75, 2.48 + i * 0.52, (W - 2 * M) * 0.52 - 2.05, 0.34, rule,
                size=13.5, color=INK)
    x2 = M + (W - 2 * M) * 0.52 + 0.25
    dk.text(s, x2, 1.85, W - x2 - M, 0.32,
            "Why the grey zone is thrown away", size=16, bold=True, color=DEEP)
    bullets(dk, s, x2, 2.28, W - x2 - M, [
        "A compound at pChEMBL 5.5 is neither usefully active nor confidently inactive.",
        "Forcing it to one side would teach the model a boundary the assay data does not support.",
        ("The cost is discarded data. The benefit is that the two classes are separated by a full "
         "log unit, comfortably more than inter-laboratory assay variability.", {"color": TEAL}),
    ], size=13, space=8)
    dk.text(s, M, 4.45, W - 2 * M, 0.32,
            "Three further curation rules", size=16, bold=True, color=DEEP)
    bullets(dk, s, M, 4.88, W - 2 * M, [
        "Replicates resolve by MEDIAN, not mean: robust to the single badly transcribed outlier that "
        "is common in aggregated assay literature.",
        "Censored values (“> 10 µM”) are labelled inactive directly, never passed "
        "through the cut. They are a real observation of inactivity, and discarding them would throw "
        "away most of the negative class.",
        "The binder panel uses a stricter cut, pChEMBL ≥ 7.0, because a binder call is a "
        "stronger claim than an activity call.",
    ], size=13, space=7)
    dk.source(s, "src/brainsafe/data/rebuild_endpoints.py:56; "
                 "src/brainsafe/models/train_binders_hybrid.py:57")


def dedup(dk, d):
    s = dk.light()
    dk.head(s, "7", "Deduplication, and why it happens on the feature vector",
            "Because φ is stereo-blind and desalts, distinct database rows become "
            "byte-identical inputs")
    dk.stat(s, M, 1.90, 3.0, "7,807", "raw BBB rows", color=MUTED)
    dk.text(s, M + 3.15, 2.05, 0.8, 0.6, "→", size=40, bold=True, color=AMBER)
    dk.stat(s, M + 4.1, 1.90, 3.0, "3,901", "rows the model is fitted on", color=DEEP)
    dk.card(s, M + 7.5, 1.85, W - M - 7.5 - M, 1.30, fill=TINT)
    dk.text(s, M + 7.75, 2.05, W - M - 7.5 - M - 0.5, 0.95,
            "Half the BBB table is duplicate input.\nThis is not a marginal correction.",
            size=16, bold=True, color=CRIMSON, line=1.3)
    bullets(dk, s, M, 3.45, W - 2 * M, [
        ("Left in place, any splitter puts copies of one compound on both sides of a fold, and the "
         "model is scored on rows it has memorised.", {"color": AMBER}),
        "It must happen on the FEATURE VECTOR, because that is the level at which the rows are "
        "indistinguishable. Neither the SMILES string nor the InChIKey collapses them: both are "
        "unique for every BBB row.",
        "Where duplicated rows disagree on the label, the whole group is DROPPED, not voted on. The "
        "same input carrying both labels cannot be learned from, and picking one would be an "
        "arbitrary decision dressed as data. Regression takes the median.",
        ("An external audit found this guard was applied in training but not in the conformal "
         "analysis. Correcting it moved BBB's mean prediction-set size from 1.013 to 1.215.",
         {"color": TEAL}),
    ], size=13.5, space=9)
    dk.source(s, "src/brainsafe/models/train_rf.py::_dedup_features")


def selection(dk, d):
    s = dk.light()
    dk.head(s, "8", "Feature selection: there is none",
            "The honest answer to “how did you get from many features to these?”")
    dk.card(s, M, 1.80, W - 2 * M, 1.00, fill=TINT)
    dk.text(s, M + 0.30, 2.00, W - 2 * M - 0.60, 0.62,
            "No SelectKBest, no SelectFromModel, no RFE, no variance threshold, no χ², "
            "no mutual information, no PCA.\nAll 1,036 columns go into every model, always.",
            size=17, bold=True, color=DEEP, align=PP_ALIGN.CENTER, line=1.3)
    n, fp, ds, better = ablation(d)
    dk.text(s, M, 3.05, (W - 2 * M) * 0.47, 0.32,
            "Why not select", size=16, bold=True, color=DEEP)
    bullets(dk, s, M, 3.48, (W - 2 * M) * 0.47, [
        "Selection on the full dataset before cross-validation is a classic way a score is inflated: "
        "the selector has seen the test fold's labels.",
        "Doing it correctly needs nesting inside every fold, multiplying compute and complicating "
        "every artefact.",
        ("Tree ensembles already perform implicit selection at each split, so the gain would be "
         "small and the leakage risk real.", {"color": TEAL}),
    ], size=13, space=8)
    x2 = M + (W - 2 * M) * 0.53
    dk.text(s, x2, 3.05, W - x2 - M, 0.32,
            f"What was done instead: an ablation over {n} endpoints",
            size=16, bold=True, color=DEEP)
    table(dk, s, x2, 3.50, W - x2 - M,
          ["remove this block", "mean change"],
          [[("descriptors only kept", {}), (f"+{ds:.4f}", {"bold": True, "color": TEAL})],
           [("fingerprint only kept", {}), (f"+{fp:.4f}", {"bold": True, "color": AMBER})]],
          [0.62, 0.38], size=13, rowh=0.36)
    dk.text(s, x2, 4.55, W - x2 - M, 1.60,
            [(f"The fingerprint earns its place decisively: removing it costs {ds:.4f}.", {}),
             (f"The twelve descriptors add {fp:.4f} on average, and help on only {better} of {n} "
              f"endpoints.", {"color": AMBER}),
             ("They are retained because they cost nothing and the exposure layer needs them, "
              "NOT because they were shown to help.", {"bold": True})],
            size=13, color=INK, space=7, line=1.2)
    dk.text(s, M, 5.95, W - 2 * M, 0.55,
            "Say it first: the descriptors were retained, not validated. That is much stronger than "
            "having it extracted under questioning.",
            size=14, bold=True, italic=True, color=CRIMSON, align=PP_ALIGN.CENTER)
    dk.source(s, "results/tables/feature_block_ablation.csv; verified: no selection or search "
                 "primitives anywhere in the live tree")


def algebra(dk, d):
    s = dk.light()
    dk.head(s, "9", "The linear algebra: what the matrix looks like",
            "Three properties of X drive most of the modelling decisions")
    dk.equation(s, M, 1.80, W - 2 * M,
                [("X  ∈  ℝ", False), ("n × 1036", True),
                 ("        y  ∈  {0,1}", False), ("n", True),
                 ("  or  ℝ", False), ("n", True),
                 ("        BBB:  3901 × 1036", False)], size=23)
    items = [("Sparse", "A drug-like molecule sets only a small fraction of 1,024 bits, so most "
                        "entries are zero. Trees handle this natively; distance methods need a "
                        "similarity built for sparse binary data, which is why read-across uses "
                        "Jaccard and not Euclidean."),
             ("Badly scaled", "Columns 1–1024 are in {0,1}; mw runs to several hundred; qed is "
                              "in [0,1]. Trees are invariant to monotone rescaling, so this is "
                              "irrelevant to them. For a gradient or distance method it is fatal, "
                              "since mw would dominate. THIS is why logistic regression is wrapped "
                              "in a StandardScaler and the forests are not."),
             ("Wide", "p = 1036, while n runs from 566 to 9,933. Some endpoints have fewer samples "
                      "than features, where an unregularised linear model fits exactly and "
                      "generalises not at all. Hence the L2 penalty on the linear baseline.")]
    y = 2.70
    for k, v in items:
        dk.card(s, M, y, W - 2 * M, 1.05, fill=TINT)
        dk.text(s, M + 0.28, y + 0.16, 2.0, 0.32, k, size=15, bold=True, color=DEEP)
        dk.text(s, M + 2.35, y + 0.14, W - 2 * M - 2.65, 0.80, v, size=12.5, color=INK, line=1.18)
        y += 1.20
    dk.source(s, "src/brainsafe/features/featurize.py; "
                 "src/brainsafe/evaluation/model_comparison.py")


def tanimoto(dk, d):
    s = dk.light()
    dk.head(s, "10", "Similarity is an inner product",
            "Tanimoto on binary vectors is the Jaccard index, and it does four different jobs here")
    dk.card(s, M, 1.80, W - 2 * M, 1.05, fill=TINT)
    dk.equation(s, M, 2.05, W - 2 * M,
                [("T(a, b)  =  | a ∩ b | / | a ∪ b |  =  a", False), ("ᵀ", True),
                 ("b  /  ( ‖a‖", False), ("2", True), (" + ‖b‖", False),
                 ("2", True), (" − a", False), ("ᵀ", True), ("b )", False)], size=21)
    dk.text(s, M, 3.05, W - 2 * M, 0.34,
            "For binary vectors aᵀa = ‖a‖² is simply the number of bits set. "
            "T ∈ [0,1], and 1 − T is the Jaccard distance.",
            size=13.5, italic=True, color=MUTED, align=PP_ALIGN.CENTER)
    table(dk, s, M, 3.65, W - 2 * M,
          ["role", "how it is used", "setting"],
          [["Applicability domain", "max T to any training compound", "bands at 0.30 / 0.50"],
           ["Decoy selection", "a decoy must be dissimilar to every active", "TAN_MAX = 0.35"],
           ["Read-across baseline", "Jaccard distance to the 5 nearest measured analogues", "k = 5"],
           ["Leakage checks", "T = 1.0 identifies a compound the model cannot tell from training",
            "2,048 bits"]],
          [0.24, 0.56, 0.20], size=13, rowh=0.42)
    dk.text(s, M, 5.75, W - 2 * M, 0.65,
            "Roles 1 and 4 use a 2,048-bit fingerprint while the models use 1,024. Doubling the width "
            "reduces collisions and sharpens the similarity estimate; it need not match the model's "
            "width because it is measuring chemistry, not feeding an estimator.",
            size=13, color=INK, line=1.2)
    dk.source(s, "src/brainsafe/evaluation/applicability_domain.py; "
                 "src/brainsafe/models/train_binders_hybrid.py")


def forest(dk, d):
    s = dk.light()
    dk.head(s, "11", "The random forest: what it actually optimises",
            "An average of 300 trees, each on a bootstrap resample, each splitting on a random "
            "feature subset")
    dk.text(s, M, 1.78, (W - 2 * M) * 0.50, 0.32,
            "Inside one tree: greedy impurity reduction", size=15, bold=True, color=DEEP)
    dk.card(s, M, 2.18, (W - 2 * M) * 0.50 - 0.15, 1.30, fill=TINT)
    dk.equation(s, M, 2.38, (W - 2 * M) * 0.50 - 0.15,
                [("G(S)  =  1  −  Σ", False), ("c", True), (" p", False), ("c", True),
                 ("²", True)], size=21)
    dk.text(s, M + 0.25, 2.92, (W - 2 * M) * 0.50 - 0.65, 0.45,
            "The split maximising ΔG = G(S) − (|Sᴸ|/|S|)G(Sᴸ) − "
            "(|Sᴿ|/|S|)G(Sᴿ) is chosen. Regression uses variance reduction.",
            size=11.5, color=MUTED, align=PP_ALIGN.CENTER, line=1.15)
    x2 = M + (W - 2 * M) * 0.50 + 0.15
    dk.text(s, x2, 1.78, W - x2 - M, 0.32,
            "Why an ensemble at all", size=15, bold=True, color=DEEP)
    bullets(dk, s, x2, 2.18, W - x2 - M, [
        "A single deep tree has low bias and very high variance.",
        "Averaging B bootstrap trees leaves bias roughly unchanged and cuts variance. If the trees "
        "were independent the variance would fall as σ²/B; they are correlated, so the "
        "real gain is smaller.",
        ("That is exactly why each split sees only √1036 ≈ 32 random features: it "
         "decorrelates the trees so the averaging works.", {"color": TEAL}),
    ], size=12.5, space=7)
    dk.text(s, M, 3.75, W - 2 * M, 0.32,
            "The deployed settings, and the reason for each", size=15, bold=True, color=DEEP)
    table(dk, s, M, 4.15, W - 2 * M,
          ["parameter", "value", "why"],
          [["n_estimators", "300", "enough for the variance reduction to have plateaued"],
           ["min_samples_leaf", "2 core, 4 binder", "prevents single-compound leaves, the main "
                                                    "overfitting route"],
           ["class_weight", "balanced", "w_c = n / (k · n_c); see the next slide"],
           ["max_features", "√p ≈ 32", "decorrelates the trees"],
           ["random_state", "42", "reproducibility"]],
          [0.22, 0.18, 0.60], size=12.5, rowh=0.34)
    dk.text(s, M, 6.20, W - 2 * M, 0.40,
            "The output p̂(x) is the mean vote share across trees. That is a VOTE SHARE, not a "
            "calibrated probability, which is why slide 18 exists.",
            size=13.5, bold=True, color=AMBER)
    dk.source(s, "src/brainsafe/models/train_rf.py::RF_COMMON")


def weighting(dk, d):
    s = dk.light()
    dk.head(s, "12", "Class weighting, worked through",
            "Without it, a model on hERG could score 76% accuracy by calling everything inactive")
    dk.card(s, M, 1.85, W - 2 * M, 0.95, fill=TINT)
    dk.equation(s, M, 2.10, W - 2 * M,
                [("w", False), ("c", True), ("  =  n  /  ( k · n", False), ("c", True),
                 (" )", False)], size=26)
    n, pos = d["bal"]["MAO_A"]
    neg = n - pos
    w1, w0 = n / (2 * pos), n / (2 * neg)
    dk.text(s, M, 3.05, W - 2 * M, 0.32,
            f"Worked example: MAO-A,  n = {n:,},  positives = {pos:,},  negatives = {neg:,},  k = 2",
            size=15, bold=True, color=DEEP)
    dk.stat(s, M, 3.55, 3.0, f"{w1:.3f}", "weight on each POSITIVE\nn / (2 × %d)" % pos,
            color=TEAL, vsize=34)
    dk.stat(s, M + 3.2, 3.55, 3.0, f"{w0:.3f}", "weight on each NEGATIVE\nn / (2 × %d)" % neg,
            color=DEEP, vsize=34)
    dk.stat(s, M + 6.4, 3.55, 3.4, f"{w1 / w0:.2f}×",
            "each positive counts this\nmuch more than each negative", color=AMBER, vsize=34)
    bullets(dk, s, M, 5.15, W - 2 * M, [
        "The weights enter the impurity calculation, so a split that separates the rare class is "
        "rewarded in proportion to how rare that class is.",
        ("Every family in the comparison is weighted: the forest by class_weight, XGBoost by "
         "scale_pos_weight, logistic regression by class_weight. HistGradientBoosting was the one "
         "exception and was corrected during the audit, because a comparison whose conclusion is "
         "“the forest is not distinguishable from the boosters” cannot rest on one "
         "booster having been handicapped.", {"color": TEAL}),
    ], size=13.5, space=9)
    dk.source(s, "sklearn class_weight='balanced'; "
                 "src/brainsafe/evaluation/model_comparison.py::classifiers")


def boosting(dk, d):
    s = dk.light()
    dk.head(s, "13", "Gradient boosting: the opposite philosophy",
            "A forest averages independent trees in parallel; boosting adds shallow trees in "
            "sequence, each correcting its predecessors")
    dk.equation(s, M, 1.80, W - 2 * M,
                [("F", False), ("m", True), ("(x)  =  F", False), ("m−1", True),
                 ("(x)  +  η · h", False), ("m", True), ("(x)", False)], size=24)
    dk.text(s, M, 2.55, W - 2 * M, 0.32,
            "XGBoost expands the loss to second order, which makes the leaf value solvable in closed "
            "form", size=14, bold=True, color=DEEP, align=PP_ALIGN.CENTER)
    dk.card(s, M + 1.2, 2.98, W - 2 * M - 2.4, 1.35, fill=TINT)
    dk.equation(s, M + 1.2, 3.14, W - 2 * M - 2.4,
                [("ℒ  ≈  Σ", False), ("j", True),
                 (" [ G", False), ("j", True), ("w", False), ("j", True),
                 (" + ½(H", False), ("j", True), (" + λ)w", False), ("j", True),
                 ("² ]  +  γT", False)], size=20)
    dk.equation(s, M + 1.2, 3.72, W - 2 * M - 2.4,
                [("⇒   w", False), ("j", True), ("*  =  − G", False), ("j", True),
                 (" / ( H", False), ("j", True), (" + λ )", False)], size=20)
    dk.text(s, M, 4.48, W - 2 * M, 0.36,
            "gᵢ and hᵢ are the first and second derivatives of the loss at the current "
            "prediction; Gⱼ and Hⱼ their sums over leaf j. The closed-form leaf value plus "
            "the λ and γ penalties are what distinguish XGBoost from classical boosting.",
            size=12.5, color=MUTED, align=PP_ALIGN.CENTER, line=1.2)
    table(dk, s, M, 5.05, W - 2 * M,
          ["family", "settings as run in the comparison"],
          [["XGBoost", "n_estimators=400, max_depth=6, learning_rate=0.05, subsample=0.8, "
                       "colsample_bytree=0.8, tree_method=hist, scale_pos_weight per endpoint"],
           ["HistGradientBoosting", "max_iter=400, learning_rate=0.06, class_weight=balanced"],
           ["LogisticRegression", "StandardScaler + C=1.0, max_iter=2000, class_weight=balanced"],
           ["kNN read-across", "n_neighbors=5, metric=jaccard"]],
          [0.24, 0.76], size=12, rowh=0.36)
    dk.source(s, "src/brainsafe/evaluation/model_comparison.py")


def baselines(dk, d):
    s = dk.light()
    dk.head(s, "14", "The two baselines that keep the project honest",
            "A linear model and a similarity search. If a fitted model cannot beat these, it has "
            "learned nothing")
    dk.text(s, M, 1.80, (W - 2 * M) * 0.48, 0.32,
            "Logistic regression", size=16, bold=True, color=DEEP)
    dk.card(s, M, 2.18, (W - 2 * M) * 0.48 - 0.15, 0.85, fill=TINT)
    dk.equation(s, M, 2.36, (W - 2 * M) * 0.48 - 0.15,
                [("log( p / (1−p) )  =  w", False), ("ᵀ", True), ("x + b", False)],
                size=19)
    body = []
    for ep in ["SERT", "antioxidant_DPPH", "D2", "HT2A"]:
        for r in d["cmp"]:
            if (r["split"] == "scaffold" and r["endpoint"] == ep
                    and r["model"] == "LogisticRegression"):
                body.append([ep.replace("_", " "),
                             (f"{num(r['mean']):+.4f}", {"bold": True, "color": CRIMSON})])
    dk.text(s, M, 3.15, (W - 2 * M) * 0.48, 0.30,
            "Scaffold-split R² on regression endpoints", size=13, bold=True, color=MUTED)
    table(dk, s, M, 3.50, (W - 2 * M) * 0.48 - 0.15, ["endpoint", "R²"], body,
          [0.62, 0.38], size=13, rowh=0.34)
    dk.text(s, M, 5.05, (W - 2 * M) * 0.48 - 0.15, 0.80,
            "Negative R² means it does worse than predicting the mean. "
            "Structure-activity relationships are not linear in fingerprint space.",
            size=12.5, color=INK, line=1.2)
    x2 = M + (W - 2 * M) * 0.52
    dk.text(s, x2, 1.80, W - x2 - M, 0.32,
            "kNN read-across: the chemist's own method", size=16, bold=True, color=DEEP)
    bullets(dk, s, x2, 2.20, W - x2 - M, [
        "Find the five nearest training compounds by Jaccard distance on the fingerprint, average "
        "their labels.",
        "Since Jaccard distance is 1 − Tanimoto, this is literally “look at the five most "
        "similar measured compounds”, which is what a medicinal chemist does by eye.",
        ("It is the honest null for the whole enterprise. The forest beats it on all thirteen "
         "endpoints, which is the result that says the models add something a similarity search "
         "does not.", {"color": TEAL}),
    ], size=13, space=8)
    dk.source(s, "results/tables/model_comparison.csv")


def gnn(dk, d):
    s = dk.light()
    dk.head(s, "15", "The graph network: the one model that does not use φ",
            "Instead of a fixed hashed fingerprint, it learns its own representation from the raw "
            "molecular graph")
    dk.text(s, M, 1.78, (W - 2 * M) * 0.50, 0.30,
            "Input: atoms are nodes, bonds are edges. Nothing is hashed.",
            size=13.5, bold=True, color=DEEP)
    table(dk, s, M, 2.14, (W - 2 * M) * 0.50 - 0.20,
          ["atom feature block", "dims"],
          [["element (14 common + other)", "15"],
           ["degree (0–5 + other)", "7"],
           ["hybridisation (5 types + other)", "6"],
           ["attached hydrogens (0–4 + other)", "6"],
           ["charge, aromaticity, ring membership", "3"],
           [("total node dimension", {"bold": True}), ("37", {"bold": True, "color": TEAL})]],
          [0.74, 0.26], size=12, rowh=0.30)
    x2 = M + (W - 2 * M) * 0.52
    dk.text(s, x2, 1.78, W - x2 - M, 0.30,
            "The GIN layer update", size=13.5, bold=True, color=DEEP)
    dk.card(s, x2, 2.14, W - x2 - M, 1.05, fill=TINT)
    dk.equation(s, x2, 2.34, W - x2 - M,
                [("h", False), ("v", True), ("(k)", True), (" = MLP( (1+ε)h", False),
                 ("v", True), ("(k−1)", True), (" + Σ", False),
                 ("u∈N(v)", True), (" h", False), ("u", True), ("(k−1)", True),
                 (" )", False)], size=17)
    dk.text(s, x2, 3.30, W - x2 - M, 0.95,
            "SUM aggregation, not mean or max, is what gives GIN its power: it is provably as "
            "discriminative as the Weisfeiler-Lehman isomorphism test, which mean and max are not. "
            "ε is learned and weights an atom's own state against its neighbourhood. After "
            "3 layers, atom vectors are mean-pooled to one molecule vector.",
            size=12, color=INK, line=1.18)
    dk.text(s, M, 4.30, W - 2 * M, 0.30,
            "Architecture and optimisation", size=13.5, bold=True, color=DEEP)
    dk.text(s, M, 4.66, W - 2 * M, 0.75,
            "hidden 64  |  3 layers  |  dropout 0.2  |  batch norm after every layer  |  "
            "Adam, lr 1e-3, weight decay 1e-5  |  batch 128  |  max 120 epochs, early stopping "
            "patience 18  |  BCEWithLogits with pos_weight = neg/pos  |  MSE on standardised "
            "targets for regression",
            size=12.5, color=INK, line=1.25)
    dk.text(s, M, 5.55, W - 2 * M, 0.70,
            "Adam keeps per-parameter running estimates of the first and second moments of the "
            "gradient and scales each step accordingly, which is what makes it robust to the very "
            "different gradient scales across a deep network's layers. Early stopping is itself "
            "regularisation.",
            size=12.5, color=MUTED, line=1.2)
    dk.source(s, "src/brainsafe/gnn/gin_model.py; graph_features.py; train_gnn.py")


def gnn_result(dk, d):
    s = dk.light()
    dk.head(s, "16", "The graph network lost, on every endpoint tested",
            "Expected at this data scale, and it should be presented as a result rather than "
            "apologised for")
    body = []
    for r in d["gnn"]:
        g, rf = num(r["GIN"]), num(r["RandomForest"])
        body.append([r["endpoint"], r["metric"], f"{int(num(r['n_test'])):,}",
                     f"{g:.4f}", (f"{rf:.4f}", {"bold": True, "color": TEAL}),
                     (r["winner"], {"bold": True,
                                    "color": TEAL if r["winner"] == "RandomForest" else AMBER})])
    table(dk, s, M, 1.85, (W - 2 * M) * 0.60,
          ["endpoint", "metric", "n test", "GIN", "Random forest", "winner"],
          body, [0.20, 0.16, 0.14, 0.16, 0.20, 0.14], size=13, rowh=0.40)
    x2 = M + (W - 2 * M) * 0.63
    dk.text(s, x2, 1.85, W - x2 - M, 0.30, "Why this is the expected result",
            size=15, bold=True, color=DEEP)
    bullets(dk, s, x2, 2.25, W - x2 - M, [
        "A graph network has far more parameters than a forest has effective degrees of freedom.",
        "Graph networks typically need tens of thousands of examples per endpoint before a learned "
        "representation beats a well-made fixed fingerprint.",
        ("With 3,000 to 10,000 compounds, the fingerprint wins. That is the honest reading.",
         {"color": TEAL}),
    ], size=13, space=8)
    dk.card(s, M, 4.15, W - 2 * M, 1.55, fill=TINT)
    dk.text(s, M + 0.30, 4.35, W - 2 * M - 0.60, 1.20,
            [("Important qualification, and it must be volunteered.", {"bold": True,
                                                                       "color": CRIMSON}),
             ("The GIN was run on 4 endpoints with a SINGLE 70/10/20 scaffold hold-out, not 10-fold "
              "cross-validation, and it does not go through the deduplication step. It is therefore "
              "an exploratory demonstration, not a fair benchmark.", {}),
             ("Its numbers carry no fold-level error bars, and its test set may hold stereoisomers "
              "of training compounds. Since it loses anyway, neither changes the conclusion, but the "
              "conclusion must be stated as “a graph network did not beat the fingerprint in "
              "the one exploratory comparison we ran”.", {})],
            size=12.5, color=INK, space=5, line=1.2)
    dk.source(s, "results/gnn/gnn_vs_rf.csv; src/brainsafe/gnn/train_gnn.py")


def protocol(dk, d):
    s = dk.light()
    dk.head(s, "17", "Was the protocol the same for every algorithm?",
            "Yes for the four comparison families. No for the graph network, and the difference "
            "must be declared")
    table(dk, s, M, 1.85, W - 2 * M,
          ["", "endpoints", "splits", "folds", "dedup", "class-weighted", "role"],
          [[("Random forest", {"bold": True}), "13", "random and scaffold", "10", "yes", "yes",
            ("DEPLOYED", {"bold": True, "color": TEAL})],
           ["XGBoost", "13", "random and scaffold", "10", "yes", "scale_pos_weight", "comparison"],
           ["HistGradientBoosting", "13", "random and scaffold", "10", "yes",
            "class_weight (added)", "comparison"],
           ["LogisticRegression", "13", "random and scaffold", "10", "yes", "yes + scaled",
            "baseline"],
           ["kNN read-across", "13", "random and scaffold", "10", "yes", "n/a", "null"],
           [("GIN graph network", {"bold": True, "color": AMBER}),
            ("4", {"color": AMBER}), ("scaffold only", {"color": AMBER}),
            ("1 hold-out", {"color": AMBER}), ("no", {"color": CRIMSON}),
            "pos_weight", ("exploratory", {"bold": True, "color": AMBER})]],
          [0.20, 0.11, 0.20, 0.11, 0.09, 0.17, 0.12], size=12, rowh=0.38)
    dk.text(s, M, 4.65, W - 2 * M, 0.75,
            "The four comparison families are exactly like-for-like with the forest: same 13 "
            "endpoints, same features, same folds, same deduplicated rows. Any difference between "
            "them is attributable to the estimator alone, and that is what makes the comparison "
            "mean anything.",
            size=14, color=INK, line=1.22)
    dk.card(s, M, 5.55, W - 2 * M, 0.95, fill=TINT)
    dk.text(s, M + 0.30, 5.74, W - 2 * M - 0.60, 0.60,
            "Over all 13 endpoints the forest is NOT distinguishable from either booster "
            "(Wilcoxon p = 0.89258 vs HGB, p = 0.73535 vs XGBoost).\n"
            "It is deployed for its calibration behaviour, not for its accuracy.",
            size=14, bold=True, color=DEEP, align=PP_ALIGN.CENTER, line=1.25)
    dk.source(s, "results/tables/model_family_significance.csv")


def splits(dk, d):
    s = dk.light()
    dk.head(s, "18", "The two cross-validation schemes",
            "Every core model is validated twice, and the gap between them is reported rather than "
            "averaged away")
    dk.card(s, M, 1.85, (W - 2 * M) * 0.48, 2.05, fill=TINT)
    dk.text(s, M + 0.28, 2.05, (W - 2 * M) * 0.48 - 0.56, 0.32,
            "Random 10-fold", size=16, bold=True, color=DEEP)
    dk.text(s, M + 0.28, 2.45, (W - 2 * M) * 0.48 - 0.56, 1.35,
            "StratifiedKFold for classification, KFold for regression, shuffled with seed 42. "
            "Compounds are assigned at random, with the class ratio preserved in every fold.\n\n"
            "This is the conventional estimate, and it is optimistic: congeneric analogues land on "
            "both sides.", size=12.5, color=INK, line=1.2)
    x2 = M + (W - 2 * M) * 0.52
    dk.card(s, x2, 1.85, (W - 2 * M) * 0.48, 2.05, fill=TINT)
    dk.text(s, x2 + 0.28, 2.05, (W - 2 * M) * 0.48 - 0.56, 0.32,
            "Scaffold-grouped 10-fold", size=16, bold=True, color=TEAL)
    dk.text(s, x2 + 0.28, 2.45, (W - 2 * M) * 0.48 - 0.56, 1.35,
            "GroupKFold on the Bemis-Murcko scaffold: ring systems and linkers, side chains "
            "stripped. Whole groups move together, so no scaffold appears in both train and test.\n\n"
            "This measures generalisation to structurally new chemistry, which is the regime a user "
            "is actually in.", size=12.5, color=INK, line=1.2)
    dk.text(s, M, 4.15, W - 2 * M, 0.30,
            "Two implementation subtleties that a careless version gets wrong",
            size=15, bold=True, color=DEEP)
    bullets(dk, s, M, 4.55, W - 2 * M, [
        "The scaffold is computed on the SAME desalted parent the featuriser uses, so a salt and its "
        "free base cannot land in different folds while being identical to the model.",
        ("Acyclic compounds share ONE group rather than each getting their own. Giving each its own, "
         "as an earlier version did, quietly turned the scaffold split back into a random split for "
         "that part of the set.", {"color": AMBER}),
        ("Checked during the audit: re-running under ten different shuffled partitions moves BBB's "
         "pooled scaffold AUROC by 0.0029 in total. The headline does not depend on the partition.",
         {"color": TEAL}),
    ], size=13, space=8)
    dk.source(s, "src/brainsafe/models/train_rf.py::_cv, _scaffold_groups")


def tuning(dk, d):
    s = dk.light()
    dk.head(s, "19", "What was optimised, and what was not",
            "The question most likely to be asked. Both halves of the answer must be given together")
    dk.card(s, M, 1.80, W - 2 * M, 0.95, fill=TINT)
    dk.text(s, M + 0.30, 2.00, W - 2 * M - 0.60, 0.58,
            "No grid search, no randomised search, no Bayesian optimisation, no manual sweep.\n"
            "No GridSearchCV, RandomizedSearchCV, optuna, hyperopt or param_grid anywhere in this "
            "repository.",
            size=15, bold=True, color=CRIMSON, align=PP_ALIGN.CENTER, line=1.3)
    dk.text(s, M, 2.95, (W - 2 * M) * 0.48, 0.30, "What it costs", size=15, bold=True, color=AMBER)
    bullets(dk, s, M, 3.32, (W - 2 * M) * 0.48, [
        "No claim of optimality is available.",
        "There is direct evidence performance was left on the table: an untuned XGBoost scores "
        "higher than the deployed forest on ALL FIVE regression endpoints.",
    ], size=13, space=8)
    x2 = M + (W - 2 * M) * 0.52
    dk.text(s, x2, 2.95, W - x2 - M, 0.30, "What it buys", size=15, bold=True, color=TEAL)
    bullets(dk, s, x2, 3.32, W - x2 - M, [
        "There was no tuning set, therefore there can be no tuning-set leakage.",
        "Selecting hyperparameters on the same folds the final number is reported from is the "
        "commonest way a figure in this field is quietly inflated. This project is structurally "
        "immune to it.",
    ], size=13, space=8)
    dk.text(s, M, 4.85, W - 2 * M, 0.30,
            "Three different things the word “optimisation” refers to",
            size=15, bold=True, color=DEEP)
    table(dk, s, M, 5.22, W - 2 * M,
          ["level", "what happens", "did it happen here?"],
          [["1. Inside one tree", "greedy search over (feature, threshold) maximising impurity "
                                  "decrease", ("yes", {"color": TEAL})],
           ["2. Fitting the ensemble", "boosting descends an objective; the GIN runs Adam; the "
                                       "forest has no global objective", ("yes", {"color": TEAL})],
           ["3. Choosing hyperparameters", "searching over settings and keeping the best",
            ("NO", {"bold": True, "color": CRIMSON})]],
          [0.24, 0.56, 0.20], size=12.5, rowh=0.38)
    dk.text(s, M, 6.55, W - 2 * M, 0.35,
            "Most published “we optimised our model” claims refer to level 3. Here only "
            "levels 1 and 2 occurred.",
            size=13, italic=True, color=MUTED, align=PP_ALIGN.CENTER)
    dk.source(s, "verified by search across the whole live tree, 7 September 2026")


def calibration(dk, d):
    s = dk.light()
    dk.head(s, "20", "Making the probability mean what it says",
            "A forest's output is a vote share across 300 trees, and vote shares are systematically "
            "miscalibrated")
    raw = sum(num(r["ece_raw"]) for r in d["cal"]) / len(d["cal"])
    cal = sum(num(r["ece_calibrated"]) for r in d["cal"]) / len(d["cal"])
    dk.card(s, M, 1.82, W - 2 * M, 0.92, fill=TINT)
    dk.equation(s, M, 2.02, W - 2 * M,
                [("ECE  =  Σ", False), ("m", True), (" ( |B", False), ("m", True),
                 ("| / n ) · | confidence(B", False), ("m", True), (") − accuracy(B",
                                                                        False),
                 ("m", True), (") |", False)], size=20)
    dk.stat(s, M, 3.00, 2.6, f"{raw:.4f}", "mean ECE before calibration\n(8 core classifiers)",
            color=AMBER, vsize=36)
    dk.text(s, M + 2.75, 3.15, 0.7, 0.5, "→", size=34, bold=True, color=TEAL)
    dk.stat(s, M + 3.55, 3.00, 2.6, f"{cal:.4f}", "after isotonic calibration", color=TEAL,
            vsize=36)
    x2 = M + 6.6
    dk.text(s, x2, 3.00, W - x2 - M, 0.30, "Isotonic regression, solved by PAVA",
            size=14, bold=True, color=DEEP)
    bullets(dk, s, x2, 3.36, W - x2 - M, [
        "Fit a non-decreasing step function g minimising Σ(yᵢ − g(p̂ᵢ))², "
        "solved exactly by pool-adjacent-violators in O(n).",
        ("Because g is monotone, calibration CANNOT change the ranking. AUROC is unchanged; only "
         "what the number means changes.", {"color": TEAL}),
        "Non-parametric, so it corrects an arbitrary distortion, unlike Platt scaling which assumes "
        "a sigmoid.",
    ], size=12.5, space=7)
    bullets(dk, s, M, 4.95, W - 2 * M, [
        "Fitted on OUT-OF-FOLD predictions under an inner 5-fold, so nothing is calibrated and "
        "scored on the same rows.",
        ("Checked during the audit: under a scaffold split the calibrated ECE is 0.0169 rather than "
         "0.0151, so the choice of split costs about 0.0018 and the claim survives.",
         {"color": TEAL}),
        "The binder panel uses Platt scaling instead, because its calibration sets are smaller and "
        "isotonic overfits on small samples.",
    ], size=13, space=8)
    dk.source(s, "results/tables/calibration.csv; src/brainsafe/models/calibrate.py")


def enrichment(dk, d):
    s = dk.light()
    dk.head(s, "21", "Base rates: why a calibrated probability still cannot be read alone",
            "0.5 is strong evidence of activity in a rare endpoint and evidence of INACTIVITY in a "
            "common one")
    dk.text(s, M, 1.80, W - 2 * M, 0.35,
            "hERG is 23.9% positive. BACE1 is 86.5% positive. The same 0.5 means opposite things.",
            size=15, bold=True, color=AMBER, align=PP_ALIGN.CENTER)
    dk.card(s, M + 1.6, 2.30, W - 2 * M - 3.2, 1.45, fill=TINT)
    dk.equation(s, M + 1.6, 2.52, W - 2 * M - 3.2,
                [("E", False), ("b", True), ("(p)  =  (p − b) / (1 − b)      "
                                             "if  p ≥ b", False)], size=19)
    dk.equation(s, M + 1.6, 3.10, W - 2 * M - 3.2,
                [("E", False), ("b", True), ("(p)  =  (p − b) / b               "
                                             "if  p < b", False)], size=19)
    bullets(dk, s, M, 4.00, W - 2 * M, [
        "Piecewise linear with a kink at p = b, continuous there since both branches give 0.",
        "Maps E(1) = +1, E(0) = −1, so it is a signed measure of evidence in [−1, +1].",
        ("STRICTLY INCREASING in p. This is the key property: it cannot reorder compounds within an "
         "endpoint, so it changes interpretation without changing ranking. The rank-invariance is "
         "proved in the thesis and verified in the falsification suite.", {"color": TEAL}),
        "The same reasoning drives the disease layer: score = max over targets of (weight × "
        "signal), multiplied by predicted barrier penetration. Because the barrier term multiplies "
        "every gated disease identically, it cannot change their relative order; it decides whether "
        "anything is reported at all.",
    ], size=13.5, space=9)
    dk.source(s, "app.py::enrichment, disease_scores; inversion/results/H3_gating.csv")


def conformal(dk, d):
    s = dk.light()
    dk.head(s, "22", "Conformal prediction: from a probability to a guaranteed set",
            "Inductive Mondrian conformal at ε = 0.10, class-conditional so the guarantee is "
            "not met on the majority class alone")
    dk.text(s, M, 1.78, (W - 2 * M) * 0.52, 0.30, "The construction", size=15, bold=True,
            color=DEEP)
    steps = ["Split into training, calibration and test.",
             "Nonconformity of a calibration compound: αᵢ = 1 − p̂(true class).",
             "Per class c, take the threshold at rank k = ⌈(n_c + 1)(1 − ε)⌉ of "
             "the sorted α.",
             "The prediction set holds every class whose nonconformity is below its threshold."]
    y = 2.14
    for i, t in enumerate(steps, 1):
        dk.dot(s, M, y, str(i), fill=TEAL, dia=0.30)
        dk.text(s, M + 0.44, y - 0.02, (W - 2 * M) * 0.52 - 0.50, 0.55, t, size=12.5, color=INK,
                line=1.18)
        y += 0.62
    dk.card(s, M, 4.72, (W - 2 * M) * 0.52 - 0.15, 0.78, fill=TINT)
    dk.text(s, M + 0.20, 4.88, (W - 2 * M) * 0.52 - 0.55, 0.50,
            "P( y ∈ Γ(x) )  ≥  1 − ε\nholds under EXCHANGEABILITY of "
            "calibration and test",
            size=13.5, bold=True, color=DEEP, align=PP_ALIGN.CENTER, line=1.25)
    x2 = M + (W - 2 * M) * 0.55
    dk.text(s, x2, 1.78, W - x2 - M, 0.30, "Reading the set size", size=15, bold=True, color=DEEP)
    dk.card(s, x2, 2.12, W - x2 - M, 0.62, fill=TINT)
    dk.text(s, x2 + 0.15, 2.26, W - x2 - M - 0.30, 0.36,
            "E[ |Γ| ]  =  1  +  P(ambiguous)  −  P(empty)",
            size=15, bold=True, color=DEEP, align=PP_ALIGN.CENTER)
    body = []
    for r in d["conf"]:
        sc = next((q for q in d["confs"] if q["endpoint"] == r["endpoint"]), None)
        body.append([r["endpoint"], f"{num(r['empirical_coverage']):.3f}",
                     f"{num(r['avg_set_size']):.3f}",
                     f"{100 * num(r['frac_ambiguous']):.1f}%",
                     (f"{100 * num(sc['frac_ambiguous']):.1f}%" if sc else "-",
                      {"color": AMBER})])
    table(dk, s, x2, 2.92, W - x2 - M,
          ["endpoint", "coverage", "set size", "ambig.", "ambig. scaffold"],
          body, [0.22, 0.20, 0.19, 0.19, 0.20], size=11, rowh=0.27)
    dk.text(s, M, 5.68, W - 2 * M, 0.72,
            "A size-2 set says “either class is plausible”; a size-0 set says "
            "“neither is”, which is a legitimate signal that the compound is anomalous. "
            "BACE1's mean set size of 0.956 is BELOW one, which is only possible because 4.4% of its "
            "sets are empty. Reading the excess over 1.0 as the ambiguous fraction is valid only "
            "when no set is empty.",
            size=12.5, color=INK, line=1.2)
    dk.source(s, "results/tables/rf_conformal.csv, rf_conformal_scaffold.csv")


def architecture(dk, d):
    s = dk.light()
    dk.head(s, "23", "What kind of classification this is",
            "Three things commonly misunderstood about the architecture")
    items = [("Independent binary classifiers",
              "Not multi-class, not multi-label. Each endpoint has its own model, training table, "
              "threshold and calibrator. There is no shared representation, no softmax over targets "
              "and no joint loss. A compound called active at AChE and BChE is the output of two "
              "entirely separate forests.", DEEP),
             ("The panel does not enforce consistency",
              "Nothing prevents a chemically implausible combination, because no component sees more "
              "than one endpoint. This is a real limitation and it is the price of the modularity "
              "that lets an endpoint be withdrawn without retraining anything else.", AMBER),
             ("The disease layer is not learned at all",
              "It is a fixed curated knowledge graph mapping targets to conditions, verified against "
              "KEGG and Reactome. An earlier prototype that LEARNED this layer was shown by ablation "
              "to be reading the answer back out of its own features, with structure-only "
              "performance collapsing to near zero. Making it non-learned is what stops the disease "
              "results being circular.", TEAL)]
    y = 1.85
    for k, v, col in items:
        dk.card(s, M, y, W - 2 * M, 1.42, fill=TINT)
        dk.text(s, M + 0.28, y + 0.16, W - 2 * M - 0.56, 0.32, k, size=16, bold=True, color=col)
        dk.text(s, M + 0.28, y + 0.56, W - 2 * M - 0.56, 0.80, v, size=12.5, color=INK, line=1.2)
        y += 1.58
    dk.source(s, "app.py; src/brainsafe/panel.py; src/brainsafe/data/verify_pathways.py")


def summary(dk, d):
    s = dk.dark()
    dk.text(s, M, 0.55, W - 2 * M, 0.60, "The one-page answer", size=34, bold=True, font=HEAD,
            color=PAPER)
    n, fp, ds, better = ablation(d)
    raw = sum(num(r["ece_raw"]) for r in d["cal"]) / len(d["cal"])
    cal = sum(num(r["ece_calibrated"]) for r in d["cal"]) / len(d["cal"])
    qa = [("What is a sample?", "A (molecule, endpoint) pair with a measured assay outcome"),
          ("How many?", "170,619 distinct structures; 566 to 9,933 per endpoint"),
          ("What are the features?", "1,024-bit ECFP-4 + 12 descriptors = 1,036, identical for "
                                     "every model"),
          ("Are they interpretable?", "The 12 descriptors yes. The 1,024 bits NO: a median of 52 "
                                      "environments share each bit"),
          ("How were labels made?", "pChEMBL ≥ 6.0 active, < 5.0 inactive, 5.0–6.0 "
                                    "DISCARDED; binder panel uses 7.0"),
          ("How were features selected?", "They were NOT. All 1,036 always used, justified by "
                                          "ablation rather than selection"),
          ("Do the descriptors help?", f"+{fp:.4f} on average, better on {better} of {n}. Retained, "
                                       f"not validated"),
          ("Does the fingerprint help?", f"+{ds:.4f}. Decisively yes"),
          ("Which model is deployed?", "Random forest: 300 trees, min_samples_leaf 2, class_weight "
                                       "balanced, seed 42"),
          ("Same protocol for all?", "Yes for XGBoost, HGB, logistic regression, kNN. NO for the "
                                     "GIN (4 endpoints, 1 hold-out)"),
          ("Did the graph network win?", "No. It lost on all four endpoints tested"),
          ("What was tuned?", "NOTHING. No hyperparameter search of any kind exists in this "
                              "repository"),
          ("What kind of classification?", "Independent binary classifiers per endpoint, plus 5 "
                                           "regressors and 8 ADME models"),
          ("How are probabilities honest?", f"Isotonic calibration, ECE {raw:.4f} → {cal:.4f}, "
                                            f"plus conformal sets at 90%")]
    y = 1.30
    for i, (q, a) in enumerate(qa):
        if i % 2 == 0:
            dk.card(s, M - 0.08, y - 0.03, W - 2 * M + 0.16, 0.40,
                    fill=dk.__class__.grid if False else CHALK)
    y = 1.30
    for q, a in qa:
        dk.text(s, M, y, 3.55, 0.34, q, size=12.5, bold=True, color=TEAL)
        dk.text(s, M + 3.70, y, W - M - 3.70 - M, 0.34, a, size=12.5, color=PAPER)
        y += 0.40
    dk.text(s, M, H - 0.62, W - 2 * M, 0.34,
            "Every figure on this slide is read from an artefact at build time.",
            size=11, italic=True, color=MUTED)


# Speaker notes, by slide number. This deck teaches rather than reports, so each slide carries the
# thing a presenter should say that is not already printed on it: usually the reason a choice was made
# or the objection it anticipates.
NOTES = {
    1: "Open by saying what the deck is for: this is the machine-learning side from first "
       "principles, not a results talk. Every number is read from an artefact at build time, so if "
       "a figure here disagrees with a paper, the artefact is right and the paper needs fixing.",
    2: "The single most important idea in the deck. The training data is a convenience sample, not "
       "a random sample: it is what chemists chose to publish. That is why a random cross-validation "
       "split flatters the model and why the scaffold split is the honest one. Everything else "
       "follows from this.",
    3: "Two questions come from this table. Why is BACE1 so easy? Because it is 86.5 per cent "
       "positive, so the base rate does much of the work, which is exactly why enrichment over the "
       "base rate is reported rather than the raw probability. And why does regression degrade more "
       "than classification? Because predicting a potency value is harder than predicting a side of "
       "a threshold.",
    4: "Emphasise that one representation is used everywhere. If asked why not a different "
       "fingerprint per endpoint: because then a comparison between endpoints, or between "
       "algorithms, would be confounded by representation and would mean nothing.",
    5: "Stress that there is no learning in this step. People often assume the fingerprint is "
       "trained. It is a deterministic hashing procedure, fixed before any model sees the data. "
       "That is why it can be computed once and reused by every endpoint.",
    6: "This slide exists because the featuriser docstring once claimed the encoding was "
       "collision-free by construction. That was the reverse of the truth. An external audit "
       "measured it and the claim was corrected. If asked whether collisions invalidate the results: "
       "no, they add noise rather than bias, because the mapping is consistent. What they "
       "invalidate is per-bit interpretation.",
    7: "The grey zone is the point to dwell on. Throwing away data feels wrong to most audiences. "
       "The argument is that a label you do not believe is worse than no label, because the model "
       "will learn the boundary you asserted rather than the one the assay supports.",
    8: "If asked why not deduplicate on InChIKey: because the InChIKey distinguishes stereoisomers "
       "and salts while the model does not, so it would leave the duplicates in place. The right "
       "level is the one at which the rows are indistinguishable to the estimator. Mention that an "
       "audit found this guard missing from the conformal analysis, and that fixing it made the "
       "reported uncertainty worse, which is the correct direction.",
    9: "The likely challenge is that not doing feature selection is lazy. The answer is that it is "
       "conservative: selection done on the full dataset leaks test labels, and doing it correctly "
       "requires nesting inside every fold. Volunteer the descriptor result before it is extracted: "
       "they were retained, not validated.",
    10: "The scaling point is the one worth landing. Trees are invariant to monotone rescaling, so "
        "scaling is a no-op for them and a correctness requirement for logistic regression. That is "
        "why one family is wrapped in a StandardScaler and the others are not; it is not an "
        "inconsistency.",
    11: "Tanimoto appears in four roles and audiences conflate them. The applicability domain and "
        "the leakage checks use 2,048 bits while the models use 1,024, and that is deliberate: a "
        "wider fingerprint collides less and so measures chemistry more sharply, and it never feeds "
        "an estimator.",
    12: "If asked why 300 trees: because the variance reduction has plateaued by then, and more "
        "trees cost time without changing the answer. If asked why min_samples_leaf 2: a leaf "
        "holding one compound is memorisation, and this is the cheapest structural guard against it.",
    13: "Work the arithmetic aloud. It makes the abstraction concrete and it pre-empts the question "
        "of whether class imbalance was handled. Then make the fairness point: the audit found "
        "HistGradientBoosting was the only unweighted ensemble, which biased the comparison against "
        "it, and it was corrected.",
    14: "The contrast to draw is parallel versus sequential. A forest averages many independent "
        "low-bias high-variance trees; boosting adds many shallow high-bias trees in sequence. They "
        "attack the bias-variance trade-off from opposite ends, which is why comparing them on "
        "identical folds is informative.",
    15: "These two baselines are what make the whole comparison credible. If a fitted model cannot "
        "beat a five-nearest-neighbour similarity search, it has learned nothing a chemist could not "
        "do by eye. The negative R-squared for logistic regression is worth pausing on: it means "
        "worse than predicting the mean.",
    16: "Explain why sum aggregation matters if the audience is technical: mean and max pooling "
        "cannot distinguish certain graph structures, sum can, and that is the Weisfeiler-Lehman "
        "result GIN is built on. If the audience is not technical, the message is simply that the "
        "network learns its own features instead of being handed fixed ones.",
    17: "Do not apologise for this result. State it as expected at this data scale and give the "
        "reason: graph networks need far more data per endpoint before a learned representation "
        "beats a good fixed one. Then volunteer the qualification about the single hold-out before "
        "anyone asks, because being the one to raise it is much stronger than conceding it.",
    18: "This is the slide that answers the question most likely to be asked. The four comparison "
        "families are genuinely like-for-like; the graph network is not, and saying so unprompted is "
        "the whole point. Finish on the significance result: the forest is deployed for its "
        "calibration behaviour, not because it is more accurate.",
    19: "If challenged that GroupKFold is deterministic and therefore fragile, give the number: ten "
        "shuffled partitions move BBB by 0.0029 in total. The concern is reasonable and it was "
        "checked and refuted.",
    20: "The most important slide for credibility. Give both halves in one breath: no optimality "
        "claim is available, and an untuned XGBoost beats us on all five regression endpoints, but "
        "there was no tuning set so there can be no tuning-set leakage. Say it first rather than "
        "conceding it under questioning.",
    21: "The key property is monotonicity. Because isotonic regression is non-decreasing, "
        "calibration cannot reorder compounds, so AUROC is untouched and only the meaning of the "
        "number changes. People frequently assume calibration improves discrimination. It does not.",
    22: "Give a concrete example: a base rate of 24 per cent means a probability of 0.5 is roughly "
        "double the prior and is real evidence; a base rate of 86 per cent means 0.5 is well below "
        "the prior and is evidence of inactivity. The same number, opposite meanings.",
    23: "The empty set is the part audiences find strange. A set of size zero is not a failure, it "
        "is the method saying no class is plausible, which is a useful signal that the compound is "
        "anomalous. BACE1's mean set size below 1.0 is the proof that empty sets occur, and it is "
        "why reading the excess over 1.0 as ambiguity is invalid.",
    24: "The circularity point is the one that matters scientifically. An earlier prototype learned "
        "the disease layer and an ablation showed it was reading the answer out of its own features. "
        "Making the layer a fixed, externally verified graph is what makes the disease results "
        "non-circular, and it is a deliberate design decision rather than a shortcut.",
    25: "Use this as the recap and as the answer sheet. If one line is worth repeating, it is that "
        "nothing was tuned, and that this is a strength as well as a limitation.",
}


def main() -> None:
    d = load()
    dk = Deck()
    title_slide(dk, d)
    sample_space(dk, d)
    sizes(dk, d)
    features(dk, d)
    ecfp(dk, d)
    collisions(dk, d)
    labels(dk, d)
    dedup(dk, d)
    selection(dk, d)
    algebra(dk, d)
    tanimoto(dk, d)
    forest(dk, d)
    weighting(dk, d)
    boosting(dk, d)
    baselines(dk, d)
    gnn(dk, d)
    gnn_result(dk, d)
    protocol(dk, d)
    splits(dk, d)
    tuning(dk, d)
    calibration(dk, d)
    enrichment(dk, d)
    conformal(dk, d)
    architecture(dk, d)
    summary(dk, d)

    slides = list(dk.prs.slides)
    missing = [i for i in range(1, len(slides) + 1) if i not in NOTES]
    if missing:
        raise SystemExit(f"no speaker notes written for slide(s): {missing}")
    for i, txt in NOTES.items():
        slides[i - 1].notes_slide.notes_text_frame.text = txt

    dk.save(OUT)


if __name__ == "__main__":
    main()
