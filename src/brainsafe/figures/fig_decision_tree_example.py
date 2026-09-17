"""One tree from one forest, shown rather than described.

The manuscript states that a random forest of 300 trees is fitted per endpoint; this figure answers
the question that leaves unanswered, which is what one of those trees actually looks like.

The deployed BBB endpoint is served by `BBB_calibrated.joblib`, a `CalibratedClassifierCV` holding
five internal random forests (one per calibration fold, see docs/ML_METHODS_AND_FORMULAS.md section
3), not by the standalone `BBB.joblib` refit. This figure therefore extracts tree 0 of fold 0 of the
object that actually serves predictions, checked directly: querying this exact tree's parent forest
through the full calibrated pipeline on donepezil reproduces the probability the app reports, to
five decimal places.

Tree 0 is shown to depth 3 (the full tree is exported as text alongside this figure, not omitted,
because a reader should be able to check that depth 3 is a genuine excerpt rather than a
simplification chosen to look clean).

Read plainly: this is one vote among 300 in one of five calibration folds, on one bootstrap resample
of that fold's training rows, considering a random subset of the 1,036 features at each split
(ceil(sqrt(1036)) = 32, per scikit-learn's default for classification). No single tree is the model;
the deployed probability is the mean, across five folds, of each fold's isotonically-calibrated
forest-vote average.

Reads models_rf/BBB_calibrated.joblib, models_rf/feature_names.json.
Writes manuscript/figures/FigureS_decision_tree_example.png (and .pdf) and
results/tables/decision_tree_example_full.txt (the complete tree, as text rules).

Run:  python src/brainsafe/figures/fig_decision_tree_example.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
from sklearn.tree import export_text, plot_tree

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))
import style as S  # noqa: E402

MODEL = "BBB"
MAX_DEPTH_SHOWN = 3


def main() -> None:
    S.use()
    cal = joblib.load(ROOT / "models_rf" / f"{MODEL}_calibrated.joblib")
    fold0_forest = cal.calibrated_classifiers_[0].estimator
    names = json.load((ROOT / "models_rf" / "feature_names.json").open())
    tree = fold0_forest.estimators_[0]

    full_n_nodes = tree.tree_.node_count
    full_depth = tree.tree_.max_depth

    # This is a standalone explanatory figure, not constrained to a journal column width, and at
    # depth 3 the bottom row holds eight boxes: narrower than this and their text collides.
    fig, ax = plt.subplots(figsize=(18, 8))
    artists = plot_tree(tree, max_depth=MAX_DEPTH_SHOWN, feature_names=names,
                        class_names=["non-penetrant", "penetrant"], filled=True,
                        impurity=False, proportion=True, rounded=True, fontsize=8,
                        ax=ax, precision=2)

    # Recolour to the house palette: sklearn's default fill is an orange/blue majority-class
    # gradient. Reads each box's own annotated class-proportion text to decide which side of the
    # split it leans rather than trusting patch RGB, since that is the actual data the box shows.
    for a in artists:
        txt = a.get_text()
        box = a.get_bbox_patch()
        if box is None:
            continue
        if "value = [" not in txt:
            continue
        try:
            frac = txt.split("value = [")[1].split("]")[0].split(",")
            non_pen, pen = float(frac[0]), float(frac[1])
        except Exception:
            continue
        lean = pen  # fraction predicted penetrant
        if lean >= 0.5:
            box.set_facecolor(S.TARGET)
            box.set_alpha(0.20 + 0.55 * (lean - 0.5) * 2)
        else:
            box.set_facecolor(S.SAFETY)
            box.set_alpha(0.20 + 0.55 * (0.5 - lean) * 2)
        box.set_edgecolor(S.MUTED)
        box.set_linewidth(0.7)

    S.note(fig,
           f"Tree 0 of 300, from fold 0 of the 5-fold calibrated ensemble that serves the deployed "
           f"{MODEL} model, shown to depth {MAX_DEPTH_SHOWN} ({full_n_nodes:,} nodes total in this "
           f"tree; full depth {full_depth}). The complete tree is exported as text in "
           f"results/tables/decision_tree_example_full.txt. Colour shows which side of the split "
           f"the node's samples lean, green towards barrier-penetrant, vermillion towards "
           f"non-penetrant; it is one vote of 1,500 (300 trees x 5 folds), not the model.")
    out = S.save(fig, "FigureS_decision_tree_example")
    print(f"  wrote {out.relative_to(ROOT)}")

    full_text = export_text(tree, feature_names=names, show_weights=True)
    out_txt = ROOT / "results" / "tables" / "decision_tree_example_full.txt"
    header = (f"Tree 0 of 300, from fold 0 of the 5-fold CalibratedClassifierCV ensemble that "
              f"serves the deployed {MODEL} model ({MODEL}_calibrated.joblib).\n"
              f"{full_n_nodes:,} nodes, {tree.get_n_leaves():,} leaves, max depth {full_depth}.\n"
              f"Feature names from models_rf/feature_names.json; weights are [non-penetrant, "
              f"penetrant] sample counts reaching that node.\n"
              f"Generated by src/brainsafe/figures/fig_decision_tree_example.py\n"
              + "=" * 90 + "\n\n")
    out_txt.write_text(header + full_text, encoding="utf-8")
    print(f"  wrote {out_txt.relative_to(ROOT)} ({len(full_text):,} chars, "
          f"{full_n_nodes:,} nodes)")


if __name__ == "__main__":
    main()
