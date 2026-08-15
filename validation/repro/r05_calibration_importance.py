"""Calibration curves and feature importance, recomputed from a re-run of the pipeline.

Calibration. The manuscript reports expected calibration error before and after isotonic regression.
ECE is a summary and can hide the shape of the miscalibration, so the reliability curve is drawn as
well: predicted probability against observed frequency, in ten equal-width bins, with the count in
each bin shown, because a bin holding nine compounds and a bin holding nine hundred should not look
alike. The calibrator is fitted on out-of-fold predictions using an inner split, so no compound
contributes to the calibrator that scores it.

Feature importance. Two methods, because impurity importance is known to favour high-cardinality
features and the fingerprint block is 1,024 binary columns against 12 continuous descriptors, which
is exactly the situation in which the two disagree:

  impurity      mean decrease in Gini, read from the fitted forest, effectively free
  permutation   drop in AUROC when one column is shuffled in the held-out fold, which measures
                what the fitted model actually uses rather than what it split on while fitting

SHAP is NOT computed. The package is not installed in this environment, and installing it would
change the environment recorded for this reproduction. It is reported as a blocker rather than
substituted for silently: permutation importance answers a related question but is not SHAP, and
presenting one as the other would be the failure this exercise exists to avoid.

Output: validation/repro/calibration_curve.csv, calibration_curve.png,
        validation/repro/feature_importance.csv

Run:  python validation/repro/r05_calibration_importance.py
      python validation/repro/r05_calibration_importance.py BBB AChE
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                                  # noqa: E402
import numpy as np                                                               # noqa: E402
import pandas as pd                                                              # noqa: E402
from sklearn.ensemble import RandomForestClassifier                              # noqa: E402
from sklearn.isotonic import IsotonicRegression                                  # noqa: E402
from sklearn.model_selection import StratifiedKFold                              # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
sys.path.insert(0, str(ROOT / "validation" / "repro"))
sys.path.insert(0, str(ROOT / "src" / "brainsafe" / "figures"))
from features.featurize import MORGAN_BITS, featurize, feature_names             # noqa: E402
from models.train_rf import (CLASSIFICATION, N_SPLITS, RF_COMMON, SEED,          # noqa: E402
                             _dedup_features, _load, _scaffold_groups)
from r02_recompute_cv import auroc                                               # noqa: E402

OUT = ROOT / "validation" / "repro"
N_BINS = 10
N_PERM_REPEATS = 3
TOP_K = 25


def reliability(y, p, n_bins=N_BINS):
    rows, ece, n = [], 0.0, len(y)
    edges = np.linspace(0, 1, n_bins + 1)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (p >= lo) & (p < hi) if hi < 1.0 else (p >= lo) & (p <= 1.0)
        if not m.sum():
            rows.append({"bin_low": round(lo, 2), "bin_high": round(hi, 2), "n": 0,
                         "mean_predicted": None, "observed_frequency": None})
            continue
        mp, of = float(p[m].mean()), float(y[m].mean())
        ece += (m.sum() / n) * abs(of - mp)
        rows.append({"bin_low": round(lo, 2), "bin_high": round(hi, 2), "n": int(m.sum()),
                     "mean_predicted": round(mp, 5), "observed_frequency": round(of, 5)})
    return rows, float(ece)


def oof_predictions(X, y, groups):
    """Out-of-fold forest votes, and isotonic-calibrated versions of the same.

    The calibrator for each fold is fitted on the other folds' out-of-fold predictions, so the
    compound being calibrated never contributed to the calibrator applied to it.
    """
    raw = np.full(len(y), np.nan)
    skf = StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED)
    folds = list(skf.split(X, y))
    for tr, te in folds:
        m = RandomForestClassifier(class_weight="balanced", **RF_COMMON)
        m.fit(X[tr], y[tr])
        raw[te] = m.predict_proba(X[te])[:, 1]

    cal = np.full(len(y), np.nan)
    for tr, te in folds:
        iso = IsotonicRegression(out_of_bounds="clip")
        iso.fit(raw[tr], y[tr])
        cal[te] = iso.predict(raw[te])
    return raw, cal


def importances(X, y, names):
    """Impurity importance from a fitted forest, and permutation importance on a held-out fold."""
    skf = StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED)
    tr, te = next(iter(skf.split(X, y)))
    m = RandomForestClassifier(class_weight="balanced", **RF_COMMON)
    m.fit(X[tr], y[tr])
    imp = m.feature_importances_

    base = auroc(y[te], m.predict_proba(X[te])[:, 1])
    rng = np.random.default_rng(SEED)
    # permuting all 1,036 columns is not affordable; the columns worth testing are those the forest
    # says it used, plus every descriptor, since the descriptor block is the interpretable one
    cand = sorted(set(np.argsort(-imp)[:TOP_K].tolist()) | set(range(MORGAN_BITS, X.shape[1])))
    perm = {}
    for j in cand:
        drops = []
        for _ in range(N_PERM_REPEATS):
            Xp = X[te].copy()
            Xp[:, j] = rng.permutation(Xp[:, j])
            drops.append(base - auroc(y[te], m.predict_proba(Xp)[:, 1]))
        perm[j] = float(np.mean(drops))
    return imp, perm, float(base), int(len(te))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Calibration curves and feature importance.")
    ap.add_argument("endpoints", nargs="*")
    args = ap.parse_args(argv)
    eps = args.endpoints or list(CLASSIFICATION)

    t0 = time.time()
    names = feature_names()
    cal_rows, imp_rows, summary = [], [], []

    for ep in eps:
        df = _load(ep).dropna(subset=["smiles", "label"]).reset_index(drop=True)
        X, mask = featurize(df["smiles"].tolist())
        df = df.loc[mask].reset_index(drop=True)
        y = df["label"].to_numpy().astype(int)
        groups = _scaffold_groups(df["smiles"].tolist())
        X, y, groups, _s, _r = _dedup_features(X, y, groups, df["smiles"].tolist(),
                                               "classification")

        raw, cal = oof_predictions(X, y, groups)
        r_rows, r_ece = reliability(y, raw)
        c_rows, c_ece = reliability(y, cal)
        for kind, rows_, e in (("raw", r_rows, r_ece), ("isotonic", c_rows, c_ece)):
            for r in rows_:
                cal_rows.append({"endpoint": ep, "calibration": kind, **r})
            summary.append({"endpoint": ep, "calibration": kind, "ece": round(e, 5),
                            "brier": round(float(np.mean(((raw if kind == 'raw' else cal) - y) ** 2)),
                                           5), "n": int(len(y))})
        print(f"[{ep:6s}] ECE raw {r_ece:.4f} -> isotonic {c_ece:.4f}", flush=True)

        imp, perm, base, n_te = importances(X, y, names)
        for j, v in sorted(perm.items(), key=lambda kv: -kv[1]):
            imp_rows.append({"endpoint": ep, "feature_index": int(j), "feature": names[j],
                             "block": "descriptor" if j >= MORGAN_BITS else "fingerprint",
                             "impurity_importance": round(float(imp[j]), 8),
                             "permutation_auroc_drop": round(v, 6),
                             "fold_auroc_baseline": round(base, 5), "n_test": n_te,
                             "n_permutation_repeats": N_PERM_REPEATS, "seed": SEED})

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(cal_rows).to_csv(OUT / "calibration_curve.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT / "calibration_summary.csv", index=False)
    pd.DataFrame(imp_rows).to_csv(OUT / "feature_importance.csv", index=False)

    # ---- the plot -----------------------------------------------------------------------------
    cal_df = pd.DataFrame(cal_rows).dropna(subset=["mean_predicted"])
    n = len(eps)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.9 * nrow), squeeze=False)
    for ax, ep in zip(axes.ravel(), eps):
        ax.plot([0, 1], [0, 1], color="#98A4AE", lw=0.8, ls=(0, (3, 2)), zorder=1)
        for kind, col in (("raw", "#98A4AE"), ("isotonic", "#127C71")):
            g = cal_df[(cal_df.endpoint == ep) & (cal_df.calibration == kind)]
            if not len(g):
                continue
            ax.plot(g.mean_predicted, g.observed_frequency, "-o", ms=3, lw=1.1, color=col,
                    label=kind, zorder=3)
            ax.scatter(g.mean_predicted, g.observed_frequency,
                       s=2 + 26 * g.n / max(g.n.max(), 1), facecolor=col, alpha=0.25,
                       edgecolor="none", zorder=2)
        e_raw = next((s["ece"] for s in summary if s["endpoint"] == ep and s["calibration"] == "raw"), None)
        e_cal = next((s["ece"] for s in summary if s["endpoint"] == ep and s["calibration"] == "isotonic"), None)
        ax.set_title(f"{ep}   ECE {e_raw:.3f} to {e_cal:.3f}", fontsize=8, loc="left")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("mean predicted probability", fontsize=7)
        ax.set_ylabel("observed frequency", fontsize=7)
        ax.tick_params(labelsize=6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    axes.ravel()[0].legend(fontsize=6, frameon=False)
    fig.suptitle("Reliability curves, out-of-fold; marker area is the number of compounds in the bin",
                 fontsize=8, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / "calibration_curve.png", dpi=300)
    plt.close(fig)

    meta = {"commit": json.loads((OUT / "environment.json").read_text())["commit"], "seed": SEED,
            "n_bins": N_BINS, "permutation_repeats": N_PERM_REPEATS, "top_k_by_impurity": TOP_K,
            "shap": "NOT COMPUTED: package 'shap' is not installed in this environment; "
                    "installing it would change the recorded environment. Reported as a blocker.",
            "wall_clock_s": round(time.time() - t0, 1)}
    (OUT / "calibration_importance_meta.json").write_text(json.dumps(meta, indent=2),
                                                          encoding="utf-8")
    print(f"\nwrote calibration_curve.csv/.png, calibration_summary.csv, feature_importance.csv "
          f"({meta['wall_clock_s']}s)")
    print("SHAP: NOT COMPUTED (package absent) - recorded as a blocker, not substituted")


if __name__ == "__main__":
    main(sys.argv[1:])
