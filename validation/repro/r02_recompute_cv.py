"""Re-run the core cross-validation from the endpoint tables and score it independently.

This reproduces the pipeline rather than reimplementing it: the data loading, featurisation,
deduplication, scaffold grouping and estimator settings are imported from the pipeline, so what is
being tested is whether running them again produces the reported numbers. What is NOT imported is
the metric computation or the cross-validation loop, both of which are written out again here. A
discrepancy therefore localises to one of those two, or to non-determinism.

Metrics computed per fold and over the pooled out-of-fold predictions:
  AUROC, AUPRC, sensitivity, specificity, balanced accuracy, MCC, Brier, ECE   (classification)
  R2, RMSE, MAE, Spearman                                                       (regression)

Sensitivity, specificity, balanced accuracy and MCC need an operating point. The pipeline's summary
reports them, so this uses the same 0.5 cut on the calibrated-free forest vote, and records the
threshold in the output so the choice is visible rather than implied.

Bootstrap 95% confidence intervals are computed on the pooled out-of-fold predictions by resampling
compounds, not folds, with a fixed seed.

Output: validation/repro/recomputed_folds.csv, recomputed_summary.csv, recomputed_bootstrap.csv

Run:  python validation/repro/r02_recompute_cv.py
      python validation/repro/r02_recompute_cv.py BBB AChE      (subset)
      python validation/repro/r02_recompute_cv.py --no-bootstrap
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import GroupKFold, KFold, StratifiedKFold

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize                                    # noqa: E402
from models.train_rf import (CLASSIFICATION, N_SPLITS, REGRESSION, RF_COMMON,  # noqa: E402
                             SEED, _dedup_features, _load, _scaffold_groups)

OUT = ROOT / "validation" / "repro"
BOOTSTRAP_N = 2000
THRESHOLD = 0.5


# ----------------------------------------------------------------------------------------------
# metrics, written here rather than imported, so this can disagree with the pipeline
# ----------------------------------------------------------------------------------------------
def auroc(y, p) -> float:
    """Rank-based AUROC with explicit tie handling, equal to the Mann-Whitney U statistic."""
    y = np.asarray(y).astype(int)
    n1, n0 = int(y.sum()), int(len(y) - y.sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(p, kind="mergesort")
    ranks = np.empty(len(p), float)
    sp = np.asarray(p)[order]
    i = 0
    while i < len(sp):
        j = i
        while j + 1 < len(sp) and sp[j + 1] == sp[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def auprc(y, p) -> float:
    """Average precision: sum over thresholds of precision weighted by the gain in recall."""
    y = np.asarray(y).astype(int)
    if y.sum() == 0:
        return float("nan")
    order = np.argsort(-np.asarray(p), kind="mergesort")
    ys = y[order]
    tp = np.cumsum(ys)
    precision = tp / np.arange(1, len(ys) + 1)
    return float((precision * ys).sum() / y.sum())


def confusion(y, p, thr=THRESHOLD):
    y = np.asarray(y).astype(int)
    pred = (np.asarray(p) >= thr).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    return tp, tn, fp, fn


def clf_metrics(y, p) -> dict:
    tp, tn, fp, fn = confusion(y, p)
    sens = tp / (tp + fn) if tp + fn else float("nan")
    spec = tn / (tn + fp) if tn + fp else float("nan")
    denom = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / denom) if denom > 0 else float("nan")
    ppv = tp / (tp + fp) if tp + fp else float("nan")
    f1 = (2 * ppv * sens / (ppv + sens)) if (ppv == ppv and sens == sens and ppv + sens) else float("nan")
    y = np.asarray(y).astype(int)
    p = np.asarray(p, dtype=float)
    brier = float(np.mean((p - y) ** 2))
    # 10-bin expected calibration error, equal-width bins on the probability axis
    ece, n = 0.0, len(y)
    for lo in np.arange(0, 1.0, 0.1):
        m = (p >= lo) & (p < lo + 0.1) if lo < 0.9 else (p >= lo) & (p <= 1.0)
        if m.sum():
            ece += (m.sum() / n) * abs(y[m].mean() - p[m].mean())
    return {"roc_auc": auroc(y, p), "pr_auc": auprc(y, p), "sensitivity": sens,
            "specificity": spec, "balanced_acc": (sens + spec) / 2 if sens == sens and spec == spec
            else float("nan"), "mcc": mcc, "f1": f1, "brier": brier, "ece": float(ece),
            "tp": tp, "tn": tn, "fp": fp, "fn": fn}


def reg_metrics(y, p) -> dict:
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    ss_res = float(((y - p) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return {"r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
            "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
            "mae": float(np.mean(np.abs(y - p))),
            "spearman": float(spearmanr(y, p).statistic) if len(y) > 2 else float("nan")}


def bootstrap_ci(y, p, fn, n=BOOTSTRAP_N, seed=SEED):
    """Percentile interval over resampled compounds. Compounds, not folds: the claim is about
    generalisation to other compounds."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    p = np.asarray(p)
    vals = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        v = fn(y[idx], p[idx])
        if v == v:
            vals.append(v)
    if not vals:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


# ----------------------------------------------------------------------------------------------
def run_endpoint(ep: str, do_bootstrap: bool):
    task = "classification" if ep in CLASSIFICATION else "regression"
    target = "label" if task == "classification" else REGRESSION[ep]
    df = _load(ep).dropna(subset=["smiles", target]).reset_index(drop=True)
    X, mask = featurize(df["smiles"].tolist())
    df = df.loc[mask].reset_index(drop=True)
    y = df[target].to_numpy()
    y = y.astype(int) if task == "classification" else y.astype(float)
    groups = _scaffold_groups(df["smiles"].tolist())
    X, y, groups, smiles, dedup = _dedup_features(X, y, groups, df["smiles"].tolist(), task)

    schemes = {
        "random": (StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED)
                   if task == "classification"
                   else KFold(N_SPLITS, shuffle=True, random_state=SEED)),
        "scaffold": GroupKFold(N_SPLITS),
    }
    fold_rows, summary_rows, boot_rows = [], [], []
    for split, splitter in schemes.items():
        it = (splitter.split(X, y, groups) if split == "scaffold" else splitter.split(X, y))
        pooled = np.full(len(y), np.nan)
        for fold, (tr, te) in enumerate(it, 1):
            if task == "classification":
                m = RandomForestClassifier(class_weight="balanced", **RF_COMMON)
                m.fit(X[tr], y[tr])
                p = m.predict_proba(X[te])[:, 1]
                met = clf_metrics(y[te], p)
            else:
                m = RandomForestRegressor(**RF_COMMON)
                m.fit(X[tr], y[tr])
                p = m.predict(X[te])
                met = reg_metrics(y[te], p)
            pooled[te] = p
            fold_rows.append({"endpoint": ep, "task": task, "split": split, "fold": fold,
                              "n_train": int(len(tr)), "n_test": int(len(te)), **met})

        per_fold = pd.DataFrame([r for r in fold_rows if r["split"] == split])
        keys = ([c for c in ("roc_auc", "pr_auc", "sensitivity", "specificity", "balanced_acc",
                             "mcc", "f1", "brier", "ece") if c in per_fold]
                if task == "classification"
                else [c for c in ("r2", "rmse", "mae", "spearman") if c in per_fold])
        row = {"endpoint": ep, "task": task, "split": split, "n": int(len(y)),
               "n_rows_before_dedup": int(dedup.get("rows_in", len(y))),
               "duplicates_removed": int(dedup.get("duplicate_rows_removed", 0)),
               "threshold": THRESHOLD if task == "classification" else None}
        for k in keys:
            row[f"{k}_mean"] = round(float(per_fold[k].mean()), 6)
            row[f"{k}_sd"] = round(float(per_fold[k].std(ddof=1)), 6)
        # pooled out-of-fold, which is what a bootstrap interval should be taken on
        pooled_met = clf_metrics(y, pooled) if task == "classification" else reg_metrics(y, pooled)
        for k, v in pooled_met.items():
            if isinstance(v, (int, float)):
                row[f"pooled_{k}"] = round(float(v), 6)
        summary_rows.append(row)

        if do_bootstrap:
            fns = ({"roc_auc": auroc, "pr_auc": auprc} if task == "classification"
                   else {"r2": lambda a, b: reg_metrics(a, b)["r2"]})
            for name, fn in fns.items():
                lo, hi = bootstrap_ci(y, pooled, fn)
                boot_rows.append({"endpoint": ep, "task": task, "split": split, "metric": name,
                                  "pooled_estimate": round(float(fn(y, pooled)), 6),
                                  "ci95_low": round(lo, 6), "ci95_high": round(hi, 6),
                                  "n_boot": BOOTSTRAP_N, "seed": SEED})
    return fold_rows, summary_rows, boot_rows


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Re-run core CV and score it independently.")
    ap.add_argument("endpoints", nargs="*", help="endpoints (default: all core)")
    ap.add_argument("--no-bootstrap", action="store_true")
    args = ap.parse_args(argv)
    eps = args.endpoints or (CLASSIFICATION + list(REGRESSION))

    t0 = time.time()
    folds, summ, boots = [], [], []
    for ep in eps:
        te = time.time()
        f, s, b = run_endpoint(ep, not args.no_bootstrap)
        folds += f; summ += s; boots += b
        for row in s:
            head = "roc_auc_mean" if row["task"] == "classification" else "r2_mean"
            print(f"[{ep:16s}] {row['split']:9s} n={row['n']:6d} {head[:-5]} "
                  f"{row[head]:.4f} +/- {row[head.replace('mean', 'sd')]:.4f}", flush=True)
        print(f"    ({time.time() - te:.0f}s)", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(folds).to_csv(OUT / "recomputed_folds.csv", index=False)
    pd.DataFrame(summ).to_csv(OUT / "recomputed_summary.csv", index=False)
    if boots:
        pd.DataFrame(boots).to_csv(OUT / "recomputed_bootstrap.csv", index=False)
    meta = {"commit": json.loads((OUT / "environment.json").read_text())["commit"],
            "seed": SEED, "n_splits": N_SPLITS, "threshold": THRESHOLD,
            "bootstrap_n": BOOTSTRAP_N, "endpoints": eps,
            "wall_clock_s": round(time.time() - t0, 1)}
    (OUT / "recomputed_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"\nwrote recomputed_folds.csv, recomputed_summary.csv"
          f"{', recomputed_bootstrap.csv' if boots else ''} ({meta['wall_clock_s']}s)")


if __name__ == "__main__":
    main(sys.argv[1:])
