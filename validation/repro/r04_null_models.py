"""Null models: what score does this pipeline produce when there is nothing to learn?

A cross-validated AUROC of 0.93 means little on its own. The question a reviewer asks is what the
same pipeline, on the same compounds, with the same folds, returns when the labels carry no
information. If a permuted-label run scores near 0.5 the pipeline is measuring signal; if it scores
appreciably above 0.5 the design leaks, and the leak is in the procedure rather than in the data.

Three nulls, each isolating a different failure:

  permuted labels     labels shuffled within the endpoint before splitting. Breaks the
                      structure-activity relationship while preserving class balance and the
                      splitter. Expected AUROC 0.5.
  permuted within     labels shuffled within each training fold only, leaving test labels intact.
                      Detects whether the model can score above chance from class balance alone.
  stratified random   a predictor that ignores the molecule and draws from the training class
                      frequency. The floor any real model must clear.

The scaffold split is included because a permuted-label run under scaffold grouping is the sharper
test: if whole scaffold classes carry class-frequency information, a null model can exceed 0.5
without any structure-activity signal at all, and that would inflate the real numbers too.

Repeats are kept modest and the seed is fixed and recorded; the point is to locate the null, not to
estimate it to three decimals.

Output: validation/repro/null_models.csv

Run:  python validation/repro/r04_null_models.py
      python validation/repro/r04_null_models.py BBB AChE --repeats 3
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold, StratifiedKFold

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
sys.path.insert(0, str(ROOT / "validation" / "repro"))
from features.featurize import featurize                                        # noqa: E402
from models.train_rf import (CLASSIFICATION, N_SPLITS, RF_COMMON, SEED,          # noqa: E402
                             _dedup_features, _load, _scaffold_groups)
from r02_recompute_cv import auroc                                               # noqa: E402

OUT = ROOT / "validation" / "repro"
REPEATS = 5
# Fewer trees for the nulls: the null level is being located, not measured to three decimals, and
# this keeps a five-repeat two-split sweep affordable. Recorded in the output so it is visible.
NULL_RF = dict(RF_COMMON, n_estimators=120)


def folds_for(split, X, y, groups):
    sp = (GroupKFold(N_SPLITS) if split == "scaffold"
          else StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED))
    return list(sp.split(X, y, groups) if split == "scaffold" else sp.split(X, y))


def run_null(kind, X, y, groups, split, rng):
    """One null run: returns pooled out-of-fold AUROC."""
    pooled = np.full(len(y), np.nan)
    if kind == "permuted_labels":
        y_use = rng.permutation(y)
    else:
        y_use = y.copy()

    for tr, te in folds_for(split, X, y_use, groups):
        y_tr = y_use[tr].copy()
        if kind == "permuted_within_train":
            y_tr = rng.permutation(y_tr)
        if kind == "stratified_random":
            p = rng.random(len(te)) * 0 + y_tr.mean()
            pooled[te] = rng.random(len(te))          # ignores the molecule entirely
            continue
        if len(np.unique(y_tr)) < 2:
            continue
        m = RandomForestClassifier(class_weight="balanced", **NULL_RF)
        m.fit(X[tr], y_tr)
        pooled[te] = m.predict_proba(X[te])[:, 1]

    ok = ~np.isnan(pooled)
    # scored against the TRUE labels: a null model must be judged on the real task
    return auroc(y[ok], pooled[ok])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Null and permutation models.")
    ap.add_argument("endpoints", nargs="*")
    ap.add_argument("--repeats", type=int, default=REPEATS)
    args = ap.parse_args(argv)
    eps = args.endpoints or list(CLASSIFICATION)

    t0, rows = time.time(), []
    for ep in eps:
        df = _load(ep).dropna(subset=["smiles", "label"]).reset_index(drop=True)
        X, mask = featurize(df["smiles"].tolist())
        df = df.loc[mask].reset_index(drop=True)
        y = df["label"].to_numpy().astype(int)
        groups = _scaffold_groups(df["smiles"].tolist())
        X, y, groups, _s, _rep = _dedup_features(X, y, groups, df["smiles"].tolist(),
                                                 "classification")
        for split in ("random", "scaffold"):
            for kind in ("permuted_labels", "permuted_within_train", "stratified_random"):
                vals = []
                for rep in range(args.repeats):
                    rng = np.random.default_rng(SEED + 1000 * rep)
                    vals.append(run_null(kind, X, y, groups, split, rng))
                vals = [v for v in vals if v == v]
                rows.append({"endpoint": ep, "split": split, "null": kind,
                             "repeats": len(vals),
                             "auroc_mean": round(float(np.mean(vals)), 4),
                             "auroc_sd": round(float(np.std(vals, ddof=1)), 4) if len(vals) > 1
                             else 0.0,
                             "auroc_min": round(float(np.min(vals)), 4),
                             "auroc_max": round(float(np.max(vals)), 4),
                             "n": int(len(y)), "n_estimators": NULL_RF["n_estimators"],
                             "seed": SEED})
                print(f"[{ep:6s}] {split:9s} {kind:22s} AUROC "
                      f"{np.mean(vals):.4f} +/- {np.std(vals, ddof=1) if len(vals) > 1 else 0:.4f}",
                      flush=True)

    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "null_models.csv", index=False)

    print("\n=== null level by split (mean over endpoints) ===")
    for split in ("random", "scaffold"):
        g = out[out.split == split]
        for kind in g.null.unique():
            v = g[g.null == kind].auroc_mean
            print(f"  {split:9s} {kind:22s} {v.mean():.4f}  (max over endpoints {v.max():.4f})")
    worst = out[out.null == "permuted_labels"].auroc_mean.max()
    print(f"\nhighest permuted-label AUROC anywhere: {worst:.4f}"
          f"  {'<- above 0.55, investigate' if worst > 0.55 else '(consistent with chance)'}")
    meta = {"commit": json.loads((OUT / 'environment.json').read_text())["commit"],
            "seed": SEED, "repeats": args.repeats, "n_estimators": NULL_RF["n_estimators"],
            "wall_clock_s": round(time.time() - t0, 1)}
    (OUT / "null_models_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"wrote null_models.csv ({meta['wall_clock_s']}s)")


if __name__ == "__main__":
    main(sys.argv[1:])
