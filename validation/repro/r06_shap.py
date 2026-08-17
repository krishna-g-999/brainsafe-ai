"""SHAP attributions for the deployed classifiers, which the first reproduction pass could not run.

The reproduction recorded SHAP as a blocker because the package was absent and installing it would
have changed the environment the run was certified against. It has since been installed, and the
install was checked to be purely additive: cloudpickle, llvmlite, numba, shap and slicer were added
and no existing package moved, so numpy, scikit-learn, pandas and scipy are at the versions every
reproduced number was computed under. Nothing already in the ledger is invalidated by this file.

TreeExplainer is exact for a random forest rather than an approximation, so this is a direct
attribution and not a sampled estimate. It is still expensive on 1,036 columns, so the explainer runs
on a background sample and the sample size is recorded next to every number.

Two things are reported that permutation importance cannot give:

  direction   permutation importance says a feature matters; SHAP says which way it pushes. The
              direction is NOT the mean signed SHAP: averaged over compounds, the contributions of
              high-value and low-value molecules cancel, and a feature with a strong effect in both
              directions averages to nearly zero. Direction is therefore taken as the Spearman
              correlation between a feature's value and its SHAP value across the explained
              compounds, which is positive when a larger value pushes towards the positive class.
              For a BBB model the sign on TPSA is then a claim about medicinal chemistry that a
              reader can check against what is known.
  agreement   SHAP and permutation importance answer related but different questions. Where they
              disagree the model is probably using a feature that correlates with a better one, so
              the rank correlation between them is reported rather than one being presented as the
              truth.

Output: validation/repro/shap_values.csv, shap_summary.png, shap_vs_permutation.csv

Run:  python validation/repro/r06_shap.py
      python validation/repro/r06_shap.py BBB hERG --sample 400
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
from scipy.stats import spearmanr                                                # noqa: E402
from sklearn.ensemble import RandomForestClassifier                              # noqa: E402
from sklearn.model_selection import StratifiedKFold                              # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import MORGAN_BITS, featurize, feature_names             # noqa: E402
from models.train_rf import (CLASSIFICATION, N_SPLITS, RF_COMMON, SEED,          # noqa: E402
                             _dedup_features, _load, _scaffold_groups)

OUT = ROOT / "validation" / "repro"
SAMPLE = 300          # compounds the explainer is evaluated on
TOP_K = 20


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="SHAP attributions for the deployed classifiers.")
    ap.add_argument("endpoints", nargs="*")
    ap.add_argument("--sample", type=int, default=SAMPLE)
    args = ap.parse_args(argv)
    eps = args.endpoints or list(CLASSIFICATION)

    import shap
    t0, rows = time.time(), []
    names = feature_names()

    for ep in eps:
        df = _load(ep).dropna(subset=["smiles", "label"]).reset_index(drop=True)
        X, mask = featurize(df["smiles"].tolist())
        df = df.loc[mask].reset_index(drop=True)
        y = df["label"].to_numpy().astype(int)
        groups = _scaffold_groups(df["smiles"].tolist())
        X, y, groups, _s, _r = _dedup_features(X, y, groups, df["smiles"].tolist(),
                                               "classification")

        # the same first fold r05 used for permutation importance, so the two are comparable
        tr, te = next(iter(StratifiedKFold(N_SPLITS, shuffle=True,
                                           random_state=SEED).split(X, y)))
        model = RandomForestClassifier(class_weight="balanced", **RF_COMMON)
        model.fit(X[tr], y[tr])

        rng = np.random.default_rng(SEED)
        idx = te if len(te) <= args.sample else rng.choice(te, args.sample, replace=False)
        expl = shap.TreeExplainer(model)
        sv = expl.shap_values(X[idx], check_additivity=False)
        # sklearn forests give one array per class, or a 3-D stack; take the positive class
        if isinstance(sv, list):
            sv = sv[1]
        elif getattr(sv, "ndim", 2) == 3:
            sv = sv[:, :, 1]

        mean_abs = np.abs(sv).mean(axis=0)
        mean_signed = sv.mean(axis=0)
        Xs = X[idx]
        for j in np.argsort(-mean_abs)[:200]:
            col = Xs[:, j]
            # direction: does a larger value of this feature push towards the positive class?
            # undefined when the column is constant across the explained sample
            rho = (float(spearmanr(col, sv[:, j]).statistic)
                   if np.unique(col).size > 1 else float("nan"))
            rows.append({"endpoint": ep, "feature_index": int(j), "feature": names[j],
                         "block": "descriptor" if j >= MORGAN_BITS else "fingerprint",
                         "mean_abs_shap": round(float(mean_abs[j]), 8),
                         "mean_signed_shap": round(float(mean_signed[j]), 8),
                         "direction_spearman_value_vs_shap": (round(rho, 4) if rho == rho
                                                              else None),
                         "n_explained": int(len(idx)), "seed": SEED})
        print(f"[{ep:6s}] explained {len(idx)} compounds; top feature "
              f"{names[int(np.argmax(mean_abs))]}", flush=True)

    sh = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    sh.to_csv(OUT / "shap_values.csv", index=False)

    # ---- agreement with permutation importance ------------------------------------------------
    perm_path = OUT / "feature_importance.csv"
    agree = []
    if perm_path.exists():
        perm = pd.read_csv(perm_path)
        for ep in eps:
            a = sh[sh.endpoint == ep].set_index("feature_index").mean_abs_shap
            b = perm[perm.endpoint == ep].set_index("feature_index").permutation_auroc_drop
            common = a.index.intersection(b.index)
            if len(common) >= 5:
                rho = spearmanr(a.loc[common], b.loc[common]).statistic
                agree.append({"endpoint": ep, "n_features_compared": int(len(common)),
                              "spearman_shap_vs_permutation": round(float(rho), 4)})
        pd.DataFrame(agree).to_csv(OUT / "shap_vs_permutation.csv", index=False)

    # ---- the plot: descriptors only, where the sign is interpretable --------------------------
    desc = sh[sh.block == "descriptor"]
    if len(desc):
        # direction, not mean signed SHAP: see the module docstring for why the latter is wrong here
        piv = desc.pivot_table(index="feature", columns="endpoint",
                               values="direction_spearman_value_vs_shap")
        order = desc.groupby("feature").mean_abs_shap.mean().sort_values(ascending=False).index
        piv = piv.reindex(order)
        lim = 1.0        # a Spearman correlation, so the scale is fixed and comparable
        fig, ax = plt.subplots(figsize=(1.2 + 0.62 * piv.shape[1], 0.34 * len(piv) + 1.5))
        im = ax.imshow(piv.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
        ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=45,
                                                               ha="right", fontsize=7)
        ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=7)
        ax.set_title("Mean signed SHAP, descriptor block\n"
                     "red pushes towards the positive class, blue away from it",
                     fontsize=8, loc="left")
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cb.ax.tick_params(labelsize=6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        fig.tight_layout()
        fig.savefig(OUT / "shap_summary.png", dpi=300)
        plt.close(fig)

    meta = {"commit": json.loads((OUT / "environment.json").read_text())["commit"],
            "seed": SEED, "sample": args.sample, "explainer": "shap.TreeExplainer (exact for trees)",
            "shap_version": shap.__version__,
            "environment_note": "shap installed after the first reproduction pass; the install was "
                                "additive only and did not move numpy, scikit-learn, pandas or "
                                "scipy, so no previously reproduced number is affected",
            "wall_clock_s": round(time.time() - t0, 1)}
    (OUT / "shap_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    if agree:
        print("\nSHAP vs permutation importance, Spearman over shared features:")
        for a in agree:
            print(f"  {a['endpoint']:6s} rho {a['spearman_shap_vs_permutation']:+.3f} "
                  f"({a['n_features_compared']} features)")
    print(f"\nwrote shap_values.csv, shap_summary.png, shap_vs_permutation.csv "
          f"({meta['wall_clock_s']}s)")


if __name__ == "__main__":
    main(sys.argv[1:])
