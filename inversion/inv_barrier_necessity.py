"""H10: does the barrier model earn its place over a descriptor rule?

The falsification suite has nine hypotheses and none of them tests the blood-brain barrier
classifier. H1, H2, H3, H6 and H9 test the disease layer, H5 tests read-across, H7 and H8 test the
target panel, and H4 tests system-level specificity. The BBB model gates every disease score and
names the architecture, and its only evidence is an external validation, which is a confirmation
rather than an attempt at refutation. A thesis audit named that gap; this closes it.

*Null: the fingerprint contributes nothing to barrier prediction. Permeability is a physicochemical
property, the twelve descriptors already in the feature vector encode it, and a rule or a small model
over those descriptors alone ranks compounds as well as a 1,036-column random forest does.*

The null is not a straw man. Barrier permeability is the endpoint in this project with the strongest
prior claim to being predictable from bulk properties: polar surface area, molecular weight and
hydrogen-bond donors are the terms of every published CNS-permeability heuristic. If the fingerprint
adds nothing here, the barrier model should be replaced by a rule that a reader can evaluate by hand,
and the project would be simpler and no less accurate.

Comparators, in increasing order of what they are allowed to know:

  single descriptor      TPSA, molecular weight, cLogP, hydrogen-bond donors and QED, each used
                         alone, with the sign that makes it a permeability predictor
  published heuristic    the common CNS rule of TPSA at most 90 and molecular weight at most 400,
                         as a binary call
  descriptor regression  logistic regression on all twelve standardised descriptors
  descriptor forest      a random forest on the twelve descriptors, same hyper-parameters as the
                         deployed model, so the only difference is the absence of the fingerprint

Two populations, because they answer different questions:

  scaffold hold-out      the same Bemis-Murcko folds the deployed cross-validation used, read from
                         the stored out-of-fold predictions so the fold assignment is identical
                         rather than merely similar. Answers: on the training distribution, does the
                         fingerprint help?
  external approved      the 241 FDA-curated drugs that are absent from B3DB and are also
                         distinguishable from training chemistry in feature space. Answers: on
                         chemistry the model has never seen, does the fingerprint still help?

A single descriptor's AUROC is reported as measured, which can fall below 0.5 when the sign is
against it; that is information about the descriptor, not an error, and it is left uncorrected.

Every comparator carries a paired bootstrap against the deployed model on the same compounds, because
the question is whether a margin of two or three AUROC points survives the sample it was measured on,
and on 241 external drugs that is not obvious.

Read-only. Writes inversion/results/H10_barrier_necessity.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize, feature_names  # noqa: E402

OUT = ROOT / "inversion" / "results"
OOF = ROOT / "data" / "processed" / "cv_predictions" / "BBB_scaffold_oof.csv"
EXT = ROOT / "data" / "external" / "processed" / "external_bbb_test.csv"
SEED = 42
RF_COMMON = dict(n_estimators=300, min_samples_leaf=2, n_jobs=-1, random_state=SEED,
                 class_weight="balanced")

# Sign that turns each descriptor into a permeability score: negative where more of it means less
# brain penetration.
SINGLE = {"tpsa": -1.0, "mw": -1.0, "clogp": +1.0, "hbd": -1.0, "qed": +1.0}


def descriptor_block(smiles: list[str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """The twelve descriptors only, and the mask of structures that featurised."""
    X, mask = featurize(smiles)
    names = feature_names()
    first = next(i for i, n in enumerate(names) if not n.startswith("ecfp4_"))
    return X[:, first:], mask, names[first:]


def wilson_delta_note(a: float, b: float) -> str:
    return f"{a - b:+.4f}"


def main() -> None:
    oof = pd.read_csv(OOF)
    D, mask, dnames = descriptor_block(oof["smiles"].tolist())
    oof = oof.loc[mask].reset_index(drop=True)
    y = oof["y_true"].to_numpy().astype(int)
    fold = oof["fold"].to_numpy()
    idx = {n: i for i, n in enumerate(dnames)}
    print(f"scaffold hold-out: {len(y):,} compounds, {int(y.sum()):,} permeable, "
          f"{len(set(fold))} folds", flush=True)

    ext = pd.read_csv(EXT)
    ext = ext[ext["novel_to_model"]].reset_index(drop=True)
    Xe_full, emask = featurize(ext["canonical_smiles"].tolist())
    ext = ext.loc[emask].reset_index(drop=True)
    De = Xe_full[:, Xe_full.shape[1] - len(dnames):]
    ye = ext["bbb_status"].to_numpy().astype(int)
    print(f"external approved drugs: {len(ye)} compounds, {int(ye.sum())} permeable", flush=True)

    # The deployed model is fitted on the whole B3DB table, so descriptor comparators for the
    # external arm are fitted on the same compounds: the out-of-fold set, which is that table
    # deduplicated. Nothing here sees an external compound before scoring it.
    rows: list[dict] = []

    scores: dict[tuple[str, str], np.ndarray] = {}

    def add(method, population, n, npos, auroc, note, score=None):
        rows.append({"method": method, "population": population, "n": int(n),
                     "n_permeable": int(npos),
                     "auroc": None if auroc is None else round(float(auroc), 4), "note": note})
        if score is not None:
            scores[(method, population)] = np.asarray(score, dtype=float)

    # ---- the deployed forest -------------------------------------------------------------------
    add("deployed forest, 1,036 features", "scaffold hold-out", len(y), y.sum(),
        roc_auc_score(y, oof["prediction"].to_numpy()),
        "stored out-of-fold predictions, not recomputed", oof["prediction"].to_numpy())
    deployed = joblib.load(ROOT / "models_rf" / "BBB.joblib")
    pe = deployed.predict_proba(Xe_full)[:, 1]
    add("deployed forest, 1,036 features", "external approved", len(ye), ye.sum(),
        roc_auc_score(ye, pe), "the served model", pe)

    # ---- single descriptors --------------------------------------------------------------------
    for name, sign in SINGLE.items():
        add(f"{name} alone", "scaffold hold-out", len(y), y.sum(),
            roc_auc_score(y, sign * D[:, idx[name]]),
            f"sign {'+' if sign > 0 else '-'}, no fitting of any kind", sign * D[:, idx[name]])
        add(f"{name} alone", "external approved", len(ye), ye.sum(),
            roc_auc_score(ye, sign * De[:, idx[name]]),
            f"sign {'+' if sign > 0 else '-'}, no fitting of any kind", sign * De[:, idx[name]])

    # ---- the published heuristic ---------------------------------------------------------------
    for label, Dm, yy in [("scaffold hold-out", D, y), ("external approved", De, ye)]:
        rule = ((Dm[:, idx["tpsa"]] <= 90.0) & (Dm[:, idx["mw"]] <= 400.0)).astype(float)
        add("CNS heuristic: TPSA <= 90 and MW <= 400", label, len(yy), yy.sum(),
            roc_auc_score(yy, rule),
            f"binary call, {rule.mean():.3f} of compounds pass it", rule)

    # ---- descriptor-only models, same folds ----------------------------------------------------
    for label, make in [("descriptor logistic regression",
                         lambda: LogisticRegression(max_iter=2000, class_weight="balanced")),
                        ("descriptor forest, 12 features",
                         lambda: RandomForestClassifier(**RF_COMMON))]:
        p = np.full(len(y), np.nan)
        for f in sorted(set(fold)):
            te = fold == f
            tr = ~te
            sc = StandardScaler().fit(D[tr])
            m = make()
            m.fit(sc.transform(D[tr]), y[tr])
            p[te] = m.predict_proba(sc.transform(D[te]))[:, 1]
        add(label, "scaffold hold-out", len(y), y.sum(), roc_auc_score(y, p),
            "same Bemis-Murcko folds as the deployed cross-validation", p)

        sc = StandardScaler().fit(D)
        m = make()
        m.fit(sc.transform(D), y)
        pex = m.predict_proba(sc.transform(De))[:, 1]
        add(label, "external approved", len(ye), ye.sum(), roc_auc_score(ye, pex),
            "fitted on the same compounds the deployed model was fitted on", pex)

    # Paired bootstrap: resample compounds, not predictions, so both models are always judged on
    # the same draw and the comparison keeps the pairing that makes it powerful.
    rng = np.random.default_rng(SEED)
    B = 2000
    labels = {"scaffold hold-out": y, "external approved": ye}
    draws = {pop: rng.integers(0, len(yy), size=(B, len(yy))) for pop, yy in labels.items()}
    for r in rows:
        pop, meth = r["population"], r["method"]
        base = scores.get(("deployed forest, 1,036 features", pop))
        this = scores.get((meth, pop))
        if base is None or this is None or meth.startswith("deployed"):
            r["delta_vs_deployed"] = 0.0 if meth.startswith("deployed") else None
            r["delta_ci95_low"] = r["delta_ci95_high"] = r["bootstrap_p_deployed_better"] = None
            continue
        yy = labels[pop]
        d = []
        for row in draws[pop]:
            yb = yy[row]
            if yb.min() == yb.max():
                continue
            d.append(roc_auc_score(yb, this[row]) - roc_auc_score(yb, base[row]))
        d = np.asarray(d)
        r["delta_vs_deployed"] = round(float(roc_auc_score(yy, this) - roc_auc_score(yy, base)), 4)
        r["delta_ci95_low"] = round(float(np.percentile(d, 2.5)), 4)
        r["delta_ci95_high"] = round(float(np.percentile(d, 97.5)), 4)
        r["bootstrap_p_deployed_better"] = round(float((d >= 0).mean()), 4)

    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.set_option("display.width", 200)

    print()
    for pop in ("scaffold hold-out", "external approved"):
        g = out[out.population == pop]
        deployed_auc = float(g[g.method.str.startswith("deployed")].auroc.iloc[0])
        best_null = g[~g.method.str.startswith("deployed")].sort_values("auroc").iloc[-1]
        print(f"{pop:20} deployed {deployed_auc:.4f}  best null "
              f"{best_null.auroc:.4f} ({best_null.method}), margin "
              f"{wilson_delta_note(deployed_auc, float(best_null.auroc))}")

    # The verdict reads the interval, not the point estimate. Chapter 9 of the thesis criticises H4
    # for deciding SUPPORTED on a point-estimate comparison whose interval contains its comparator,
    # and a hypothesis added in response to that criticism must not repeat it. On this rule the
    # deployed model has to beat the strongest null on both populations AND the paired bootstrap
    # interval has to exclude zero on both, or the verdict is only WEAKENED.
    decided = []
    for pop in ("scaffold hold-out", "external approved"):
        g = out[(out.population == pop) & (~out.method.str.startswith("deployed"))]
        best = g.sort_values("auroc").iloc[-1]
        decided.append((pop, float(best.delta_vs_deployed), float(best.delta_ci95_high),
                        str(best.method)))
    beats = all(d < 0 for _, d, _, _ in decided)
    separated = all(hi < 0 for _, _, hi, _ in decided)
    verdict = ("SUPPORTED" if beats and separated else
               "WEAKENED" if beats else
               "REFUTED: the fingerprint adds nothing a descriptor model does not")
    for pop, d, hi, meth in decided:
        print(f"   {pop:20} best null {meth}: delta {d:+.4f}, upper bound {hi:+.4f}"
              + ("" if hi < 0 else "   <-- interval includes zero"))

    # One rule, in the place that holds the data. Chapter 9 of the thesis found that H1, H2 and H5
    # each have their verdict computed twice, once in the test script and once in summarise.py, and
    # that for H2 the two copies disagree. The fix it recommends is one rule per hypothesis, so this
    # hypothesis writes its verdict into its own artefact and summarise.py reads it rather than
    # deriving it a second time.
    worst = max(decided, key=lambda t: t[2])
    out["verdict"] = verdict
    out["verdict_basis"] = (
        f"deployed model beats the strongest null on both populations; on {worst[0]} the paired "
        f"bootstrap interval for that margin has an upper bound of {worst[2]:+.4f}")
    out.to_csv(OUT / "H10_barrier_necessity.csv", index=False)
    print()
    print(out.drop(columns=["verdict", "verdict_basis"]).to_string(index=False))
    print(f"\n   VERDICT H10: {verdict}")
    print(f"\nwrote {(OUT / 'H10_barrier_necessity.csv').relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
