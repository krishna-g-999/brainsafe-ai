"""H7: what does silence mean? Deployed sensitivity of every target in the panel.

H6 found that 47 of 58 approved antiepileptics score exactly zero, even though the panel contains
Nav1.1, GluA2, GABA-A, GluN2B, A1 and mGluR5, all of them epilepsy mechanisms. Probing those drugs
individually suggested a broken model: the Nav1.1 classifier returns 0.65 to 0.82 for nearly every
antiepileptic, from carbamazepine to acetazolamide to levetiracetam, against a deployed threshold of
0.796, and perampanel, a genuine AMPA antagonist, reaches 0.785 against a GluA2 threshold of 0.789.

That explanation is wrong, and this test is what shows it. Two failure modes were distinguished:

  the model cannot rank         held-out actives are not separated from random chemistry, AUROC
                                near 0.5, and no threshold can rescue it
  the operating point is strict the model ranks well but the deployed cut sits above a large share
                                of true actives, so sensitivity in deployment is low by design

Measured against 600 random PubChem structures with models that never saw the compounds scored,
every one of the 37 testable targets ranks at AUROC 0.91 or better and none is non-discriminative.
The first explanation is refuted. What the numbers show instead is the price of the second: the
thresholds hold the false-positive rate on random chemistry near 0.2 percent, and the cost is that
deployed sensitivity ranges from 0.26 to 0.98 with a wide spread across the panel.

The consequence for users is the one that matters. A silent result is not evidence that a compound
is inactive. For the strictest endpoints the tool misses roughly three true actives in four, and old
antiepileptics, being small, simple and low-affinity, are exactly the chemistry that sits just under
the cut.

Caveat stated rather than hidden: sensitivity applies the DEPLOYED threshold to the HOLD-OUT model.
The two are trained and calibrated identically on nearly the same data, so their probability scales
are close, but they are not the same estimator. The AUROC column does not depend on this. Endpoints
with fewer than about 30 held-out actives carry wide uncertainty on sensitivity and are flagged.

That caveat turned out to matter more than it reads. The false-positive rates below describe the
hold-out twins, and for Nav1.1 the twin and the deployed model diverge by a factor of fifty: 0.0017
here against 0.098 measured on the same 600 structures with the model the web server actually
served. The deployed Nav1.1 called glucose, urea and acetate binders while its twin did not. Nothing
in this file was wrong, but nothing in it licensed a claim about the deployed panel either, and the
distinction was not obvious enough. results/deployed_specificity_audit.csv now measures the served
models directly and is the authority for deployed specificity; this file is the authority for
whether an endpoint can rank at all. Nav1.1 was withdrawn from the panel on that evidence.

Read-only. Writes inversion/results/H7_target_discrimination.csv
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize  # noqa: E402

OUT = ROOT / "inversion" / "results"
HOLD = ROOT / "models_rf" / "holdout"
OUT.mkdir(parents=True, exist_ok=True)
MAX_ACTIVES = 400


def main():
    import app

    bg = pd.read_csv(OUT / "H4_distant_predictions.csv").smiles.tolist()
    Xbg, _ = featurize(bg)
    print(f"background: {Xbg.shape[0]} random PubChem structures", flush=True)

    held = json.loads((HOLD / "heldout_actives.json").read_text())
    rows = []
    for ep in sorted(held):
        p = HOLD / f"{ep}_binder.joblib"
        if not p.exists() or app.TARGET_KIND.get(ep) != "binder":
            continue
        smis = list(held[ep])[:MAX_ACTIVES]
        Xa, _ = featurize(smis)
        if Xa.shape[0] < 20:
            continue
        mdl = joblib.load(p)
        pa = mdl.predict_proba(Xa)[:, 1]
        pb = mdl.predict_proba(Xbg)[:, 1]
        thr = float(app.binder_threshold(ep))
        y = np.r_[np.ones(len(pa)), np.zeros(len(pb))]
        auroc = float(roc_auc_score(y, np.r_[pa, pb]))
        sens = float((pa >= thr).mean())          # deployed sensitivity
        fpr = float((pb >= thr).mean())           # deployed false-positive rate on random chemistry
        # how much of the probability scale the model actually uses on random chemistry: a model
        # that returns a narrow band cannot be thresholded usefully wherever that band sits
        spread = float(np.percentile(pb, 95) - np.percentile(pb, 5))
        rows.append({
            "target": ep, "n_heldout_actives": int(Xa.shape[0]),
            "auroc_vs_random": round(auroc, 4),
            "deployed_threshold": round(thr, 4),
            "deployed_sensitivity": round(sens, 4),
            "deployed_fpr_random": round(fpr, 4),
            "background_p5_p95_spread": round(spread, 4),
            "median_active_p": round(float(np.median(pa)), 4),
            "median_background_p": round(float(np.median(pb)), 4),
        })
        print(f"  {ep:9} AUROC {auroc:.3f}  thr {thr:.3f}  sensitivity {sens:.3f}  "
              f"FPR {fpr:.3f}", flush=True)

    df = pd.DataFrame(rows).sort_values("deployed_sensitivity")
    df.to_csv(OUT / "H7_target_discrimination.csv", index=False)

    weak = df[df.auroc_vs_random < 0.70]        # a model that cannot rank
    strict = df[df.deployed_sensitivity < 0.50]  # a model that ranks but rarely fires
    thin = df[df.n_heldout_actives < 30]         # too few actives to estimate sensitivity precisely
    print()
    print(df.to_string(index=False))
    print(f"\n{len(weak)} of {len(df)} targets fail to rank held-out actives against random chemistry "
          f"(AUROC below 0.70)")
    print(f"median deployed sensitivity {df.deployed_sensitivity.median():.3f}, "
          f"range {df.deployed_sensitivity.min():.3f} to {df.deployed_sensitivity.max():.3f}")
    print(f"median deployed false-positive rate on random chemistry "
          f"{df.deployed_fpr_random.median():.4f}")
    print(f"\n{len(strict)} targets fire for under half of their own held-out actives, so silence at "
          f"these endpoints carries little information:")
    for _, r in strict.iterrows():
        flag = "  (few held-out actives, wide uncertainty)" if r.n_heldout_actives < 30 else ""
        print(f"   {r.target:9} sensitivity {r.deployed_sensitivity:.3f} at threshold "
              f"{r.deployed_threshold:.3f}, AUROC {r.auroc_vs_random:.3f}{flag}")
    if len(thin):
        print(f"\n{len(thin)} endpoints have fewer than 30 held-out actives: "
              f"{', '.join(sorted(thin.target))}")

    # targets named in the knowledge graph that this test could not reach at all
    tested = set(df.target)
    graph_binders = {t for t in app.KNOWLEDGE_GRAPH if app.TARGET_KIND.get(t) == "binder"}
    untested = sorted(graph_binders - tested)
    if untested:
        print(f"\n{len(untested)} graph targets could not be tested for want of a hold-out model or "
              f"enough held-out actives: {', '.join(untested)}")

    # per-disease exposure to low-sensitivity endpoints
    rows2 = []
    sens = dict(zip(df.target, df.deployed_sensitivity))
    for d in app.DISEASE_ORDER:
        ts = {t for t, _ in app.DISEASE_CONTRIB.get(d, [])}
        vals = [sens[t] for t in ts if t in sens]
        rows2.append({"disease": d, "targets": len(ts), "tested": len(vals),
                      "untestable": len([t for t in ts if t in untested]),
                      "best_target_sensitivity": round(max(vals), 3) if vals else None,
                      "median_target_sensitivity": round(float(np.median(vals)), 3) if vals else None})
    pd.DataFrame(rows2).to_csv(OUT / "H7_disease_sensitivity.csv", index=False)
    print()
    print(pd.DataFrame(rows2).to_string(index=False))
    print("\nwrote", OUT / "H7_target_discrimination.csv")


if __name__ == "__main__":
    main()
