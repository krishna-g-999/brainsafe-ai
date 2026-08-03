"""Set every binder decision threshold from two independent requirements, computed from scratch.

Requirement A, selectivity against tested non-binders: at most TARGET_FPR of the target's measured
inactives may be called a binder. For hybrid models a portion of those inactives was seen in
training, so this component is mildly optimistic and is not used to make performance claims.

Requirement B, specificity against unrelated chemistry: at most BACKGROUND_FPR of a random sample of
the measured compound library may be called a binder. This sample is independent of every target's
training set, so it is the component that bounds false discovery when screening arbitrary compounds.
It exists because a prospective test found Nav1.1 assigning glucose a binder probability of 0.806
while scoring 0.979 against its own assay's negatives: a model judged only on its own chemistry can
be confidently wrong about everything else.

The operating threshold is the stricter of the two, recomputed from scratch each run so repeated
calibration cannot ratchet thresholds upward. BACKGROUND_FPR is 0.05; 0.02 was tested and rejected
because it drove the mu-opioid threshold to 0.999 and caused morphine to be missed.

Updates models_rf/binder_modes.json and writes results/tables/final_thresholds.csv
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize  # noqa: E402

M = ROOT / "models_rf"
TARGET_FPR = 0.10
BACKGROUND_FPR = 0.05
N_BACKGROUND = 3000
MIN_SENS = 0.50
rng = np.random.default_rng(7)


def main():
    with (M / "ad_reference.pkl").open("rb") as fh:
        bg_smiles, _ = pickle.load(fh)
    idx = rng.choice(len(bg_smiles), size=min(N_BACKGROUND, len(bg_smiles)), replace=False)
    Xbg, _ = featurize([str(bg_smiles[i]) for i in idx])
    print(f"independent background sample: {Xbg.shape[0]} compounds", flush=True)

    modes = json.loads((M / "binder_modes.json").read_text())
    rows = []
    for ep, v in modes.items():
        mp = M / f"{ep}_binder.joblib"
        if not mp.exists():
            continue
        mdl = joblib.load(mp)
        pbg = mdl.predict_proba(Xbg)[:, 1]
        thr_bg = float(np.quantile(pbg, 1.0 - BACKGROUND_FPR))

        thr_in, sens, n_in = None, None, 0
        f = ROOT / "data" / "endpoints" / f"{ep}.csv"
        pa = None
        if f.exists():
            df = pd.read_csv(f).dropna(subset=["smiles"])
            p = pd.to_numeric(df.get("pchembl"), errors="coerce")
            act = df.loc[p >= 7, "smiles"].astype(str).tolist()
            if len(act) < 30:
                act = df.loc[df["label"] == 1, "smiles"].astype(str).tolist()
            ina = df.loc[df["label"] == 0, "smiles"].astype(str).tolist()
            if len(ina) >= 15:
                Xi, _ = featurize(ina)
                pi = mdl.predict_proba(Xi)[:, 1]
                thr_in = float(np.quantile(pi, 1.0 - TARGET_FPR))
                n_in = len(pi)
            if act:
                Xa, _ = featurize(act)
                pa = mdl.predict_proba(Xa)[:, 1]

        thr = max([t for t in (thr_in, thr_bg) if t is not None])
        thr = float(np.clip(thr, 0.05, 0.999))
        if pa is not None and len(pa):
            sens = float((pa >= thr).mean())
        bgfpr = float((pbg >= thr).mean())

        v["threshold"] = round(thr, 4)
        v["threshold_basis"] = ("measured_inactives_and_background" if thr_in is not None
                                else "background_only")
        v["background_fpr_at_threshold"] = round(bgfpr, 4)
        v["n_measured_inactive_for_threshold"] = int(n_in)
        if sens is not None:
            v["sensitivity_at_threshold"] = round(sens, 3)
            v["reliable_call"] = bool(sens >= MIN_SENS
                                      and (v.get("auroc_vs_measured_inactives") or 1.0) >= 0.75)
        modes[ep] = v
        rows.append({"target": ep, "threshold": round(thr, 4),
                     "from_measured_inactives": round(thr_in, 4) if thr_in else None,
                     "from_background": round(thr_bg, 4),
                     "binding_constraint": ("background" if thr_in is None or thr_bg >= thr_in
                                            else "measured inactives"),
                     "background_fpr": round(bgfpr, 4),
                     "sensitivity": round(sens, 3) if sens is not None else None,
                     "reliable": v.get("reliable_call")})
        print(f"[{ep:8}] thr {thr:.3f} ({rows[-1]['binding_constraint']}) | bgFPR {bgfpr:.3f} | "
              f"sens {sens if sens is None else round(sens, 3)}", flush=True)

    (M / "binder_modes.json").write_text(json.dumps(modes, indent=2))
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results" / "tables" / "final_thresholds.csv", index=False)
    print(f"\nmedian background FPR {out.background_fpr.median():.3f} | "
          f"median sensitivity {out.sensitivity.median():.3f}")
    print("unreliable:", out.loc[out.reliable == False, 'target'].tolist())  # noqa: E712
    print("wrote", ROOT / "results" / "tables" / "final_thresholds.csv")


if __name__ == "__main__":
    main()
