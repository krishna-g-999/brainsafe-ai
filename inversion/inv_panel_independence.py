"""H8: does a larger panel add information, or only more chances to fire?

Two claims were made when the panel was expanded and neither was tested.

The first was a prediction I made to argue against adding targets indiscriminately: every endpoint
added is another opportunity to fire spuriously, so the probability that at least one target fires on
a compound engaging nothing should rise with panel size. That is a falsifiable prediction about this
specific change, which added four endpoints and withdrew one, and it is measured here rather than
asserted.

The second is implicit in how the results are presented. The interface reports how many targets a
compound engages, and the network draws an edge for each. Both invite the reading that three engaged
targets are three independent pieces of evidence. Nav1.6, Nav1.7 and Nav1.8 are homologous sodium
channels; alpha4beta2, alpha3beta4 and alpha7 are homologous nicotinic receptors. If those fire
together on the same compounds, the count is one observation reported several times, and a user
weighing "engages three sodium channels" is being shown false corroboration.

Both are measured on the 600 random PubChem structures from H4, which no model in this project has
seen, plus the approved drugs from H6 for a set where real engagement exists.

Note on what this can and cannot show. Correlated firing among homologous targets is not by itself an
error: a genuinely promiscuous sodium-channel blocker should light up several sodium channels, and
that is a true fact about the compound. What would be an error is presenting that as independent
corroboration. The test therefore reports the correlation and leaves the interpretation explicit
rather than scoring it pass or fail.

A first version of this test defined firing as a target signal above zero and reported that at least
one target fires on 60% of random compounds. That figure was an artefact of the definition and is not
what a user sees. The six base-rate-enrichment endpoints return a positive signal whenever the
predicted probability exceeds the training base rate, which happens for about 12% of random compounds
each, whereas the 43 binder endpoints require a calibrated threshold to be crossed and fire for about
0.6%. Pooling the two counts an ordinary above-average probability as an event. The quantity that
reaches the user is whether any condition crosses the reporting threshold, and that is now the
headline; the signal-level rates are retained separately and labelled as internal.

Read-only. Writes inversion/results/H8_panel_independence.csv
"""
from __future__ import annotations

import json
import sys
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))

OUT = ROOT / "inversion" / "results"
OUT.mkdir(parents=True, exist_ok=True)

# families of homologous targets whose independence is in question
FAMILIES = {
    "voltage-gated sodium": ["Nav1_5", "Nav1_6", "Nav1_7", "Nav1_8"],
    "neuronal nicotinic": ["a7nAChR", "a4b2nAChR", "a3b4nAChR"],
    "monoamine transporters": ["SERT", "DAT", "NET"],
    "opioid": ["OPRM1", "OPRK1"],
    "dopamine": ["D2", "D3"],
    "serotonin GPCR": ["HT1A", "HT2A", "HT6", "HT7"],
}


def main():
    import app

    bg = pd.read_csv(OUT / "H4_distant_predictions.csv").smiles.tolist()
    drugs = pd.read_csv(OUT / "H6_clinical_indication_predictions.csv").smiles.tolist()[:400]
    models = app.load_models()
    live = [t for t in app.TARGET_KIND if t != "NEURO"]
    print(f"panel: {len(live)} scorable targets; random {len(bg)}, approved drugs {len(drugs)}",
          flush=True)

    def matrix(smis, label):
        rows, actionable = [], []
        for i, s in enumerate(smis, 1):
            r = app.predict_all(s, models)
            if r is None:
                continue
            neuro = app.neuro_signal(r)
            rows.append({t: float(app.target_signal(r, neuro, t) > 0) for t in live})
            _b, _n, dz = app.disease_scores(r)
            actionable.append(float(dz[0]["gated"] >= app.MIN_ACTIONABLE_SCORE))
            if i % 200 == 0:
                print(f"  {label}: {i}/{len(smis)}", flush=True)
        return pd.DataFrame(rows), np.array(actionable)

    print("scoring random chemistry ...", flush=True)
    B, act_bg = matrix(bg, "random")
    print("scoring approved drugs ...", flush=True)
    D, act_dr = matrix(drugs, "drugs")

    # ---- 1. how often does the server report a finding for a compound that should engage nothing?
    # This is the decision-relevant rate: it is what the user sees and what costs an experiment.
    report_bg = float(act_bg.mean())
    report_dr = float(act_dr.mean())
    binder_cols = [t for t in live if app.TARGET_KIND.get(t) == "binder"]
    enrich_cols = [t for t in live if app.TARGET_KIND.get(t) != "binder"]
    any_binder = float((B[binder_cols].sum(axis=1) > 0).mean()) if binder_cols else float("nan")
    per_binder = B[binder_cols].mean(axis=0) if binder_cols else pd.Series(dtype=float)
    indep = 1.0 - float(np.prod(1.0 - per_binder.values)) if len(per_binder) else float("nan")
    print(f"\nreported finding on random chemistry: {report_bg:.3f}   (approved drugs {report_dr:.3f})")
    print(f"   any binder endpoint above its threshold: {any_binder:.3f}")
    print(f"   independent-panel expectation for the binder endpoints: {indep:.3f}")
    print(f"   internal, not user visible: mean signal-positive rate "
          f"{B[enrich_cols].values.mean():.3f} for the {len(enrich_cols)} enrichment endpoints "
          f"against {B[binder_cols].values.mean():.3f} for the {len(binder_cols)} binder endpoints")

    # ---- 2. do homologous targets fire together? ----
    fam_rows = []
    for fam, members in FAMILIES.items():
        ms = [m for m in members if m in B.columns]
        if len(ms) < 2:
            continue
        for a, b in combinations(ms, 2):
            # phi coefficient on the drug set, where real engagement exists; the random set is too
            # sparse in positives for a correlation to be estimable
            x, y = D[a].values, D[b].values
            if x.std() == 0 or y.std() == 0:
                phi, joint = np.nan, float((x * y).mean())
            else:
                phi = float(np.corrcoef(x, y)[0, 1])
                joint = float((x * y).mean())
            # of compounds firing a, what fraction also fire b
            cond = float((x * y).sum() / x.sum()) if x.sum() else np.nan
            fam_rows.append({"family": fam, "target_a": a, "target_b": b,
                             "rate_a": round(float(x.mean()), 4),
                             "rate_b": round(float(y.mean()), 4),
                             "joint_rate": round(joint, 4),
                             "phi_correlation": round(phi, 3) if phi == phi else None,
                             "p_b_given_a": round(cond, 3) if cond == cond else None})
    fam = pd.DataFrame(fam_rows)
    fam.to_csv(OUT / "H8_family_correlation.csv", index=False)

    # ---- 3. effective panel size ----
    # A panel of N perfectly correlated targets carries the information of one. Kaiser's rule on the
    # correlation matrix of the firing indicators gives a blunt but honest count of how many
    # independent directions the panel actually spans on real drugs.
    used = [c for c in D.columns if D[c].std() > 0]
    eff = np.nan
    if len(used) > 2:
        C = np.corrcoef(D[used].values, rowvar=False)
        C = np.nan_to_num(C, nan=0.0)
        ev = np.linalg.eigvalsh(C)
        eff = float((ev > 1.0).sum())
        print(f"\ntargets firing at all on the drug set: {len(used)} of {len(live)}")
        print(f"   independent directions among them (eigenvalues above 1): {eff:.0f}")

    rows = [
        {"metric": "targets in panel", "value": len(live)},
        {"metric": "reported finding on random chemistry (user visible)", "value": round(report_bg, 4)},
        {"metric": "reported finding on approved drugs", "value": round(report_dr, 4)},
        {"metric": "any binder endpoint above threshold, random chemistry", "value": round(any_binder, 4)},
        {"metric": "independent-panel expectation, binder endpoints", "value": round(indep, 4)},
        {"metric": "mean binder endpoints fired per random compound",
         "value": round(float(B[binder_cols].sum(axis=1).mean()), 3)},
        {"metric": "mean binder endpoints fired per approved drug",
         "value": round(float(D[binder_cols].sum(axis=1).mean()), 3)},
        {"metric": "targets that ever fire on the drug set", "value": len(used)},
        {"metric": "independent directions in the firing pattern", "value": eff},
    ]
    pd.DataFrame(rows).to_csv(OUT / "H8_panel_independence.csv", index=False)
    pd.set_option("display.width", 200)
    print()
    print(pd.DataFrame(rows).to_string(index=False))
    if len(fam):
        print()
        print(fam.sort_values("phi_correlation", ascending=False).head(12).to_string(index=False))
        hi = fam[fam.phi_correlation.fillna(0) > 0.5]
        if len(hi):
            print(f"\n{len(hi)} homologous pairs fire together with correlation above 0.5. Their "
                  f"joint engagement is one observation, not several, and the interface should not "
                  f"present the count as independent corroboration.")
        else:
            print("\nNo homologous pair exceeds a correlation of 0.5: within these families the "
                  "targets fire largely independently, so an engaged count is not double counting.")
    print("\nwrote", OUT / "H8_panel_independence.csv")


if __name__ == "__main__":
    main()
