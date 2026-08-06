"""Act on the deployed-specificity audit.

The audit (results/deployed_specificity_audit.csv) tested every deployed binder model at its deployed
threshold against 600 random PubChem structures and against molecules no CNS target plausibly binds.
Four endpoints failed, and the four failures need three different remedies, so each is recorded here
with its evidence rather than applied by hand to a JSON file.

  Nav1_1     WITHDRAWN. Every trivial control lands in a band 0.801 to 0.816 against a threshold of
             0.796: glucose 0.806, urea 0.809, acetate 0.802, glycine 0.816. The model returns
             approximately 0.80 for anything small and polar, so no threshold separates a metabolite
             from a sodium-channel ligand; raising the cut to the value that would hold 5 percent on
             random chemistry still calls glycine a binder. This is the compressed-band pathology
             that also explains its deployed sensitivity of 0.32 in inversion H7 and its provenance,
             78 percent of measurements from a single assay. It is withdrawn from the panel rather
             than re-thresholded, because the defect is in the model and not in the cut. Nav1.6,
             added in this expansion, covers the axonal sodium-channel mechanism with a
             random-chemistry false-positive rate of 0.000.

  Cav3_2     REJECTED before deployment, same pathology in the opposite direction: the active band is
             compressed near zero, the calibrated threshold falls to 0.065, and atenolol scores 0.084.
             A T-type calcium endpoint that calls a beta-blocker a binder would add false leads to
             exactly the epilepsy queries this expansion exists to improve.

  a3b4nAChR  RETHRESHOLDED 0.392 to 0.450. Ethanol scored 0.437. Ethanol does modulate neuronal
             nicotinic receptors, but it is not a binder in the affinity sense this panel reports,
             and a tool that answers "ethanol engages alpha3beta4" will not be believed about
             anything else. The cost is sensitivity 0.985 to 0.975 and the random-chemistry
             false-positive rate improves from 0.028 to 0.018.

  SIRT1      RETHRESHOLDED 0.620 to 0.650. It sat at 0.058 on random chemistry, above the 5 percent
             its own calibration intends. No trivial control fired. Sensitivity is unchanged at
             0.933 and the false-positive rate falls to 0.033.

Writes models_rf/binder_modes.json and results/specificity_decisions.csv
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
M = ROOT / "models_rf"
OUT = ROOT / "results"

WITHDRAW = {
    "Nav1_1": "compressed probability band; every trivial metabolite scores 0.80 to 0.82 against a "
              "0.796 threshold, so no cut separates them from real ligands",
    "Cav3_2": "active band compressed near zero, calibrated threshold 0.065, atenolol scores 0.084",
}
RETHRESHOLD = {
    "a3b4nAChR": (0.450, "ethanol scored 0.437 at the calibrated threshold of 0.392"),
    "SIRT1": (0.650, "random-chemistry false-positive rate 0.058, above the 5 percent intended"),
    # RIPK1 shows the opposite of the Nav1.1 pathology and is therefore rescuable. Its actives sit at
    # a median probability of 1.000 and its measured inactives at 0.002, so discrimination is
    # excellent; the calibrated threshold fell to the permitted floor of 0.050 only because the 90th
    # percentile of those easy kinase-screening negatives is 0.006, and that floor happens to sit
    # inside the band where acetate, glycine and lactate score. Moving the cut to 0.50 removes every
    # trivial call at a cost of 0.013 in sensitivity, where no cut at all could save Nav1.1.
    "RIPK1": (0.500, "fired on acetate, glycine and lactate at the floor threshold of 0.050"),
}


def main():
    p = M / "binder_modes.json"
    shutil.copy(p, M / "binder_modes.prespecificity.json")
    modes = json.loads(p.read_text())
    rows = []

    for ep, why in WITHDRAW.items():
        if ep not in modes:
            continue
        old = dict(modes[ep])
        modes[ep].update({"deployed": False, "reliable_call": False,
                          "withdrawn_reason": why,
                          "withdrawn_evidence": "results/deployed_specificity_audit.csv"})
        rows.append({"target": ep, "action": "withdrawn",
                     "old_threshold": old.get("threshold"), "new_threshold": None, "reason": why})
        print(f"[{ep}] WITHDRAWN: {why}", flush=True)

    for ep, (thr, why) in RETHRESHOLD.items():
        if ep not in modes:
            continue
        old = modes[ep].get("threshold")
        modes[ep]["threshold"] = thr
        modes[ep]["threshold_basis"] = "random_chemistry_specificity"
        modes[ep]["rethreshold_reason"] = why
        # the high-precision mode must never be looser than the standard one
        s = modes[ep].get("screening_threshold")
        if s is None or s < thr:
            modes[ep]["screening_threshold"] = thr
        modes[ep]["deployed"] = True
        rows.append({"target": ep, "action": "rethresholded", "old_threshold": old,
                     "new_threshold": thr, "reason": why})
        print(f"[{ep}] threshold {old} -> {thr}: {why}", flush=True)

    for ep, v in modes.items():
        v.setdefault("deployed", True)

    p.write_text(json.dumps(modes, indent=2))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "specificity_decisions.csv", index=False)
    live = [k for k, v in modes.items() if v.get("deployed", True)]
    print(f"\ndeployed binder endpoints after these decisions: {len(live)} of {len(modes)}")
    print("withdrawn:", ", ".join(k for k, v in modes.items() if not v.get("deployed", True)))
    print("wrote", p, "and", OUT / "specificity_decisions.csv")


if __name__ == "__main__":
    main()
