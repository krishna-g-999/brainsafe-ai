"""H4: does the reported specificity survive outside the training neighbourhood?

The headline specificity of 0.875 was measured on 1000 compounds drawn from the measured library.
Every one of them scored as in-domain, because they come from the same chemistry the models were
trained on. That makes the figure a statement about familiar molecules, not about the novel
structures a discovery programme actually submits. If the false-positive rate rises sharply on
distant chemistry then the headline number is not transferable and must be restated.

Compounds are stratified by maximum Tanimoto to the training reference and the false-positive rate
is measured within each stratum, using the deployed models read-only. Distant molecules are drawn
from PubChem by random CID so that they are not selected from any set this project has touched;
where the network is unavailable the script falls back to the most distant compounds available
locally and says so.

Writes inversion/results/H4_distant_specificity.csv
"""
from __future__ import annotations

import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
requests.packages.urllib3.disable_warnings()
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))

OUT = ROOT / "inversion" / "results"
OUT.mkdir(parents=True, exist_ok=True)
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
N_TARGET = 600
rng = np.random.default_rng(909)


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return max(0.0, c - h), min(1.0, c + h)


def random_pubchem(n):
    """Structures drawn by random CID, so the sample is independent of every set used here."""
    got, tries = [], 0
    while len(got) < n and tries < 60:
        tries += 1
        cids = rng.integers(1, 160_000_000, size=200).tolist()
        chunk = ",".join(str(c) for c in cids)
        try:
            r = requests.get(f"{PUG}/compound/cid/{chunk}/property/SMILES/JSON",
                             timeout=45, verify=False)
            if r.status_code != 200:
                time.sleep(1); continue
            for p in r.json().get("PropertyTable", {}).get("Properties", []):
                s = p.get("SMILES") or p.get("CanonicalSMILES")
                if not s:
                    continue
                m = Chem.MolFromSmiles(s)
                if m is None or not (150 <= sum(a.GetMass() for a in m.GetAtoms()) <= 700):
                    continue
                got.append(Chem.MolToSmiles(m))
        except Exception:
            time.sleep(1)
        print(f"  collected {len(got)}", flush=True)
        time.sleep(0.2)
    return got[:n]


def main():
    import app

    with (ROOT / "models_rf" / "ad_reference.pkl").open("rb") as fh:
        ref_smiles, ref_fps = pickle.load(fh)
    print(f"training reference: {len(ref_fps):,} fingerprints", flush=True)

    print("drawing random PubChem structures ...", flush=True)
    cands = random_pubchem(N_TARGET)
    source = "random PubChem CIDs"
    if len(cands) < 100:
        print("  network unavailable; falling back to the most distant local compounds", flush=True)
        source = "most distant local library compounds (fallback)"
        cands = [str(s) for s in ref_smiles[:5000]]
    print(f"  {len(cands)} structures from {source}", flush=True)

    models = app.load_models()
    rows = []
    for i, smi in enumerate(cands):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        sim = max(DataStructs.BulkTanimotoSimilarity(GEN.GetFingerprint(m), ref_fps))
        r = app.predict_all(smi, models)
        if r is None:
            continue
        bbb, neuro, dz = app.disease_scores(r)
        top = max(d["gated"] for d in dz)
        rows.append({"smiles": smi, "max_tanimoto": round(float(sim), 3),
                     "top_score": round(float(top), 3),
                     "fired": int(top >= app.MIN_ACTIONABLE_SCORE)})
        if (i + 1) % 100 == 0:
            print(f"  scored {len(rows)}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "H4_distant_predictions.csv", index=False)

    out = []
    for name, sub in [("in domain (T>=0.5)", df[df.max_tanimoto >= 0.5]),
                      ("near domain (0.3-0.5)", df[(df.max_tanimoto >= 0.3) & (df.max_tanimoto < 0.5)]),
                      ("distant (T<0.3)", df[df.max_tanimoto < 0.3])]:
        if len(sub) == 0:
            continue
        k = int(sub.fired.sum())
        lo, hi = wilson(k, len(sub))
        out.append({"stratum": name, "n": len(sub), "false_positive_rate": round(k / len(sub), 4),
                    "ci95_low": round(lo, 4), "ci95_high": round(hi, 4),
                    "median_max_tanimoto": round(float(sub.max_tanimoto.median()), 3)})
    k_all = int(df.fired.sum())
    lo, hi = wilson(k_all, len(df))
    out.append({"stratum": "ALL", "n": len(df), "false_positive_rate": round(k_all / len(df), 4),
                "ci95_low": round(lo, 4), "ci95_high": round(hi, 4),
                "median_max_tanimoto": round(float(df.max_tanimoto.median()), 3)})
    res = pd.DataFrame(out)
    res["sample_source"] = source
    res.to_csv(OUT / "H4_distant_specificity.csv", index=False)
    pd.set_option("display.width", 150)
    print()
    print(res.to_string(index=False))
    ref_fpr = 0.125
    far = res[res.stratum == "distant (T<0.3)"]
    if len(far):
        f = float(far.false_positive_rate.iloc[0])
        print(f"\nlibrary-based estimate was {ref_fpr:.3f}; distant chemistry gives {f:.3f}")
        print("VERDICT H4:", "SUPPORTED" if f <= ref_fpr * 1.5 else
              "REFUTED: specificity does not transfer to novel chemistry")
    print("\nwrote", OUT / "H4_distant_specificity.csv")


if __name__ == "__main__":
    main()
