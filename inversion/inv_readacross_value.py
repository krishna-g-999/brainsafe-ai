"""H5: does read-across carry information beyond what a frequency baseline would give?

Read-across is presented to users as evidence when no model fires. That is only defensible if the
targets of a compound's nearest neighbours are genuinely predictive of its own target, rather than
simply reflecting which targets are common in the index. A large target with many actives will
appear among anyone's neighbours by sheer abundance.

The test uses the scaffold hold-out compounds, whose true target is known and whose scaffold is
absent from the index build only in the sense that the models never trained on it; to avoid a
trivially correct answer the query compound itself and anything identical to it are excluded from
the neighbour set. Performance is compared against a frequency baseline that always answers with the
commonest targets in the index.

Writes inversion/results/H5_readacross_value.csv
"""
from __future__ import annotations

import json
import pickle
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))

OUT = ROOT / "inversion" / "results"
OUT.mkdir(parents=True, exist_ok=True)
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
K = 5
TOPN = 3
MAX_PER_TARGET = 120
rng = np.random.default_rng(77)


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return max(0.0, c - h), min(1.0, c + h)


def main():
    idx = pickle.loads((ROOT / "models_rf" / "readacross_index.pkl").read_bytes())
    smiles, fps, targets = idx["smiles"], idx["fps"], idx["targets"]
    print(f"read-across index: {len(smiles):,} compounds", flush=True)

    freq = Counter(t for tg in targets for t in tg)
    common = {t for t, _ in freq.most_common(TOPN)}
    print(f"commonest targets in the index: {sorted(common)}", flush=True)

    held = json.loads((ROOT / "models_rf" / "holdout" / "heldout_actives.json").read_text())
    canon_index = {s: i for i, s in enumerate(smiles)}

    hit_ra = hit_freq = tot = 0
    per_target = []
    for true_t, smis in held.items():
        smis = list(smis)[:MAX_PER_TARGET]
        h = f = n = 0
        for smi in smis:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            q = GEN.GetFingerprint(m)
            sims = np.array(DataStructs.BulkTanimotoSimilarity(q, fps))
            # exclude the query itself and any identical structure, which would be circular
            if smi in canon_index:
                sims[canon_index[smi]] = -1.0
            sims[sims >= 0.999] = -1.0
            order = np.argsort(-sims)[:K]
            # weight each neighbour's targets by similarity, then take the top-scoring targets
            score = {}
            for i in order:
                if sims[i] < 0.35:
                    continue
                for t in targets[i]:
                    score[t] = score.get(t, 0.0) + float(sims[i])
            pred = {t for t, _ in sorted(score.items(), key=lambda kv: -kv[1])[:TOPN]}
            n += 1
            h += int(true_t in pred)
            f += int(true_t in common)
        if n:
            hit_ra += h; hit_freq += f; tot += n
            per_target.append({"target": true_t, "n": n,
                               "readacross_recall": round(h / n, 3),
                               "frequency_baseline": round(f / n, 3)})
            print(f"  {true_t:9} n={n:4}  read-across {h/n:.3f}  frequency {f/n:.3f}", flush=True)

    lo, hi = wilson(hit_ra, tot)
    flo, fhi = wilson(hit_freq, tot)
    rows = [
        {"method": "read-across (top-3 of 5 neighbours)", "recall": round(hit_ra / tot, 4),
         "ci95_low": round(lo, 4), "ci95_high": round(hi, 4), "n": tot},
        {"method": "frequency baseline (commonest targets)", "recall": round(hit_freq / tot, 4),
         "ci95_low": round(flo, 4), "ci95_high": round(fhi, 4), "n": tot},
    ]
    pd.DataFrame(rows).to_csv(OUT / "H5_readacross_value.csv", index=False)
    pd.DataFrame(per_target).to_csv(OUT / "H5_readacross_per_target.csv", index=False)
    print()
    for r in rows:
        print(f"   {r['method']:42} {r['recall']:.4f}  [{r['ci95_low']}, {r['ci95_high']}]")
    gain = hit_ra / tot - hit_freq / tot
    print(f"\n   gain over frequency baseline: {gain:+.4f}")
    print("   VERDICT H5:", "SUPPORTED" if gain > 0.10 else
          "WEAKENED" if gain > 0.02 else "REFUTED: no better than answering with common targets")
    print("\nwrote", OUT / "H5_readacross_value.csv")


if __name__ == "__main__":
    main()
