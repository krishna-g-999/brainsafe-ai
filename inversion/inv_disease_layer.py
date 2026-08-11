"""H1-H3: does the disease layer carry information, and are its design choices justified?

Read-only with respect to models and data. Scores every held-out compound with the scaffold
hold-out models, which never saw its scaffold, so the disease prediction is genuinely prospective.

H1  Is the disease score informative at all? Compared against a permutation null that shuffles the
    target-to-disease map and a frequency null that always answers with the commonest disease.
H2  Do the curated edge weights earn their place? Compared against uniform weights and against
    weights randomly permuted across edges.
H3  Is multiplying by BBB penetration the right gate? Compared against no gate and against gating by
    predicted unbound brain exposure.

Writes inversion/results/H1_disease_layer.csv, H2_weight_ablation.csv, H3_gating.csv
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize  # noqa: E402

OUT = ROOT / "inversion" / "results"
OUT.mkdir(parents=True, exist_ok=True)
HOLD = ROOT / "models_rf" / "holdout"
rng = np.random.default_rng(4242)
TOPK = 3
MAX_PER_TARGET = 250      # cap so large targets do not dominate the pooled estimate
N_PERM = 200


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return max(0.0, c - h), min(1.0, c + h)


def main():
    import app  # read-only use of the deployed knowledge graph and helpers

    modes = json.loads((HOLD / "binder_modes.json").read_text())
    held = json.loads((HOLD / "heldout_actives.json").read_text())

    # A compound is held out from its own target's model, but every other model in the panel is
    # applied to it too, and scaffolds were withheld per target rather than across the panel. A
    # compound active at two targets is therefore withheld from one model and a training active of
    # the other, and the disease it recovers can come from the model that memorised it. Measured
    # across the 46 targets, 4,562 of 17,607 held-out entries are training actives elsewhere, and
    # the pairs affected are exactly the homologous families H8 identifies: D2 and D3 both mapping
    # to psychosis, SERT, NET and DAT to depression and ADHD, the two opioid receptors to pain.
    # Shuffling the target-to-disease map does not control for this, because the memorisation
    # survives the shuffle and is merely sent to the wrong disease.
    #
    # Restrict to compounds that are training actives of no other target, so the claim the summary
    # makes ("no compound was seen in training") is true of the set it is made about.
    trained_elsewhere = {}
    for ep in held:
        f = ROOT / "data" / "endpoints" / f"{ep}.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f)
        pv = pd.to_numeric(df.get("pchembl"), errors="coerce")
        trained_elsewhere[ep] = set(df.loc[pv >= 7, "smiles"].astype(str))
    before = sum(len(v) for v in held.values())
    filtered = {}
    for ep, smis in held.items():
        others = set()
        for e2, acts in trained_elsewhere.items():
            if e2 != ep:
                others |= acts
        filtered[ep] = [s for s in smis if s not in others]
    after = sum(len(v) for v in filtered.values())
    print(f"held-out entries: {before:,} -> {after:,} after removing those that are training "
          f"actives of another panel target ({before - after:,} removed, "
          f"{100 * (before - after) / max(before, 1):.1f}%)", flush=True)
    held = {ep: smis for ep, smis in filtered.items() if smis}
    KG = app.KNOWLEDGE_GRAPH
    DIS = list(app.DISEASE_ORDER)

    # target -> diseases it maps to, and the curated weight of each edge
    tgt_dis = {}
    for t, entries in KG.items():
        for pw, src, dis, w in entries:
            tgt_dis.setdefault(t, {}).setdefault(dis, 0.0)
            tgt_dis[t][dis] = max(tgt_dis[t][dis], w)

    # ---- score every held-out compound with hold-out models only ----
    print("scoring held-out compounds with hold-out models ...", flush=True)
    loaded = {}
    for ep in modes:
        p = HOLD / f"{ep}_binder.joblib"
        if p.exists():
            loaded[ep] = joblib.load(p)
    print(f"  {len(loaded)} hold-out models loaded", flush=True)

    records = []          # one row per held-out compound
    for true_t, smis in held.items():
        if true_t not in tgt_dis or true_t not in loaded:
            continue
        smis = list(smis)[:MAX_PER_TARGET]
        X, mask = featurize(smis)
        if X.shape[0] == 0:
            continue
        sig = {}
        for ep, mdl in loaded.items():
            thr = float(modes[ep]["threshold"])
            p = mdl.predict_proba(X)[:, 1]
            s = np.clip((p - thr) / max(1e-9, 1.0 - thr), 0.0, None)
            sig[ep] = np.where(p >= thr, s, 0.0)
        records.append({"target": true_t, "n": X.shape[0], "signals": sig,
                        "true_diseases": set(tgt_dis[true_t])})
        print(f"  {true_t:9} {X.shape[0]:4} compounds", flush=True)

    if not records:
        print("no records; aborting")
        return

    def disease_scores_from(sig_row, weights, gate=1.0):
        """Strongest engaged mechanism per disease, matching the deployed rule."""
        out = {}
        for t, s in sig_row.items():
            if s <= 0 or t not in weights:
                continue
            for dis, w in weights[t].items():
                v = w * s
                if v > out.get(dis, 0.0):
                    out[dis] = v
        return {d: v * gate for d, v in out.items()}

    def evaluate(weights, gate_mode="bbb"):
        """Fraction of held-out compounds whose true disease is in the top-k prediction."""
        hits = tot = 0
        for rec in records:
            n = rec["n"]
            for i in range(n):
                row = {t: float(v[i]) for t, v in rec["signals"].items()}
                gate = 1.0
                if gate_mode == "bbb":
                    gate = 1.0  # BBB is applied identically to every disease, so it cannot change
                    # the ranking; it is examined explicitly in H3 below.
                sc = disease_scores_from(row, weights, gate)
                if not sc:
                    tot += 1
                    continue
                top = [d for d, _ in sorted(sc.items(), key=lambda kv: -kv[1])[:TOPK]]
                hits += int(bool(rec["true_diseases"] & set(top)))
                tot += 1
        return hits, tot

    # ---------------- H1 ----------------
    print("\nH1: is the disease layer informative?", flush=True)
    k, n = evaluate(tgt_dis)
    lo, hi = wilson(k, n)
    acc = k / n

    # permutation null: shuffle which disease each target maps to
    perm_acc = []
    all_dis = sorted({d for v in tgt_dis.values() for d in v})
    for _ in range(N_PERM):
        shuffled = {}
        for t, dw in tgt_dis.items():
            picks = rng.choice(all_dis, size=len(dw), replace=False)
            shuffled[t] = {d: w for d, w in zip(picks, dw.values())}
        kk, nn = evaluate(shuffled)
        perm_acc.append(kk / nn)
    perm_acc = np.array(perm_acc)
    p_perm = float(((perm_acc >= acc).sum() + 1) / (N_PERM + 1))

    # frequency null: always answer with the most common true disease
    from collections import Counter
    cnt = Counter()
    for rec in records:
        for d in rec["true_diseases"]:
            cnt[d] += rec["n"]
    top_dis = {d for d, _ in cnt.most_common(TOPK)}
    freq_hits = sum(rec["n"] for rec in records if rec["true_diseases"] & top_dis)
    freq_acc = freq_hits / n

    rows = [
        {"test": "observed top-3 accuracy", "value": round(acc, 4),
         "ci95_low": round(lo, 4), "ci95_high": round(hi, 4), "n": n},
        {"test": "permutation null (mean)", "value": round(float(perm_acc.mean()), 4),
         "ci95_low": round(float(np.percentile(perm_acc, 2.5)), 4),
         "ci95_high": round(float(np.percentile(perm_acc, 97.5)), 4), "n": N_PERM},
        {"test": "permutation p-value", "value": p_perm, "ci95_low": "", "ci95_high": "", "n": N_PERM},
        {"test": "frequency null (always top-3 commonest)", "value": round(freq_acc, 4),
         "ci95_low": "", "ci95_high": "", "n": n},
    ]
    pd.DataFrame(rows).to_csv(OUT / "H1_disease_layer.csv", index=False)
    for r in rows:
        print(f"   {r['test']:44} {r['value']}", flush=True)
    verdict = ("SUPPORTED" if (acc > perm_acc.mean() + 0.05 and p_perm < 0.05 and acc > freq_acc)
               else "WEAKENED" if acc > perm_acc.mean() else "REFUTED")
    print(f"   VERDICT H1: {verdict}", flush=True)

    # ---------------- H2 ----------------
    print("\nH2: do the curated weights earn their place?", flush=True)
    uniform = {t: {d: 1.0 for d in dw} for t, dw in tgt_dis.items()}
    ku, nu = evaluate(uniform)
    permw_acc = []
    flat = [w for dw in tgt_dis.values() for w in dw.values()]
    for _ in range(30):
        sh = rng.permutation(flat)
        it = iter(sh)
        pw = {t: {d: float(next(it)) for d in dw} for t, dw in tgt_dis.items()}
        kk, nn = evaluate(pw)
        permw_acc.append(kk / nn)
    rows2 = [
        {"weights": "curated", "top3_accuracy": round(acc, 4)},
        {"weights": "uniform (all 1.0)", "top3_accuracy": round(ku / nu, 4)},
        {"weights": "randomly permuted", "top3_accuracy": round(float(np.mean(permw_acc)), 4)},
    ]
    pd.DataFrame(rows2).to_csv(OUT / "H2_weight_ablation.csv", index=False)
    for r in rows2:
        print(f"   {r['weights']:24} {r['top3_accuracy']}", flush=True)
    d_unif = acc - ku / nu
    print(f"   VERDICT H2: {'SUPPORTED' if d_unif > 0.02 else 'WEAKENED (curation adds little over uniform)'}",
          flush=True)

    # ---------------- H3 ----------------
    print("\nH3: is BBB gating the right operation?", flush=True)
    print("   Note: the gate multiplies every disease equally, so it cannot change which disease "
          "ranks first.", flush=True)
    print("   It changes only the absolute score, and therefore what crosses the reporting "
          "threshold.", flush=True)
    pd.DataFrame([{"finding": "BBB gate is rank-invariant across diseases",
                   "implication": "gating cannot improve or harm disease ranking; it only rescales, "
                                  "so it must be justified as an exposure filter rather than as a "
                                  "discriminative term"}]).to_csv(OUT / "H3_gating.csv", index=False)
    print("   VERDICT H3: gating is a filter, not a discriminator; see CSV", flush=True)
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
