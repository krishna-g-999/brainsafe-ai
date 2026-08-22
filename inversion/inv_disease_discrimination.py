"""H9: does the disease layer discriminate between compounds, or only echo which indications are common?

H6 reported top-3 accuracy of 0.352 on drugs never seen in training against a frequency null of
0.654, and concluded the layer is weakened. That comparison is unfair to the layer in one direction
and far too kind to the null in another, and this file measures the thing the comparison should have
measured.

The frequency null answers "chronic pain, depression, psychosis" for every compound it is shown. It
is right 65 per cent of the time because 40 per cent of approved CNS drugs treat chronic pain, not
because it knows anything. A constant predictor cannot rank one molecule against another, cannot
tell a user which of two candidates to make, and has no use at a bench. Top-k accuracy against it
rewards guessing the base rate, which is precisely what a triage tool must not do.

Two metrics are computed here that a constant predictor cannot pass:

  per-indication AUROC   for one condition, do the drugs that treat it score higher than the drugs
                         that do not? A constant predictor scores 0.5 by construction, whatever its
                         top-k accuracy, because it assigns every compound the same value.

  macro-averaged recall  the mean of per-indication recall rather than the pooled figure. A constant
                         predictor naming three conditions gets recall 1.0 on those and 0.0 on the
                         other thirteen, so it collapses towards 3/16 however common those three are.

Neither replaces H6. H6 asks "is the top-3 list right", which is a fair question about the output a
user reads, and its answer stands. These ask "does the layer respond to the compound", which is the
question that decides whether the layer is doing anything at all, and it is the one a frequency null
is built to fail.

Read-only. Writes inversion/results/H9_disease_discrimination.csv

Run:  python inversion/inv_disease_discrimination.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "inversion" / "results"
MIN_POSITIVES = 5          # below this an AUROC over a handful of drugs is not evidence


def auroc(y, s) -> float:
    """Rank-based AUROC, written here so the number does not come from the pipeline's own code."""
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    n1, n0 = int(y.sum()), int(len(y) - y.sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), float)
    sp, i = s[order], 0
    while i < len(sp):
        j = i
        while j + 1 < len(sp) and sp[j + 1] == sp[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def main() -> None:
    import app as A

    pred = pd.read_csv(OUT / "H6_clinical_indication_predictions.csv")
    pred = pred[pred.in_training == 0].reset_index(drop=True)
    print(f"scoring {len(pred)} approved drugs whose structures are absent from training", flush=True)

    models = A.load_models()
    diseases = list(A.DISEASE_ORDER)
    scores, kept = [], []
    for i, r in enumerate(pred.itertuples(), 1):
        res = A.predict_all(str(r.smiles), models)
        if res is None:
            continue
        _bbb, _neuro, dz = A.disease_scores(res)
        by = {d["disease"]: float(d["gated"]) for d in dz}
        scores.append([by.get(d, 0.0) for d in diseases])
        kept.append(r)
        if i % 40 == 0:
            print(f"  {i}/{len(pred)}", flush=True)

    S = np.asarray(scores)
    truth = np.zeros_like(S, dtype=int)
    for k, r in enumerate(kept):
        for lab in str(r.true_indications).split(" | "):
            lab = lab.strip()
            if lab in diseases:
                truth[k, diseases.index(lab)] = 1

    rows = []
    for j, d in enumerate(diseases):
        n_pos = int(truth[:, j].sum())
        if n_pos < MIN_POSITIVES:
            continue
        rows.append({"indication": d, "n_drugs_with_it": n_pos,
                     "auroc_model": round(auroc(truth[:, j], S[:, j]), 4),
                     "auroc_frequency_null": 0.5})
    per = pd.DataFrame(rows).sort_values("auroc_model", ascending=False)

    # Macro-averaged top-3 recall. The frequency null is evaluated exactly as it behaves: the same
    # three commonest indications for every compound.
    freq_order = [diseases[j] for j in np.argsort(-truth.sum(0))]
    top3_null = set(freq_order[:3])
    rec_model, rec_null = [], []
    for j, d in enumerate(diseases):
        pos = truth[:, j] == 1
        if pos.sum() < MIN_POSITIVES:
            continue
        top3 = np.argsort(-S, axis=1)[:, :3]
        hit = np.array([j in top3[k] for k in range(len(S))])
        rec_model.append(hit[pos].mean())
        rec_null.append(1.0 if d in top3_null else 0.0)

    mean_auroc = float(per.auroc_model.mean())
    print()
    print(per.to_string(index=False))
    print()
    print(f"  indications with at least {MIN_POSITIVES} drugs : {len(per)}")
    print(f"  mean per-indication AUROC, model               : {mean_auroc:.4f}")
    print(f"  mean per-indication AUROC, frequency null      : 0.5000  (by construction)")
    print(f"  indications where the model beats chance       : {int((per.auroc_model > 0.5).sum())}/{len(per)}")
    print()
    print(f"  macro-averaged top-3 recall, model             : {np.mean(rec_model):.4f}")
    print(f"  macro-averaged top-3 recall, frequency null    : {np.mean(rec_null):.4f}")

    OUT.mkdir(parents=True, exist_ok=True)
    per.to_csv(OUT / "H9_disease_discrimination.csv", index=False)
    summary = pd.DataFrame([{
        "metric": "mean per-indication AUROC", "model": round(mean_auroc, 4),
        "frequency_null": 0.5,
        "note": "a constant predictor scores 0.5 whatever its top-k accuracy"},
        {"metric": "macro-averaged top-3 recall", "model": round(float(np.mean(rec_model)), 4),
         "frequency_null": round(float(np.mean(rec_null)), 4),
         "note": "per-indication mean, so naming only the common conditions cannot carry it"}])
    summary.to_csv(OUT / "H9_disease_discrimination_summary.csv", index=False)
    verdict = ("SUPPORTED" if mean_auroc > 0.6 else
               "WEAKENED" if mean_auroc > 0.55 else "REFUTED")
    print(f"\n   VERDICT H9: {verdict}")
    print(f"\nwrote {OUT / 'H9_disease_discrimination.csv'}")


if __name__ == "__main__":
    main()
