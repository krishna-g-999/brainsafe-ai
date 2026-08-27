"""Paired significance test on the model-family comparison.

`results/tables/model_comparison.csv` reports five families over thirteen endpoints under two split
regimes, and every document that cites it quotes point deltas alone. The decisions log makes exactly
this criticism of an earlier comparison ("point deltas alone do not establish significance") and
answers it with DeLong tests, but those were run on the superseded ensemble and were never re-run for
the deployed random forests. This script closes that gap for the current panel.

What it does and does not measure, stated plainly because the two are easy to confuse:

  It pairs the families by endpoint and asks, over the thirteen endpoints, whether the forest sits
  above an alternative more often than chance would give. The statistic is Wilcoxon signed-rank on
  the per-endpoint means. Because the panel mixes AUROC and R-squared, which are not the same
  quantity, the test is rank-based over "is the forest higher" and the median delta it reports is a
  summary of that ordering rather than an effect size with units.

  It is NOT a DeLong test. DeLong compares two ROC curves on identical compounds and accounts for
  their correlation, which is the stronger instrument where it applies. It applies per endpoint, on
  classification only, and would need the fold-level predictions rather than the summary table. That
  remains the better test and is not done here.

Run:  brainsafe_env/Scripts/python.exe src/brainsafe/evaluation/model_family_significance.py
Out:  results/tables/model_family_significance.csv
"""
from __future__ import annotations

import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results" / "tables" / "model_comparison.csv"
OUT = ROOT / "results" / "tables" / "model_family_significance.csv"

REFERENCE = "RandomForest"
CLASSIFIERS = {"BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_A", "MAO_B", "hERG"}


def load() -> dict[str, dict[str, dict[str, float]]]:
    out: dict[str, dict[str, dict[str, float]]] = defaultdict(lambda: defaultdict(dict))
    with open(SRC, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            out[r["split"]][r["endpoint"]][r["model"]] = float(r["mean"])
    return out


def compare(scores: dict[str, dict[str, float]], endpoints: list[str], alt: str) -> dict:
    ref = [scores[e][REFERENCE] for e in endpoints]
    other = [scores[e][alt] for e in endpoints]
    delta = [a - b for a, b in zip(ref, other)]
    # wilcoxon is undefined when every difference is zero; report it rather than crashing
    if all(d == 0 for d in delta):
        p = 1.0
    else:
        p = float(wilcoxon(ref, other).pvalue)
    return {
        "alternative": alt,
        "n_endpoints": len(endpoints),
        "forest_higher_on": sum(1 for d in delta if d > 0),
        "median_delta": round(st.median(delta), 4),
        "mean_delta": round(st.mean(delta), 4),
        "min_delta": round(min(delta), 4),
        "max_delta": round(max(delta), 4),
        "wilcoxon_p": round(p, 5),
        "distinguishable_at_0.05": p < 0.05,
    }


def main() -> None:
    data = load()
    rows = []
    for split in ("scaffold", "random"):
        scores = data[split]
        subsets = {
            "all 13 endpoints (AUROC and R-squared pooled by rank)": sorted(scores),
            "8 classification endpoints (AUROC only)": sorted(e for e in scores if e in CLASSIFIERS),
            "5 regression endpoints (R-squared only)": sorted(e for e in scores if e not in CLASSIFIERS),
        }
        for label, eps in subsets.items():
            alts = sorted({m for e in eps for m in scores[e]} - {REFERENCE})
            for alt in alts:
                rows.append({"split": split, "subset": label, **compare(scores, eps, alt)})

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {OUT.relative_to(ROOT).as_posix()}  ({len(rows)} comparisons)")
    print()
    print("scaffold split, all 13 endpoints:")
    for r in rows:
        if r["split"] == "scaffold" and r["subset"].startswith("all 13"):
            verdict = "distinguishable" if r["distinguishable_at_0.05"] else "NOT distinguishable"
            print(f"  forest vs {r['alternative']:<22} higher on {r['forest_higher_on']:>2}/"
                  f"{r['n_endpoints']}  median {r['median_delta']:+.4f}  "
                  f"p={r['wilcoxon_p']:.4f}  {verdict}")


if __name__ == "__main__":
    main()
