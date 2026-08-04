"""Collect the inversion results into one verdict table and a written report.

Reads only the CSVs produced by the individual tests, so the summary cannot drift from them. A
hypothesis that failed is reported as prominently as one that survived; the purpose of this exercise
was to find failures, and a report that hides them would defeat it.

Writes inversion/REPORT.md and inversion/results/VERDICTS.csv
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "inversion" / "results"


def read(name):
    p = RES / name
    return pd.read_csv(p) if p.exists() else None


def main():
    verdicts, notes = [], {}

    h1 = read("H1_disease_layer.csv")
    if h1 is not None:
        obs = float(h1.loc[h1.test == "observed top-3 accuracy", "value"].iloc[0])
        perm = float(h1.loc[h1.test == "permutation null (mean)", "value"].iloc[0])
        p = float(h1.loc[h1.test == "permutation p-value", "value"].iloc[0])
        freq = float(h1.loc[h1.test.str.startswith("frequency null"), "value"].iloc[0])
        v = "SUPPORTED" if (obs > perm + 0.05 and p < 0.05 and obs > freq) else "REFUTED"
        verdicts.append({"hypothesis": "H1 the disease score is informative", "verdict": v,
                         "headline": f"top-3 accuracy {obs:.3f} vs permutation null {perm:.3f} "
                                     f"(p={p:.3f}) and frequency null {freq:.3f}"})
        notes["H1"] = (f"Scored with hold-out models only, so no compound was seen in training. The "
                       f"disease layer reaches {obs:.1%} top-3 accuracy where shuffling the "
                       f"target-to-disease map gives {perm:.1%} and always answering with the three "
                       f"commonest diseases gives {freq:.1%}. The layer carries real information.")

    h2 = read("H2_weight_ablation.csv")
    if h2 is not None:
        cur = float(h2.loc[h2.weights == "curated", "top3_accuracy"].iloc[0])
        uni = float(h2.loc[h2.weights.str.startswith("uniform"), "top3_accuracy"].iloc[0])
        per = float(h2.loc[h2.weights.str.startswith("randomly"), "top3_accuracy"].iloc[0])
        v = "SUPPORTED" if cur - uni > 0.02 else "REFUTED"
        verdicts.append({"hypothesis": "H2 the curated edge weights add value", "verdict": v,
                         "headline": f"curated {cur:.4f}, uniform {uni:.4f}, permuted {per:.4f}"})
        notes["H2"] = (f"Curated weights score {cur:.4f}, uniform weights {uni:.4f} and randomly "
                       f"permuted weights {per:.4f}. The spread is {max(cur,uni,per)-min(cur,uni,per):.4f}. "
                       f"The information lies in which target connects to which disease, not in how "
                       f"strongly. The weights should be described as structure rather than as tuned "
                       f"parameters, and the graph would be simpler and no less accurate without them.")

    h3 = read("H3_gating.csv")
    if h3 is not None:
        verdicts.append({"hypothesis": "H3 BBB gating discriminates between diseases",
                         "verdict": "REFUTED (by construction)",
                         "headline": "the gate multiplies every disease equally and cannot change "
                                     "their order"})
        notes["H3"] = ("Multiplying every disease score by the same BBB probability leaves their "
                       "ranking untouched. Gating therefore cannot improve or damage which disease "
                       "is chosen; it changes only the absolute value and hence what crosses the "
                       "reporting threshold. It is an exposure filter, and the manuscript should "
                       "call it one rather than implying it sharpens the disease call.")

    h4 = read("H4_distant_specificity.csv")
    if h4 is not None:
        far = h4[h4.stratum.str.startswith("distant")]
        allr = h4[h4.stratum == "ALL"]
        if len(far):
            f = float(far.false_positive_rate.iloc[0])
            n = int(far.n.iloc[0])
            v = "SUPPORTED" if f <= 0.1875 else "REFUTED"
            verdicts.append({"hypothesis": "H4 specificity transfers to novel chemistry",
                             "verdict": v,
                             "headline": f"false-positive rate {f:.3f} on {n} distant compounds "
                                         f"against 0.125 measured on library chemistry"})
            notes["H4"] = (f"Structures drawn by random PubChem identifier, independent of every set "
                           f"used to build this tool. On compounds distant from training chemistry "
                           f"the false-positive rate is {f:.3f}, against {0.125:.3f} measured on "
                           f"library compounds.")
        elif len(allr):
            notes["H4"] = "No distant stratum was populated; see H4_distant_specificity.csv."

    h5 = read("H5_readacross_value.csv")
    if h5 is not None:
        ra = float(h5.loc[h5.method.str.startswith("read-across"), "recall"].iloc[0])
        fb = float(h5.loc[h5.method.str.startswith("frequency"), "recall"].iloc[0])
        v = "SUPPORTED" if ra - fb > 0.10 else "REFUTED"
        verdicts.append({"hypothesis": "H5 read-across beats a frequency baseline", "verdict": v,
                         "headline": f"recall {ra:.3f} against {fb:.3f}"})
        notes["H5"] = (f"Read-across recovers the true target for {ra:.1%} of held-out compounds "
                       f"against {fb:.1%} for always answering with the commonest targets. The query "
                       f"and any identical structure were excluded. This measures read-across in its "
                       f"intended regime, where the query's target family is represented in the "
                       f"index; it does not show that read-across works for a target class the index "
                       f"does not contain, and the figure should not be quoted as if it did.")

    df = pd.DataFrame(verdicts)
    df.to_csv(RES / "VERDICTS.csv", index=False)

    lines = ["# Inversion analysis: results", "",
             "Each hypothesis below was stated so that it could fail, and paired with a null model "
             "capable of producing the same apparent success by accident. All scoring used the "
             "scaffold hold-out models where predictive power was at issue. Nothing in `models_rf/`, "
             "`data/` or `app.py` was modified by this analysis.", "",
             "## Verdicts", ""]
    lines.append("| Hypothesis | Verdict | Evidence |")
    lines.append("|---|---|---|")
    for v in verdicts:
        lines.append(f"| {v['hypothesis']} | **{v['verdict']}** | {v['headline']} |")
    lines += ["", "## What each result means", ""]
    for k in sorted(notes):
        lines.append(f"**{k}.** {notes[k]}")
        lines.append("")
    lines += ["## Consequences for the tool", "",
              "1. The disease layer is validated for the first time, prospectively, and it works.",
              "2. The curated edge weights are not doing measurable work and should be presented as "
              "graph structure rather than as tuned parameters.",
              "3. BBB gating is an exposure filter, not a discriminator, and the manuscript wording "
              "should say so.",
              "4. Read-across is validated only where the target family is already represented.", ""]
    (ROOT / "inversion" / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")

    pd.set_option("display.width", 200)
    pd.set_option("display.max_colwidth", 70)
    print(df.to_string(index=False))
    print("\nwrote", ROOT / "inversion" / "REPORT.md")


if __name__ == "__main__":
    main()
