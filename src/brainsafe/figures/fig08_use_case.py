"""Figure 4 of the manuscript. A worked profile, and the behaviour that matters as much: silence.

NAR asks for a use case showing what the server returns for a real query. The useful demonstration
here is not a single spectacular hit but the pair of behaviours the design is built around. For an
approved CNS drug the server names the driving mechanism, not only the condition. For a compound that
acts peripherally it returns nothing above the reporting threshold, because a target score is
admitted only in proportion to predicted barrier penetration and a compound that does not arrive
cannot generate a disease call.

Both halves are scored by the deployed pipeline at run time, through the same entry point the web
interface uses, so this figure cannot drift from the server. The compounds are fixed here; their
scores are not.

Output: manuscript/figures/Figure8_use_case.png (and .pdf)

Run:  python src/brainsafe/figures/fig08_use_case.py
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

logging.getLogger("streamlit").setLevel(logging.ERROR)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))
import style as S  # noqa: E402

# Four approved CNS drugs whose primary mechanism is not in dispute, and four compounds that act
# outside the brain. Naming them here rather than sampling keeps the figure stable between runs.
CNS = [
    ("donepezil", "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2", "AChE"),
    ("haloperidol", "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1", "D2"),
    ("morphine", "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "OPRM1"),
    ("fluoxetine", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1", "SERT"),
]
PERIPHERAL = [
    ("atorvastatin", "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"),
    ("metformin", "CN(C)C(=N)NC(N)=N"),
    ("losartan", "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1"),
    ("hydrochlorothiazide", "NS(=O)(=O)c1cc2c(cc1Cl)NCNS2(=O)=O"),
]
REPORT_THRESHOLD = 0.30      # the score below which the server reports nothing


def main() -> None:
    S.use()
    import app

    models = app.load_models()
    rows = []
    for label, smi, expected in ([(n, s, d) for n, s, d in CNS]
                                 + [(n, s, None) for n, s in PERIPHERAL]):
        r = app.predict_all(smi, models)
        if r is None:
            raise SystemExit(f"{label} did not featurise; the figure would be describing nothing")
        bbb, _neuro, dz = app.disease_scores(r)
        top = dz[0] if dz else None
        # the driver reported by the server, not the one this file expected
        driver = (top["driver"][0] if top and top.get("driver") else None)
        rows.append({"compound": label, "expected": expected, "driver": driver,
                     "bbb": float(bbb),
                     "top": (top["disease"] if top else "none"),
                     "score": float(top["gated"]) if top else 0.0,
                     "cns": expected is not None})
        print(f"  {label:20s} BBB {bbb:.2f}  top {rows[-1]['top']} {rows[-1]['score']:.2f} "
              f"via {driver}", flush=True)

    # A CNS probe whose reported driver is not the pharmacologically expected one is a finding, not
    # something to draw over: say so on stdout so it cannot pass unnoticed into the figure.
    for r in rows:
        if r["cns"] and r["expected"] and r["driver"] != r["expected"]:
            print(f"  NOTE {r['compound']}: expected driver {r['expected']}, "
                  f"server reports {r['driver']}", flush=True)

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(S.DOUBLE, 2.9),
                                 gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.42})

    # ---- left: the disease call, gated by exposure ------------------------------------------
    ys = np.arange(len(rows))[::-1]
    for y, r in zip(ys, rows):
        col = S.TARGET if r["cns"] else S.FAINT
        w = max(r["score"], 0.0)
        ax.barh(y, w, height=0.62, color=col, alpha=0.85 if r["cns"] else 0.55, edgecolor="none")
        # below threshold the server reports nothing, so naming a driver there would describe an
        # internal value the user never sees
        called = r["score"] >= REPORT_THRESHOLD
        lab = (f"{r['top']}" + (f"  via {r['driver']}" if r["driver"] else "")) if called \
            else "no call above threshold"
        ax.text(max(w, 0.02) + 0.02, y, lab, va="center", fontsize=5.8,
                color=S.INK if called else S.MUTED)
    ax.axvline(REPORT_THRESHOLD, color=S.WARN, lw=0.9, ls=(0, (3, 2)), zorder=3)
    ax.text(REPORT_THRESHOLD, len(rows) - 0.35, " reporting threshold", fontsize=5.4,
            color=S.WARN, va="bottom")
    ax.set_yticks(ys); ax.set_yticklabels([r["compound"] for r in rows], fontsize=6.4)
    ax.set_xlim(0, 1.45); ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("top disease score, after exposure gating")
    S.strip(ax, x=True, y=False)

    # ---- right: exposure is what separates them ----------------------------------------------
    # Compounds cluster tightly in both corners, so labels are stepped apart within a cluster and
    # the high-exposure group is labelled to the left to stay inside the axes.
    placed = []
    for r in sorted(rows, key=lambda q: (-q["score"], -q["bbb"])):
        bx.plot(r["bbb"], r["score"], "o", ms=6 if r["cns"] else 5,
                mfc=S.TARGET if r["cns"] else S.FAINT, mec="white", mew=0.7, zorder=3)
        dy = 0.0
        while any(abs(r["bbb"] - px) < 0.18 and abs(r["score"] + dy / 200 - py) < 0.045
                  for px, py in placed):
            dy -= 9.0
        right = r["bbb"] < 0.6
        bx.annotate(r["compound"], (r["bbb"], r["score"]), textcoords="offset points",
                    xytext=(7 if right else -7, dy - 1.5), fontsize=5.2,
                    ha="left" if right else "right",
                    color=S.INK if r["cns"] else S.MUTED)
        placed.append((r["bbb"], r["score"] + dy / 200))
    bx.axhline(REPORT_THRESHOLD, color=S.WARN, lw=0.9, ls=(0, (3, 2)))
    bx.set_xlabel("predicted barrier penetration")
    bx.set_ylabel("top disease score", linespacing=1.6)
    bx.set_xlim(0, 1.05); bx.set_ylim(-0.05, 1.08)
    S.strip(bx, x=True, y=True)

    S.panel(ax, "A", "what the server returns", dx=-0.30, dy=1.06, gap=0.055)
    S.panel(bx, "B", "exposure gates the call", dx=-0.26, dy=1.06, gap=0.070)
    S.save(fig, "Figure8_use_case")


if __name__ == "__main__":
    main()
