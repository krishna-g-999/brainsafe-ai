"""The graphical abstract NAR requires as a separate submission file.

A graphical abstract is read in seconds, at thumbnail size, by someone deciding whether to open the
paper at all. It earns that attention by showing the one idea the paper adds, not by summarising the
whole results section. The idea here is the coupling in the title: a target score only counts if the
compound predicted to engage it is also predicted to reach it. Two panels carry that: a submitted
structure passing through the two questions the server answers, and the same real query (donepezil,
the paper's own worked example) shown reaching every value a reader would want to check against the
manuscript's own numbers, both read here from the same artefacts the manuscript cites rather than
retyped.

This replaces figures/graphical_abstract.png, which quoted 64,474 training records against the
228,200 in the current tables and predates the exposure-gating architecture the paper is actually
about. A stale graphical abstract is a worse first impression than none, since it is the one figure a
reader compares against nothing else in the paper before deciding whether to read further.

Output: manuscript/figures/GraphicalAbstract.png (and .pdf)

Run:  python src/brainsafe/figures/fig_graphical_abstract.py
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

logging.getLogger("streamlit").setLevel(logging.ERROR)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))
import style as S  # noqa: E402


def box(ax, x, y, w, h, text, fc, tc="white", fs=9.5, weight="bold"):
    b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.03",
                       linewidth=0, facecolor=fc, transform=ax.transAxes, zorder=2)
    ax.add_patch(b)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, fontweight=weight,
            color=tc, transform=ax.transAxes, linespacing=1.35, zorder=3)


def arrow(ax, x0, y, x1):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), transform=ax.transAxes,
                                 arrowstyle="-|>", mutation_scale=11,
                                 linewidth=1.4, color=S.MUTED, zorder=2))


def main() -> None:
    S.use()
    import app

    models = app.load_models()
    r = app.predict_all("COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2", models)  # donepezil
    bbb, _neuro, dz = app.disease_scores(r)
    ache = r["targets"]["AChE"]
    top = dz[0]

    facts = app.panel_facts()
    shape = app.panel_shape()
    n_records = facts["n_records"]
    n_targets = shape["targets"]
    n_estimators, n_deployed = shape["trained"], shape["deployed"]

    fig = plt.figure(figsize=(7.2, 4.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.955, "BrainSafe AI", fontsize=17, fontweight="bold", ha="center",
            color=S.INK, transform=ax.transAxes)
    ax.text(0.5, 0.89, "Does it reach the brain, and what does it do when it gets there?",
            fontsize=10.5, ha="center", color=S.MUTED, style="italic", transform=ax.transAxes)

    # ---- top row: the coupled question, worked on the paper's own example --------------------
    row_y, row_h = 0.60, 0.19
    box(ax, 0.02, row_y, 0.15, row_h, "submitted\nstructure\n(SMILES or name)", S.INK, fs=8.5)
    arrow(ax, 0.175, row_y + row_h / 2, 0.225)
    box(ax, 0.225, row_y, 0.20, row_h,
        f"predicted brain\nexposure\n{bbb:.2f}", S.EXPOSURE, fs=9)
    arrow(ax, 0.435, row_y + row_h / 2, 0.485)
    box(ax, 0.485, row_y, 0.24, row_h,
        f"target engagement,\ngated by exposure\nAChE {ache:.2f}", S.TARGET, fs=9)
    arrow(ax, 0.735, row_y + row_h / 2, 0.785)
    box(ax, 0.785, row_y, 0.21, row_h,
        f"disease relevance\n{top['disease']}\n{top['gated']:.2f}", S.BINDER, fs=8.7)

    ax.text(0.5, row_y - 0.045,
            "A target score is admitted only in proportion to predicted exposure, "
            "so potency a compound cannot reach contributes nothing.",
            fontsize=8.3, ha="center", color=S.MUTED, style="italic", transform=ax.transAxes)

    # ---- middle: what the panel is built on -----------------------------------------------
    stat_y = 0.345
    stats = [(f"{n_records:,}", "measured\ncompound-endpoint\nrecords"),
             (f"{n_targets}", "molecular targets,\nplus exposure, ADME\nand safety"),
             (f"{n_deployed}", f"deployed estimators\nof {n_estimators} trained"),
             ("scaffold-\nheld-out", "cross-validation,\ncalibrated and\nconformal")]
    w = 0.235
    for i, (val, lab) in enumerate(stats):
        x = 0.02 + i * (w + 0.007)
        ax.text(x + w / 2, stat_y + 0.10, val, ha="center", va="center", fontsize=14.5,
                fontweight="bold", color=S.INK, transform=ax.transAxes)
        ax.text(x + w / 2, stat_y - 0.02, lab, ha="center", va="top", fontsize=7.6,
                color=S.MUTED, transform=ax.transAxes, linespacing=1.3)

    ax.plot([0.02, 0.98], [0.245, 0.245], color=S.HAIR, lw=1.0, transform=ax.transAxes)

    # ---- bottom: the discipline, not a single headline metric -----------------------------
    ax.text(0.5, 0.185,
            "Validated, calibrated, and honest about its limits:", fontsize=9, fontweight="bold",
            ha="center", color=S.INK, transform=ax.transAxes)
    line2 = ("negative class recovered from measurement, not decoys   ·   thresholds set on "
             "a pool disjoint from the one that measures them   ·   every validation reported "
             "whichever way it falls")
    ax.text(0.5, 0.125, line2, fontsize=8.2, ha="center", color=S.MUTED,
            transform=ax.transAxes, wrap=True)

    ax.text(0.5, 0.045,
            "freely available, no registration · huggingface.co/spaces/Krishnag999/brainsafe-ai",
            fontsize=8.4, ha="center", color=S.EXPOSURE, fontweight="bold",
            transform=ax.transAxes)

    S.save(fig, "GraphicalAbstract")


if __name__ == "__main__":
    main()
