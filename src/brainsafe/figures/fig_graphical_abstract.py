"""The graphical abstract NAR requires as a separate submission file.

A graphical abstract is read in seconds, at thumbnail size, by someone deciding whether to open the
paper at all. The idea it has to carry is the coupling in the title: a target score only counts if
the compound predicted to engage it is also predicted to reach it. That is drawn here as a mechanism
diagram, not a product screenshot: a real submitted structure (donepezil, the paper's own worked
example, drawn by RDKit rather than described) crossing a schematic blood-brain barrier, engaging a
schematic receptor on the other side, and resolving to the disease that mechanism drives. Every
number on the figure is read live from the same artefacts the manuscript cites.

The previous version of this figure used flat coloured cards with white text and four oversized
"hero" statistics, a layout borrowed from software marketing rather than from how a mechanism is
actually drawn in this literature. This version keeps the same house style as the other ten
figures (style.py): thin strokes, restrained colour used to carry meaning rather than to decorate,
and validation figures set as a caption-sized line rather than as headline numbers.

Output: manuscript/figures/GraphicalAbstract.png (and .pdf)

Run:  python src/brainsafe/figures/fig_graphical_abstract.py
"""
from __future__ import annotations

import io
import logging
import sys
from pathlib import Path

logging.getLogger("streamlit").setLevel(logging.ERROR)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, PathPatch, Wedge
from matplotlib.path import Path as MplPath

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))
import style as S  # noqa: E402

DONEPEZIL = "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2"


def molecule_image(smiles: str, size=(560, 420)):
    """The literal submitted structure, drawn rather than named."""
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    mol = Chem.MolFromSmiles(smiles)
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = d.drawOptions()
    opts.bondLineWidth = 2
    opts.clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return plt.imread(io.BytesIO(d.GetDrawingText()), format="png")


def membrane(ax, x, y0, y1, n=7):
    """Two parallel lines with a phospholipid-bilayer hatch: the barrier, drawn schematically."""
    for side in (-1, 1):
        xline = x + side * 0.006
        ax.plot([xline, xline], [y0, y1], color=S.MUTED, lw=1.1, zorder=2,
                transform=ax.transAxes)
    ys = np.linspace(y0 + 0.01, y1 - 0.01, n)
    for yy in ys:
        for side in (-1, 1):
            ax.plot([x + side * 0.006, x + side * 0.022], [yy, yy], color=S.HAIR, lw=2.2,
                    solid_capstyle="round", zorder=1, transform=ax.transAxes)


def receptor(ax, x, y, color, r=0.028):
    """A binding pocket, drawn as a notched circle rather than borrowed clip-art."""
    ax.add_patch(Wedge((x, y), r, 200, 160, width=r * 0.55, facecolor=color, edgecolor="none",
                       transform=ax.transAxes, zorder=2, alpha=0.9))
    ax.add_patch(Circle((x, y), r * 0.42, facecolor=color, edgecolor="white", linewidth=1.0,
                        transform=ax.transAxes, zorder=3))


def curved_arrow(ax, xy0, xy1, color, rad=-0.25, lw=1.3):
    ax.add_patch(FancyArrowPatch(xy0, xy1, transform=ax.transAxes, connectionstyle=f"arc3,rad={rad}",
                                 arrowstyle="-|>", mutation_scale=9, linewidth=lw, color=color,
                                 zorder=2))


def main() -> None:
    S.use()
    import app

    models = app.load_models()
    r = app.predict_all(DONEPEZIL, models)
    bbb, _neuro, dz = app.disease_scores(r)
    ache = r["targets"]["AChE"]
    top = dz[0]

    facts = app.panel_facts()
    shape = app.panel_shape()

    fig = plt.figure(figsize=(7.2, 4.55))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Explicit vertical grid. Every element is assigned one named level rather than tuned in
    # isolation, because tuning y-values one at a time is exactly how two lines end up sharing
    # the same 0.02 of figure height without either edit looking wrong by itself.
    Y_TITLE, Y_SUBTITLE, Y_RULE1 = 0.965, 0.922, 0.888
    Y_COLHEAD = 0.845
    DIAG_TOP, DIAG_BOT = 0.80, 0.46
    DIAG_MID = (DIAG_TOP + DIAG_BOT) / 2
    Y_VALUE = 0.375                 # predicted-exposure / AChE / disease-score labels
    Y_CAPTION = 0.305                # "submitted structure (SMILES or name)"
    Y_NOTE = 0.24                    # the gating sentence
    Y_RULE2 = 0.205
    Y_STATS1, Y_STATS2, Y_URL = 0.165, 0.10, 0.035

    ax.text(0.0, Y_TITLE, "BrainSafe AI", fontsize=15, fontweight="bold", ha="left", color=S.INK,
            transform=ax.transAxes)
    ax.text(0.0, Y_SUBTITLE, "does a molecule reach the brain, and what does it do once it is there?",
            fontsize=9, ha="left", color=S.MUTED, style="italic", transform=ax.transAxes)
    ax.plot([0.0, 1.0], [Y_RULE1, Y_RULE1], color=S.HAIR, lw=1.0, transform=ax.transAxes)

    # ---- the mechanism, drawn left to right --------------------------------------------------
    mol_ax = fig.add_axes([0.015, DIAG_BOT, 0.235, DIAG_TOP - DIAG_BOT])
    mol_ax.imshow(molecule_image(DONEPEZIL))
    mol_ax.axis("off")
    ax.text(0.13, Y_CAPTION, "submitted structure (SMILES or name)", ha="center", va="center",
            fontsize=7.4, color=S.MUTED, transform=ax.transAxes)

    curved_arrow(ax, (0.275, DIAG_MID), (0.345, DIAG_MID), S.MUTED, rad=-0.15)

    bx = 0.40
    membrane(ax, bx, DIAG_BOT + 0.03, DIAG_TOP - 0.03)
    ax.text(bx, Y_COLHEAD, "blood-brain barrier", ha="center", fontsize=7.6, color=S.MUTED,
            transform=ax.transAxes)
    ax.plot([0.345, bx - 0.03], [DIAG_MID, DIAG_MID], color=S.EXPOSURE, lw=1.6,
            transform=ax.transAxes, zorder=3)
    ax.add_patch(FancyArrowPatch((bx - 0.03, DIAG_MID), (bx + 0.05, DIAG_MID),
                                 transform=ax.transAxes, arrowstyle="-|>", mutation_scale=9,
                                 linewidth=1.6, color=S.EXPOSURE, zorder=3))
    ax.text(bx, Y_VALUE, f"predicted exposure {bbb:.2f}", ha="center", va="center", fontsize=8.3,
            color=S.EXPOSURE, fontweight="bold", transform=ax.transAxes)

    rx = 0.575
    receptor(ax, rx, DIAG_MID, S.TARGET)
    ax.text(rx, Y_COLHEAD, "target engagement", ha="center", fontsize=7.6, color=S.MUTED,
            transform=ax.transAxes)
    ax.text(rx, Y_VALUE, f"AChE {ache:.2f}", ha="center", va="center", fontsize=8.3,
            color=S.TARGET, fontweight="bold", transform=ax.transAxes)

    curved_arrow(ax, (rx + 0.035, DIAG_MID), (0.755, DIAG_MID), S.BINDER, rad=-0.15)
    ax.text(0.885, Y_COLHEAD, "disease relevance", ha="center", fontsize=7.6, color=S.MUTED,
            transform=ax.transAxes)
    ax.text(0.885, DIAG_MID + 0.03, top["disease"], ha="center", va="center", fontsize=8.6,
            color=S.INK, fontweight="bold", transform=ax.transAxes, linespacing=1.3, wrap=True)
    ax.text(0.885, Y_VALUE, f"{top['gated']:.2f}", ha="center", va="center", fontsize=10,
            color=S.BINDER, fontweight="bold", transform=ax.transAxes)

    ax.text(0.5, Y_NOTE,
            "a target score is admitted only in proportion to predicted exposure, so potency\n"
            "the compound cannot reach contributes nothing",
            fontsize=7.8, ha="center", va="center", color=S.MUTED, style="italic",
            transform=ax.transAxes, linespacing=1.5)

    ax.plot([0.0, 1.0], [Y_RULE2, Y_RULE2], color=S.HAIR, lw=1.0, transform=ax.transAxes)

    # ---- the evidence, set as a caption rather than as headline numbers ----------------------
    ax.text(0.0, Y_STATS1,
            f"{facts['n_records']:,} measured compound-endpoint records"
            f"  ·  {shape['targets']} molecular targets plus exposure, ADME and safety"
            f"  ·  {shape['deployed']} of {shape['trained']} trained estimators deployed"
            f"  ·  scaffold-held-out, calibrated and conformal",
            fontsize=7.6, ha="left", va="top", color=S.INK, transform=ax.transAxes, wrap=True)
    ax.text(0.0, Y_STATS2,
            "negative class recovered from measurement, not decoys  ·  thresholds set on a "
            "pool disjoint from the one that measures them  ·  every validation reported "
            "whichever way it falls",
            fontsize=7.4, ha="left", va="top", color=S.MUTED, transform=ax.transAxes, wrap=True)

    ax.text(0.0, Y_URL,
            "freely available, no registration · huggingface.co/spaces/Krishnag999/brainsafe-ai",
            fontsize=7.8, ha="left", va="center", color=S.EXPOSURE, fontweight="bold",
            transform=ax.transAxes)

    S.save(fig, "GraphicalAbstract")


if __name__ == "__main__":
    main()
