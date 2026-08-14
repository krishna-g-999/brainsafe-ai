"""Figure 2. What the estimator is actually given, computed for a real molecule.

Reviewers asked, reasonably, what the input to these models is. The answer is a fixed-length numeric
vector, and this figure shows it for donepezil rather than describing it: the structure, the
substructure environments that set bits, the 1,024-bit fingerprint drawn as a grid with the set bits
marked, and the twelve descriptors with their computed values.

The figure also states the two properties of this representation that a reader should hold on to,
because both bound what the models can do. Folding means a bit is not a unique substructure. Leaving
chirality out means two enantiomers are one row, which is why the deduplication step exists at all.

Everything is computed here by the same featuriser the models use, imported rather than reimplemented.

Run:  python src/brainsafe/figures/fig02_feature_vector.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
import style as S  # noqa: E402
from features.featurize import (MORGAN_BITS, MORGAN_RADIUS, _DESCRIPTORS,  # noqa: E402
                                featurize_one)

# Donepezil: an approved AChE inhibitor, in the training data, and a molecule a CNS reader knows.
SMILES = "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2"
NAME = "donepezil"

DESC_LABEL = {
    "mw": "molecular weight", "clogp": "Crippen logP", "tpsa": "TPSA (A^2)",
    "hbd": "H-bond donors", "hba": "H-bond acceptors", "rotatable_bonds": "rotatable bonds",
    "aromatic_rings": "aromatic rings", "fraction_csp3": "fraction sp3 C",
    "ring_count": "rings", "heavy_atoms": "heavy atoms", "formal_charge": "formal charge",
    "qed": "QED drug-likeness",
}


def structure_image():
    """Draw the molecule with RDKit, or report that it could not be drawn."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        mol = Chem.MolFromSmiles(SMILES)
        d = rdMolDraw2D.MolDraw2DCairo(760, 470)
        opts = d.drawOptions()
        opts.bondLineWidth = 2
        opts.clearBackground = False
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
        d.FinishDrawing()
        import io
        return plt.imread(io.BytesIO(d.GetDrawingText()), format="png")
    except Exception as exc:                                  # pragma: no cover - drawing only
        print(f"  structure not drawn: {exc}")
        return None


def main() -> None:
    S.use()
    vec = featurize_one(SMILES)
    if vec is None:
        raise SystemExit(f"{NAME} did not parse; the figure would be describing nothing")
    vec = np.asarray(vec, dtype=float)
    bits, desc = vec[:MORGAN_BITS], vec[MORGAN_BITS:]
    on = np.flatnonzero(bits)

    fig = plt.figure(figsize=(S.DOUBLE, 4.05))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.06, 1.30, 0.92], wspace=0.20,
                          left=0.035, right=0.985, top=0.855, bottom=0.075)

    # ---- A: the molecule -------------------------------------------------------------------
    a = fig.add_subplot(gs[0]); a.axis("off")
    S.panel(a, "A", f"the molecule: {NAME}", dx=-0.02, dy=1.035, gap=0.075)
    img = structure_image()
    if img is not None:
        a.imshow(img)
    a.text(0.5, -0.045, SMILES, transform=a.transAxes, ha="center", fontsize=5.2, color=S.MUTED,
           family="monospace", wrap=True)
    a.text(0.5, -0.135, "standardised first: largest organic fragment, salts stripped,\n"
                        "sanitised, then keyed by the InChIKey of that parent",
           transform=a.transAxes, ha="center", va="top", fontsize=5.6, color=S.MUTED,
           linespacing=1.6)

    # ---- B: the fingerprint, all 1,024 bits ------------------------------------------------
    b = fig.add_subplot(gs[1]); b.axis("off")
    S.panel(b, "B", f"the fingerprint: {MORGAN_BITS} bits, {int(bits.sum())} set",
            dx=-0.03, dy=1.035, gap=0.052)
    side = int(np.sqrt(MORGAN_BITS))                       # 32 x 32
    grid = bits.reshape(side, side)
    b.imshow(grid, cmap=plt.matplotlib.colors.ListedColormap(["#EDF1F4", S.EXPOSURE]),
             vmin=0, vmax=1, interpolation="nearest", aspect="equal", extent=(0, side, side, 0))
    for k in range(0, side + 1, 8):
        b.axhline(k, color="white", lw=0.7); b.axvline(k, color="white", lw=0.7)
    b.set_xlim(0, side); b.set_ylim(side, 0)
    b.text(0.0, 1.008, "each cell is one bit, read left to right, top to bottom",
           transform=b.transAxes, fontsize=5.6, color=S.MUTED, va="bottom")
    b.text(0.0, -0.045, f"Morgan / ECFP-4, radius {MORGAN_RADIUS}, folded to {MORGAN_BITS} "
                        "bits, chirality NOT included.\n"
                        "Folding means a set bit reports that some environment hashing\n"
                        "to that index is present, not which one: several substructures\n"
                        "share a bit. Excluding chirality means two enantiomers produce\n"
                        "identical rows, which is why identical rows are collapsed before\n"
                        "any split rather than left to fall on both sides of one.",
           transform=b.transAxes, fontsize=5.6, color=S.MUTED, va="top", linespacing=1.7)

    # ---- C: the descriptors, with the values this molecule has -----------------------------
    c = fig.add_subplot(gs[2]); c.axis("off")
    S.panel(c, "C", f"the descriptors: {len(_DESCRIPTORS)} values", dx=-0.06, dy=1.035,
            gap=0.085)
    names = list(_DESCRIPTORS)          # dict, ordered as the featuriser concatenates them
    for i, (nm, val) in enumerate(zip(names, desc)):
        y = 0.945 - i * 0.0645
        c.add_patch(Rectangle((0.0, y - 0.020), 1.0, 0.052, transform=c.transAxes,
                              facecolor="#F4F7F9" if i % 2 == 0 else "white", edgecolor="none"))
        c.text(0.03, y, DESC_LABEL.get(nm, nm), transform=c.transAxes, fontsize=5.9,
               color=S.MUTED, va="center")
        c.text(0.97, y, f"{val:,.2f}" if abs(val) < 1e4 else f"{val:,.0f}",
               transform=c.transAxes, fontsize=5.9, color=S.TARGET, va="center", ha="right",
               fontweight="bold")
    c.text(0.0, 0.108, "Unscaled. A random forest splits on thresholds\nand is unchanged by any "
                       "monotone rescaling, so\nno scaler is fitted and none can leak across a\n"
                       "split.",
           transform=c.transAxes, fontsize=5.6, color=S.MUTED, va="top", linespacing=1.7)
    c.add_patch(FancyBboxPatch((0.0, -0.115), 1.0, 0.078, transform=c.transAxes,
                               boxstyle="round,pad=0,rounding_size=0.02", clip_on=False,
                               facecolor="#F4F7F9", edgecolor=S.HAIR, lw=0.6))
    c.text(0.5, -0.076, f"{MORGAN_BITS} + {len(names)} = "
                        f"{MORGAN_BITS + len(names):,} columns, every endpoint",
           transform=c.transAxes, ha="center", va="center", fontsize=6.6, color=S.INK,
           fontweight="bold")

    S.save(fig, "Figure2_feature_vector")


if __name__ == "__main__":
    main()
