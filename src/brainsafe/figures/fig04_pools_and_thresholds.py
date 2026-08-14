"""Figure 4. Where the negatives come from, and why the threshold is set on one pool and measured
on another.

A binder classifier needs negatives, and for most targets almost none are published: a paper reports
what bound, not what did not. The panel therefore trains against property-matched decoys drawn from
a background library, and that creates the failure this figure is about. If the decision threshold
is chosen as a quantile of a sample and the false-positive rate is then computed on that same sample,
the rate cannot exceed the quantile. It is a restatement of the target, not a measurement of it.

The background library is split into three disjoint pools by a stable hash of the structure, so a
compound's pool never depends on run order and no compound appears in two. Decoys are drawn from the
first, thresholds are set on the second, and the false-positive rate is measured on the third.

Panel C is the evidence the fix took: measured on a pool it was not set on, the background
false-positive rate now exceeds its 0.05 target for some endpoints, which under the old procedure
was arithmetically impossible.

Reads models_rf/binder_modes.json and the pool sizes from models.pools.

Run:  python src/brainsafe/figures/fig04_pools_and_thresholds.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
import style as S  # noqa: E402
from models.pools import SHARES  # noqa: E402

BM = ROOT / "models_rf" / "binder_modes.json"
BACKGROUND_TARGET = 0.05          # the intended background false-positive rate


def pool_sizes() -> dict[str, int]:
    """Actual pool sizes, taken from a binder record rather than recomputed.

    Every endpoint scores the same evaluation pool, so any record carries its size; the other two
    follow from the fixed share and are reported as such rather than as separate measurements.
    """
    bm = json.loads(BM.read_text(encoding="utf-8"))
    n_eval = next((v["n_background_evaluation"] for v in bm.values()
                   if v.get("n_background_evaluation")), None)
    if n_eval is None:
        raise SystemExit("no endpoint records n_background_evaluation; nothing to draw")
    total = round(n_eval / (SHARES["evaluation"] / sum(SHARES.values())))
    return {k: round(total * v / sum(SHARES.values())) for k, v in SHARES.items()} | {
        "evaluation": n_eval, "total": total}


def panel_a(ax, sizes) -> None:
    """The hash partition, and what each pool is used for."""
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.add_patch(FancyBboxPatch((0.02, 0.80), 0.96, 0.155,
                                boxstyle="round,pad=0,rounding_size=0.02",
                                facecolor="#F4F7F9", edgecolor=S.HAIR, lw=0.7))
    ax.text(0.5, 0.915, f"background library, {sizes['total']:,} compounds",
            ha="center", fontsize=7.2, color=S.INK, fontweight="bold")
    ax.text(0.5, 0.845, "assigned by blake2b(salt + canonical SMILES) mod 100, so a compound's pool "
                        "is a property of the\nstructure: it never depends on run order, and "
                        "re-running the split cannot move one",
            ha="center", va="center", fontsize=5.7, color=S.MUTED, linespacing=1.7)

    pools = [
        ("DECOY POOL", "decoy", S.BINDER, 0.02,
         "property-matched negatives\nfor training, Tanimoto <= 0.35\nto any active"),
        ("THRESHOLD POOL", "threshold", S.WITHHELD, 0.353,
         "the decision threshold is\nchosen here, together with\nheld-out measured inactives"),
        ("EVALUATION POOL", "evaluation", S.TARGET, 0.686,
         "the false-positive rate is\nmeasured here, on compounds\nneither step has seen"),
    ]
    for name, key, col, x, sub in pools:
        w = 0.294
        ax.add_patch(FancyBboxPatch((x, 0.30), w, 0.40,
                                    boxstyle="round,pad=0,rounding_size=0.02",
                                    facecolor=col, alpha=0.10, edgecolor=col, lw=0.8))
        ax.add_patch(Rectangle((x, 0.30), w, 0.014, facecolor=col, edgecolor="none"))
        ax.text(x + w / 2, 0.640, name, ha="center", fontsize=6.4, color=col, fontweight="bold")
        ax.text(x + w / 2, 0.545, f"{sizes[key]:,}", ha="center", fontsize=13, color=col,
                fontweight="bold")
        ax.text(x + w / 2, 0.495, f"{SHARES[key]} per cent", ha="center", fontsize=5.6,
                color=S.MUTED)
        ax.text(x + w / 2, 0.435, sub, ha="center", va="top", fontsize=5.6, color=S.MUTED,
                linespacing=1.7)
        ax.add_patch(FancyArrowPatch((0.5, 0.795), (x + w / 2, 0.706), arrowstyle="-|>",
                                     mutation_scale=7, color=S.FAINT, lw=0.9, shrinkA=0, shrinkB=0))

    ax.text(0.5, 0.215, "no compound appears in more than one pool, and the overlap was measured "
                        "rather than assumed",
            ha="center", fontsize=5.9, color=S.INK, fontweight="bold")
    ax.text(0.5, 0.145, "A threshold chosen on a sample and scored on that same sample returns the "
                        "quantile it was given.\nSeparating the two is what allows the measurement "
                        "to disagree with the target.",
            ha="center", va="top", fontsize=5.7, color=S.MUTED, linespacing=1.7)


def panel_b(ax, bm) -> None:
    """In-sample against held-out false-positive rate, per endpoint."""
    xs, ys, names = [], [], []
    for ep, v in bm.items():
        a, b = v.get("background_fpr_in_sample"), v.get("background_fpr_held_out")
        if a is None or b is None:
            continue
        xs.append(a); ys.append(b); names.append(ep)
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)

    hi = max(xs.max(), ys.max()) * 1.12
    ax.plot([0, hi], [0, hi], color=S.HAIR, lw=0.9, zorder=1)
    ax.axhline(BACKGROUND_TARGET, color=S.WARN, lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.axvline(BACKGROUND_TARGET, color=S.FAINT, lw=0.8, ls=(0, (3, 2)), zorder=2)

    over = ys > BACKGROUND_TARGET
    ax.plot(xs[~over], ys[~over], "o", ms=3.4, mfc=S.TARGET, mec="white", mew=0.5, alpha=0.85,
            zorder=3)
    ax.plot(xs[over], ys[over], "o", ms=4.6, mfc=S.WARN, mec="white", mew=0.6, zorder=4)
    for i in np.flatnonzero(over):
        ax.annotate(names[i], (xs[i], ys[i]), textcoords="offset points", xytext=(4.5, -1.2),
                    fontsize=5.3, color=S.WARN)

    ax.set_xlabel("false-positive rate on the pool the threshold was set on")
    ax.set_ylabel("measured on the held-out\nevaluation pool", linespacing=1.6)
    ax.set_xlim(0, hi); ax.set_ylim(0, hi)
    S.strip(ax, x=True, y=True)
    ax.text(0.03, 0.955, f"{int(over.sum())} of {len(xs)} endpoints exceed the "
                         f"{BACKGROUND_TARGET:.2f} target\nwhen measured on a pool they were "
                         "not tuned on",
            transform=ax.transAxes, fontsize=5.7, color=S.WARN, va="top", linespacing=1.7)
    ax.text(0.97, 0.06, "points above the diagonal are\nendpoints the in-sample rate flattered",
            transform=ax.transAxes, fontsize=5.5, color=S.MUTED, ha="right", va="bottom",
            linespacing=1.7)


def panel_c(ax, bm) -> None:
    """The two operating points every deployed endpoint carries."""
    rows = [(ep, v) for ep, v in bm.items() if v.get("deployed", True)
            and v.get("sensitivity_at_threshold") is not None
            and v.get("screening_sensitivity") is not None]
    rows.sort(key=lambda r: r[1]["sensitivity_at_threshold"])
    ys = np.arange(len(rows))
    tri = np.array([v["sensitivity_at_threshold"] for _, v in rows])
    scr = np.array([v["screening_sensitivity"] for _, v in rows])

    for y, a, b in zip(ys, tri, scr):
        ax.plot([a, b], [y, y], color=S.HAIR, lw=1.6, solid_capstyle="round", zorder=1)
    ax.plot(tri, ys, "o", ms=2.6, mfc=S.BINDER, mec="none", zorder=3, label="triage threshold")
    ax.plot(scr, ys, "o", ms=2.6, mfc=S.TARGET, mec="none", zorder=3, label="screening threshold")

    ax.set_yticks([]); ax.set_xlabel("sensitivity, held-out actives by scaffold")
    ax.set_ylabel(f"{len(rows)} deployed endpoints,\nordered by triage sensitivity",
                  linespacing=1.6)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-1, len(rows))
    S.strip(ax, x=True, y=False)
    ax.legend(loc="upper left", handletextpad=0.3, borderpad=0.2)
    ax.text(0.97, 0.06, f"median {np.median(tri):.2f} at triage,\n"
                        f"{np.median(scr):.2f} at screening",
            transform=ax.transAxes, fontsize=5.6, color=S.MUTED, ha="right", va="bottom",
            linespacing=1.7)


def main() -> None:
    S.use()
    bm = json.loads(BM.read_text(encoding="utf-8"))
    sizes = pool_sizes()

    fig = plt.figure(figsize=(S.DOUBLE, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.02, 1.0], width_ratios=[1.0, 1.0],
                          hspace=0.42, wspace=0.30,
                          left=0.075, right=0.985, top=0.925, bottom=0.085)
    a = fig.add_subplot(gs[0, :]); b = fig.add_subplot(gs[1, 0]); c = fig.add_subplot(gs[1, 1])
    S.panel(a, "A", "one background library, three pools that never overlap", dx=-0.012, dy=1.02,
            gap=0.030)
    S.panel(b, "B", "a rate that can now disagree with its target", dx=-0.20, dy=1.045, gap=0.048)
    S.panel(c, "C", "two operating points per endpoint", dx=-0.20, dy=1.045, gap=0.048)
    panel_a(a, sizes); panel_b(b, bm); panel_c(c, bm)
    S.save(fig, "Figure4_pools_and_thresholds")


if __name__ == "__main__":
    main()
