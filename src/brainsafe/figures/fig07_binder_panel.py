"""Figure 7. The binder panel, all 49 endpoints, including the ones that do not work.

The panel is not uniform and a single mean over it would hide that. These classifiers answer "does
this compound bind this target" for targets whose evidence ranges from thousands of measured actives
with hundreds of measured non-binders, to a few dozen actives and almost no measured negatives at
all. Performance tracks that, and the figure is ordered so it can be read off.

Panel A places every endpoint by what it discriminates and what it recovers. Panel B lists all of
them, so a reader can look up any target rather than take a panel average on trust, with the two
withdrawn endpoints and the ones that fail the reliability gate marked rather than omitted.

An endpoint is withdrawn when its probability band is too compressed for any threshold to separate
real ligands from trivial metabolites. That is a property of the fitted model, found by testing it
against chemistry it should reject, and it is reported here rather than being quietly dropped.

Reads models_rf/binder_modes.json.

Run:  python src/brainsafe/figures/fig07_binder_panel.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import style as S  # noqa: E402

BM = ROOT / "models_rf" / "binder_modes.json"


def table() -> pd.DataFrame:
    bm = json.loads(BM.read_text(encoding="utf-8"))
    rows = []
    for ep, v in bm.items():
        rows.append({
            "endpoint": ep,
            "auroc": v.get("auroc_vs_measured_inactives"),
            "sens": v.get("sensitivity_at_threshold"),
            "n_pos": v.get("n_positive"),
            "n_inact": v.get("n_measured_inactive") or v.get("n_measured_inactive_holdout"),
            "bg_fpr": v.get("background_fpr_held_out"),
            "reliable": bool(v.get("reliable_call", True)),
            "deployed": bool(v.get("deployed", True)),
            "measured_labels": v.get("mode") == "measured_labels_holdout",
        })
    return pd.DataFrame(rows).dropna(subset=["auroc", "sens"])


def style_of(r):
    if not r.deployed:
        return S.WARN, "withdrawn after specificity testing"
    if not r.reliable:
        return S.WITHHELD, "below the reliability gate"
    return S.TARGET, "deployed"


def panel_a(ax, d) -> None:
    # No shaded "bad" quadrant: an L-shaped wash covered most of the panel and asserted a cut that
    # the reliability gate already draws explicitly, in the marker colour.
    sizes = 4 + 22 * np.sqrt(d.n_pos / d.n_pos.max())
    seen = set()
    for (_, r), s in zip(d.iterrows(), sizes):
        col, lab = style_of(r)
        ax.scatter(r.auroc, r.sens, s=s, facecolor=col, edgecolor="white", linewidth=0.5,
                   alpha=0.9, zorder=3, label=lab if lab not in seen else None)
        seen.add(lab)
    # Label the endpoints a reader needs to identify, stepping the offset for any that would
    # otherwise print on top of a neighbour.
    weak = d[(d.sens < 0.8) | (d.auroc < 0.85) | (~d.deployed)].sort_values("sens")
    placed = []
    for _, r in weak.iterrows():
        dy = -1.5
        while any(abs(r.auroc - px) < 0.045 and abs(r.sens + dy / 260 - py) < 0.028
                  for px, py in placed):
            dy -= 7.0
        ax.annotate(r.endpoint, (r.auroc, r.sens), textcoords="offset points", xytext=(5.5, dy),
                    fontsize=6.5, color=style_of(r)[0])
        placed.append((r.auroc, r.sens + dy / 260))
    ax.set_xlabel("AUROC against measured non-binders")
    ax.set_ylabel("sensitivity at the triage threshold,\nheld-out actives by scaffold",
                  linespacing=1.6)
    # The lower bound follows the data. Fixed at 0.55 this panel silently dropped the three
    # endpoints withdrawn at AUROC 0.392, 0.479 and 0.539 off the left of the axis, so a panel
    # captioned as showing every endpoint showed 49 of 52, and the three it hid were the failures.
    lo = min(0.55, float(d.auroc.min()) - 0.03)
    ax.set_xlim(lo, 1.01); ax.set_ylim(-0.04, 1.04)
    S.strip(ax, x=True, y=True)
    ax.legend(loc="lower right", handletextpad=0.2, borderpad=0.3, labelspacing=0.35,
              scatterpoints=1)
    ax.text(0.02, 0.97, "marker area is the number of measured actives",
            transform=ax.transAxes, fontsize=6.5, color=S.MUTED, va="top")


def panel_b(ax, d) -> None:
    d = d.sort_values("auroc", ascending=True).reset_index(drop=True)
    ys = np.arange(len(d))
    for y, r in zip(ys, d.itertuples()):
        col, _ = style_of(r)
        ax.plot([r.sens, r.auroc], [y, y], color=S.HAIR, lw=1.3, solid_capstyle="round", zorder=1)
        ax.plot(r.auroc, y, "o", ms=2.9, mfc=col, mec="none", zorder=3)
        ax.plot(r.sens, y, "o", ms=2.9, mfc=S.BINDER, mec="none", alpha=0.75, zorder=3)

    labels = []
    for r in d.itertuples():
        mark = "" if r.deployed else "  (withdrawn)"
        if r.deployed and not r.reliable:
            mark = "  (below gate)"
        labels.append(f"{r.endpoint}{mark}")
    ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=6.5)
    for tick, r in zip(ax.get_yticklabels(), d.itertuples()):
        tick.set_color(style_of(r)[0] if (not r.deployed or not r.reliable) else S.MUTED)
    ax.set_xlabel("AUROC (coloured) and sensitivity (violet)")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.8, len(d) - 0.2)
    S.strip(ax, x=True, y=False)
    # Below the axis, not inside it: at the foot of the data area this sat on top of the four
    # weakest endpoints, which are the rows a reader most needs to be able to read.
    # The mean is quoted over the deployed set, because that is the panel the server offers and the
    # set the manuscript describes. Averaging in the withdrawn endpoints reports a number for a
    # panel nobody can query, and lowers it by exactly the failures that caused the withdrawal.
    dep = d[d.deployed]
    n_dep = len(dep)
    ax.text(0.0, -0.075,
            f"{len(d)} endpoints trained, {n_dep} deployed.   "
            f"Over the deployed set: mean AUROC {dep.auroc.mean():.3f}, "
            f"mean sensitivity {dep.sens.mean():.3f}.   "
            f"Median background false-positive rate {dep.bg_fpr.median():.4f}.",
            transform=ax.transAxes, fontsize=6.5, color=S.INK, va="top")


def main() -> None:
    S.use()
    d = table()
    fig = plt.figure(figsize=(S.DOUBLE, 6.83))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.16, 1.0], wspace=0.30,
                          left=0.085, right=0.985, top=0.905, bottom=0.115)
    a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1])
    S.panel(a, "A", "what each endpoint discriminates, and what it recovers", dx=-0.16, dy=1.045,
            gap=0.038)
    S.panel(b, "B", "every endpoint, named", dx=-0.30, dy=1.045, gap=0.045)
    panel_a(a, d); panel_b(b, d)
    S.save(fig, "Figure7_binder_panel")


if __name__ == "__main__":
    main()
