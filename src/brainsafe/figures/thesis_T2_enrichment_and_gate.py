"""Figure T2. The two transforms that are easiest to misread: enrichment, and the exposure gate.

Chapter 7 argues both in prose. Neither has ever been drawn, and both are geometric statements that a
picture settles faster than a paragraph.

A. The enrichment map. A calibrated probability is converted to a signed distance from the endpoint's
   own base rate, positive above and negative below, so that 0.60 at an endpoint where 24 per cent of
   compounds are active and 0.60 at one where 86 per cent are mean opposite things. The map is
   piecewise linear with a kink at p = b, and the three curves are drawn at three real base rates
   from the deployed panel rather than at round numbers.
B. Why that matters here. The eight core base rates span 0.2363 to 0.8621, so a single probability
   threshold applied across the panel would be a different question at every endpoint.
C. The gate cannot change a disease ranking. Multiplying the fourteen non-peripheral conditions by
   the same barrier probability is a common positive scaling, and a common positive scaling is
   rank-invariant. Donepezil's conditions are drawn at three values of the gate; the bars shrink
   together and the order never changes. This is H3, refuted by construction.
D. The exception no document in the project states. Migraine and multiple sclerosis are exempt from
   the gate, because their mechanisms act outside the barrier. For those two the multiplier is 1
   while the other fourteen are scaled, so the gate does change the ordering across that boundary.
   Teriflunomide is the case: a DHODH inhibitor used in multiple sclerosis, with a barrier
   probability of 0.560 that its own condition never sees.

Every number is computed by importing app and scoring the structures.

Run:  brainsafe_env/Scripts/python.exe src/brainsafe/figures/thesis_T2_enrichment_and_gate.py
Out:  thesis/figures/FigureT2_enrichment_and_gate.png and .pdf
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))
import style as S  # noqa: E402

OUTDIR = ROOT / "thesis" / "figures"
DONEPEZIL = "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2"
TERIFLUNOMIDE = r"C/C(O)=C(\C#N)C(=O)Nc1ccc(C(F)(F)F)cc1"
CORE = ["MAO_A", "hERG", "BChE", "AChE", "BBB", "GSK3B", "BACE1", "MAO_B"]
LABEL = {"AChE": "AChE", "BChE": "BChE", "BACE1": "BACE1", "GSK3B": "GSK-3β",
         "MAO_A": "MAO-A", "MAO_B": "MAO-B", "hERG": "hERG", "BBB": "BBB"}
SHORT = {"Neuroprotection / oxidative stress": "Neuroprotection",
         "Cognition (cholinergic)": "Cognition",
         "Alzheimer's disease": "Alzheimer's",
         "Depression / anxiety": "Depression"}
GAMMAS = [1.0, 0.6, 0.3]


def facts() -> dict:
    import app
    m = app.load_models()
    rates = {ep: app.base_rate(ep) for ep in CORE}

    rd = app.predict_all(DONEPEZIL, m)
    _, _, rows_d = app.disease_scores(rd)
    central = [r for r in rows_d if not r["peripheral"] and r["signal"] > 0]

    rt = app.predict_all(TERIFLUNOMIDE, m)
    bbb_t, _, rows_t = app.disease_scores(rt)
    ms = next(r for r in rows_t if r["disease"] == "Multiple sclerosis")

    return {"rates": rates, "central": central, "bbb_t": bbb_t, "ms": ms,
            "enrich": app.enrichment}


def panel_a(ax, F) -> None:
    p = np.linspace(0, 1, 601)
    picks = [("MAO_A", S.EXPOSURE), ("AChE", S.TARGET), ("BACE1", S.BINDER)]
    for ep, col in picks:
        b = F["rates"][ep]
        e = np.where(p >= b, (p - b) / (1 - b), (p - b) / b)
        ax.plot(p, e, color=col, linewidth=1.5, zorder=3,
                label=f"{LABEL[ep]},  b = {b:.3f}")
        ax.plot([b], [0], marker="o", markersize=3.6, color=col,
                markeredgecolor="white", markeredgewidth=0.5, zorder=4)
    ax.axhline(0, color=S.MUTED, linewidth=0.7, zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel("calibrated probability  $\\hat{q}_t$", fontsize=S.pt(7))
    ax.set_ylabel("enrichment  $E_t$", fontsize=S.pt(7))
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.legend(loc="upper left", fontsize=S.pt(6.5), handlelength=1.3, borderpad=0.2,
              labelspacing=0.35)
    S.strip(ax, x=True, y=True)
    ax.text(0.98, -0.94, "dot marks the kink at $p = b$", ha="right", va="bottom",
            fontsize=S.pt(6.5), color=S.MUTED)
    ax.text(0.985, 0.06, "engaged", ha="right", va="bottom", fontsize=S.pt(6.5),
            color=S.GOOD, style="italic")
    ax.text(0.985, -0.06, "clipped to zero", ha="right", va="top", fontsize=S.pt(6.5),
            color=S.FAINT, style="italic")


def panel_b(ax, F) -> None:
    order = sorted(F["rates"].items(), key=lambda kv: kv[1])
    y = list(range(len(order)))
    for i, (ep, b) in enumerate(order):
        ax.plot([0, b], [i, i], color=S.HAIR, linewidth=1.2, solid_capstyle="round", zorder=1)
        ax.plot([b], [i], marker="o", markersize=4.4, color=S.TARGET,
                markeredgecolor="white", markeredgewidth=0.5, zorder=3)
        ax.text(b + 0.028, i, f"{b:.3f}", va="center", ha="left", fontsize=S.pt(6.5),
                color=S.INK)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[ep] for ep, _ in order], fontsize=S.pt(7))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.95, len(order) - 0.35)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xlabel("base rate of the endpoint's own table", fontsize=S.pt(7))
    lo, hi = order[0][1], order[-1][1]
    ax.annotate("", xy=(lo, -0.68), xytext=(hi, -0.68),
                arrowprops=dict(arrowstyle="<->", linewidth=0.8, color=S.WITHHELD))
    ax.text((lo + hi) / 2, -0.38, f"spread {hi - lo:.3f}", ha="center", va="center",
            fontsize=S.pt(6.5), color=S.WITHHELD, fontweight="bold")
    S.strip(ax, x=True, y=False)


def panel_c(ax, F) -> None:
    central = sorted(F["central"], key=lambda r: r["signal"], reverse=True)
    n = len(central)
    x = np.arange(n)
    width = 0.26
    shades = [S.EXPOSURE, "#5A93C4", "#A9C6DE"]
    for k, (g, col) in enumerate(zip(GAMMAS, shades)):
        ax.bar(x + (k - 1) * width, [r["signal"] * g for r in central], width=width * 0.92,
               color=col, edgecolor="none", zorder=2,
               label=f"$\\gamma$ = {g:.1f}")
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT.get(r["disease"], r["disease"]) for r in central],
                       fontsize=S.pt(6.5))
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("gated condition score", fontsize=S.pt(7))
    ax.legend(loc="upper right", fontsize=S.pt(6.5), ncol=3, handlelength=1.0,
              columnspacing=0.9, borderpad=0.2)
    S.strip(ax)
    ax.text(0.0, -0.26, "the order never changes: a common positive scaling is rank-invariant",
            transform=ax.transAxes, fontsize=S.pt(6.8), color=S.BINDER, fontweight="bold")


def panel_d(ax, F) -> None:
    ms, g = F["ms"], F["bbb_t"]
    exempt, hypothetical = ms["gated"], ms["gated"] * g
    bars = [("exempt from the gate,\nas deployed", exempt, S.TARGET),
            ("if it were gated\nlike the other fourteen", hypothetical, S.HAIR)]
    y = [1, 0]
    for i, (lab, v, col) in zip(y, bars):
        ax.barh(i, v, height=0.44, color=col, edgecolor="none", zorder=2)
        ax.text(v + 0.02, i, f"{v:.3f}", va="center", ha="left", fontsize=S.pt(7.5),
                color=S.INK, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([b[0] for b in bars], fontsize=S.pt(6.8))
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.75, 1.75)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xlabel("multiple sclerosis score, teriflunomide", fontsize=S.pt(7))
    ax.axvline(0.30, color=S.WITHHELD, linewidth=0.9, zorder=4)
    ax.text(0.315, -0.62, "reporting threshold", fontsize=S.pt(6.5), color=S.WITHHELD,
            ha="left", va="center", fontweight="bold")
    S.strip(ax, x=True, y=False)
    ax.text(0.0, -0.26,
            f"barrier probability {g:.3f}, so the exemption is worth a factor of "
            f"{exempt / hypothetical:.2f} here.\nThe gate does reorder a peripheral condition "
            f"against a central one. No document states this.",
            transform=ax.transAxes, fontsize=S.pt(6.8), color=S.WARN, linespacing=1.45,
            va="top")


def main() -> None:
    S.use()
    F = facts()
    fig = plt.figure(figsize=(S.DOUBLE, 6.35))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.90], width_ratios=[1.0, 1.0],
                          hspace=0.58, wspace=0.34,
                          left=0.078, right=0.972, top=0.945, bottom=0.215)
    a = fig.add_subplot(gs[0, 0])
    b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[1, 0])
    d = fig.add_subplot(gs[1, 1])

    panel_a(a, F)
    panel_b(b, F)
    panel_c(c, F)
    panel_d(d, F)

    S.panel(a, "A", "the enrichment map", dx=-0.145, dy=1.045, gap=0.062)
    S.panel(b, "B", "why it is needed here", dx=-0.185, dy=1.045, gap=0.062)
    S.panel(c, "C", "the gate cannot reorder conditions", dx=-0.135, dy=1.055, gap=0.058)
    S.panel(d, "D", "except across the exemption", dx=-0.42, dy=1.055, gap=0.098)

    S.note(fig, "A and B use the deployed base rates from models_rf/endpoint_context.json. C uses "
                "donepezil, D teriflunomide, both scored through app. The scaling in C is applied "
                "to a real compound's scores; the gate's actual value for donepezil is 0.991.",
           y=0.020)
    S.save(fig, "FigureT2_enrichment_and_gate", outdir=OUTDIR)


if __name__ == "__main__":
    main()
