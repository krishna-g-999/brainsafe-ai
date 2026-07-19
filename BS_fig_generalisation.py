"""
BS_fig_generalisation.py
Main-text similarity-binned generalisation curve (reviewer point 1): AUROC as a
function of a query's nearest-neighbour Tanimoto to the training set, per endpoint,
with the n-weighted mean. Reads supplementary/STable5_similarity_binned_auroc.csv
(produced by BS_external_validation.py). Nothing hand-entered.
Output: figures/fig8_generalisation.png
"""
import os; os.chdir(os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
OI = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.titlesize": 12,
    "axes.titleweight": "bold", "figure.dpi": 300, "savefig.dpi": 300,
    "axes.spines.top": False, "axes.spines.right": False})
df = pd.read_csv("supplementary/STable5_similarity_binned_auroc.csv")
order = ["T[0.0-0.4)", "T[0.4-0.6)", "T[0.6-0.8)", "T[0.8-1.0)"]
xlab = ["<0.4", "0.4–0.6", "0.6–0.8", "≥0.8"]
EPS = ["BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_B", "MAO_A", "hERG"]
disp = {"MAO_A": "MAO-A", "MAO_B": "MAO-B", "GSK3B": "GSK-3β"}
fig, ax = plt.subplots(figsize=(8.2, 5.2))
for i, ep in enumerate(EPS):
    sub = df[df.endpoint == ep].set_index("tanimoto_bin").reindex(order)
    ax.plot(range(4), sub["AUROC"].values, marker="o", color=OI[i % 8], lw=1.8, ms=6, label=disp.get(ep, ep))
w = df[df.endpoint.isin(EPS)]; wm = []
for b in order:
    bb = w[w.tanimoto_bin == b]; wm.append((bb["AUROC"] * bb["n"]).sum() / bb["n"].sum())
ax.plot(range(4), wm, color="black", lw=3, ls="--", marker="s", ms=8, label="n-weighted mean", zorder=10)
ax.axhline(0.5, color="#999", lw=1, ls=":"); ax.text(2.6, 0.505, "chance", color="#777", fontsize=9)
ax.set_xticks(range(4)); ax.set_xticklabels(xlab); ax.set_xlim(-0.15, 3.15)
ax.set_xlabel("Nearest-neighbour Tanimoto to training set\n(← more novel chemistry            more familiar chemistry →)")
ax.set_ylabel("AUROC"); ax.set_ylim(0.5, 1.0)
ax.set_title("Generalisation versus chemical novelty (similarity-binned)")
ax.legend(fontsize=8.5, ncol=3, loc="lower right")
plt.tight_layout(); plt.savefig("figures/fig8_generalisation.png", bbox_inches="tight"); plt.close()
print("n-weighted mean AUROC by bin (<0.4 .. >=0.8):", [round(x, 3) for x in wm])
print("saved figures/fig8_generalisation.png")
