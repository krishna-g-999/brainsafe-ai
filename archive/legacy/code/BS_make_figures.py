"""Generate publication figures strictly from the validated report files.
Every value is read from models_brain/*_meta.json (no hardcoded/fabricated numbers)."""
import os, glob, json, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

os.makedirs("figures", exist_ok=True)
M = {}
for mp in glob.glob("models_brain/*_meta.json"):
    m = json.load(open(mp)); M[m["endpoint"]] = m

DEPLOY = ["BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_B", "MAO_A", "hERG"]
ALL = sorted(M, key=lambda k: -M[k].get("mcc", 0))
plt.rcParams.update({"font.size": 10, "axes.spixne.top": False} if False else {"font.size": 10})

# ---- Fig 1: AUROC across validation regimes (honest generalisation) ----
fig, ax = plt.subplots(figsize=(8.5, 4.2))
x = np.arange(len(DEPLOY)); w = 0.26
sc = [M[e]["auroc"] for e in DEPLOY]
cl = [M[e].get("cluster_split_auroc") or np.nan for e in DEPLOY]
tm = [M[e].get("temporal_auroc") or np.nan for e in DEPLOY]
ax.bar(x - w, sc, w, label="Scaffold-CV", color="#2C7FB8")
ax.bar(x, cl, w, label="Leave-cluster-out", color="#7FCDBB")
ax.bar(x + w, tm, w, label="Temporal (future cpds)", color="#EDA63A")
ax.axhline(0.5, ls="--", c="grey", lw=0.8); ax.set_ylim(0.5, 1.0)
ax.set_xticks(x); ax.set_xticklabels(DEPLOY, rotation=0)
ax.set_ylabel("AUROC"); ax.legend(loc="lower left", fontsize=8, ncol=3)
ax.set_title("Figure 1. AUROC under increasingly strict validation (measured-data endpoints)")
for i, v in enumerate(tm):
    if not np.isnan(v): ax.text(x[i] + w, v + 0.01, f"{v:.2f}", ha="center", fontsize=7)
plt.tight_layout(); plt.savefig("figures/fig1_validation_regimes.png", dpi=300); plt.close()

# ---- Fig 2: Conformal calibration (empirical coverage vs 0.90 target) ----
fig, ax = plt.subplots(figsize=(7.5, 4.0))
cov = [M[e].get("conformal_coverage") or np.nan for e in DEPLOY]
bars = ax.bar(DEPLOY, cov, color="#41B6C4")
ax.axhline(0.90, ls="--", c="crimson", lw=1.2, label="Target coverage 0.90")
ax.set_ylim(0.80, 0.95); ax.set_ylabel("Empirical coverage (Mondrian conformal)")
ax.legend(fontsize=8); ax.set_title("Figure 2. Conformal prediction coverage (held-out)")
for b, v in zip(bars, cov):
    if not np.isnan(v): ax.text(b.get_x()+b.get_width()/2, v+0.002, f"{v:.3f}", ha="center", fontsize=7)
plt.tight_layout(); plt.savefig("figures/fig2_conformal_coverage.png", dpi=300); plt.close()

# ---- Fig 3: dataset size + class balance ----
fig, ax1 = plt.subplots(figsize=(8.0, 4.0))
n = [M[e]["n"] for e in DEPLOY]; pos = [M[e]["pos_rate"]*100 for e in DEPLOY]
ax1.bar(DEPLOY, n, color="#BDBDBD"); ax1.set_ylabel("n compounds")
ax2 = ax1.twinx(); ax2.plot(DEPLOY, pos, "o-", color="#D7301F"); ax2.set_ylabel("% positive (active/penetrant)", color="#D7301F")
ax2.set_ylim(0, 100); ax2.axhline(50, ls=":", c="grey", lw=0.7)
ax1.set_title("Figure 3. Training set size and class balance (deployed endpoints)")
for i, v in enumerate(n): ax1.text(i, v+80, str(v), ha="center", fontsize=7)
plt.tight_layout(); plt.savefig("figures/fig3_datasets.png", dpi=300); plt.close()

# ---- Fig 4: MCC quality gate (why receptor targets were excluded) ----
fig, ax = plt.subplots(figsize=(8.5, 4.0))
mcc = [M[e].get("mcc", 0) for e in ALL]
cols = ["#238B45" if (M[e].get("mcc",0) >= 0.45 or e in ("BBB","hERG")) else "#CB181D" for e in ALL]
ax.bar(ALL, mcc, color=cols)
ax.axhline(0.45, ls="--", c="black", lw=1.0, label="Deployment gate (MCC = 0.45)")
ax.set_ylabel("Matthews correlation coefficient (scaffold-CV)")
ax.legend(fontsize=8); ax.set_title("Figure 4. Fact-based quality gate: deployed (green) vs excluded receptor targets (red)")
for i, v in enumerate(mcc): ax.text(i, v+0.01, f"{v:.2f}", ha="center", fontsize=7)
plt.xticks(rotation=0); plt.tight_layout(); plt.savefig("figures/fig4_mcc_gate.png", dpi=300); plt.close()

print("Figures written:", os.listdir("figures"))
