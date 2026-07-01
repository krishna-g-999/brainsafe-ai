"""
BS_figures.py — publication-grade figures generated from ACTUAL predictions + saved metrics.
Classifier ROC/reliability and regression scatter are recomputed via scaffold GroupKFold(5)
out-of-fold predictions (same protocol/seed as deployment); comparison panels use the saved
validated JSON metrics. Output: figures/*.png (300 dpi).
"""
import os, json, glob, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_curve, roc_auc_score, r2_score
from sklearn.isotonic import IsotonicRegression
from BS_predictive_model import morgan, descriptors, scaffold
from BS_train_endpoints import models as clf_models, canon as clf_canon
from BS_train_regression import reg_ensemble, canon as reg_canon

FIG = "figures"; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
                     "axes.spinetop": False} if False else {"font.size": 11})
DEPLOYED = ["BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_B", "MAO_A", "hERG"]
PALETTE = plt.cm.tab10(np.linspace(0, 1, 10))

rs = json.load(open("BS_randomsplit_benchmark.json"))
ext = json.load(open("BS_external_validation_report.json"))
metas = {os.path.basename(f).replace("_meta.json", ""): json.load(open(f))
         for f in glob.glob("models_brain/*_meta.json")}


def feats(smi):
    return np.hstack([morgan(smi), descriptors(smi)])


def clf_oof(ep):
    d = clf_canon(pd.read_csv(f"data/endpoints/{ep}.csv"))
    smi = d["smiles"].tolist(); y = d["label"].values.astype(int)
    X = feats(smi); g = np.array([scaffold(s) for s in smi]); oof = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, groups=g):
        oof[te] = np.mean([m.fit(X[tr], y[tr]).predict_proba(X[te])[:, 1] for m in clf_models().values()], axis=0)
    return y, oof


def reg_oof(ep, anti=False):
    if anti:
        df = pd.read_csv("data/endpoints_reg/antioxidant_dpph.csv").rename(columns={"y": "pchembl"})
        d = reg_canon(df.assign(label=0))
    else:
        d = reg_canon(pd.read_csv(f"data/endpoints/{ep}.csv"))
    smi = d["smiles"].tolist(); y = d["y"].values
    X = feats(smi); g = np.array([scaffold(s) for s in smi]); oof = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, groups=g):
        oof[te] = np.mean([m.fit(X[tr], y[tr]).predict(X[te]) for m in reg_ensemble()], axis=0)
    return y, oof


print("Recomputing classifier OOF (ROC/reliability)...")
CLF = {ep: clf_oof(ep) for ep in DEPLOYED}
print("Recomputing regression OOF (scatter)...")
REG = {"antioxidant": reg_oof("antioxidant", anti=True)}
for ep in ["A2A", "HT2A", "D2", "SERT"]:
    REG[ep] = reg_oof(ep)

# ---------- FIG 1: validation hierarchy ----------
fig, ax = plt.subplots(figsize=(11, 5.5))
splits = ["AUROC_random", "auroc", "cluster_split_auroc", "temporal_auroc"]
labels = ["Random", "Scaffold", "Cluster (leave-out)", "Temporal (future)"]
x = np.arange(len(DEPLOYED)); w = 0.2
for i, (s, lab) in enumerate(zip(splits, labels)):
    vals = [rs.get(ep) if s == "AUROC_random" else metas[ep].get(s) for ep in DEPLOYED]
    vals = [v if v is not None else np.nan for v in vals]
    ax.bar(x + (i - 1.5) * w, vals, w, label=lab, color=PALETTE[i])
ax.axhline(0.5, ls="--", c="grey", lw=1, label="random chance")
ax.set_xticks(x); ax.set_xticklabels([e.replace("_", "-") for e in DEPLOYED]); ax.set_ylim(0.5, 1.0)
ax.set_ylabel("AUROC"); ax.set_title("Figure 1. Validation hierarchy: AUROC across split rigour (BrainSafe AI)")
ax.legend(ncol=5, fontsize=9, loc="lower center", bbox_to_anchor=(0.5, -0.22)); plt.tight_layout()
plt.savefig(f"{FIG}/Fig1_validation_hierarchy.png", dpi=300, bbox_inches="tight"); plt.close()

# ---------- FIG 2: ROC curves ----------
fig, ax = plt.subplots(figsize=(7, 7))
for i, ep in enumerate(DEPLOYED):
    y, p = CLF[ep]; fpr, tpr, _ = roc_curve(y, p)
    ax.plot(fpr, tpr, lw=2, color=PALETTE[i], label=f"{ep.replace('_','-')} (AUROC {roc_auc_score(y,p):.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=1)
ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
ax.set_title("Figure 2. Scaffold-CV ROC curves (out-of-fold)"); ax.legend(fontsize=9, loc="lower right")
plt.tight_layout(); plt.savefig(f"{FIG}/Fig2_roc_curves.png", dpi=300); plt.close()

# ---------- FIG 3: reliability / calibration ----------
fig, ax = plt.subplots(figsize=(7, 7))
for i, ep in enumerate(["BBB", "AChE", "BACE1", "hERG"]):
    y, p = CLF[ep]; pc = IsotonicRegression(out_of_bounds="clip").fit(p, y).predict(p)
    bins = np.linspace(0, 1, 11); idx = np.digitize(pc, bins) - 1
    xs, ys = [], []
    for b in range(10):
        m = idx == b
        if m.sum() > 5: xs.append(pc[m].mean()); ys.append(y[m].mean())
    ax.plot(xs, ys, "o-", color=PALETTE[i], label=ep.replace("_", "-"))
ax.plot([0, 1], [0, 1], "k--", lw=1, label="perfect calibration")
ax.set_xlabel("Predicted probability (calibrated)"); ax.set_ylabel("Observed frequency")
ax.set_title("Figure 3. Calibration reliability (isotonic, scaffold-CV)"); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig3_calibration.png", dpi=300); plt.close()

# ---------- FIG 4: conformal coverage ----------
fig, ax = plt.subplots(figsize=(9, 4.8))
cov = [metas[ep].get("conformal_coverage") for ep in DEPLOYED]
ax.bar([e.replace("_", "-") for e in DEPLOYED], cov, color=PALETTE[3])
ax.axhline(0.90, ls="--", c="red", label="target 90% coverage")
ax.set_ylim(0.8, 1.0); ax.set_ylabel("Empirical coverage")
ax.set_title("Figure 4. Conformal prediction: empirical coverage vs 90% target"); ax.legend()
plt.tight_layout(); plt.savefig(f"{FIG}/Fig4_conformal_coverage.png", dpi=300); plt.close()

# ---------- FIG 5: benchmark vs literature ----------
fig, ax = plt.subplots(figsize=(9, 5))
lit = {"BBB": (0.88, 0.96), "hERG": (0.86, 0.93), "AChE": (0.90, 0.97), "BACE1": (0.90, 0.96),
       "MAO_B": (0.88, 0.96), "MAO_A": (0.85, 0.95), "BChE": (0.90, 0.97), "GSK3B": (0.88, 0.95)}
for i, ep in enumerate(DEPLOYED):
    lo, hi = lit[ep]; ax.plot([i, i], [lo, hi], lw=8, color="lightsteelblue", solid_capstyle="round")
    ax.plot(i, rs[ep], "D", color="crimson", ms=9, zorder=5)
ax.plot([], [], lw=8, color="lightsteelblue", label="published random-split range")
ax.plot([], [], "D", color="crimson", label="BrainSafe (random split)")
ax.set_xticks(range(len(DEPLOYED))); ax.set_xticklabels([e.replace("_", "-") for e in DEPLOYED])
ax.set_ylim(0.82, 1.0); ax.set_ylabel("AUROC")
ax.set_title("Figure 5. BrainSafe vs published state of the art (like-for-like random splits)"); ax.legend()
plt.tight_layout(); plt.savefig(f"{FIG}/Fig5_benchmark_vs_sota.png", dpi=300); plt.close()

# ---------- FIG 6: regression scatter ----------
fig = plt.figure(figsize=(15, 6)); gs = GridSpec(2, 3, figure=fig)
order = ["antioxidant", "A2A", "HT2A", "D2", "SERT"]
titles = {"antioxidant": "Antioxidant (DPPH pIC50)", "A2A": "A2A (pKi)", "HT2A": "5-HT2A (pKi)",
          "D2": "D2 (pKi)", "SERT": "SERT (pKi)"}
for i, ep in enumerate(order):
    axx = fig.add_subplot(gs[i // 3, i % 3]); y, p = REG[ep]
    axx.scatter(y, p, s=6, alpha=0.3, color=PALETTE[i])
    lo, hi = min(y.min(), p.min()), max(y.max(), p.max())
    axx.plot([lo, hi], [lo, hi], "k--", lw=1)
    axx.set_title(f"{titles[ep]}  R²={r2_score(y,p):.2f}"); axx.set_xlabel("measured"); axx.set_ylabel("predicted (OOF)")
fig.suptitle("Figure 6. Regression endpoints: predicted vs measured (scaffold-CV OOF)", y=1.02)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig6_regression_scatter.png", dpi=300, bbox_inches="tight"); plt.close()

# ---------- FIG 7: dataset overview ----------
fig, ax = plt.subplots(figsize=(9, 5))
ns = [metas[ep]["n"] for ep in DEPLOYED]; pos = [metas[ep]["pos_rate"] * 100 for ep in DEPLOYED]
ax.bar([e.replace("_", "-") for e in DEPLOYED], ns, color="steelblue", label="n compounds")
ax.set_ylabel("n compounds"); ax2 = ax.twinx(); ax2.plot([e.replace("_", "-") for e in DEPLOYED], pos, "ro-", label="% active")
ax2.set_ylabel("% active"); ax2.set_ylim(0, 100)
ax.set_title("Figure 7. Dataset size and class balance per endpoint (measured ChEMBL/B3DB)")
plt.tight_layout(); plt.savefig(f"{FIG}/Fig7_dataset_overview.png", dpi=300); plt.close()

print("Figures written to", FIG + "/")
for f in sorted(glob.glob(f"{FIG}/*.png")):
    print("  ", os.path.basename(f))
