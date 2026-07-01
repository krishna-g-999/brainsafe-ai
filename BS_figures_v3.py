# -*- coding: utf-8 -*-
"""Submission-standard figure set: unified style, Okabe–Ito colorblind palette, consistent
fonts/sizes, panel labels (A/B), bootstrap 95% CIs on AUROC, and a graphical abstract.
Classifier/regression OOF recomputed via scaffold GroupKFold(5) and cached."""
import os, json, glob, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_curve, roc_auc_score, r2_score
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from rdkit import DataStructs
from BS_predictive_model import morgan, descriptors, scaffold, bvs
from BS_train_endpoints import models as clf_models, canon as clf_canon
from BS_train_regression import reg_ensemble, canon as reg_canon

FIG = "figures"; os.makedirs(FIG, exist_ok=True)
# Okabe–Ito colorblind-safe palette
OI = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]
NAVY = "#0D2137"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11, "axes.titlesize": 12, "axes.titleweight": "bold",
    "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 9,
    "figure.dpi": 300, "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False,
})
DEPLOYED = ["BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_B", "MAO_A", "hERG"]
LAB = [e.replace("_", "-") for e in DEPLOYED]
rs = json.load(open("BS_randomsplit_benchmark.json"))
metas = {os.path.basename(f).replace("_meta.json", ""): json.load(open(f)) for f in glob.glob("models_brain/*_meta.json")}

def feats(smi): return np.hstack([morgan(smi), descriptors(smi)])
def panel(ax, letter): ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")
def boot_auc(y, p, B=1000):
    rng = np.random.default_rng(42); n = len(y); a = []
    for _ in range(B):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) == 2: a.append(roc_auc_score(y[idx], p[idx]))
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))

print("Recomputing OOF + bootstrap CIs (resumable cache; baselines reused)...", flush=True)
CACHE = f"{FIG}/_oof_cache.npz"
cache = {k: np.load(CACHE, allow_pickle=True)[k] for k in np.load(CACHE, allow_pickle=True).files} if os.path.exists(CACHE) else {}
CLF, CIs, REG = {}, {}, {}
BASE = json.load(open("BS_baseline_comparison.json"))   # reuse baselines from prior run
for ep in DEPLOYED:
    ky = f"clf_{ep}"
    if ky + "_p" in cache:
        y, oof = cache[ky + "_y"], cache[ky + "_p"]
    else:
        d = clf_canon(pd.read_csv(f"data/endpoints/{ep}.csv"))
        smi = d["smiles"].tolist(); y = d["label"].values.astype(int); X = feats(smi)
        g = np.array([scaffold(s) for s in smi]); oof = np.zeros(len(y))
        for tr, te in GroupKFold(5).split(X, groups=g):
            oof[te] = np.mean([m.fit(X[tr], y[tr]).predict_proba(X[te])[:, 1] for m in clf_models().values()], axis=0)
        cache[ky + "_y"] = y; cache[ky + "_p"] = oof; np.savez(CACHE, **cache); print("  cached", ep, flush=True)
    CLF[ep] = (y, oof); CIs[ep] = boot_auc(y, oof)
for ep, anti in [("antioxidant", True), ("A2A", False), ("HT2A", False), ("D2", False), ("SERT", False)]:
    ky = f"reg_{ep}"
    if ky + "_p" in cache:
        y, oof = cache[ky + "_y"], cache[ky + "_p"]
    else:
        if anti:
            df = pd.read_csv("data/endpoints_reg/antioxidant_dpph.csv").rename(columns={"y": "pchembl"}); dd = reg_canon(df.assign(label=0))
        else:
            dd = reg_canon(pd.read_csv(f"data/endpoints/{ep}.csv"))
        smi = dd["smiles"].tolist(); y = dd["y"].values; X = feats(smi); g = np.array([scaffold(s) for s in smi]); oof = np.zeros(len(y))
        for tr, te in GroupKFold(5).split(X, groups=g):
            oof[te] = np.mean([m.fit(X[tr], y[tr]).predict(X[te]) for m in reg_ensemble()], axis=0)
        cache[ky + "_y"] = y; cache[ky + "_p"] = oof; np.savez(CACHE, **cache); print("  cached reg", ep, flush=True)
    REG[ep] = (y, oof)
json.dump({e: {"auroc": round(roc_auc_score(*CLF[e]), 3), "ci95": [round(CIs[e][0], 3), round(CIs[e][1], 3)]} for e in DEPLOYED},
          open("BS_auroc_cis.json", "w"), indent=2)
print("OOF + CIs done.", flush=True)

# ===== GRAPHICAL ABSTRACT =====
fig, ax = plt.subplots(figsize=(12, 5)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 5)
def gbox(x, y, w, h, t, fc, tc="white", fs=9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06,rounding_size=0.18", fc=fc, ec="none"))
    ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", color=tc, fontsize=fs, fontweight="bold")
ax.text(6, 4.7, "BrainSafe AI: structure-to-brain-effect prediction from measured data", ha="center", fontsize=13, fontweight="bold", color=NAVY)
gbox(0.3, 2.6, 2.0, 1.3, "Any compound\n(SMILES / name)\ne.g. flavonoid", OI[0])
gbox(2.7, 2.6, 2.2, 1.3, "Measured data\nChEMBL_37 + B3DB\n64,474 records", OI[2])
gbox(5.3, 2.6, 2.3, 1.3, "Calibrated ensemble\n+ conformal\n+ evidence", OI[1], tc="black")
for x in (2.3, 4.9): ax.add_patch(FancyArrowPatch((x, 3.25), (x + 0.4, 3.25), arrowstyle="-|>", mutation_scale=15, color="#555"))
ax.add_patch(FancyArrowPatch((7.6, 3.25), (8.0, 3.25), arrowstyle="-|>", mutation_scale=15, color="#555"))
# example output mini-bars
ax.text(10.0, 4.05, "Brain-effect profile", ha="center", fontsize=10, fontweight="bold", color=NAVY)
ex = [("BBB", 0.92, OI[0]), ("Alzheimer", 0.84, OI[2]), ("Parkinson", 0.30, OI[5]), ("hERG risk", 0.20, OI[3])]
for i, (lab, v, c) in enumerate(ex):
    yy = 3.5 - i * 0.42
    ax.add_patch(FancyBboxPatch((8.3, yy), 2.6 * v, 0.3, boxstyle="round,pad=0.01", fc=c, ec="none"))
    ax.text(8.25, yy + 0.15, lab, ha="right", va="center", fontsize=8)
    ax.text(8.3 + 2.6 * v + 0.05, yy + 0.15, f"{v:.2f}", va="center", fontsize=7.5)
gbox(1.4, 0.62, 9.2, 1.0, "Validation:  random AUROC 0.94–0.98   ·   scaffold / cluster 0.87–0.95   ·   "
     "temporal 0.61–0.92   ·   conformal coverage ≈ 0.90", OI[4], fs=8.8)
ax.text(6, 0.22, "Research use · pending peer review · predicts target engagement, not clinical efficacy", ha="center", fontsize=8, style="italic", color="#777")
plt.tight_layout(); plt.savefig(f"{FIG}/graphical_abstract.png", bbox_inches="tight"); plt.close()

# ===== FIG1 WORKFLOW =====
fig, ax = plt.subplots(figsize=(13, 5.5)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 6)
def wbox(x, y, w, h, t, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.2", fc=fc, ec="#333", lw=1))
    ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=9)
def war(x1, y1, x2, y2): ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=13, color="#555", lw=1.3))
r1 = 4.3
wbox(0.2, r1, 2.3, 1.4, "Measured data\nChEMBL_37 pChEMBL,\nB3DB, DPPH", "#CDE6F5")
wbox(2.9, r1, 2.0, 1.4, "Curation\ncanonicalise,\ndedupe, label", "#CDEFE0")
wbox(5.3, r1, 2.0, 1.4, "Features\nECFP-1024 +\n24 descriptors", "#FBE7C2")
wbox(7.7, r1, 2.1, 1.4, "Ensemble\nRF+ExtraTrees\n+HistGB", "#F5D6C6")
wbox(10.2, r1, 2.5, 1.4, "Calibration\nisotonic +\nconformal", "#E5D4EC")
for x in (2.5, 4.9, 7.3, 9.8): war(x, r1 + 0.7, x + 0.4, r1 + 0.7)
r2 = 1.7
wbox(1.8, r2, 3.4, 1.4, "Integration\nBBB-gated disease scores,\nnearest-analogue evidence,\ndruggability, clinical precedent", "#CDE6F5")
wbox(6.0, r2, 3.0, 1.4, "Outputs\ndisease engagement,\nsafety, antioxidant,\nconformal confidence", "#CDEFE0")
wbox(9.7, r2, 3.0, 1.4, "Validation\nrandom / scaffold /\ncluster / temporal", "#FBE7C2")
war(11.4, r1, 3.5, r2 + 1.4); war(5.2, r2 + 0.7, 6.0, r2 + 0.7); war(9.0, r2 + 0.7, 9.7, r2 + 0.7)
ax.set_title("BrainSafe AI pipeline: data, training, integration and validation", fontsize=12, pad=8)
plt.tight_layout(); plt.savefig(f"{FIG}/fig1_workflow.png", bbox_inches="tight"); plt.close()

# ===== FIG2 DATASET =====
fig, ax = plt.subplots(figsize=(9, 5))
ax.bar(LAB, [metas[e]["n"] for e in DEPLOYED], color=OI[0]); ax.set_ylabel("Number of compounds")
ax2 = ax.twinx(); ax2.plot(LAB, [metas[e]["pos_rate"] * 100 for e in DEPLOYED], "o-", color=OI[3]); ax2.set_ylabel("Active (%)"); ax2.set_ylim(0, 100); ax2.spines["top"].set_visible(False)
ax.set_title("Training-set size and class balance per endpoint"); plt.tight_layout(); plt.savefig(f"{FIG}/fig2_dataset.png"); plt.close()

# ===== FIG3 VALIDATION (with CI on scaffold) =====
fig, ax = plt.subplots(figsize=(11, 5.5)); x = np.arange(len(DEPLOYED)); w = 0.2
regs = [("AUROC_random", "Random", OI[0]), ("auroc", "Scaffold", OI[1]), ("cluster_split_auroc", "Cluster", OI[2]), ("temporal_auroc", "Temporal", OI[3])]
for i, (s, lab, c) in enumerate(regs):
    vals = [rs.get(e) if s == "AUROC_random" else metas[e].get(s) for e in DEPLOYED]; vals = [v if v is not None else np.nan for v in vals]
    if s == "auroc":
        err = np.array([[metas[e]["auroc"] - CIs[e][0] for e in DEPLOYED], [CIs[e][1] - metas[e]["auroc"] for e in DEPLOYED]])
        ax.bar(x + (i - 1.5) * w, vals, w, label=lab, color=c, yerr=err, capsize=2, error_kw={"lw": 0.8})
    else:
        ax.bar(x + (i - 1.5) * w, vals, w, label=lab, color=c)
ax.set_xticks(x); ax.set_xticklabels(LAB); ax.set_ylim(0.5, 1.0); ax.set_ylabel("AUROC")
ax.set_title("AUROC across validation regimes (scaffold bars: 95% bootstrap CI)")
ax.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.2)); plt.tight_layout(); plt.savefig(f"{FIG}/fig3_validation.png", bbox_inches="tight"); plt.close()

# ===== FIG4 ROC + CALIBRATION =====
fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for i, ep in enumerate(DEPLOYED):
    y, p = CLF[ep]; fpr, tpr, _ = roc_curve(y, p)
    axes[0].plot(fpr, tpr, lw=1.8, color=OI[i % 8], label=f"{LAB[i]} ({roc_auc_score(y,p):.3f})")
axes[0].plot([0, 1], [0, 1], "k--", lw=1); axes[0].set_xlabel("False positive rate"); axes[0].set_ylabel("True positive rate")
axes[0].set_title("Scaffold cross-validation ROC"); axes[0].legend(fontsize=8, loc="lower right", title="AUROC"); panel(axes[0], "A")
for i, ep in enumerate(["BBB", "AChE", "BACE1", "hERG"]):
    y, p = CLF[ep]; pc = IsotonicRegression(out_of_bounds="clip").fit(p, y).predict(p)
    b = np.linspace(0, 1, 11); idx = np.digitize(pc, b) - 1; xs, ys = [], []
    for k in range(10):
        m = idx == k
        if m.sum() > 5: xs.append(pc[m].mean()); ys.append(y[m].mean())
    axes[1].plot(xs, ys, "o-", color=OI[i], label=ep.replace("_", "-"))
axes[1].plot([0, 1], [0, 1], "k--", lw=1); axes[1].set_xlabel("Predicted probability"); axes[1].set_ylabel("Observed frequency")
axes[1].set_title("Probability calibration (isotonic)"); axes[1].legend(fontsize=9); panel(axes[1], "B")
plt.tight_layout(); plt.savefig(f"{FIG}/fig4_roc_calibration.png", bbox_inches="tight"); plt.close()

# ===== FIG5 CONFORMAL + MODEL COMPARISON =====
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
axes[0].bar(LAB, [metas[e].get("conformal_coverage") for e in DEPLOYED], color=OI[4])
axes[0].axhline(0.90, ls="--", c=OI[3], label="0.90 target"); axes[0].set_ylim(0.8, 1.0); axes[0].set_ylabel("Empirical coverage")
axes[0].set_title("Conformal prediction coverage"); axes[0].legend(); axes[0].tick_params(axis="x", rotation=30); panel(axes[0], "A")
xx = np.arange(len(DEPLOYED)); ww = 0.26
_enskey = "Ensemble (deployed)" if "Ensemble (deployed)" in BASE["BBB"] else "Ensemble"
for i, (key, lab) in enumerate([("kNN-Tanimoto", "kNN-Tanimoto"), ("Logistic regression", "Logistic regression"), (_enskey, "Ensemble (this work)")]):
    axes[1].bar(xx + (i - 1) * ww, [BASE[e][key] for e in DEPLOYED], ww, label=lab, color=OI[i])
axes[1].set_xticks(xx); axes[1].set_xticklabels(LAB, rotation=30); axes[1].set_ylim(0.5, 1.0); axes[1].set_ylabel("AUROC (scaffold CV)")
axes[1].set_title("Ensemble versus baselines"); axes[1].legend(fontsize=8); panel(axes[1], "B")
plt.tight_layout(); plt.savefig(f"{FIG}/fig5_conformal_comparison.png", bbox_inches="tight"); plt.close()

# ===== FIG6 BENCHMARK =====
fig, ax = plt.subplots(figsize=(9, 5))
lit = {"BBB": (0.88, 0.96), "hERG": (0.86, 0.93), "AChE": (0.90, 0.97), "BACE1": (0.90, 0.96), "MAO_B": (0.88, 0.96), "MAO_A": (0.85, 0.95), "BChE": (0.90, 0.97), "GSK3B": (0.88, 0.95)}
for i, ep in enumerate(DEPLOYED):
    lo, hi = lit[ep]; ax.plot([i, i], [lo, hi], lw=8, color="#9ECAE1", solid_capstyle="round"); ax.plot(i, rs[ep], "D", color=OI[3], ms=9, zorder=5)
ax.plot([], [], lw=8, color="#9ECAE1", label="published random-split range"); ax.plot([], [], "D", color=OI[3], label="this work (random split)")
ax.set_xticks(range(len(DEPLOYED))); ax.set_xticklabels(LAB); ax.set_ylim(0.82, 1.0); ax.set_ylabel("AUROC")
ax.set_title("Per-endpoint AUROC versus published random-split ranges"); ax.legend(); plt.tight_layout(); plt.savefig(f"{FIG}/fig6_benchmark.png"); plt.close()

# ===== FIG7 REGRESSION =====
fig, axes = plt.subplots(2, 3, figsize=(15, 9)); axes = axes.ravel()
order = ["antioxidant", "A2A", "HT2A", "D2", "SERT"]
tt = {"antioxidant": "Antioxidant (DPPH pIC50)", "A2A": "A2A (pKi)", "HT2A": "5-HT2A (pKi)", "D2": "D2 (pKi)", "SERT": "SERT (pKi)"}
for i, ep in enumerate(order):
    y, p = REG[ep]; axes[i].scatter(y, p, s=6, alpha=0.3, color=OI[i]); lo, hi = min(y.min(), p.min()), max(y.max(), p.max())
    axes[i].plot([lo, hi], [lo, hi], "k--", lw=1); axes[i].set_title(f"{tt[ep]}  R² = {r2_score(y,p):.2f}")
    axes[i].set_xlabel("Measured"); axes[i].set_ylabel("Predicted (OOF)"); panel(axes[i], "ABCDE"[i])
axes[5].axis("off")
fig.suptitle("Regression endpoints: predicted versus measured potency (scaffold cross-validation)", fontsize=13, fontweight="bold")
plt.tight_layout(); plt.savefig(f"{FIG}/fig7_regression.png", bbox_inches="tight"); plt.close()

for old in glob.glob(f"{FIG}/fig_*.png"): os.remove(old)
print("v3 figures written:")
for f in sorted(glob.glob(f"{FIG}/*.png")): print("  ", os.path.basename(f))
