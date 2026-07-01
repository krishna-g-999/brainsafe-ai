# -*- coding: utf-8 -*-
"""Publication figure set (v2): professional titles only (no editorial wording), adds a
workflow/architecture schematic and a model-comparison panel. Classifier/regression OOF are
recomputed via scaffold GroupKFold(5) (same protocol/seed as deployment) and cached."""
import os, json, glob, warnings
os.chdir(os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.gridspec import GridSpec
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
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans", "axes.titlesize": 12})
DEPLOYED = ["BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_B", "MAO_A", "hERG"]
LAB = [e.replace("_", "-") for e in DEPLOYED]
PAL = plt.cm.tab10(np.linspace(0, 1, 10))
rs = json.load(open("BS_randomsplit_benchmark.json"))
metas = {os.path.basename(f).replace("_meta.json", ""): json.load(open(f)) for f in glob.glob("models_brain/*_meta.json")}

def feats(smi): return np.hstack([morgan(smi), descriptors(smi)])

# ---------- recompute OOF (cached) ----------
CACHE = f"{FIG}/_oof_cache.npz"
CLF, BASE, REG = {}, {}, {}
print("Computing OOF (ensemble + baselines + regression)...")
for ep in DEPLOYED:
    d = clf_canon(pd.read_csv(f"data/endpoints/{ep}.csv"))
    smi = d["smiles"].tolist(); y = d["label"].values.astype(int)
    X = feats(smi); g = np.array([scaffold(s) for s in smi]); bv = bvs(smi)
    oof = np.zeros(len(y)); knn = np.zeros(len(y)); lr = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, groups=g):
        oof[te] = np.mean([m.fit(X[tr], y[tr]).predict_proba(X[te])[:, 1] for m in clf_models().values()], axis=0)
        # baseline 1: Tanimoto kNN (k=5)
        bt = [bv[i] for i in tr]
        for i in te:
            sims = np.array(DataStructs.BulkTanimotoSimilarity(bv[i], bt))
            idx = np.argsort(sims)[::-1][:5]; w = sims[idx]
            knn[i] = (w * y[tr][idx]).sum() / w.sum() if w.sum() > 0 else y[tr].mean()
        # baseline 2: logistic regression
        sc = StandardScaler().fit(X[tr])
        lr[te] = LogisticRegression(max_iter=2000, class_weight="balanced").fit(sc.transform(X[tr]), y[tr]).predict_proba(sc.transform(X[te]))[:, 1]
    CLF[ep] = (y, oof)
    BASE[ep] = {"kNN-Tanimoto": roc_auc_score(y, knn), "Logistic regression": roc_auc_score(y, lr),
                "Ensemble (deployed)": roc_auc_score(y, oof)}
for ep, anti in [("antioxidant", True), ("A2A", False), ("HT2A", False), ("D2", False), ("SERT", False)]:
    if anti:
        df = pd.read_csv("data/endpoints_reg/antioxidant_dpph.csv").rename(columns={"y": "pchembl"})
        dd = reg_canon(df.assign(label=0))
    else:
        dd = reg_canon(pd.read_csv(f"data/endpoints/{ep}.csv"))
    smi = dd["smiles"].tolist(); y = dd["y"].values; X = feats(smi); g = np.array([scaffold(s) for s in smi])
    oof = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, groups=g):
        oof[te] = np.mean([m.fit(X[tr], y[tr]).predict(X[te]) for m in reg_ensemble()], axis=0)
    REG[ep] = (y, oof)
print("OOF done.")

# ============ FIGURE: WORKFLOW / ARCHITECTURE ============
fig, ax = plt.subplots(figsize=(13, 6)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 6)
def box(x, y, w, h, text, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec="#333", lw=1.2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9.5, wrap=True)
def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14, color="#555", lw=1.4))
row = 4.2
box(0.2, row, 2.3, 1.4, "Measured data\nChEMBL pChEMBL,\nB3DB (BBB), DPPH", "#D6EAF8")
box(2.9, row, 2.1, 1.4, "Curation\ncanonicalise,\nInChIKey dedupe,\nlabel", "#D5F5E3")
box(5.4, row, 2.1, 1.4, "Features\nECFP-1024 +\n24 RDKit\ndescriptors", "#FCF3CF")
box(7.9, row, 2.2, 1.4, "Ensemble\nRF + ExtraTrees\n+ HistGB", "#FADBD8")
box(10.5, row, 2.3, 1.4, "Calibration\nisotonic +\nMondrian conformal", "#E8DAEF")
for x in (2.5, 5.0, 7.5, 10.1): arrow(x, row + 0.7, x + 0.4, row + 0.7)
row2 = 1.6
box(2.0, row2, 3.4, 1.4, "Integration\nBBB-gated disease scores,\nnearest-analogue evidence,\ndruggability / CNS-MPO,\nclinical precedent", "#D6EAF8")
box(6.2, row2, 3.0, 1.4, "Outputs\nper-disease engagement,\nsafety (hERG), antioxidant,\nconformal confidence", "#D1F2EB")
box(9.9, row2, 2.9, 1.4, "Validation\nrandom / scaffold /\ncluster / temporal", "#FDEBD0")
arrow(11.6, row, 4.0, row2 + 1.4)   # calibration -> integration
arrow(5.4, row2 + 0.7, 6.2, row2 + 0.7)
arrow(9.2, row2 + 0.7, 9.9, row2 + 0.7)
ax.set_title("BrainSafe AI pipeline: data, model training and integration", fontsize=13, pad=10)
plt.tight_layout(); plt.savefig(f"{FIG}/fig_workflow.png", dpi=300, bbox_inches="tight"); plt.close()

# ============ FIGURE: DATASET OVERVIEW ============
fig, ax = plt.subplots(figsize=(9, 5))
ns = [metas[e]["n"] for e in DEPLOYED]; pos = [metas[e]["pos_rate"] * 100 for e in DEPLOYED]
ax.bar(LAB, ns, color="#4C72B0"); ax.set_ylabel("Number of compounds")
ax2 = ax.twinx(); ax2.plot(LAB, pos, "o-", color="#C44E52"); ax2.set_ylabel("Active (%)"); ax2.set_ylim(0, 100)
ax.set_title("Training set size and class balance per endpoint"); plt.tight_layout()
plt.savefig(f"{FIG}/fig_dataset.png", dpi=300); plt.close()

# ============ FIGURE: VALIDATION HIERARCHY ============
fig, ax = plt.subplots(figsize=(11, 5.5)); x = np.arange(len(DEPLOYED)); w = 0.2
for i, (s, lab) in enumerate(zip(["AUROC_random", "auroc", "cluster_split_auroc", "temporal_auroc"],
                                 ["Random", "Scaffold", "Cluster", "Temporal"])):
    vals = [rs.get(e) if s == "AUROC_random" else metas[e].get(s) for e in DEPLOYED]
    vals = [v if v is not None else np.nan for v in vals]
    ax.bar(x + (i - 1.5) * w, vals, w, label=lab, color=PAL[i])
ax.set_xticks(x); ax.set_xticklabels(LAB); ax.set_ylim(0.5, 1.0); ax.set_ylabel("AUROC")
ax.set_title("AUROC across validation regimes"); ax.legend(ncol=4, fontsize=9, loc="lower center", bbox_to_anchor=(0.5, -0.2))
plt.tight_layout(); plt.savefig(f"{FIG}/fig_validation_hierarchy.png", dpi=300, bbox_inches="tight"); plt.close()

# ============ FIGURE: ROC ============
fig, ax = plt.subplots(figsize=(7, 7))
for i, ep in enumerate(DEPLOYED):
    y, p = CLF[ep]; fpr, tpr, _ = roc_curve(y, p)
    ax.plot(fpr, tpr, lw=2, color=PAL[i], label=f"{LAB[i]} ({roc_auc_score(y,p):.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=1); ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
ax.set_title("Scaffold cross-validation ROC curves"); ax.legend(fontsize=9, loc="lower right", title="AUROC")
plt.tight_layout(); plt.savefig(f"{FIG}/fig_roc.png", dpi=300); plt.close()

# ============ FIGURE: CALIBRATION ============
fig, ax = plt.subplots(figsize=(7, 7))
for i, ep in enumerate(["BBB", "AChE", "BACE1", "hERG"]):
    y, p = CLF[ep]; pc = IsotonicRegression(out_of_bounds="clip").fit(p, y).predict(p)
    bins = np.linspace(0, 1, 11); idx = np.digitize(pc, bins) - 1; xs, ys = [], []
    for b in range(10):
        m = idx == b
        if m.sum() > 5: xs.append(pc[m].mean()); ys.append(y[m].mean())
    ax.plot(xs, ys, "o-", color=PAL[i], label=ep.replace("_", "-"))
ax.plot([0, 1], [0, 1], "k--", lw=1, label="ideal"); ax.set_xlabel("Predicted probability")
ax.set_ylabel("Observed frequency"); ax.set_title("Probability calibration (isotonic)"); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(f"{FIG}/fig_calibration.png", dpi=300); plt.close()

# ============ FIGURE: CONFORMAL ============
fig, ax = plt.subplots(figsize=(9, 4.8))
ax.bar(LAB, [metas[e].get("conformal_coverage") for e in DEPLOYED], color="#8172B3")
ax.axhline(0.90, ls="--", c="#C44E52", label="0.90 target"); ax.set_ylim(0.8, 1.0)
ax.set_ylabel("Empirical coverage"); ax.set_title("Conformal prediction coverage"); ax.legend()
plt.tight_layout(); plt.savefig(f"{FIG}/fig_conformal.png", dpi=300); plt.close()

# ============ FIGURE: MODEL COMPARISON (vs baselines) ============
fig, ax = plt.subplots(figsize=(11, 5.5)); x = np.arange(len(DEPLOYED)); w = 0.26
for i, name in enumerate(["kNN-Tanimoto", "Logistic regression", "Ensemble (deployed)"]):
    ax.bar(x + (i - 1) * w, [BASE[e][name] for e in DEPLOYED], w, label=name, color=PAL[i])
ax.set_xticks(x); ax.set_xticklabels(LAB); ax.set_ylim(0.5, 1.0); ax.set_ylabel("AUROC (scaffold CV)")
ax.set_title("Model comparison: ensemble versus baselines"); ax.legend(fontsize=9, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.2))
plt.tight_layout(); plt.savefig(f"{FIG}/fig_model_comparison.png", dpi=300, bbox_inches="tight"); plt.close()

# ============ FIGURE: BENCHMARK vs LITERATURE ============
fig, ax = plt.subplots(figsize=(9, 5))
lit = {"BBB": (0.88, 0.96), "hERG": (0.86, 0.93), "AChE": (0.90, 0.97), "BACE1": (0.90, 0.96),
       "MAO_B": (0.88, 0.96), "MAO_A": (0.85, 0.95), "BChE": (0.90, 0.97), "GSK3B": (0.88, 0.95)}
for i, ep in enumerate(DEPLOYED):
    lo, hi = lit[ep]; ax.plot([i, i], [lo, hi], lw=8, color="#AEC7E8", solid_capstyle="round")
    ax.plot(i, rs[ep], "D", color="#C44E52", ms=9, zorder=5)
ax.plot([], [], lw=8, color="#AEC7E8", label="published random-split range")
ax.plot([], [], "D", color="#C44E52", label="this work (random split)")
ax.set_xticks(range(len(DEPLOYED))); ax.set_xticklabels(LAB); ax.set_ylim(0.82, 1.0); ax.set_ylabel("AUROC")
ax.set_title("Per-endpoint AUROC versus published random-split ranges"); ax.legend()
plt.tight_layout(); plt.savefig(f"{FIG}/fig_benchmark.png", dpi=300); plt.close()

# ============ FIGURE: REGRESSION SCATTER ============
fig = plt.figure(figsize=(15, 6)); gs = GridSpec(2, 3, figure=fig)
order = ["antioxidant", "A2A", "HT2A", "D2", "SERT"]
tt = {"antioxidant": "Antioxidant (DPPH pIC50)", "A2A": "A2A (pKi)", "HT2A": "5-HT2A (pKi)", "D2": "D2 (pKi)", "SERT": "SERT (pKi)"}
for i, ep in enumerate(order):
    axx = fig.add_subplot(gs[i // 3, i % 3]); y, p = REG[ep]
    axx.scatter(y, p, s=6, alpha=0.3, color=PAL[i]); lo, hi = min(y.min(), p.min()), max(y.max(), p.max())
    axx.plot([lo, hi], [lo, hi], "k--", lw=1)
    axx.set_title(f"{tt[ep]}  R2={r2_score(y,p):.2f}"); axx.set_xlabel("Measured"); axx.set_ylabel("Predicted")
fig.suptitle("Regression endpoints: predicted versus measured potency (scaffold cross-validation)", y=1.02)
plt.tight_layout(); plt.savefig(f"{FIG}/fig_regression_scatter.png", dpi=300, bbox_inches="tight"); plt.close()

# remove old capitalised figures to avoid duplicates
for old in glob.glob(f"{FIG}/Fig[0-9]_*.png"): os.remove(old)
json.dump({e: BASE[e] for e in DEPLOYED}, open("BS_baseline_comparison.json", "w"), indent=2)
print("Figures (v2) written:")
for f in sorted(glob.glob(f"{FIG}/fig_*.png")): print("  ", os.path.basename(f))
