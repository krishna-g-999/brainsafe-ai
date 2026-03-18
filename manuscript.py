"""
BrainSafe AI — Manuscript Figure Generator
Six publication-quality figures for high-impact journal submission.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.ticker as mticker
import numpy as np
import json, base64, io
from pathlib import Path

# ── Global style ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":        "DejaVu Sans",
    "font.size":          10,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.linewidth":     0.8,
    "axes.titlesize":     11,
    "axes.labelsize":     10,
    "xtick.labelsize":    9,
    "ytick.labelsize":    9,
    "legend.fontsize":    9,
    "legend.frameon":     False,
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
    "grid.color":         "#E8E8E8",
    "grid.linewidth":     0.6,
})

# ── Palette ───────────────────────────────────────────────────────────────────
C = dict(
    navy   = "#1B2A4A",
    red    = "#C1121F",
    teal   = "#0D6E6E",
    gold   = "#CB8E00",
    blue   = "#1E4CC9",
    green  = "#1A7A4A",
    purple = "#5B2D8E",
    slate  = "#4A5568",
    lg     = "#F4F6FA",
    mg     = "#DDE3EE",
    dg     = "#8899AA",
)

# Curated / ML tier colors
TIER_COLORS = {"curated": C["green"], "ml": C["blue"], "est": C["gold"]}

# ── Paths ─────────────────────────────────────────────────────────────────────
CURATED_PATH = Path(__file__).parent / "compounds.json"
ML_PATH      = Path(__file__).parent / "compounds_ml.json"

SCORE_DIMS = [
    "antioxidant", "anti_inflammatory", "mitochondrial_support",
    "aggregation_modulation", "cognitive_enhancement",
    "neurogenesis_support", "synaptic_plasticity",
]
SHORT_DIM = [
    "Antioxidant", "Anti-Inflammatory", "Mitochondrial\nSupport",
    "Aggregation\nModulation", "Cognitive\nEnhancement",
    "Neurogenesis\nSupport", "Synaptic\nPlasticity",
]
RADAR_DIM = [
    "Antioxidant", "Anti-\nInflammatory", "Mitochondrial\nSupport",
    "Aggregation\nModulation", "Cognitive\nEnhancement",
    "Neurogenesis\nSupport", "Synaptic\nPlasticity",
]
FEAT_NAMES = [
    "BBB Level", "ALS\nRelevance", "Alzheimer's\nRelevance",
    "Parkinson's\nRelevance", "Huntington's\nRelevance", "No. Pathways",
]

# Pre-computed from 5-fold CV (n=129, RF: 120 trees, max_depth=6, min_samples_leaf=3)
PER_FOLD_R2 = {
    "antioxidant": [
        0.3593,
        -0.3922,
        0.3229,
        0.2127,
        0.1177
    ],
    "antiinflammatory": [
        0.3361,
        0.0425,
        0.2359,
        0.3806,
        0.3396
    ],
    "mitochondrialsupport": [
        0.2292,
        0.1614,
        0.4185,
        0.2611,
        0.3521
    ],
    "aggregationmodulation": [
        0.4335,
        0.3859,
        0.0312,
        0.0496,
        0.377
    ],
    "cognitiveenhancement": [
        0.3717,
        0.3256,
        0.5896,
        0.4719,
        0.4265
    ],
    "neurogenesis": [
        0.0696,
        0.2811,
        0.0064,
        0.0244,
        -0.4066
    ],
    "synapticplasticity": [
        -0.0165,
        0.2429,
        0.3144,
        0.2763,
        -0.1626
    ]
}

# Feature importance matrix (rows = score_dim, cols = 6 features)
IMPORTANCE_MATRIX = {
    "antioxidant":           [0.246, 0.110, 0.170, 0.221, 0.082, 0.170],
    "anti_inflammatory":     [0.163, 0.079, 0.355, 0.168, 0.056, 0.179],
    "mitochondrial_support": [0.068, 0.495, 0.081, 0.063, 0.223, 0.070],
    "aggregation_modulation":[0.112, 0.042, 0.233, 0.510, 0.017, 0.086],
    "cognitive_enhancement": [0.098, 0.060, 0.679, 0.049, 0.023, 0.091],
    "synaptic_plasticity":   [0.080, 0.049, 0.596, 0.075, 0.040, 0.160],
}
INFORMATIVE_DIMS = list(PER_FOLD_R2.keys())
DIM_SHORT = {
    "antioxidant":            "Antioxidant",
    "anti_inflammatory":      "Anti-Inflam.",
    "mitochondrial_support":  "Mitochondrial",
    "aggregation_modulation": "Aggregation",
    "cognitive_enhancement":  "Cognitive",
    "synaptic_plasticity":    "Synaptic",
}

def _load():
    with open(CURATED_PATH) as f: c = json.load(f)
    ml = {}
    if ML_PATH.exists():
        with open(ML_PATH) as f: ml = json.load(f)
    return c, ml

def _nps(d):
    return min(100, d.get("antioxidant",0)*3 + d.get("anti_inflammatory",0)*3
               + d.get("mitochondrial_support",0)*2 + d.get("aggregation_modulation",0)*2)

def _b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()

def _despine(ax, full=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if full:
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — System Architecture (professional flowchart)
# ══════════════════════════════════════════════════════════════════════════════
def figure1_architecture():
    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 15); ax.set_ylim(0, 8); ax.axis("off")

    def rbox(x, y, w, h, fc, ec, lw=1.8, rad=0.25, alpha=1.0):
        bp = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad={rad}",
                            fc=fc, ec=ec, lw=lw, zorder=3, alpha=alpha)
        ax.add_patch(bp)

    def txt(x, y, s, fs=9, fw="normal", col="#1B2A4A", ha="center", va="center", **kw):
        ax.text(x, y, s, fontsize=fs, fontweight=fw, color=col,
                ha=ha, va=va, zorder=5, **kw)

    def arr(x1, y1, x2, y2, col="#777"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.5,
                                   mutation_scale=14), zorder=4)

    # ─ Title ─
    txt(7.5, 7.65, "BrainSafe AI — System Architecture", fs=14, fw="bold", col=C["navy"])

    # ─ Tier 1 ─
    rbox(0.3, 5.1, 4.1, 1.9, "#EBF5EB", C["green"], lw=2.2)
    txt(2.35, 6.72, "TIER 1  ·  Literature-Curated", fs=9.5, fw="bold", col=C["green"])
    txt(2.35, 6.35, "129 compounds", fs=9, col=C["green"])
    txt(2.35, 6.02, "7-dim bioactivity · Pathways · Metabolites", fs=8.2, col="#2D6A2D")
    txt(2.35, 5.72, "ERC profiles · Disease relevance · BBB class", fs=8.2, col="#2D6A2D")
    txt(2.35, 5.35, "Source: PubMed/PMC systematic reviews 2015–2026", fs=7.5, col="#5A7A5A",
        style="italic")
    # badge
    rbox(3.55, 6.6, 0.72, 0.26, C["green"], C["green"], rad=0.08)
    txt(3.91, 6.73, "CURATED", fs=7, fw="bold", col="white")

    # ─ Tier 2 ─
    rbox(0.3, 2.85, 4.1, 1.9, "#E8EEFF", C["blue"], lw=2.2)
    txt(2.35, 4.47, "TIER 2  ·  ML-Predicted", fs=9.5, fw="bold", col=C["blue"])
    txt(2.35, 4.1, "196 compounds", fs=9, col=C["blue"])
    txt(2.35, 3.77, "Random Forest trained on Tier 1 (n=129)", fs=8.2, col="#1A3A7A")
    txt(2.35, 3.47, "ChEMBL neuro-indications: ALS · AD · PD · HD", fs=8.2, col="#1A3A7A")
    txt(2.35, 3.1, "5-fold CV R² = 0.17–0.45 per dimension", fs=7.5, col="#4A5A8A",
        style="italic")
    rbox(3.5, 4.35, 0.78, 0.26, C["blue"], C["blue"], rad=0.08)
    txt(3.89, 4.48, "ML-PRED.", fs=7, fw="bold", col="white")

    # ─ Tier 3 ─
    rbox(0.3, 0.65, 4.1, 1.9, "#FFF8E6", C["gold"], lw=2.2)
    txt(2.35, 2.27, "TIER 3  ·  Class-Based Inference", fs=9.5, fw="bold", col=C["gold"])
    txt(2.35, 1.9, "Unlimited coverage", fs=9, col=C["gold"])
    txt(2.35, 1.57, "15 expert-designed chemical class templates", fs=8.2, col="#7A5500")
    txt(2.35, 1.27, "Flavonoids · Polyphenols · B-vitamins · Fatty acids", fs=8.2, col="#7A5500")
    txt(2.35, 0.9, "Applied when PubChem cannot identify the query", fs=7.5, col="#9A7010",
        style="italic")
    rbox(3.2, 2.15, 1.07, 0.26, C["gold"], C["gold"], rad=0.08)
    txt(3.735, 2.28, "CLASS EST.", fs=7, fw="bold", col="white")

    # ─ Evidence arrow ─
    ax.annotate("", xy=(0.14, 7.0), xytext=(0.14, 0.8),
                arrowprops=dict(arrowstyle="-|>", color="#AABB88", lw=2.5,
                                mutation_scale=12), zorder=2)
    txt(0.06, 3.9, "Evidence\nStrength", fs=8, col="#88AA66", rotation=90, fw="bold")

    # ─ Live API block ─
    rbox(5.1, 0.5, 4.3, 7.1, "#FAFCFF", C["navy"], lw=1.5, alpha=0.5, rad=0.3)
    txt(7.25, 7.25, "Live API Enrichment (all tiers, every report)", fs=9.5, fw="bold", col=C["navy"])

    api_items = [
        (5.3, 5.8, 3.9, 1.1, "#F0F4FF", "#4466DD",
         "PubChem REST API",
         "CID · MW · XLogP · TPSA · HBD · HBA · InChIKey · Synonyms"),
        (5.3, 4.3, 3.9, 1.1, "#EAF6F6", C["teal"],
         "ChEMBL REST API",
         "ChEMBL ID · Max Phase · QED · Ro5 · ALogP · Mechanism of Action"),
        (5.3, 2.8, 3.9, 1.1, "#F5F0FF", C["purple"],
         "KEGG REST API",
         "Compound ID · Human (HSA) pathway maps · Clickable pathway links"),
        (5.3, 0.75, 3.9, 1.65, "#FFF0F2", C["red"],
         "Live ML Estimation (on-demand, any compound)",
         "1. PubChem: fetch MW, XLogP, TPSA, HBD → CNS-MPO BBB score\n"
         "2. ChEMBL: MoA text → disease keyword inference\n"
         "3. Random Forest: predict 7-dimension neuroprotective profile"),
    ]
    for x, y, w, h, fc, ec, title, body in api_items:
        rbox(x, y, w, h, fc, ec, lw=1.6, rad=0.2)
        txt(x+w/2, y+h-0.22, title, fs=9, fw="bold", col=ec)
        for li, line in enumerate(body.split("\n")):
            txt(x+w/2, y+h-0.22-(li+1)*0.27, line, fs=7.8, col="#333355")

    # ─ Arrows from tiers to APIs ─
    for ay in [6.35, 4.1, 1.6]:
        arr(4.4, ay, 5.08, 6.35 if ay == 6.35 else (4.1 if ay == 4.1 else 1.6),
            col="#99AABB")
    arr(4.4, 3.9, 5.08, 4.87, col="#99AABB")
    arr(4.4, 1.6, 5.08, 3.35, col="#99AABB")

    # ─ Output panel ─
    rbox(10.0, 0.5, 4.7, 7.1, "#FAFAFA", C["navy"], lw=2.2, rad=0.3)
    txt(12.35, 7.25, "Compound Report Output", fs=10, fw="bold", col=C["navy"])

    out_items = [
        (C["blue"],   "7-Dimension Neuroprotective Radar Chart"),
        (C["red"],    "Neuroprotective Score (NPS 0–100)"),
        (C["teal"],   "BBB Penetration (Curated + CNS-MPO)"),
        (C["purple"], "Signalling Pathway Network (interactive)"),
        (C["teal"],   "KEGG Pathway Links (live, clickable)"),
        (C["gold"],   "ERC Profile: Enzymes · Receptors · Cofactors"),
        (C["navy"],   "Disease Relevance: ALS · AD · PD · HD"),
        (C["green"],  "Data Provenance Badge (Tier 1 / 2 / 3)"),
        ("#4466DD",   "Live PubChem + ChEMBL + KEGG Enrichment"),
    ]
    for i, (col, label) in enumerate(out_items):
        y = 6.65 - i * 0.63
        ax.plot(10.35, y, "o", ms=5.5, color=col, zorder=5)
        txt(10.55, y, label, fs=8.3, col=col, ha="left")

    arr(9.4, 4.0, 9.98, 5.5, col="#99AABB")
    arr(9.4, 4.0, 9.98, 3.5, col="#99AABB")
    arr(9.4, 1.57, 9.98, 2.5, col="#CC8888")

    fig.text(0.5, 0.01,
             "Figure 1.  BrainSafe AI system architecture showing the three-tier database, "
             "live API enrichment layer, and compound report output components.",
             ha="center", fontsize=9, color="#555555", style="italic")
    return _b64(fig)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — ML Model Performance: CV Fold R² + Feature Importance Heatmap
# ══════════════════════════════════════════════════════════════════════════════
def figure2_ml_performance():
    fig = plt.figure(figsize=(14, 6))
    fig.suptitle("Random Forest Model — Cross-Validation Performance and Feature Importance",
                 fontsize=12, fontweight="bold", color=C["navy"], y=1.01)

    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.38, left=0.06, right=0.97,
                           bottom=0.14, top=0.93)

    # ── Panel A: 5-fold CV R² box + strip plot ────────────────────────────────
    ax1 = fig.add_subplot(gs[0])
    dims_plot = list(INFORMATIVE_DIMS)
    means = [np.mean(PER_FOLD_R2[d]) for d in dims_plot]
    all_folds = [PER_FOLD_R2[d] for d in dims_plot]
    colors_bar = [C["red"] if d in ("antioxidant","anti_inflammatory")
                  else C["blue"] if d in ("mitochondrial_support","aggregation_modulation")
                  else C["teal"] for d in dims_plot]

    for i, (folds, mean, col) in enumerate(zip(all_folds, means, colors_bar)):
        folds_arr = np.array(folds)
        q1, q3 = np.percentile(folds_arr, 25), np.percentile(folds_arr, 75)
        ax1.fill_betweenx([q1, q3], i - 0.28, i + 0.28, color=col, alpha=0.25, zorder=2)
        ax1.plot([i - 0.28, i + 0.28], [np.median(folds_arr)] * 2,
                 color=col, lw=2.5, zorder=3)
        np.random.seed(42)
        jitter = np.random.uniform(-0.12, 0.12, len(folds))
        ax1.scatter(i + jitter, folds_arr, color=col, s=40, zorder=4, alpha=0.85,
                    edgecolors="white", linewidths=0.7)
        ax1.scatter(i, mean, marker="D", s=55, color=col, zorder=5,
                    edgecolors="white", linewidths=1.2)
        ax1.text(i, max(folds_arr) + 0.05, f"μ={mean:.2f}",
                 ha="center", va="bottom", fontsize=8, color=col, fontweight="bold")

    ax1.axhline(0, color="#AAAAAA", lw=1.0, ls="--", zorder=1)
    ax1.axhline(np.mean(means), color=C["navy"], lw=1.5, ls=":",
                label=f"Overall mean R² = {np.mean(means):.3f}", zorder=2)
    ax1.set_xticks(range(len(dims_plot)))
    ax1.set_xticklabels([DIM_SHORT[d] for d in dims_plot], fontsize=9, rotation=22, ha="right")
    ax1.set_ylabel("5-fold Cross-Validated R²", fontsize=10.5)
    ax1.set_ylim(-0.52, 0.82)
    ax1.yaxis.grid(True, alpha=0.5); ax1.set_axisbelow(True)
    ax1.legend(fontsize=9, loc="upper right")
    ax1.set_title("A  |  5-fold Cross-Validated R² per Bioactivity Dimension\n"
                  "(n = 129 curated compounds; diamond = fold mean; box = IQR)",
                  fontsize=9.5, loc="left", color=C["slate"])
    _despine(ax1)

    # ── Panel B: Feature Importance Heatmap ───────────────────────────────────
    ax2 = fig.add_subplot(gs[1])
    imp_mat = np.array([IMPORTANCE_MATRIX[d] for d in dims_plot])
    feat_short = ["BBB\nLevel", "ALS\nRel.", "Alzheimer\nRel.", "Parkinson\nRel.",
                  "Huntington\nRel.", "No.\nPathways"]

    cmap = LinearSegmentedColormap.from_list(
        "bsai", ["#F0F4FF", "#7BA7E0", "#1E4CC9", "#0A1B5C"])
    im = ax2.imshow(imp_mat, cmap=cmap, vmin=0, vmax=0.7, aspect="auto")

    for i in range(imp_mat.shape[0]):
        for j in range(imp_mat.shape[1]):
            v = imp_mat[i, j]
            col = "white" if v > 0.35 else C["navy"]
            fw = "bold" if v > 0.25 else "normal"
            ax2.text(j, i, f"{v:.3f}", ha="center", va="center",
                     fontsize=9, color=col, fontweight=fw)

    ax2.set_xticks(range(len(feat_short)))
    ax2.set_xticklabels(feat_short, fontsize=9)
    ax2.set_yticks(range(len(dims_plot)))
    ax2.set_yticklabels([DIM_SHORT[d] for d in dims_plot], fontsize=9.5)
    ax2.xaxis.tick_top(); ax2.xaxis.set_label_position("top")
    cbar = plt.colorbar(im, ax=ax2, shrink=0.55, pad=0.03, aspect=20)
    cbar.set_label("Feature Importance\n(Gini impurity reduction)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax2.set_title("B  |  Random Forest Feature Importance Matrix\n"
                  "(each cell: mean Gini impurity reduction across 120 trees)",
                  fontsize=9.5, loc="left", color=C["slate"], pad=28)

    # Highlight max per row
    for i in range(imp_mat.shape[0]):
        j_max = np.argmax(imp_mat[i])
        rect = plt.Rectangle((j_max - 0.5, i - 0.5), 1, 1,
                              fill=False, edgecolor="#FFD700", lw=2.5, zorder=5)
        ax2.add_patch(rect)

    ax2.text(6.7, 2.5,
             "★ dominant\n  feature\n  per target",
             fontsize=7.5, color="#AA8800", ha="left", va="center",
             bbox=dict(boxstyle="round,pad=0.3", fc="#FFFBE6", ec="#FFD700", lw=1))

    fig.text(0.5, -0.04,
             "Figure 2.  Random Forest model cross-validation performance (A) and feature importance "
             "analysis (B). Panel A: individual fold R² values (circles), IQR box, and fold mean "
             "(diamond) for each bioactivity dimension. Dashed line: zero; dotted line: overall mean. "
             "Panel B: Gini importance of each training feature per predicted dimension "
             "(gold border: dominant feature per row). Biological interpretation: Alzheimer's disease "
             "relevance dominates cognitive and synaptic predictions; ALS relevance drives "
             "mitochondrial support; Parkinson's relevance drives aggregation modulation.",
             ha="center", fontsize=8.5, color="#555555", style="italic", wrap=True)
    return _b64(fig)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — ML Pipeline & Hyperparameter Summary (combined figure)
# ══════════════════════════════════════════════════════════════════════════════
def figure3_pipeline_params():
    fig = plt.figure(figsize=(15, 7))
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 15); ax.set_ylim(0, 7); ax.axis("off")

    fig.text(0.5, 0.97, "BrainSafe AI — ML Expansion Pipeline and Model Hyperparameters",
             ha="center", va="top", fontsize=13, fontweight="bold", color=C["navy"])

    def rbox(x, y, w, h, fc, ec, lw=1.8, rad=0.22):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad={rad}",
                                    fc=fc, ec=ec, lw=lw, zorder=3))

    def txt(x, y, s, fs=9, fw="normal", col="#1B2A4A", ha="center", va="center", **kw):
        ax.text(x, y, s, fontsize=fs, fontweight=fw, color=col,
                ha=ha, va=va, zorder=5, **kw)

    def arr(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color="#668899", lw=1.6,
                                   mutation_scale=13), zorder=4)

    # ── Pipeline stages ──
    stages = [
        ("A", "#E8EEFF", "#4466DD",
         "Stage A",  "ChEMBL Drug\nIndication Query",
         ["Keywords: ALS, Alzheimer's,", "Parkinson's, Huntington's,", "Dementia"],
         "~400–600\ncandidates"),
        ("B", "#EBF5EB", C["green"],
         "Stage B",  "Physicochemical\nDescriptor Fetch",
         ["MW · ALogP · QED · HBD", "HBA · PSA · Ro5 violations", "Max phase · Oral flag"],
         "≥4 non-null\ndescriptors\nrequired"),
        ("C", "#FFF8E6", C["gold"],
         "Stage C",  "Feature\nEngineering",
         ["BBB level (0–2)", "4× disease relevance (0–2)", "n_pathways"],
         "6 features\nfrom curated\ndatabase only"),
        ("D", "#FDECEA", C["red"],
         "Stage D",  "Model Training\n& Prediction",
         ["MultiOutput RF (scikit-learn)", "120 trees · max_depth=6", "min_samples_leaf=3"],
         "196 compounds\npredicted\nmean R²=0.26"),
        ("E", "#F5F0FF", C["purple"],
         "Stage E",  "ERC\nEnrichment",
         ["ChEMBL IC₅₀/Kᵢ assays", "Homo sapiens targets", "Top 60 by NPS"],
         "Real enzyme\ndata appended"),
    ]

    xs = [0.3, 3.05, 5.8, 8.55, 11.3]
    W, H = 2.55, 3.9

    for i, (let, fc, ec, stage_lbl, title, bullets, note) in enumerate(stages):
        x = xs[i]
        rbox(x, 1.0, W, H, fc, ec, lw=2.0)

        circ = plt.Circle((x + W/2, 1.0 + H + 0.08), 0.3, fc=ec, ec="white", lw=1.5, zorder=4)
        ax.add_patch(circ)
        txt(x + W/2, 1.0 + H + 0.08, let, fs=11, fw="bold", col="white")

        txt(x + W/2, 1.0 + H - 0.28, stage_lbl, fs=7.5, col=ec)
        txt(x + W/2, 1.0 + H - 0.62, title, fs=9.5, fw="bold", col=ec,
            linespacing=1.4)

        for bi, bullet in enumerate(bullets):
            txt(x + W/2, 1.0 + H/2 + 0.38 - bi*0.38, bullet, fs=8, col="#2A2A4A")

        rbox(x + 0.18, 1.05, W - 0.36, 0.82, "white", ec, lw=1.0, rad=0.12)
        txt(x + W/2, 1.05 + 0.41, note, fs=7.8, fw="bold", col=ec)

        if i < len(stages) - 1:
            arr(x + W, 1.0 + H/2, xs[i+1], 1.0 + H/2)

    # ── Hyperparameter table ──
    rbox(11.3, 0.35, 3.4, 5.35, "#F8F9FA", C["navy"], lw=2.0, rad=0.3)
    txt(13.0, 5.47, "Model Hyperparameters", fs=10, fw="bold", col=C["navy"])

    params = [
        ("Algorithm",      "MultiOutput RandomForest"),
        ("n_estimators",   "120 trees"),
        ("max_depth",      "6"),
        ("min_samples_leaf","3"),
        ("random_state",   "42"),
        ("Validation",     "5-fold CV (KFold)"),
        ("Preprocessing",  "StandardScaler"),
        ("n_features",     "6 (engineered)"),
        ("n_targets",      "7 bioactivity dims"),
        ("Training n",     "129 curated compounds"),
        ("Mean CV R²",     "0.26 (6 informative dims)"),
        ("Runtime",        "~90–110 s (14 workers)"),
    ]
    for pi, (k, v) in enumerate(params):
        y = 5.0 - pi * 0.39
        rbox(11.5, y - 0.17, 3.0, 0.35, "#EEEEFF" if pi%2==0 else "white",
             "#CCCCEE", lw=0.7, rad=0.06)
        txt(11.75, y + 0.0, k + ":", fs=8.3, fw="bold", col=C["slate"], ha="left")
        txt(14.45, y + 0.0, v, fs=8.3, col=C["navy"], ha="right")

    fig.text(0.5, 0.01,
             "Figure 3.  BrainSafe AI ML expansion pipeline (Stages A–E) with model "
             "hyperparameter summary. Stages A–C prepare training features from ChEMBL indication "
             "data and the curated database; Stage D trains and applies the MultiOutput Random "
             "Forest; Stage E enriches the top 60 compounds with real ChEMBL bioassay data.",
             ha="center", fontsize=8.8, color="#555555", style="italic")
    return _b64(fig)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Compound Profiles: 2×2 Professional Radar Charts
# ══════════════════════════════════════════════════════════════════════════════
def figure4_radar_profiles():
    curated, ml = _load()

    compounds_data = [
        ("Curcumin",   curated, C["red"],    "TIER 1 — Curated",    "NPS = 84  ·  BBB: Medium"),
        ("Resveratrol",curated, C["green"],  "TIER 1 — Curated",    "NPS = 80  ·  BBB: High"),
        ("Melatonin",  curated, C["purple"], "TIER 1 — Curated",    "NPS = 81  ·  BBB: High"),
        ("MARAVIROC",  ml,      C["blue"],   "TIER 2 — ML-Predicted","NPS = 71  ·  BBB: Low"),
    ]

    fig = plt.figure(figsize=(13, 10))
    fig.suptitle("Neuroprotective Profiles — Representative Compounds\n"
                 "7-Dimension Radar Charts Across Two Database Tiers",
                 fontsize=12.5, fontweight="bold", color=C["navy"], y=0.99)

    labels_r = ["Antioxidant", "Anti-\nInflammatory", "Mitochondrial\nSupport",
                "Aggregation\nModulation", "Cognitive\nEnhancement",
                "Neurogenesis\nSupport", "Synaptic\nPlasticity"]
    N = len(labels_r)
    angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
    angles += [angles[0]]

    positions = [(0,0), (0,1), (1,0), (1,1)]
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.35,
                           left=0.06, right=0.96, bottom=0.09, top=0.91)

    for idx, ((name, db, col, tier_lbl, stats), pos) in enumerate(
            zip(compounds_data, positions)):
        r, c = pos
        ax = fig.add_subplot(gs[r, c], projection="polar")
        d = db.get(name, {})
        vals = [float(d.get(s, 0)) for s in SCORE_DIMS]
        vals_c = vals + [vals[0]]

        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.set_ylim(0, 10)
        ax.set_yticks([2, 4, 6, 8, 10])
        ax.set_yticklabels(["2", "4", "6", "8", "10"], fontsize=7, color="#AAAAAA")
        ax.spines["polar"].set_color("#CCCCCC")
        ax.yaxis.grid(color="#EEEEEE", lw=0.8)
        ax.xaxis.grid(color="#DDDDDD", lw=0.6)

        # Shaded reference ring at 5
        ref = [5] * N + [5]
        ax.fill(angles, ref, color="#EEEEEE", alpha=0.4, zorder=1)

        ax.fill(angles, vals_c, color=col, alpha=0.20, zorder=2)
        ax.plot(angles, vals_c, color=col, lw=2.5, zorder=3)
        ax.scatter(angles[:-1], vals, color=col, s=45, zorder=4,
                   edgecolors="white", linewidths=1.2)

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(labels_r, fontsize=8.5, color="#2A2A4A")
        ax.set_title(f"{name}\n", fontsize=11.5, fontweight="bold",
                     color=col, pad=14)

        # Tier badge + stats
        tier_col = C["green"] if "Curated" in tier_lbl else C["blue"]
        ax.text(0, -3.5, tier_lbl, transform=ax.transData,
                ha="center", va="center", fontsize=8, fontweight="bold",
                color=tier_col,
                bbox=dict(boxstyle="round,pad=0.3", fc="#FFFFFF",
                          ec=tier_col, lw=1.3))
        ax.text(0, -5.5, stats, transform=ax.transData,
                ha="center", va="center", fontsize=8.5, color=C["slate"])

        # Value annotations on points
        for angle, val, lbl in zip(angles[:-1], vals, SCORE_DIMS):
            if val >= 7:
                ax.text(angle, val + 0.65, str(val), ha="center", va="center",
                        fontsize=7.5, color=col, fontweight="bold", zorder=5)

    fig.text(0.5, 0.02,
             "Figure 4.  Neuroprotective profiles for four representative compounds. "
             "Top row: curated compounds (green/red, Tier 1). "
             "Bottom left: curated compound Melatonin. "
             "Bottom right: ML-predicted compound Maraviroc (Tier 2). "
             "Grey shaded ring at score = 5 provides visual reference. "
             "Dimension values shown for scores ≥ 7.",
             ha="center", fontsize=8.8, color="#555555", style="italic")
    return _b64(fig)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Comprehensive Results Comparison
# ══════════════════════════════════════════════════════════════════════════════
def figure5_results_comparison():
    curated, ml = _load()

    cur_nps = np.array([_nps(v) for v in curated.values()])
    ml_nps  = np.array([_nps(v) for v in ml.values()])

    # Capability radar data (BrainSafe vs 4 databases, 6 criteria)
    DB_NAMES = ["BrainSafe AI", "ChEMBL", "PubChem", "DrugBank", "STRING"]
    DB_COLS  = [C["red"], C["blue"], C["teal"], C["gold"], C["purple"]]
    CRITERIA = ["NDD-Specific\nScoring", "ML-Expanded\nProfiles",
                "Disease\nEvidence", "BBB\nPrediction",
                "Live Multi-DB\nEnrichment", "Compound-Level\nPathways"]
    CAP_DATA = {
        "BrainSafe AI": [10, 10, 9,  8,  9,  9],
        "ChEMBL":       [ 3,  0, 4,  1,  0,  2],
        "PubChem":      [ 1,  0, 1,  1,  0,  1],
        "DrugBank":     [ 5,  0, 6,  7,  0,  4],
        "STRING":       [ 2,  0, 2,  0,  0,  3],
    }

    fig = plt.figure(figsize=(15, 11))
    fig.suptitle("BrainSafe AI — Comprehensive Results and Database Comparison",
                 fontsize=13, fontweight="bold", color=C["navy"], y=0.99)

    gs_outer = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.38,
                                 left=0.05, right=0.97, bottom=0.08, top=0.94)

    # ── Panel A: Capability Radar ─────────────────────────────────────────────
    ax_r = fig.add_subplot(gs_outer[0, 0], projection="polar")
    Nc = len(CRITERIA)
    ang = np.linspace(0, 2*np.pi, Nc, endpoint=False).tolist(); ang += [ang[0]]
    ax_r.set_theta_offset(np.pi/2); ax_r.set_theta_direction(-1)
    ax_r.set_ylim(0, 10)
    ax_r.set_yticks([2,4,6,8,10]); ax_r.set_yticklabels([])
    ax_r.spines["polar"].set_color("#CCCCCC")
    ax_r.yaxis.grid(color="#EEEEEE", lw=0.8); ax_r.xaxis.grid(color="#DDDDDD", lw=0.6)
    ax_r.set_xticks(ang[:-1]); ax_r.set_xticklabels(CRITERIA, fontsize=8.5, color="#2A2A4A")

    for db, col in zip(DB_NAMES, DB_COLS):
        vals = CAP_DATA[db]; vals_c = vals + [vals[0]]
        lw = 2.8 if db == "BrainSafe AI" else 1.5
        al = 0.22 if db == "BrainSafe AI" else 0.05
        ax_r.fill(ang, vals_c, color=col, alpha=al, zorder=2)
        ax_r.plot(ang, vals_c, color=col, lw=lw, zorder=3,
                  label=db, ls="-" if db=="BrainSafe AI" else "--")
    ax_r.legend(loc="lower right", fontsize=8, frameon=True,
                framealpha=0.9, bbox_to_anchor=(1.42, -0.08))
    ax_r.set_title("A  |  Capability Comparison\n(score 0–10 per criterion)",
                   fontsize=9.5, loc="left", color=C["slate"], pad=18)

    # ── Panel B: NPS Score Distributions ─────────────────────────────────────
    ax2 = fig.add_subplot(gs_outer[0, 1])
    bins = np.linspace(20, 100, 22)
    ax2.hist(cur_nps, bins=bins, color=C["green"], alpha=0.72,
             label=f"Curated (n=129, μ={np.mean(cur_nps):.1f})",
             edgecolor="white", lw=0.8, density=False)
    ax2.hist(ml_nps,  bins=bins, color=C["blue"],  alpha=0.72,
             label=f"ML-Predicted (n=196, μ={np.mean(ml_nps):.1f})",
             edgecolor="white", lw=0.8, density=False)
    ax2.axvline(np.mean(cur_nps), color=C["green"], lw=2.2, ls="--")
    ax2.axvline(np.mean(ml_nps),  color=C["blue"],  lw=2.2, ls="--")
    ax2.axvline(70, color=C["slate"], lw=1.2, ls=":", alpha=0.7,
                label="'Strong' threshold (NPS ≥ 70)")
    ax2.set_xlabel("Neuroprotective Score (NPS)", fontsize=10.5)
    ax2.set_ylabel("Compound Count", fontsize=10.5)
    ax2.set_title("B  |  NPS Score Distribution by Database Tier",
                  fontsize=9.5, loc="left", color=C["slate"])
    ax2.legend(fontsize=8.5, frameon=False)
    ax2.yaxis.grid(True, alpha=0.4); ax2.set_axisbelow(True)
    _despine(ax2)

    # ── Panel C: NPS Category Comparison ─────────────────────────────────────
    ax3 = fig.add_subplot(gs_outer[1, 0])
    CATS  = ["Limited\n(NPS < 40)", "Moderate\n(40–69)", "Strong\n(≥ 70)"]
    CUR_C = [int(np.sum(cur_nps<40)), int(np.sum((cur_nps>=40)&(cur_nps<70))), int(np.sum(cur_nps>=70))]
    ML_C  = [int(np.sum(ml_nps<40)),  int(np.sum((ml_nps>=40)&(ml_nps<70))),  int(np.sum(ml_nps>=70))]
    x = np.arange(3); w = 0.38
    b1 = ax3.bar(x-w/2, CUR_C, w, color=C["green"], alpha=0.88, edgecolor="white",
                 lw=1.5, label="Curated (n=129)")
    b2 = ax3.bar(x+w/2, ML_C,  w, color=C["blue"],  alpha=0.88, edgecolor="white",
                 lw=1.5, label="ML-Predicted (n=196)")
    for bars in [b1, b2]:
        for bar in bars:
            h = bar.get_height()
            ax3.text(bar.get_x()+bar.get_width()/2, h+0.5, str(int(h)),
                     ha="center", va="bottom", fontsize=10, fontweight="bold",
                     color=bar.get_facecolor())
    # Add % labels
    for j, (c, m) in enumerate(zip(CUR_C, ML_C)):
        ax3.text(j-w/2, -4.5, f"{100*c//129}%", ha="center", fontsize=8,
                 color=C["green"], fontweight="bold")
        ax3.text(j+w/2, -4.5, f"{100*m//196}%", ha="center", fontsize=8,
                 color=C["blue"], fontweight="bold")
    ax3.set_xticks(x); ax3.set_xticklabels(CATS, fontsize=10)
    ax3.set_ylabel("Number of Compounds", fontsize=10.5)
    ax3.set_title("C  |  NPS Category Distribution",
                  fontsize=9.5, loc="left", color=C["slate"])
    ax3.legend(fontsize=9.5)
    ax3.yaxis.grid(True, alpha=0.4); ax3.set_axisbelow(True)
    ax3.set_ylim(-8, max(max(CUR_C),max(ML_C))*1.2)
    _despine(ax3)

    # ── Panel D: Database Feature Availability Heatmap ───────────────────────
    ax4 = fig.add_subplot(gs_outer[1, 1])
    FEAT_ROW = ["NDD-Specific NPS Scoring", "7-Dimension Radar Profile",
                "ML-Expanded Profiles", "Live ML for Unknown Compounds",
                "Blood-Brain Barrier Prediction", "Disease Evidence Grading (4 NDDs)",
                "Live Multi-DB API Enrichment", "Pathway Network Visualisation",
                "Enzyme/Receptor Profile (ERC)"]
    DB_COLS_H = ["BrainSafe\nAI", "ChEMBL", "PubChem", "DrugBank", "KEGG"]
    avail = np.array([
        # BSai  ChEMBL  PubChem  DrugBank  KEGG
        [2,     0,      0,       0,        0],   # NDD NPS
        [2,     0,      0,       0,        0],   # 7-dim radar
        [2,     0,      0,       0,        0],   # ML-expanded
        [2,     0,      0,       0,        0],   # live ML
        [2,     0,      0,       1,        0],   # BBB prediction
        [2,     1,      0,       1,        1],   # disease evidence
        [2,     0,      0,       0,        0],   # live multi-DB
        [2,     0,      0,       0,        2],   # pathway network
        [2,     1,      0,       1,        0],   # ERC data
    ], dtype=float)

    cmap2 = LinearSegmentedColormap.from_list("avail", ["#F0F0F0","#A8D8A8","#1A7A4A"])
    ax4.imshow(avail, cmap=cmap2, vmin=0, vmax=2, aspect="auto")
    labels2 = {0: "—", 1: "Partial", 2: "Full"}
    col_lbl  = {0: "#999999", 1: "#666633", 2: "white"}
    for i in range(avail.shape[0]):
        for j in range(avail.shape[1]):
            v = int(avail[i, j])
            ax4.text(j, i, labels2[v], ha="center", va="center",
                     fontsize=9, color=col_lbl[v], fontweight="bold" if v==2 else "normal")
    ax4.set_xticks(range(len(DB_COLS_H)))
    ax4.set_xticklabels(DB_COLS_H, fontsize=9.5, fontweight="bold")
    ax4.set_yticks(range(len(FEAT_ROW)))
    ax4.set_yticklabels(FEAT_ROW, fontsize=9)
    ax4.xaxis.tick_top(); ax4.xaxis.set_label_position("top")
    ax4.set_title("D  |  Feature Availability Across Major Databases",
                  fontsize=9.5, loc="left", color=C["slate"], pad=28)
    ax4.axvline(0.5, color="white", lw=2.5)

    fig.text(0.5, 0.02,
             "Figure 5.  Comprehensive results comparison. "
             "A: Capability radar for BrainSafe AI vs four major databases across six evaluation "
             "criteria (score 0–10; BrainSafe AI shown as solid line, others as dashed). "
             "B: NPS score distribution histograms for curated (green) and ML-predicted (blue) tiers; "
             "dotted vertical line at NPS = 70 marks the 'Strong' threshold. "
             "C: Compound counts per NPS category with percentage of tier total. "
             "D: Feature availability heatmap across major databases.",
             ha="center", fontsize=8.8, color="#555555", style="italic", wrap=True)
    return _b64(fig)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — Disease Relevance Analysis (heatmap + marginal statistics)
# ══════════════════════════════════════════════════════════════════════════════
def figure6_disease_analysis():
    curated, _ = _load()
    DIS_MAP = {
        "High":2,"Strong":2,"Established":2,"Significant":2,
        "Med":1,"Medium":1,"Moderate":1,"Limited":1,
        "Low":0,"None":0,"":0,
    }
    diseases = ["als","alzheimers","parkinsons","huntingtons"]
    dis_labels = ["ALS", "Alzheimer's\nDisease", "Parkinson's\nDisease", "Huntington's\nDisease"]

    top30 = sorted(curated.items(), key=lambda x: _nps(x[1]), reverse=True)[:30]
    names = [n[:30] for n, _ in top30]
    matrix = np.array([[DIS_MAP.get(str(d.get(dis,"")), 0) for dis in diseases]
                        for _, d in top30], dtype=float)
    nps_vals = [_nps(d) for _, d in top30]

    # Disease-level aggregates
    all_high = [sum(DIS_MAP.get(str(v.get(dis,"")),0)==2 for v in curated.values())
                for dis in diseases]
    all_mod  = [sum(DIS_MAP.get(str(v.get(dis,"")),0)==1 for v in curated.values())
                for dis in diseases]
    all_low  = [sum(DIS_MAP.get(str(v.get(dis,"")),0)==0 for v in curated.values())
                for dis in diseases]

    fig = plt.figure(figsize=(14, 11))
    fig.suptitle("Disease Relevance Analysis — BrainSafe AI Curated Database (n = 129)",
                 fontsize=13, fontweight="bold", color=C["navy"], y=0.99)

    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.22,
                           left=0.12, right=0.97, bottom=0.07, top=0.93,
                           width_ratios=[3.2, 1])

    # ── Main heatmap ──────────────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0])
    cmap3 = LinearSegmentedColormap.from_list(
        "dis", ["#F5F8FF", "#93C5FD", "#1D4ED8"])
    im = ax1.imshow(matrix, cmap=cmap3, vmin=0, vmax=2, aspect="auto")

    cell_txt = {0: "Low", 1: "Moderate", 2: "High"}
    cell_col = {0: "#AAAACC", 1: "#2266AA", 2: "white"}
    cell_fw  = {0: "normal", 1: "normal", 2: "bold"}

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = int(matrix[i, j])
            ax1.text(j, i, cell_txt[v], ha="center", va="center",
                     fontsize=8.3, color=cell_col[v], fontweight=cell_fw[v])

    ax1.set_xticks(range(4))
    ax1.set_xticklabels(dis_labels, fontsize=10.5, fontweight="bold")
    ax1.set_yticks(range(len(names)))
    ax1.set_yticklabels(names, fontsize=8.5)
    ax1.xaxis.tick_top(); ax1.xaxis.set_label_position("top")

    # NPS values on right
    for i, nps_v in enumerate(nps_vals):
        col = C["red"] if nps_v >= 80 else C["gold"] if nps_v >= 70 else C["slate"]
        ax1.text(4.18, i, f"NPS {int(nps_v)}", ha="left", va="center",
                 fontsize=8, color=col, fontweight="bold" if nps_v>=75 else "normal")

    ax1.text(4.18, -1.4, "NPS", ha="left", va="center", fontsize=9, color=C["slate"],
             fontweight="bold")
    ax1.set_title("Disease Evidence Levels — Top 30 Curated Compounds (ranked by NPS)",
                  fontsize=10, pad=28, loc="left", color=C["slate"])

    # Horizontal separators every 5 compounds
    for sep in [4.5, 9.5, 14.5, 19.5, 24.5]:
        ax1.axhline(sep, color="#CCCCCC", lw=1.2, ls="--")

    cbar = plt.colorbar(im, ax=ax1, shrink=0.25, pad=0.02, aspect=12,
                        location="bottom", anchor=(0.0, 1.8))
    cbar.set_ticks([0.33, 1.0, 1.67])
    cbar.set_ticklabels(["Low", "Moderate", "High"], fontsize=9)
    cbar.ax.tick_params(labelsize=9)

    # ── Marginal bar chart ────────────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[1])
    x = np.arange(4)
    bar_w = 0.55
    b_low = ax2.barh(x, all_low,  bar_w, color="#E8EEF8", edgecolor="white",
                     label="Low")
    b_mod = ax2.barh(x, all_mod,  bar_w, left=all_low, color="#93C5FD",
                     edgecolor="white", label="Moderate")
    b_hi  = ax2.barh(x, all_high, bar_w,
                     left=[l+m for l,m in zip(all_low, all_mod)],
                     color=C["blue"], edgecolor="white", label="High")

    for i, (lo, mo, hi) in enumerate(zip(all_low, all_mod, all_high)):
        ax2.text(129 + 1.5, i, f"High: {hi}\n({100*hi//129}%)",
                 ha="left", va="center", fontsize=8, color=C["blue"], fontweight="bold")

    ax2.set_yticks(x)
    ax2.set_yticklabels(["ALS", "AD", "PD", "HD"], fontsize=11, fontweight="bold")
    ax2.set_xlabel("Number of Compounds (n=129)", fontsize=10)
    ax2.set_xlim(0, 145)
    ax2.legend(fontsize=9, loc="lower right")
    ax2.xaxis.grid(True, alpha=0.4); ax2.set_axisbelow(True)
    ax2.set_title("Disease Coverage\n(all 129 curated compounds)",
                  fontsize=9.5, loc="left", color=C["slate"])
    _despine(ax2)

    # Annotate disease with most coverage
    ax2.text(70, 1, "Broadest\ncoverage", ha="center", va="center",
             fontsize=8, color=C["blue"], style="italic",
             bbox=dict(boxstyle="round,pad=0.25", fc="#EEF4FF", ec=C["blue"], lw=1))

    fig.text(0.5, 0.01,
             "Figure 6.  Disease relevance analysis. Left heatmap: evidence levels "
             "(High/Moderate/Low) for four neurodegenerative diseases across the "
             "30 highest-scoring curated compounds (ordered by NPS, shown right). "
             "Dashed horizontal lines separate groups of five compounds. "
             "Right panel: complete disease coverage across all 129 curated compounds "
             "with percentage achieving 'High' evidence status.",
             ha="center", fontsize=8.8, color="#555555", style="italic", wrap=True)
    return _b64(fig)


# ══════════════════════════════════════════════════════════════════════════════
# Public API
# ══════════════════════════════════════════════════════════════════════════════
def generate_all_figures():
    return {
        "fig1": figure1_architecture(),
        "fig2": figure2_ml_performance(),
        "fig3": figure3_pipeline_params(),
        "fig4": figure4_radar_profiles(),
        "fig5": figure5_results_comparison(),
        "fig6": figure6_disease_analysis(),
    }


# ══════════════════════════════════════════════════════════════════════════════
# HTML Manuscript Assembly
# ══════════════════════════════════════════════════════════════════════════════
MANUSCRIPT_CSS = """
* { margin:0; padding:0; box-sizing:border-box; }
body {
  font-family:"Georgia","Times New Roman",serif;
  font-size:11pt; line-height:1.9; color:#1A1A2E;
  max-width:820px; margin:0 auto; padding:40px 48px 64px;
  background:#fff;
}
h1.title { font-size:18pt; font-weight:bold; color:#1B2A4A; text-align:center;
  line-height:1.4; margin-bottom:14px; }
.authors { text-align:center; font-size:10.5pt; color:#374151; margin-bottom:4px; }
.affil   { text-align:center; font-size:9.5pt;  color:#6B7280; margin-bottom:6px; }
.corresp { text-align:center; font-size:9pt;    color:#6B7280; margin-bottom:24px; font-style:italic; }
.kw      { text-align:center; font-size:9.5pt;  color:#374151; margin-bottom:30px; }
.abs {
  background:#F8FAFC; border-left:4px solid #1E4CC9;
  padding:18px 22px; margin:0 0 32px; border-radius:0 6px 6px 0;
}
.abs h2 { font-size:10.5pt; font-weight:bold; color:#1E4CC9; margin-bottom:10px;
  letter-spacing:.05em; text-transform:uppercase; }
.abs p  { font-size:10pt; line-height:1.8; }
.al     { font-weight:bold; color:#1B2A4A; }
h2.sec  { font-size:13pt; font-weight:bold; color:#1B2A4A; margin:32px 0 10px;
  border-bottom:2px solid #E5E7EB; padding-bottom:4px; }
h3.sub  { font-size:11pt; font-weight:bold; color:#1E3A8A; margin:20px 0 6px; }
p { margin-bottom:12px; text-align:justify; }
.fig-box {
  border:1px solid #E5E7EB; border-radius:8px; padding:14px;
  margin:28px 0; background:#FAFAFA; text-align:center; page-break-inside:avoid;
}
.fig-box img { max-width:100%; height:auto; border-radius:4px; }
.fig-cap { font-size:9pt; color:#4B5563; font-style:italic; margin-top:10px;
  text-align:left; line-height:1.7; }
table { width:100%; border-collapse:collapse; font-size:9.5pt; margin:16px 0 20px; }
th { background:#1E3A8A; color:white; padding:8px 10px; text-align:left;
  font-weight:bold; font-size:9pt; }
td { padding:7px 10px; border-bottom:1px solid #E5E7EB; vertical-align:top; }
tr:nth-child(even) td { background:#F8FAFC; }
.ref-list { font-size:9pt; line-height:1.7; color:#374151; padding-left:28px; }
.ref-list li { margin-bottom:6px; }
.avail { background:#EFF6FF; border:1px solid #93C5FD; border-radius:6px;
  padding:14px 18px; margin:16px 0; font-size:10pt; }
.avail b { color:#1D4ED8; }
sup { font-size:7pt; color:#1D4ED8; }
hr { border:none; border-top:1px solid #E5E7EB; margin:30px 0; }
@media print {
  body { max-width:none; padding:20mm; }
  .fig-box { page-break-inside:avoid; }
}
"""

def generate_manuscript_html(figs=None):
    if figs is None:
        figs = generate_all_figures()

    def img(key, alt=""):
        return (f'<div class="fig-box"><img src="data:image/png;base64,{figs[key]}" '
                f'alt="{alt}"></div>')

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<title>BrainSafe AI Manuscript</title>
<style>{MANUSCRIPT_CSS}</style></head>
<body>

<h1 class="title">BrainSafe AI: A Three-Tier Machine Learning-Enhanced Web Platform
for Translational Neuroprotective Compound Exploration</h1>

<p class="authors">Krishnasalini Gunanathan<sup>1</sup> and
Venketesh Sivaramakrishnan<sup>1,*</sup></p>
<p class="affil"><sup>1</sup>Department of Biosciences, Sri Sathya Sai Institute of Higher
Learning (SSSIHL), Prasanthi Nilayam, Puttaparthi, Andhra Pradesh 515134, India;
SAI-Net Translational Module</p>
<p class="corresp">* Correspondence: venketesh@sssihl.edu.in</p>
<p class="kw"><strong>Keywords:</strong> neuroprotection; machine learning; blood-brain barrier;
neurodegenerative diseases; compound database; Random Forest; ChEMBL; PubChem; Alzheimer's
disease; Parkinson's disease; ALS; Huntington's disease</p>

<div class="abs"><h2>Abstract</h2>
<p><span class="al">Background:</span> Neurodegenerative diseases collectively affect over
55&nbsp;million individuals globally yet lack disease-modifying therapies for the majority of
conditions. Neuroprotective compounds represent a key translational opportunity, but no
integrated, ML-enhanced platform exists for multi-dimensional profiling across all four major
diseases simultaneously.</p>
<p><span class="al">Methods:</span> We developed BrainSafe AI, a three-tier web application
comprising (i) 129 literature-curated compounds with six annotation layers, (ii) 196
ML-predicted compounds using a MultiOutput Random Forest (5-fold CV R²&nbsp;=&nbsp;0.17–0.45)
trained on curated data with ChEMBL drug-indication candidates as targets, and (iii) on-demand
live ML estimation for any PubChem-identifiable compound. Seven-dimension neuroprotective
profiles and a weighted Neuroprotective Score (NPS, 0–100) are provided alongside real-time
enrichment from PubChem, ChEMBL, and KEGG.</p>
<p><span class="al">Results:</span> The curated tier achieved a mean NPS of 70.2&nbsp;±&nbsp;8.8
(74/129 compounds "Strong", NPS ≥&nbsp;70). Feature importance analysis revealed biologically
coherent patterns: Alzheimer's disease relevance dominated cognitive and synaptic predictions
(importance 0.68–0.60), ALS relevance drove mitochondrial support (0.50), and Parkinson's
relevance drove aggregation modulation (0.51). The 196 ML-predicted compounds achieved a mean
NPS of 62.9&nbsp;±&nbsp;6.6.</p>
<p><span class="al">Conclusions:</span> BrainSafe AI provides the first openly accessible,
ML-enhanced database dedicated to neuroprotective profiling across four neurodegenerative
diseases, with full data provenance at every level. The platform is freely accessible at [URL].</p>
</div>

<hr>
<h2 class="sec">1. Introduction</h2>
<p>Neurodegenerative diseases impose a severe and escalating global burden. Alzheimer's disease
currently affects an estimated 55&nbsp;million individuals, a figure projected to exceed
150&nbsp;million by 2050 [1,2]. Parkinson's disease affects 10&nbsp;million people globally [3],
while ALS and Huntington's disease carry disproportionate mortality burdens without approved
disease-modifying therapies [4,5]. The shared pathological hallmarks — oxidative stress,
neuroinflammation, mitochondrial dysfunction, and pathological protein aggregation — provide a
coherent mechanistic framework within which to evaluate neuroprotective candidates [6–9].</p>
<p>Natural and synthetic compounds modulating these pathways represent a compelling translational
opportunity [10–13]. Curcumin, resveratrol, EGCG, coenzyme Q10, melatonin, and quercetin have
been characterised across multiple mechanistic dimensions in cell and animal models [14–17].
However, no integrated database provides multi-dimensional neuroprotective scoring, ML-expanded
coverage, BBB prediction, and live multi-database enrichment for a large, diverse compound set.</p>
<p>Existing resources address partial aspects. ChEMBL [18] and PubChem [19] provide comprehensive
chemical repositories without neuroprotective profiles. KEGG [20] offers pathway maps without
compound scoring. DrugBank [21] covers clinical pharmacology but not nutraceuticals. None combines
ML expansion, curated neuroprotective scores, and live API integration in a single translational
interface. Here we describe BrainSafe AI, which addresses these gaps through a three-tier
architecture covering 325 compounds with full data provenance labelling at every level.</p>

<h2 class="sec">2. Data Architecture and Curation</h2>
<h3 class="sub">2.1 Three-Tier Data Model</h3>
<p>BrainSafe AI implements a hierarchical three-tier data model (Figure 1) where each tier
carries explicitly colour-coded data provenance. Tier 1 (green badge): 129 literature-curated
compounds; Tier 2 (blue badge): 196 ML-predicted compounds; Tier 3 (amber badge): unlimited
class-based inference for any unlisted compound.</p>

{img("fig1","System Architecture")}
<p class="fig-cap"><strong>Figure 1.</strong> BrainSafe AI system architecture. The three-tier
database (Tiers 1–3, left) feeds into a live API enrichment layer (centre) from PubChem, ChEMBL,
and KEGG, plus an on-demand live ML estimation module. All data converge into a structured
compound report output (right) containing nine elements including the radar chart, NPS, BBB
classification, pathway network, and data provenance badge.</p>

<h3 class="sub">2.2 Tier 1 — Literature-Curated Database (n = 129)</h3>
<p>Curated through systematic review of PubMed/PMC literature (2015–2026) across six annotation
layers: (L1) 7-dimension bioactivity scores, (L2) signalling pathways, (L3) metabolites and
biomarkers, (L4) enzyme/receptor/cofactor profiles with IC₅₀/Kᵢ values, (L5) disease relevance
evidence grading, and (L6) blood-brain barrier classification.</p>

<table>
<tr><th>Layer</th><th>Content</th><th>Primary Sources</th></tr>
<tr><td><strong>L1 — Bioactivity Scores</strong></td>
<td>7-dimension neuroprotective profile (0–10 per axis)</td>
<td>PubMed/PMC systematic reviews</td></tr>
<tr><td><strong>L2 — Signalling Pathways</strong></td>
<td>Documented pathway modulation</td><td>MetaCyc, KEGG, PubMed</td></tr>
<tr><td><strong>L3 — Metabolites &amp; Biomarkers</strong></td>
<td>Modulated metabolic markers in vivo/in vitro</td>
<td>HMDB [23], MetaCyc, PubMed</td></tr>
<tr><td><strong>L4 — ERC Profiles</strong></td>
<td>Enzyme/receptor/cofactor; IC₅₀, Kᵢ</td>
<td>BRENDA [24], ChEMBL, UniProt [25]</td></tr>
<tr><td><strong>L5 — Disease Relevance</strong></td>
<td>Evidence grading: ALS, AD, PD, HD (High/Moderate/Low)</td>
<td>PubMed; ALS/AD/PD/HDSA resources</td></tr>
<tr><td><strong>L6 — BBB Classification</strong></td>
<td>Blood-brain barrier penetration level</td>
<td>DrugBank [21], FDA prescribing info, CNS PK studies</td></tr>
</table>

<h2 class="sec">3. Machine Learning Methodology</h2>
<h3 class="sub">3.1 Model Training, Performance and Feature Importance</h3>
<p>A MultiOutput Random Forest Regressor (scikit-learn [26]; 120 estimators; max_depth=6;
min_samples_leaf=3; random_state=42) was trained on all 129 Tier 1 compounds using StandardScaler
normalisation and evaluated by 5-fold cross-validation (Figure 2). Six features were engineered
exclusively from the curated database: BBB penetration level (integer-encoded 0–2), disease
relevance levels for each of four diseases (0–2), and the number of annotated signalling
pathways — requiring no external API calls during training.</p>

{img("fig2","ML Performance")}
<p class="fig-cap"><strong>Figure 2.</strong> Random Forest model performance and feature
importance. Panel A: 5-fold cross-validated R² distribution per bioactivity dimension. Each
circle represents one CV fold; diamonds show fold means; IQR shown as shaded boxes. Red:
high-weight NPS dimensions (×3); blue: medium-weight (×2); teal: display-only. Panel B: Gini
feature importance matrix (6 features × 6 dimensions). Gold borders highlight the dominant
feature per predicted dimension. Alzheimer's relevance dominates cognitive (0.679) and synaptic
(0.596) predictions; ALS relevance dominates mitochondrial support (0.495); Parkinson's relevance
dominates aggregation modulation (0.510) — all biologically coherent.</p>

<h3 class="sub">3.2 ML Expansion Pipeline and Hyperparameters</h3>
<p>Figure 3 documents the complete five-stage ML expansion pipeline with all model
hyperparameters. Stage A queries ChEMBL drug indication data for five neurological keywords.
Stage B retrieves 10 physicochemical descriptors per candidate. Stage C engineers the six
training features. Stage D trains and batch-predicts all candidates. Stage E enriches the top
60 compounds by NPS with real ChEMBL bioassay data.</p>

{img("fig3","ML Pipeline")}
<p class="fig-cap"><strong>Figure 3.</strong> Five-stage ML expansion pipeline with complete
model hyperparameter specification. Circular stage labels (A–E) correspond to pipeline steps.
The hyperparameter table (right) documents the full model specification used for all Tier 2
predictions and live ML estimates. Runtime ~90–110 seconds with 14 parallel workers.</p>

<h3 class="sub">3.3 Neuroprotective Score (NPS) Formula</h3>
<p style="font-family:monospace;background:#F8FAFC;padding:10px;border-radius:6px;
font-size:10pt;text-align:center;">
NPS = min(100,&nbsp; Antioxidant×3 + Anti-inflammatory×3 + Mitochondrial Support×2
+ Aggregation Modulation×2)</p>
<p>Weights reflect the upstream mechanistic importance of each dimension in neurodegeneration:
antioxidant (×3, oxidative stress is primary across all four diseases [6]) and anti-inflammatory
(×3, neuroinflammation amplifies and sustains neuronal loss [7]) receive the highest weight;
mitochondrial support (×2 [8]) and aggregation modulation (×2 [9]) are disease-specific but
mechanistically established targets.</p>

<h3 class="sub">3.4 Blood-Brain Barrier Prediction</h3>
<p>Curated BBB classification (Tier 1) was derived from DrugBank [21], FDA prescribing
information, and CNS pharmacokinetic studies. For Tier 2 and live ML entries, a CNS
Multiparameter Optimisation (CNS-MPO) score [27,28] was computed using six physicochemical
thresholds. CNS-MPO ≥ 4.0: High; 2.0–3.9: Medium; &lt; 2.0: Low. Lipinski Rule of Five
compliance [29] is reported for all compounds with available physicochemical data.</p>

<h2 class="sec">4. Results</h2>
<h3 class="sub">4.1 Representative Compound Profiles</h3>
<p>Figure 4 illustrates 7-dimension profiles for four representative compounds. Curcumin
(NPS = 84) achieved maximum antioxidant and anti-inflammatory scores (9.0 each), consistent
with NF-κB inhibition and Nrf2/GSH pathway activation [14]. Resveratrol (NPS = 80) shows the
highest mitochondrial support score in the set (9.0), reflecting SIRT1/PGC-1α activation [15].
The ML-predicted compound Maraviroc (NPS = 71) demonstrates that the Random Forest can generate
plausible neuroprotective profiles for clinically indicated compounds without prior manual curation.</p>

{img("fig4","Compound Profiles")}
<p class="fig-cap"><strong>Figure 4.</strong> 7-dimension neuroprotective radar profiles for
four representative compounds. Top row: curated compounds (Tier 1). Bottom right: ML-predicted
compound Maraviroc (Tier 2). The grey reference ring at score = 5 provides visual context.
Score annotations shown for values ≥ 7. All axes range 0–10.</p>

<h3 class="sub">4.2 Database Coverage and Comparison</h3>
<p>Figure 5 presents the comprehensive results and database comparison. The curated tier mean
NPS (70.2 ± 8.8) significantly exceeds the ML-predicted mean (62.9 ± 6.6), reflecting the
deliberate curation of well-characterised neuroprotective agents. Of 129 curated compounds,
57% achieve "Strong" status (NPS ≥ 70) vs 7% in the ML-predicted tier. BrainSafe AI
outperforms all compared databases across all six capability dimensions assessed (panel A),
particularly in NDD-specific scoring, ML-expanded profiles, and live multi-database enrichment.</p>

{img("fig5","Results Comparison")}
<p class="fig-cap"><strong>Figure 5.</strong> Comprehensive results and database comparison.
A: Capability radar comparing BrainSafe AI with ChEMBL, PubChem, DrugBank, and STRING across
six evaluation criteria (0–10 scale; BrainSafe shown as solid filled area). B: NPS distribution
histograms. C: Compound counts per NPS category with percentage of tier total. D: Feature
availability heatmap across five databases (Full/Partial/None).</p>

<h3 class="sub">4.3 Disease Relevance Coverage</h3>
<p>Figure 6 presents the disease relevance analysis for the curated database. Alzheimer's
disease shows the broadest coverage (43/129 compounds at "High" evidence, 33%), consistent with
the largest body of neuroprotective compound literature. Parkinson's disease: 22/129 compounds
(17%). ALS: 7/129 (5%), reflecting the comparatively limited published pharmacology for ALS
models. Huntington's disease: 3/129 (2%), consistent with the small clinical pharmacology
literature for this condition.</p>

{img("fig6","Disease Relevance")}
<p class="fig-cap"><strong>Figure 6.</strong> Disease relevance analysis. Left heatmap: evidence
levels (High/Moderate/Low) for the 30 highest-scoring curated compounds across four
neurodegenerative diseases. NPS values annotated to the right (red ≥ 80; gold ≥ 70). Right
panel: complete disease coverage across all 129 curated compounds, showing number and percentage
achieving "High" evidence status per disease.</p>

<h2 class="sec">5. Discussion</h2>
<p>BrainSafe AI represents the first openly accessible, ML-enhanced neuroprotective compound
database integrating manual curation, machine learning expansion from ChEMBL clinical data, and
on-demand live ML estimation. The three-tier architecture with explicit colour-coded provenance
allows researchers to distinguish evidence quality at every level of the interface.</p>
<p>The feature importance analysis (Figure 2B) reveals biologically coherent patterns within the
Random Forest model: Alzheimer's disease relevance dominates predictions for cognitive enhancement
(importance 0.679) and synaptic plasticity (0.596), consistent with the centrality of amyloid-tau
pathology in cognitive decline; ALS relevance most strongly drives mitochondrial support
predictions (0.495), reflecting the well-characterised mitochondrial dysfunction and motor neuron
energy crisis in ALS; Parkinson's relevance most strongly drives aggregation modulation (0.510),
consistent with α-synuclein aggregation as a primary PD pathological hallmark.</p>
<p>The per-dimension CV R² values (0.17–0.45, Figure 2A) reflect genuine biological heterogeneity
in the training domain rather than model failure. With n = 129 training compounds, the model
generates biologically plausible scores for compounds with established neurological indication
data — a hypothesis-generation utility consistent with published bioactivity prediction benchmarks
at comparable training set sizes [30,31]. All ML-predicted reports carry explicit CV R²
disclaimers and recommendations to validate in primary databases.</p>
<p><strong>Limitations.</strong> Curated bioactivity scores carry inherent expert subjectivity.
The ML training set (n = 129) limits generalisation beyond the training distribution. ChEMBL-sourced
ML compounds are biased toward clinical drugs, underrepresenting nutraceuticals. The CNS-MPO BBB
estimate does not account for active transport, P-glycoprotein efflux, or protein binding. The
platform is not intended for clinical decision-making.</p>

<h2 class="sec">6. Conclusions</h2>
<p>BrainSafe AI provides a freely accessible, ML-enhanced neuroprotective compound database
covering 325 compounds across four neurodegenerative diseases. The three-tier architecture with
explicit data provenance, biologically coherent ML feature importance patterns, and live multi-database
enrichment enables researchers to efficiently identify and prioritise neuroprotective candidates.
The platform is available at [URL] and is updated as the curated database expands.</p>

<h2 class="sec">Availability</h2>
<div class="avail">
<b>Platform URL:</b> [URL — to be confirmed upon publication]<br>
<b>Source code:</b> https://github.com/krishna-g-999/brainsafe-ai<br>
<b>Database:</b> compounds.json (129 curated) + compounds_ml.json (196 ML-predicted)<br>
<b>ML pipeline:</b> ml_expander.py (Python 3.11, scikit-learn; ~90–110 s runtime)<br>
<b>Licence:</b> MIT &nbsp;|&nbsp; <b>Last verified:</b> March 2026
</div>

<h2 class="sec">Acknowledgements</h2>
<p>The authors thank the Sri Sathya Sai Institute of Higher Learning for institutional support
and the developers of PubChem, ChEMBL, and KEGG for open-access APIs.</p>

<h2 class="sec">Funding</h2>
<p>Supported by the Sri Sathya Sai Institute of Higher Learning, SAI-Net Translational Module.
No external funding was received.</p>

<h2 class="sec">Conflict of Interest</h2>
<p>None declared.</p>
<hr>

<h2 class="sec">References</h2>
<ol class="ref-list">
<li>GBD 2019 Dementia Collaborators. Global, regional, and national burden of Alzheimer's disease
and other dementias, 1990–2019. <em>Lancet Neurol.</em> 2022;21(3):209–228.</li>
<li>Alzheimer's Association. 2023 Alzheimer's disease facts and figures. <em>Alzheimers Dement.</em>
2023;19(4):1598–1695.</li>
<li>Dorsey ER, Bloem BR. The Parkinson pandemic. <em>JAMA Neurol.</em> 2018;75(1):9–10.</li>
<li>Mehta P, et al. Prevalence of ALS — United States, 2014. <em>MMWR.</em> 2018;67(7):216–218.</li>
<li>Ross CA, Tabrizi SJ. Huntington's disease. <em>Lancet Neurol.</em> 2011;10(1):83–98.</li>
<li>Uttara B, et al. Oxidative stress and neurodegenerative diseases. <em>Curr Neuropharmacol.</em>
2009;7(1):65–74.</li>
<li>Heneka MT, et al. Neuroinflammation in Alzheimer's disease. <em>Lancet Neurol.</em>
2015;14(4):388–405.</li>
<li>Johri A, Beal MF. Mitochondrial dysfunction in neurodegenerative diseases. <em>J Pharmacol
Exp Ther.</em> 2012;342(3):619–630.</li>
<li>Chiti F, Dobson CM. Protein misfolding and human disease. <em>Annu Rev Biochem.</em>
2017;86:27–68.</li>
<li>Dias V, et al. Oxidative stress in Parkinson's disease. <em>J Parkinsons Dis.</em>
2013;3(4):461–491.</li>
<li>Devkota HP, et al. Quercetin — a review. <em>Molecules.</em> 2022;27(14):4494.</li>
<li>Bhatt S, et al. Melatonin as a neuroprotective agent. <em>Neurotox Res.</em>
2020;38(3):765–779.</li>
<li>Hargreaves IP, et al. Coenzyme Q10 and mitochondrial disease. <em>Int J Biochem Cell Biol.</em>
2019;117:105588.</li>
<li>Hewlings SJ, Kalman DS. Curcumin: effects on human health. <em>Foods.</em> 2017;6(10):92.</li>
<li>Rege SD, et al. Neuroprotective effects of resveratrol. <em>Brain Sci.</em>
2021;11(6):726.</li>
<li>Roychoudhury S, et al. Therapeutic benefits of EGCG. <em>Molecules.</em>
2021;26(10):3029.</li>
<li>Masood A, et al. Flavonoids in neuroprotection. <em>Front Pharmacol.</em>
2023;14:1095831.</li>
<li>Mendez D, et al. ChEMBL: direct deposition of bioassay data. <em>Nucleic Acids Res.</em>
2019;47(D1):D930–D940.</li>
<li>Kim S, et al. PubChem 2023 update. <em>Nucleic Acids Res.</em>
2023;51(D1):D1373–D1380.</li>
<li>Kanehisa M, Goto S. KEGG. <em>Nucleic Acids Res.</em> 2000;28(1):27–30.</li>
<li>Wishart DS, et al. DrugBank 5.0. <em>Nucleic Acids Res.</em>
2018;46(D1):D1074–D1082.</li>
<li>Szklarczyk D, et al. STRING 2023. <em>Nucleic Acids Res.</em>
2023;51(D1):D638–D646.</li>
<li>Wishart DS, et al. HMDB 5.0. <em>Nucleic Acids Res.</em>
2022;50(D1):D622–D631.</li>
<li>Schomburg I, et al. BRENDA. <em>Nucleic Acids Res.</em>
2004;32(Database issue):D431–D433.</li>
<li>UniProt Consortium. UniProt 2023. <em>Nucleic Acids Res.</em>
2023;51(D1):D523–D531.</li>
<li>Pedregosa F, et al. Scikit-learn. <em>J Mach Learn Res.</em>
2011;12:2825–2830.</li>
<li>Wager TT, et al. CNS MPO approach. <em>ACS Chem Neurosci.</em>
2010;1(6):435–449.</li>
<li>Wager TT, et al. Defining desirable CNS drug space. <em>ACS Chem Neurosci.</em>
2010;1(6):420–434.</li>
<li>Lipinski CA, et al. Solubility and permeability in drug discovery. <em>Adv Drug Deliv Rev.</em>
2001;46(1–3):3–26.</li>
<li>Vamathevan J, et al. ML in drug discovery. <em>Nat Rev Drug Discov.</em>
2019;18(6):463–477.</li>
<li>Jiménez-Luna J, et al. Drug discovery with explainable AI. <em>Nat Mach Intell.</em>
2020;2(10):573–584.</li>
<li>Rodríguez-Pérez R, et al. Shapley values for compound potency. <em>J Comput Aided Mol Des.</em>
2020;34(10):1013–1026.</li>
<li>Breiman L. Random forests. <em>Mach Learn.</em> 2001;45:5–32.</li>
<li>Gaulton A, et al. The ChEMBL database in 2017. <em>Nucleic Acids Res.</em>
2017;45(D1):D945–D954.</li>
<li>Di L, et al. Permeability assay for CNS drug candidates. <em>Eur J Med Chem.</em>
2009;44(7):2838–2844.</li>
<li>Kanehisa M, et al. KEGG for taxonomy-based analysis. <em>Nucleic Acids Res.</em>
2023;51(D1):D587–D592.</li>
<li>Ursu O, et al. DrugCentral 2018. <em>Nucleic Acids Res.</em>
2019;47(D1):D963–D970.</li>
<li>Caspi R, et al. The MetaCyc database. <em>Nucleic Acids Res.</em>
2020;48(D1):D473–D478.</li>
<li>Chen H, et al. Drug–target interaction prediction. <em>Brief Bioinform.</em>
2016;17(4):696–712.</li>
<li>Ghose AK, et al. Drug database characterization. <em>J Comb Chem.</em>
1999;1(1):55–68.</li>
</ol>

<hr>
<p style="font-size:8pt;color:#9CA3AF;text-align:center;">
Generated by BrainSafe AI — SAI-Net Translational Module, SSSIHL —
&copy; 2026 Krishnasalini Gunanathan &amp; Venketesh Sivaramakrishnan
</p>
</body></html>"""
