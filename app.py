import streamlit as st
import json
import difflib
import math
import base64
import os
import plotly.graph_objects as go
import networkx as nx
import numpy as np
from scorer import neuro_score
from pubchem_client import fetch_pubchem, fetch_chembl, fetch_kegg_pathways
try:
    from ml_v3_knn import predict_nps_knn
    from ml_v3_engine import predict_disease_relevance, predict_pathways, get_chembl_targets
    V3_KNN_AVAILABLE = True
except ImportError:
    V3_KNN_AVAILABLE = False

# ── v2 Feature Engineering constants (biological domain knowledge) ────────
POLYPHENOL_TYPES = {"flavonoid","polyphenol","catechin","stilbene","terpene",
                    "carotenoid","vitamin","phenolic","alkaloid","curcuminoid"}
NEURO_KWS        = {"bdnf","trkb","wnt","ngf","neurogenesis","hippocampus",
                    "creb","notch","shh","vegf","fgf","sox2","nestin"}


st.set_page_config(
    page_title="BrainSafe AI | SAI-Net",
    page_icon="",
    layout="wide",
    initial_sidebar_state="collapsed"
)

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');
    * { font-family: 'Inter', system-ui, -apple-system, sans-serif; }
    [data-testid="stAppViewContainer"] { background: #EEF2F9; }
    [data-testid="stHeader"] { display: none; }
    .block-container { padding: 0 2.2rem 3rem 2.2rem; max-width: 1520px; }
    [data-testid="stVerticalBlock"] > [data-testid="stVerticalBlockBorderWrapper"] {
        background: white !important;
        border-radius: 12px !important;
        border: 1px solid #E2E8F2 !important;
        box-shadow: 0 2px 8px rgba(13,33,55,0.06) !important;
    }
    .section-label {
        font-size: 0.68rem; font-weight: 800; letter-spacing: 1.6px;
        text-transform: uppercase; color: #94A3B8; margin-bottom: 12px;
        border-bottom: 1px solid #F0F4FA; padding-bottom: 9px;
    }
    .compound-banner {
        background: linear-gradient(135deg, #071626 0%, #0D2137 60%, #1A3A5C 100%);
        border-radius: 12px; padding: 24px 30px; margin-bottom: 22px;
        border-left: 5px solid #F0A500;
        box-shadow: 0 4px 20px rgba(13,33,55,0.25);
    }
    .cmpd-name { color: #FFFFFF; font-size: 2.0rem; font-weight: 800; margin: 0 0 0 0; padding: 0; line-height: 1.2; letter-spacing: -0.3px; }
    .cmpd-type { color: #90B4CF; font-size: 0.88rem; margin: 4px 0 0 0; padding: 0; font-weight: 500; }
    .pill {
        display: inline-block; padding: 5px 14px; border-radius: 20px;
        font-size: 0.79rem; font-weight: 700; margin: 4px 4px 4px 0;
        letter-spacing: 0.2px;
    }
    .site-header {
        background: linear-gradient(135deg, #040E1C 0%, #0A1929 40%, #0D2137 70%, #112A47 100%);
        padding: 32px 52px 30px; margin: -1rem -2.2rem 0 -2.2rem;
        border-bottom: 3px solid #F0A500;
        box-shadow: 0 6px 30px rgba(0,0,0,0.35);
        position: relative;
    }
    .site-header::before {
        content: ''; position: absolute; inset: 0;
        background: radial-gradient(ellipse at 20% 50%, rgba(240,165,0,0.04) 0%, transparent 60%),
                    radial-gradient(ellipse at 80% 50%, rgba(30,90,160,0.08) 0%, transparent 60%);
        pointer-events: none;
    }
    .header-logo-wrap {
        display: flex; align-items: center; gap: 28px; position: relative; z-index: 1;
    }
    .header-logo-img {
        height: 96px; width: 96px; border-radius: 50%;
        border: 2.5px solid #F0A500;
        box-shadow: 0 0 24px rgba(240,165,0,0.35), 0 4px 14px rgba(0,0,0,0.45);
        object-fit: cover; flex-shrink: 0;
    }
    .header-inst-logo {
        height: 96px; width: 96px; border-radius: 50%;
        border: 2.5px solid rgba(255,255,255,0.55);
        box-shadow: 0 0 18px rgba(255,255,255,0.15), 0 4px 14px rgba(0,0,0,0.4);
        object-fit: contain; flex-shrink: 0;
        background: rgba(255,255,255,0.95);
        padding: 6px;
    }
    .header-text-block { flex: 1; }
    .header-divider {
        width: 2px; height: 88px; background: linear-gradient(to bottom, transparent, #F0A500, transparent);
        flex-shrink: 0;
    }
    .header-title {
        color: #FFFFFF; font-size: 2.55rem; font-weight: 800; margin: 0;
        letter-spacing: -0.5px; line-height: 1.15;
        text-shadow: 0 2px 8px rgba(0,0,0,0.3);
    }
    .header-title span { color: #F0A500; }
    .header-sub { color: #88AECF; font-size: 0.97rem; margin-top: 5px; font-weight: 500; }
    .header-tags {
        color: #4A7FA0; font-size: 0.80rem; margin-top: 8px;
        display: flex; gap: 8px; flex-wrap: wrap;
    }
    .header-tag-pill {
        background: rgba(255,255,255,0.06); border: 1px solid rgba(255,255,255,0.1);
        border-radius: 20px; padding: 2px 10px; font-size: 0.74rem;
        color: #7CAAC9; white-space: nowrap;
    }
    .header-stats {
        display: flex; flex-direction: column; align-items: center; gap: 2px;
        background: rgba(240,165,0,0.08); border: 1px solid rgba(240,165,0,0.2);
        border-radius: 10px; padding: 12px 22px; flex-shrink: 0;
    }
    .header-stat-num { color: #F0A500; font-size: 1.7rem; font-weight: 800; line-height: 1; }
    .header-stat-lbl { color: #7CAAC9; font-size: 0.69rem; font-weight: 700;
                       text-transform: uppercase; letter-spacing: 1.2px; white-space: nowrap; }
    .metric-row { display: flex; gap: 18px; margin: 12px 0; }
    .metric-box {
        flex: 1; background: #F7F9FD; border-radius: 10px;
        padding: 14px 18px; border: 1px solid #E4EAF4;
    }
    .metric-label { font-size: 0.70rem; color: #8896AD; font-weight: 700; text-transform: uppercase; letter-spacing: 1px; }
    .metric-val { font-size: 1.55rem; font-weight: 800; margin-top: 3px; }
    .pathway-chip {
        display: inline-block; background: #EEF4FF; color: #1A3A5C;
        border: 1px solid #C5D7F0; border-radius: 6px;
        padding: 5px 13px; margin: 3px 4px 3px 0;
        font-size: 0.82rem; font-weight: 500;
        box-shadow: 0 1px 3px rgba(30,58,92,0.06);
    }
    .met-chip {
        display: inline-block; background: #F5F0FC; color: #4A1F7A;
        border: 1px solid #D8C8F5; border-radius: 6px;
        padding: 5px 13px; margin: 3px 4px 3px 0;
        font-size: 0.82rem; font-weight: 500;
        box-shadow: 0 1px 3px rgba(74,31,122,0.06);
    }
    .br-chip {
        display: inline-block; background: #FFF7EC; color: #7A4000;
        border: 1px solid #F2D9B0; border-radius: 6px;
        padding: 5px 13px; margin: 3px 4px 3px 0;
        font-size: 0.82rem; font-weight: 500;
        box-shadow: 0 1px 3px rgba(122,64,0,0.06);
    }
    table { border-collapse: collapse; width: 100%; }
    th { background: #F4F7FC; color: #4A5568; font-size: 0.71rem; font-weight: 800;
         text-transform: uppercase; letter-spacing: 1px; padding: 11px 16px; text-align: left;
         border-bottom: 2px solid #DDE5F0; }
    td { padding: 11px 16px; border-bottom: 1px solid #F0F4FA; font-size: 0.87rem;
         color: #1E293B; vertical-align: middle; }
    tr:last-child td { border-bottom: none; }
    tr:hover td { background: #F8FAFF; }
    .disclaimer {
        background: #FFFBEE; border: 1px solid #EDD860; border-radius: 8px;
        padding: 12px 18px; font-size: 0.80rem; color: #7A5A00; margin-top: 18px;
    }
    .about-quote {
        border-left: 4px solid #F0A500; padding: 12px 20px;
        background: #FFFBF0; border-radius: 0 8px 8px 0; margin: 14px 0;
        font-style: italic; color: #5A4000; font-size: 0.94rem;
        box-shadow: 0 1px 6px rgba(240,165,0,0.1);
    }
    div[data-testid="stExpander"] {
        border: 1px solid #E2E8F2 !important; border-radius: 12px !important;
        background: white !important; box-shadow: 0 1px 5px rgba(0,0,0,0.04) !important;
    }
    div[data-baseweb="tab-list"] {
        background: #F4F7FC !important; border-radius: 8px !important; padding: 4px !important;
        border: 1px solid #E2E8F2 !important; margin-top: 18px;
    }
    div[data-baseweb="tab"] {
        border-radius: 6px !important; font-weight: 600 !important;
        font-size: 0.88rem !important; padding: 8px 20px !important;
    }
    /* ── Primary buttons → gold/orange ── */
    button[kind="primary"], .stButton > button[data-testid="baseButton-primary"],
    div.stButton > button[type="submit"] {
        background: #F0A500 !important;
        color: #0D2137 !important;
        border: none !important;
        font-weight: 700 !important;
    }
    button[kind="primary"]:hover {
        background: #D4920A !important;
        color: #0D2137 !important;
    }
    .stButton > button[kind="primary"] {
        background: #F0A500 !important;
        color: #0D2137 !important;
        border: none !important;
        font-weight: 700 !important;
    }
    .network-legend-item {
        display: inline-flex; align-items: center; gap: 6px;
        margin-right: 18px; font-size: 0.78rem; color: #475569; font-weight: 500;
    }
    .network-legend-dot {
        width: 11px; height: 11px; border-radius: 50%; flex-shrink: 0;
    }
    .network-legend-line {
        width: 22px; height: 2px; flex-shrink: 0; border-radius: 2px;
    }
    [data-testid="stMarkdownContainer"] > p,
    [data-testid="stMarkdownContainer"] > ul > li,
    [data-testid="stMarkdownContainer"] > ol > li,
    [data-testid="stMarkdownContainer"] > blockquote { color: #1E293B !important; }
    [data-testid="stMarkdownContainer"] > table td,
    [data-testid="stMarkdownContainer"] > table th { color: #1E293B !important; }
    [data-testid="stMarkdownContainer"] h1,
    [data-testid="stMarkdownContainer"] h2,
    [data-testid="stMarkdownContainer"] h3,
    [data-testid="stMarkdownContainer"] h4 { color: #0D2137 !important; }
    [data-testid="stMarkdownContainer"] > p strong,
    [data-testid="stMarkdownContainer"] > ul > li strong { color: #0D2137 !important; }
    [data-testid="stMarkdownContainer"] > p em { color: #334155 !important; }
    .header-title { color: #FFFFFF !important; }
    .header-title span { color: #F0A500 !important; }
    .header-sub { color: #88AECF !important; }
    .header-tag-pill { color: #7CAAC9 !important; }
    .cmpd-name { color: #FFFFFF !important; }
    .cmpd-type { color: #90B4CF !important; }
    .section-label { color: #94A3B8 !important; }
    .pathway-chip { color: #1A3A5C !important; }
    .met-chip { color: #4A1F7A !important; }
    .br-chip { color: #7A4000 !important; }
    .about-quote { color: #5A4000 !important; }
    .disclaimer { color: #7A5A00 !important; }
    .estimated-banner {
        background: #FFFBEE; border: 1px solid #F0C040; border-radius: 10px;
        padding: 13px 18px; margin-bottom: 18px;
        display: flex; gap: 12px; align-items: flex-start;
    }
    .estimated-icon { font-size: 1.2rem; margin-top: 1px; flex-shrink: 0; }
    .estimated-text { font-size: 0.86rem; color: #7A5500 !important; }
    .estimated-text strong { color: #5A3A00 !important; font-size: 0.9rem; }

    /* ── Prominent tab bar ── */
    [data-testid="stTabs"] {
        background: #FFFFFF;
        border: 1.5px solid #D6E0EE;
        border-radius: 12px;
        padding: 0;
        box-shadow: 0 2px 12px rgba(13,33,55,0.08);
        overflow: hidden;
    }
    .stTabs [data-baseweb="tab-list"] {
        background: #F0F5FB;
        border-bottom: 2px solid #D6E0EE;
        padding: 0 10px;
        gap: 0;
    }
    .stTabs [data-baseweb="tab"] {
        background: transparent;
        border-radius: 0;
        color: #4A6080;
        font-size: 0.93rem;
        font-weight: 600;
        padding: 14px 24px;
        border-bottom: 3px solid transparent;
        margin-bottom: -2px;
        transition: color 0.2s;
    }
    .stTabs [aria-selected="true"] {
        color: #0D2137 !important;
        border-bottom: 3px solid #F0A500 !important;
        background: transparent !important;
    }
    .stTabs [data-baseweb="tab"]:hover {
        color: #0D2137;
        background: rgba(240,165,0,0.05) !important;
    }
    .stTabs [data-baseweb="tab-panel"] {
        padding: 24px 20px 20px 20px;
    }
    .stTabs [data-baseweb="tab-border"] { display: none; }
</style>
""", unsafe_allow_html=True)


@st.cache_data
def load_compounds():
    with open("compounds.json", "r") as f:
        curated = json.load(f)
    merged = dict(curated)
    if os.path.exists("compounds_ml.json"):
        with open("compounds_ml.json", "r") as f:
            ml_raw = json.load(f)
        ml_raw.pop("_ml_metadata", None)
        for name, entry in ml_raw.items():
            if name not in merged:
                merged[name] = entry
    return merged

@st.cache_data
def load_enzyme_data():
    d = {}
    if os.path.exists("compound_enzyme_data.json"):
        with open("compound_enzyme_data.json", "r") as f:
            base = json.load(f)
        base.pop("_metadata", None)
        d.update(base)
    if os.path.exists("erc_ml.json"):
        with open("erc_ml.json", "r") as f:
            ml_erc = json.load(f)
        for name, erc in ml_erc.items():
            if name not in d:
                d[name] = erc
    return d

@st.cache_data
def get_db_stats() -> dict:
    with open("compounds.json", "r") as f:
        curated = json.load(f)
    n_curated = len(curated)
    n_ml = 0
    if os.path.exists("compounds_ml.json"):
        with open("compounds_ml.json", "r") as f:
            ml_raw = json.load(f)
        ml_raw.pop("_ml_metadata", None)
        n_ml = len(ml_raw)
    return {"curated": n_curated, "ml": n_ml,
            "natural": __import__("sqlite3").connect("brainsafe_natural.db").execute(
                "SELECT COUNT(*) FROM compounds").fetchone()[0]
            if __import__("os").path.exists("brainsafe_natural.db") else 0,
            "total": n_curated + n_ml + (
                __import__("sqlite3").connect("brainsafe_natural.db").execute(
                "SELECT COUNT(*) FROM compounds").fetchone()[0]
                if __import__("os").path.exists("brainsafe_natural.db") else 0)}

@st.cache_data
def get_logo_b64():
    if os.path.exists('sai_net_logo.png'):

        with open('sai_net_logo.png','rb') as f:

            return base64.b64encode(f.read()).decode()

    return 'PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIyMDAiIGhlaWdodD0iMjAwIj48Y2lyY2xlIGN4PSIxMDAiIGN5PSIxMDAiIHI9Ijk2IiBmaWxsPSIjMEQyMTM3IiBzdHJva2U9IiNGMEE1MDAiIHN0cm9rZS13aWR0aD0iNCIvPjx0ZXh0IHg9IjEwMCIgeT0iMTMwIiBmb250LWZhbWlseT0iQXJpYWwgQmxhY2ssc2Fucy1zZXJpZiIgZm9udC1zaXplPSIxMTAiIGZvbnQtd2VpZ2h0PSI5MDAiIGZpbGw9IiNGMEE1MDAiIHRleHQtYW5jaG9yPSJtaWRkbGUiPkI8L3RleHQ+PC9zdmc+'

def get_sssihl_logo_b64():
    if os.path.exists('sssihl_logo.png'):

        with open('sssihl_logo.png','rb') as f:

            return base64.b64encode(f.read()).decode()

    return 'PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIyMDAiIGhlaWdodD0iMjAwIj48Y2lyY2xlIGN4PSIxMDAiIGN5PSIxMDAiIHI9Ijk2IiBmaWxsPSIjMUEzQTVDIiBzdHJva2U9IiNGRkZGRkYiIHN0cm9rZS13aWR0aD0iMyIvPjx0ZXh0IHg9IjEwMCIgeT0iMTA1IiBmb250LWZhbWlseT0iQXJpYWwsc2Fucy1zZXJpZiIgZm9udC1zaXplPSIyNiIgZm9udC13ZWlnaHQ9IjcwMCIgZmlsbD0iI0ZGRkZGRiIgdGV4dC1hbmNob3I9Im1pZGRsZSI+U1NTSUhMPC90ZXh0Pjx0ZXh0IHg9IjEwMCIgeT0iMTM1IiBmb250LWZhbWlseT0iQXJpYWwsc2Fucy1zZXJpZiIgZm9udC1zaXplPSIxMyIgZm9udC13ZWlnaHQ9IjQwMCIgZmlsbD0iI0FBQ0NFRSIgdGV4dC1hbmNob3I9Im1pZGRsZSI+UHJhc2FudGhpIE5pbGF5YW08L3RleHQ+PC9zdmc+'

COMPOUNDS = load_compounds()

import sqlite3 as _sqlite3

@st.cache_data
def load_natural_products_db():
    """Tier 2B: 400K+ COCONUT/PubChem natural products — SQLite fast lookup"""
    if not os.path.exists("brainsafe_natural.db"):
        return {}
    try:
        conn = _sqlite3.connect("brainsafe_natural.db")
        rows = conn.execute("SELECT * FROM compounds").fetchall()
        conn.close()
        cols = ["name","compound_type","bbb","antioxidant","anti_inflammatory",
                "mitochondrial_support","aggregation_modulation","cognitive_enhancement",
                "neurogenesis","synaptic_plasticity","nps","confidence","source","pubchem_cid"]
        return {r[0]: dict(zip(cols, r)) for r in rows}
    except Exception:
        return {}

NATURAL_DB = load_natural_products_db()

COMPOUND_NAMES = sorted(list(COMPOUNDS.keys()))
DB_STATS = get_db_stats()
LOGO_B64 = get_logo_b64()
SSSIHL_LOGO_B64 = get_sssihl_logo_b64()
ENZYME_DATA = load_enzyme_data()

_ML_SCORE_COLS = [
    "antioxidant", "anti_inflammatory", "mitochondrial_support",
    "aggregation_modulation", "cognitive_enhancement",
    "neurogenesis", "synaptic_plasticity",
]
_BBB_MAP  = {"Low": 0, "Low-Med": 1, "Medium": 2, "High": 3}
_DIS_MAP  = {"Low": 0, "Med": 1, "High": 2}
_DIS_KW: dict[str, list[str]] = {
    "alzheimers": ["alzheimer", "cholinesterase", "acetylcholinesterase", "bace", "amyloid", "donepezil", "memantine", "galantamine", "rivastigmine", "aricept"],
    "parkinsons":  ["parkinson", "dopaminerg", "levodopa", "alpha-synuclein", "mao-b", "selegiline", "rasagiline", "pramipexole"],
    "als":        ["amyotrophic", "motor neuron", "sod1", "tdp-43", "riluzole", "edaravone"],
    "huntingtons": ["huntington", "striatum", "hdac", "tetrabenazine"],
}

@st.cache_resource
def build_ml_predictor():
    """
    Train a MultiOutput RandomForest on the 129 curated compounds using disease-relevance features.
    Used at runtime to predict scores for user-typed compounds PubChem can identify.
    Cached once per session via @st.cache_resource (non-serializable object).
    """
    try:
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.multioutput import MultiOutputRegressor
        from sklearn.preprocessing import StandardScaler
        with open("compounds.json") as f:
            curated = json.load(f)
        X_rows, y_rows = [], []
        for entry in curated.values():
            _ct   = (entry.get('compound_type','') or '').lower()
            _poly = 1.0 if any(p_ in _ct for p_ in POLYPHENOL_TYPES) else 0.0
            _pt   = ' '.join(entry.get('pathways',[])).lower()
            _nn   = float(sum(1 for kw in NEURO_KWS if kw in _pt))
            _nm   = float(len(entry.get('metabolites',[])))
            X_rows.append([
                float(BBB_MAP.get(entry.get('bbb', 'Low'), 0)),
                float(DIS_MAP.get(entry.get('als', 'Low'), 0)),
                float(DIS_MAP.get(entry.get('alzheimers', 'Low'), 0)),
                float(DIS_MAP.get(entry.get('parkinsons', 'Low'), 0)),
                float(DIS_MAP.get(entry.get('huntingtons', 'Low'), 0)),
                float(len(entry.get('pathways', []))),
                _poly, _nn, _nm,      # v2.2: polyphenol, neuro_kws, n_mets
            ])
            y_rows.append([float(entry.get(c, 5.0)) for c in _ML_SCORE_COLS])
        X = np.array(X_rows, dtype=float)
        y = np.array(y_rows, dtype=float)
        scaler = StandardScaler()
        X_s = scaler.fit_transform(X)
        base = RandomForestRegressor(n_estimators=150, max_depth=6, min_samples_leaf=2, random_state=42, n_jobs=1)
        model = MultiOutputRegressor(base)
        model.fit(X_s, y)
        return model, scaler
    except Exception:
        return None, None


def predict_unknown_via_ml(compound_name: str) -> dict | None:
    """
    v3: PubChem SMILES + KNN NPS + disease relevance + pathway annotation.
    Returns None if PubChem cannot identify the compound.
    """
    model, scaler = build_ml_predictor()
    pc = fetch_pubchem(compound_name)
    if pc.get("error"):
        return None
    try:
        mw   = float(pc.get("mw",   999))
        logp = float(pc.get("xlogp", 5.0))  if pc.get("xlogp") not in (None, "N/A") else 5.0
        tpsa = float(pc.get("tpsa", 999))   if pc.get("tpsa")  not in (None, "N/A") else 999.0
        hbd  = float(pc.get("hbd",   99))   if pc.get("hbd")   not in (None, "N/A") else 99.0
    except (TypeError, ValueError):
        return None

    if   mw <= 360 and 1.0 <= logp <= 3.0 and tpsa <= 60 and hbd <= 1:
        bbb_str, bbb_num = "High",    3
    elif mw <= 450 and 0.0 <= logp <= 4.0 and tpsa <= 90 and hbd <= 3:
        bbb_str, bbb_num = "Medium",  2
    elif mw <= 500 and tpsa <= 120:
        bbb_str, bbb_num = "Low-Med", 1
    else:
        bbb_str, bbb_num = "Low",     0

    v3_nps, v3_neighbours, v3_max_sim, v3_confidence = None, [], 0.0, "Low"
    v3_diseases, v3_pathways = {}, []
    smiles = (pc.get("isomeric_smiles") or pc.get("smiles") or "")
    if V3_KNN_AVAILABLE and smiles:
        try:
            v3_nps, v3_neighbours, v3_max_sim, v3_confidence = predict_nps_knn(smiles)
            v3_targets  = get_chembl_targets(compound_name)
            v3_diseases = predict_disease_relevance(compound_name, smiles, v3_targets)
            v3_pathways = predict_pathways(compound_name, v3_targets)
        except Exception:
            pass

    scores = {col: 5.0 for col in _ML_SCORE_COLS}
    ch = fetch_chembl(compound_name)
    mech_text = " ".join(
        (m.get("target", "") + " " + m.get("action", "")).lower()
        for m in (ch.get("mechanisms") or [])
    )
    dis_nums: dict[str, int] = {}
    for dis, kws in _DIS_KW.items():
        dis_nums[dis] = 2 if any(kw in mech_text for kw in kws) else 0
    if model is not None and scaler is not None:
        try:
            _mt   = " ".join(m.get("target","") or "" for m in ch.get("mechanisms",[]))
            _poly = 1.0 if any(p_ in compound_name.lower() for p_ in POLYPHENOL_TYPES) else 0.0
            _nn   = float(sum(1 for kw in NEURO_KWS if kw in _mt.lower()))
            _nm   = float(min(len(ch.get("mechanisms",[])), 8))
            feat  = np.array([[
                float(bbb_num),
                float(dis_nums.get("als", 0)),
                float(dis_nums.get("alzheimers", 0)),
                float(dis_nums.get("parkinsons", 0)),
                float(dis_nums.get("huntingtons", 0)),
                float(max(2, min(sum(
                    1 for _m in ch.get("mechanisms", [])
                    for _kw in TARGETTOPATHWAY
                    if _kw in (_m.get("target", "") or "").lower()
                ), 8))),
                _poly, _nn, _nm,
            ]], dtype=float)
            preds  = model.predict(scaler.transform(feat))[0]
            scores = {col: round(float(np.clip(v, 1.0, 10.0)), 1)
                      for col, v in zip(_ML_SCORE_COLS, preds)}
        except Exception:
            pass

    dis_levels = v3_diseases if v3_diseases else {
        d: ("High" if v == 2 else "Low") for d, v in dis_nums.items()
    }
    nps_final = v3_nps if v3_nps is not None else round(
        float(np.mean([scores.get(c, 5.0) for c in _ML_SCORE_COLS])) * 10, 1)

    return {
        "compound_type":    ch.get("molecule_type", "Unknown"),
        "bbb":              bbb_str,
        **scores,
        **dis_levels,
        "pathways":         v3_pathways if v3_pathways else ["NF-kB", "Nrf2/GSH"],
        "metabolites":      ["ROS", "IL-6"],
        "brain_regions":    ["Cortex", "Hippocampus"],
        "data_source":      "pubchem_ml_predicted",
        "confidence":       v3_confidence.lower() if v3_confidence else "low",
        "ml_predicted":     True,
        "chembl_id":        ch.get("chembl_id", ""),
        "chembl_max_phase": ch.get("max_phase"),
        "indication_diseases": [d for d, v in dis_nums.items() if v > 0],
        "model_cv_r2":      0.20,
        "pubchem_cid":      pc.get("cid", ""),
        "mw": mw, "xlogp": logp, "tpsa": tpsa,
        "_live_ml":         True,
        "nps":              nps_final,
        "nps_knn":          v3_nps,
        "nps_confidence":   v3_confidence,
        "nps_neighbours":   v3_neighbours,
        "nps_max_sim":      round(v3_max_sim, 3),
        "v3_diseases":      v3_diseases,
        "v3_pathways":      v3_pathways,
        "smiles":           smiles,
    }


ENZYME_ALIASES = {
    "N-Acetylcysteine": "NAC",
    "Omega-3 (DHA/EPA)": "DHA",
    "Coenzyme Q10": "Coenzyme Q10",
    "Vitamin D3": "Vitamin D",
    "Nicotinamide Riboside": "Nicotinamide",
    "Alpha Lipoic Acid": "Alpha-Lipoic Acid",
    "α-Lipoic Acid": "Alpha-Lipoic Acid",
    "R-Lipoic Acid": "Alpha-Lipoic Acid",
    "EGCG (Epigallocatechin Gallate)": "EGCG",
    "Green Tea Catechin": "EGCG",
    "Bacopa": "Bacopa Monnieri",
    "Lion's Mane Mushroom": "Lion's Mane",
    "Hericium erinaceus": "Lion's Mane",
    "Pterostilbene (trans)": "Pterostilbene",
    "CoQ10": "Coenzyme Q10",
    "Ubiquinol": "Coenzyme Q10",
    "DHA (Docosahexaenoic Acid)": "DHA",
    "Methylcobalamin": "Vitamin B12",
    "Cyanocobalamin": "Vitamin B12",
}

_ESTI_NOTE = "Class-based estimate — consult BRENDA or PubMed for compound-specific kinetic data."

CLASS_ENZYME_TEMPLATES = {
    "flavonoid": {
        "is_estimated": True,
        "enzymes": [
            {"name": "COX-2", "action": "Inhibition", "strength": "Moderate", "ec": "EC 1.14.99.1",
             "note": "Anti-inflammatory prostaglandin suppression; typical for polyphenolic flavonoids. " + _ESTI_NOTE, "source": "BRENDA; PMID:22178138"},
            {"name": "MAO-B", "action": "Inhibition", "strength": "Moderate", "ec": "EC 1.14.14.1",
             "note": "Reduces dopamine catabolism; common flavonoid class activity. " + _ESTI_NOTE, "source": "BRENDA; PMID:25975702"},
            {"name": "Nrf2", "action": "Activation", "strength": "Moderate", "ec": "Transcription factor",
             "note": "Antioxidant response element (ARE) induction; upregulates HO-1, NQO1, GSH synthesis. " + _ESTI_NOTE, "source": "MetaCyc; PMID:22178138"},
            {"name": "AChE", "action": "Inhibition", "strength": "Weak", "ec": "EC 3.1.1.7",
             "note": "Variable across flavonoid subclasses; typically weak. " + _ESTI_NOTE, "source": "PMID:28695743"},
        ],
        "transporters": [
            {"name": "P-glycoprotein / ABCB1", "role": "Substrate / Inhibitor", "location": "BBB efflux (luminal)",
             "note": "Common for polyphenolic flavonoids; may limit CNS accumulation. " + _ESTI_NOTE, "source": "PMID:14642381"},
            {"name": "BCRP / ABCG2", "role": "Substrate", "location": "BBB efflux (luminal)",
             "note": "Frequent efflux substrate; class-based estimate. " + _ESTI_NOTE, "source": "PMID:15901794"},
        ],
        "cofactors": ["NADPH", "GSH", "Mg²⁺"],
        "sources": ["BRENDA (flavonoid class)", "MetaCyc", "PMID:22178138", "PMID:28695743"],
    },
    "vitamin_b_thiamine": {
        "is_estimated": True,
        "enzymes": [
            {"name": "Pyruvate dehydrogenase complex (PDC)", "action": "Cofactor", "strength": "Strong", "ec": "EC 1.2.4.1",
             "note": "Thiamine pyrophosphate (ThPP) is covalently bound cofactor of E1 subunit; converts pyruvate to acetyl-CoA. Critical in brain energy metabolism.", "source": "MetaCyc PWY-5481; BRENDA"},
            {"name": "α-Ketoglutarate dehydrogenase", "action": "Cofactor", "strength": "Strong", "ec": "EC 1.2.4.2",
             "note": "ThPP essential in TCA cycle; deficiency causes Wernicke encephalopathy and neuronal death.", "source": "MetaCyc; BRENDA"},
            {"name": "Transketolase (TK)", "action": "Cofactor", "strength": "Strong", "ec": "EC 2.2.1.1",
             "note": "ThPP-dependent enzyme of pentose phosphate pathway; maintains NADPH for antioxidant defence.", "source": "BRENDA"},
        ],
        "transporters": [
            {"name": "ThT1 / SLC19A2", "role": "Substrate", "location": "BBB and intestinal epithelium",
             "note": "High-affinity thiamine transporter; mutations cause thiamine-responsive megaloblastic anaemia.", "source": "UniProt O60779"},
            {"name": "ThT2 / SLC19A3", "role": "Substrate", "location": "BBB and neuronal membranes",
             "note": "Second thiamine transporter; expressed in brain; mutations cause biotin-thiamine-responsive basal ganglia disease.", "source": "UniProt Q9BZV2"},
        ],
        "cofactors": ["ThPP (thiamine pyrophosphate — active form)", "Mg²⁺ (enzyme activation)", "ATP"],
        "sources": ["MetaCyc PWY-5481", "BRENDA EC 2.2.1.1", "UniProt O60779", "DrugBank DB00648"],
    },
    "vitamin_b_riboflavin": {
        "is_estimated": True,
        "enzymes": [
            {"name": "Riboflavin kinase", "action": "Substrate", "strength": "Strong", "ec": "EC 2.7.1.26",
             "note": "Phosphorylates riboflavin to FMN; first step of FAD biosynthesis.", "source": "BRENDA; MetaCyc"},
            {"name": "FMN adenylyltransferase (FADS)", "action": "Substrate", "strength": "Strong", "ec": "EC 2.7.7.2",
             "note": "Converts FMN to FAD; FAD serves as cofactor for Complex I and many oxidoreductases.", "source": "MetaCyc; BRENDA"},
            {"name": "Mitochondrial Complex I / II", "action": "Cofactor (as FAD/FMN)", "strength": "Strong", "ec": "EC 7.1.1.2",
             "note": "FAD and FMN are essential electron carriers in the mitochondrial electron transport chain.", "source": "MetaCyc PWY-3781"},
        ],
        "transporters": [
            {"name": "RFVT1 / SLC52A1", "role": "Substrate", "location": "Intestinal absorption",
             "note": "Riboflavin transporter; mutations cause Brown-Vialetto-Van Laere syndrome.", "source": "UniProt Q9HAB3"},
            {"name": "RFVT3 / SLC52A3", "role": "Substrate", "location": "BBB and brain",
             "note": "Primary riboflavin transporter in CNS; expressed in neurons and choroid plexus.", "source": "UniProt Q9Y2Y1"},
        ],
        "cofactors": ["FMN (flavin mononucleotide — active form)", "FAD (flavin adenine dinucleotide — active form)", "ATP"],
        "sources": ["MetaCyc", "BRENDA EC 2.7.1.26", "UniProt Q9HAB3", "DrugBank DB00140"],
    },
    "vitamin_b_niacin": {
        "is_estimated": True,
        "enzymes": [
            {"name": "NAPRT (Nicotinate phosphoribosyltransferase)", "action": "Substrate", "strength": "Strong", "ec": "EC 2.4.2.11",
             "note": "Converts niacin (nicotinic acid) to NaMN → NMN → NAD+ via Preiss-Handler pathway.", "source": "MetaCyc PWY-3981; BRENDA"},
            {"name": "GPR109A (niacin receptor)", "action": "Agonist", "strength": "Strong", "ec": "GPCR / Receptor",
             "note": "Nicotinic acid activates GPR109A; mediates vasodilation (flushing) and lipid-modifying effects.", "source": "UniProt Q8TDS4"},
            {"name": "SIRT1 / SIRT3 (indirect, via NAD+)", "action": "Activation", "strength": "Moderate", "ec": "EC 3.5.1.98",
             "note": "As NAD+ precursor, niacin elevates SIRT activity; promotes neuroprotection and mitochondrial quality control.", "source": "PMID:27524370"},
        ],
        "transporters": [
            {"name": "SLC5A8 (SMCT1)", "role": "Substrate", "location": "Intestinal and CNS transport",
             "note": "Sodium-coupled monocarboxylate transporter; facilitates niacin uptake across gut and possibly BBB.", "source": "PMID:16091360"},
            {"name": "Monocarboxylate transporter (MCT)", "role": "Substrate", "location": "BBB",
             "note": "Niacin can enter the brain via MCT-mediated transport due to its carboxylate structure.", "source": "Literature"},
        ],
        "cofactors": ["NAD+", "NADP+", "PRPP (phosphoribosyl pyrophosphate)"],
        "sources": ["MetaCyc PWY-3981", "BRENDA EC 2.4.2.11", "PMID:27524370", "DrugBank DB00627"],
    },
    "vitamin_b5": {
        "is_estimated": True,
        "enzymes": [
            {"name": "Pantothenate kinase (PANK)", "action": "Substrate", "strength": "Strong", "ec": "EC 2.7.1.33",
             "note": "Rate-limiting step converting pantothenate to 4'-phosphopantothenate → CoA. PANK2 mutations cause NBIA (neurodegeneration with brain iron accumulation).", "source": "BRENDA; MetaCyc"},
            {"name": "CoA-dependent enzymes (TCA cycle, fatty acid metabolism)", "action": "Cofactor (as CoA)", "strength": "Strong", "ec": "Multiple EC classes",
             "note": "CoA is the active form; essential for acetyl-CoA production, fatty acid synthesis and beta-oxidation.", "source": "MetaCyc"},
        ],
        "transporters": [
            {"name": "SMVT / SLC5A6", "role": "Substrate", "location": "Intestinal absorption and BBB",
             "note": "Sodium-dependent multivitamin transporter; co-transports pantothenate with biotin and lipoate.", "source": "UniProt Q9Y289"},
        ],
        "cofactors": ["CoA (active form)", "Acetyl-CoA", "ATP"],
        "sources": ["MetaCyc", "BRENDA EC 2.7.1.33", "UniProt Q9Y289", "DrugBank DB00205"],
    },
    "vitamin_b6": {
        "is_estimated": True,
        "enzymes": [
            {"name": "Pyridoxal kinase (PDXK)", "action": "Substrate", "strength": "Strong", "ec": "EC 2.7.1.35",
             "note": "Activates pyridoxal to PLP (pyridoxal 5'-phosphate); active coenzyme form.", "source": "BRENDA; MetaCyc"},
            {"name": "Aminotransferases (ALT, AST, GABA-T)", "action": "Cofactor (as PLP)", "strength": "Strong", "ec": "EC 2.6.1.x",
             "note": "PLP is Schiff base cofactor for >150 enzymatic reactions; includes GABA synthesis (GAD), serotonin and dopamine decarboxylation.", "source": "MetaCyc; BRENDA"},
            {"name": "DOPA decarboxylase (AADC)", "action": "Cofactor (as PLP)", "strength": "Strong", "ec": "EC 4.1.1.28",
             "note": "PLP-dependent enzyme converting L-DOPA to dopamine and 5-HTP to serotonin; deficiency causes movement disorders.", "source": "BRENDA; PMID:9427359"},
        ],
        "transporters": [
            {"name": "Pyridoxine transporter (SLC25A39)", "role": "Substrate", "location": "Mitochondrial transport",
             "note": "Mitochondrial carrier for PLP, supporting intramitochondrial PLP-dependent reactions.", "source": "UniProt Q9Y601"},
            {"name": "Carrier-mediated transport", "role": "Substrate", "location": "BBB",
             "note": "Pyridoxal crosses BBB preferentially; phosphorylated in brain endothelium.", "source": "Literature"},
        ],
        "cofactors": ["PLP (pyridoxal 5'-phosphate — active form)", "ATP (phosphorylation)"],
        "sources": ["MetaCyc", "BRENDA EC 2.7.1.35", "PMID:9427359", "DrugBank DB00165"],
    },
    "vitamin_b9": {
        "is_estimated": True,
        "enzymes": [
            {"name": "DHFR (Dihydrofolate reductase)", "action": "Substrate", "strength": "Strong", "ec": "EC 1.5.1.3",
             "note": "Converts dietary folate to THF (tetrahydrofolate); inhibited by methotrexate.", "source": "BRENDA; MetaCyc"},
            {"name": "MTHFR (Methylenetetrahydrofolate reductase)", "action": "Substrate", "strength": "Strong", "ec": "EC 1.5.1.20",
             "note": "Converts 5,10-MTHF to 5-MTHF (methyl donor for methionine synthase). MTHFR C677T polymorphism increases dementia risk.", "source": "MetaCyc; PMID:12519925"},
            {"name": "Methionine synthase (MS)", "action": "Cofactor (as 5-MTHF)", "strength": "Strong", "ec": "EC 2.1.1.13",
             "note": "5-MTHF donates methyl group to convert homocysteine to methionine; links one-carbon metabolism to methylation and B12 pathway.", "source": "MetaCyc PWY-3221"},
        ],
        "transporters": [
            {"name": "PCFT / SLC46A1", "role": "Substrate", "location": "Intestinal absorption and choroid plexus",
             "note": "Proton-coupled folate transporter; primary route for dietary folate absorption and CNS folate delivery.", "source": "UniProt Q96NT5"},
            {"name": "RFC1 / SLC19A1", "role": "Substrate", "location": "BBB and cellular uptake",
             "note": "Reduced folate carrier; ubiquitous transporter for reduced folates and antifolate drugs.", "source": "UniProt P41440"},
            {"name": "Folate receptor alpha (FOLR1)", "role": "Receptor-mediated endocytosis", "location": "Choroid plexus (CNS)",
             "note": "High-affinity folate receptor; critical for CSF folate delivery; defects cause cerebral folate deficiency.", "source": "UniProt P15328"},
        ],
        "cofactors": ["NADPH (DHFR reaction)", "SAM (methylation product)", "Vitamin B12 (interdependent)"],
        "sources": ["MetaCyc PWY-3221", "BRENDA EC 1.5.1.3", "UniProt Q96NT5", "PMID:12519925", "DrugBank DB00158"],
    },
    "vitamin_c": {
        "is_estimated": True,
        "enzymes": [
            {"name": "Prolyl hydroxylases (PHD1/2/3)", "action": "Cofactor", "strength": "Strong", "ec": "EC 1.14.11.2",
             "note": "Vitamin C is obligate cofactor for HIF prolyl hydroxylases; controls HIF-1α degradation and oxygen sensing in neurons.", "source": "MetaCyc; BRENDA"},
            {"name": "Dopamine beta-hydroxylase (DBH)", "action": "Cofactor", "strength": "Strong", "ec": "EC 1.14.17.1",
             "note": "Ascorbate-dependent conversion of dopamine to norepinephrine in catecholaminergic neurons; deficiency causes autonomic failure.", "source": "BRENDA; MetaCyc"},
            {"name": "Glutathione reductase (GR)", "action": "Recycling partner", "strength": "Moderate", "ec": "EC 1.8.1.7",
             "note": "Ascorbate regenerates oxidised glutathione (GSSG → GSH) and vitamin E via the antioxidant network.", "source": "MetaCyc"},
            {"name": "Nrf2", "action": "Activation", "strength": "Moderate", "ec": "Transcription factor",
             "note": "High-dose ascorbate activates Nrf2/ARE pathway; upregulates antioxidant enzyme expression.", "source": "PMID:22178138"},
        ],
        "transporters": [
            {"name": "SVCT2 / SLC23A2", "role": "Substrate", "location": "BBB and neuronal membranes",
             "note": "Sodium-dependent vitamin C transporter 2; primary route for ascorbate uptake in neurons; CNS-specific.", "source": "UniProt I6L9M3"},
            {"name": "GLUT1 (for DHAA form)", "role": "Substrate", "location": "BBB",
             "note": "Dehydroascorbic acid (DHAA, oxidised form) crosses BBB via GLUT1; reduced back to ascorbate inside cells.", "source": "PMID:9878065"},
        ],
        "cofactors": ["Fe²⁺ (enzyme activation)", "Cu²⁺ (DBH cofactor)", "O₂", "Dehydroascorbate (redox couple)"],
        "sources": ["MetaCyc", "BRENDA EC 1.14.17.1", "UniProt I6L9M3", "PMID:9878065", "DrugBank DB00126"],
    },
    "vitamin_e": {
        "is_estimated": True,
        "enzymes": [
            {"name": "COX-2", "action": "Inhibition", "strength": "Moderate", "ec": "EC 1.14.99.1",
             "note": "Tocotrienols show stronger COX-2 inhibition than tocopherols; reduces neuroinflammatory prostaglandin synthesis.", "source": "BRENDA; PMID:23537148"},
            {"name": "5-LOX", "action": "Inhibition", "strength": "Moderate", "ec": "EC 1.13.11.34",
             "note": "Reduces leukotriene production; anti-neuroinflammatory mechanism.", "source": "BRENDA"},
            {"name": "12-LOX (neuronal isoform)", "action": "Inhibition", "strength": "Moderate", "ec": "EC 1.13.11.31",
             "note": "Tocotrienols specifically inhibit glutamate-induced 12-LOX activation in neurons; prevents oxidative neuronal death.", "source": "PMID:23537148"},
            {"name": "PKC (Protein Kinase C)", "action": "Activation", "strength": "Moderate", "ec": "EC 2.7.11.13",
             "note": "Alpha-tocopherol activates PKC in smooth muscle and neuronal cells; influences cell signalling.", "source": "PMID:9109418"},
        ],
        "transporters": [
            {"name": "ABCA1 (ABC transporter)", "role": "Substrate / Regulatory", "location": "Cellular efflux and BBB",
             "note": "Mediates vitamin E and cholesterol efflux; modulates neuronal membrane composition.", "source": "UniProt O95477"},
            {"name": "SR-B1 (Scavenger receptor B1)", "role": "Receptor-mediated uptake", "location": "BBB, liver",
             "note": "Facilitates selective uptake of alpha-tocopherol from HDL particles.", "source": "UniProt Q8WTV0"},
        ],
        "cofactors": ["CoQ10 / Ubiquinol (recycling partner — regenerates tocopheroxyl radical)", "GSH", "Ascorbate (vitamin C, recycling)"],
        "sources": ["BRENDA", "MetaCyc", "PMID:23537148", "DrugBank DB00163"],
    },
    "carotenoid": {
        "is_estimated": True,
        "enzymes": [
            {"name": "BCMO1 (Beta-carotene monooxygenase)", "action": "Substrate (cleavage)", "strength": "Strong", "ec": "EC 1.14.99.36",
             "note": "Symmetrically cleaves beta-carotene to retinal (vitamin A). Does not act on xanthophyll carotenoids (lutein, zeaxanthin).", "source": "BRENDA; MetaCyc"},
            {"name": "COX-2", "action": "Inhibition", "strength": "Moderate", "ec": "EC 1.14.99.1",
             "note": "Carotenoids reduce neuroinflammatory eicosanoid production; particularly astaxanthin and lycopene.", "source": "BRENDA; PMID:23537148"},
            {"name": "NF-κB (IKK)", "action": "Inhibition", "strength": "Moderate", "ec": "EC 2.7.11.10",
             "note": "Anti-inflammatory transcription factor suppression; lycopene and astaxanthin well-documented.", "source": "Literature"},
            {"name": "Nrf2", "action": "Activation", "strength": "Moderate", "ec": "Transcription factor",
             "note": "Antioxidant gene upregulation via ARE; astaxanthin is the most potent carotenoid Nrf2 activator.", "source": "PMID:22178138"},
        ],
        "transporters": [
            {"name": "SR-B1 (Scavenger receptor B1)", "role": "Receptor-mediated uptake", "location": "Intestine, liver, BBB",
             "note": "Micelle-associated carotenoid uptake by SR-B1 and CD36 at intestinal brush border.", "source": "UniProt Q8WTV0"},
            {"name": "CD36 (Scavenger receptor)", "role": "Uptake facilitator", "location": "Intestinal epithelium",
             "note": "Facilitates absorption of carotenoid-containing lipid micelles.", "source": "UniProt P16671"},
        ],
        "cofactors": ["Retinol (vitamin A — product of provitamin A carotenoids)", "O₂ (cleavage)", "Fe²⁺ (BCMO1 activity)"],
        "sources": ["BRENDA EC 1.14.99.36", "MetaCyc", "UniProt Q8WTV0", "PMID:23537148"],
    },
    "amino_acid": {
        "is_estimated": True,
        "enzymes": [
            {"name": "Glutamate-cysteine ligase (GCL)", "action": "Substrate / Cofactor", "strength": "Strong", "ec": "EC 6.3.2.2",
             "note": "Rate-limiting GSH synthesis enzyme; amino acid substrates (cysteine, glutamate, glycine) are essential for neuronal antioxidant defence.", "source": "MetaCyc PWY-7295; BRENDA"},
            {"name": "Aminotransferases (general)", "action": "Substrate", "strength": "Strong", "ec": "EC 2.6.1.x",
             "note": "PLP-dependent transamination; amino acid pool is substrate for neurotransmitter and metabolite synthesis.", "source": "MetaCyc"},
            {"name": "NF-κB (IKK)", "action": "Inhibition", "strength": "Moderate", "ec": "EC 2.7.11.10",
             "note": "Many amino acid derivatives reduce neuroinflammation via NF-κB pathway modulation; class-based estimate. " + _ESTI_NOTE, "source": "Literature"},
        ],
        "transporters": [
            {"name": "System L amino acid transporter (LAT1/SLC7A5)", "role": "Substrate", "location": "BBB",
             "note": "Large neutral amino acid transporter; primary route for brain entry of aromatic and large amino acids.", "source": "UniProt Q01650"},
            {"name": "xCT / SLC7A11 (for cysteine derivatives)", "role": "Substrate / Modulator", "location": "Neuronal membrane",
             "note": "Cystine/glutamate antiporter; mediates cysteine uptake for GSH synthesis; modulates extracellular glutamate.", "source": "UniProt Q9UPY5"},
        ],
        "cofactors": ["ATP (activation)", "GSH (synthesis product)", "PLP (transamination cofactor)"],
        "sources": ["MetaCyc PWY-7295", "BRENDA EC 6.3.2.2", "UniProt Q01650", "Literature"],
    },
    "alkaloid": {
        "is_estimated": True,
        "enzymes": [
            {"name": "MAO-A / MAO-B", "action": "Inhibition", "strength": "Moderate", "ec": "EC 1.14.14.1",
             "note": "Many plant alkaloids inhibit monoamine oxidase; increases dopamine, serotonin, and norepinephrine availability. " + _ESTI_NOTE, "source": "BRENDA; PMID:25975702"},
            {"name": "AChE", "action": "Inhibition", "strength": "Weak", "ec": "EC 3.1.1.7",
             "note": "Alkaloid class often shows cholinesterase inhibition; strength varies widely. " + _ESTI_NOTE, "source": "PMID:28695743"},
            {"name": "CYP450 enzymes (CYP1A2, CYP3A4)", "action": "Substrate / Inhibitor", "strength": "Moderate", "ec": "EC 1.14.14.1",
             "note": "Alkaloids are frequently metabolised by hepatic CYP enzymes; drug interaction potential. " + _ESTI_NOTE, "source": "BRENDA"},
        ],
        "transporters": [
            {"name": "P-glycoprotein / ABCB1", "role": "Substrate", "location": "BBB efflux",
             "note": "Alkaloids are frequently P-gp substrates, limiting CNS bioavailability. " + _ESTI_NOTE, "source": "Literature"},
            {"name": "Organic cation transporters (OCT1/2)", "role": "Substrate", "location": "Hepatic and renal",
             "note": "Many cationic alkaloids are OCT substrates affecting distribution and elimination.", "source": "Literature"},
        ],
        "cofactors": ["FAD (MAO cofactor, target of interaction)", "NAD+"],
        "sources": ["BRENDA", "MetaCyc", "PMID:25975702", "Literature"],
    },
    "fatty_acid": {
        "is_estimated": True,
        "enzymes": [
            {"name": "COX-1 / COX-2", "action": "Substrate (competitive)", "strength": "Strong", "ec": "EC 1.14.99.1",
             "note": "Polyunsaturated fatty acids compete with arachidonic acid at COX; generate anti-inflammatory lipid mediators (resolvins, protectins). " + _ESTI_NOTE, "source": "MetaCyc; BRENDA"},
            {"name": "Phospholipase A2 (cPLA2)", "action": "Substrate", "strength": "Strong", "ec": "EC 3.1.1.4",
             "note": "PLA2 liberates fatty acids from membrane phospholipids on neuronal activation, initiating lipid mediator synthesis.", "source": "MetaCyc; PMID:22021740"},
            {"name": "Acyl-CoA synthetase (ACSL)", "action": "Substrate (activation)", "strength": "Strong", "ec": "EC 6.2.1.3",
             "note": "Activates fatty acids to acyl-CoA for phospholipid incorporation and beta-oxidation.", "source": "BRENDA"},
        ],
        "transporters": [
            {"name": "FATP4 / SLC27A4", "role": "Substrate", "location": "BBB and intestine",
             "note": "Fatty acid transport protein; facilitates brain uptake of long-chain fatty acids.", "source": "PMID:19497945"},
            {"name": "Serum albumin", "role": "Binding / transport protein", "location": "Blood plasma",
             "note": "Non-esterified fatty acids bind albumin for systemic transport.", "source": "HMDB"},
        ],
        "cofactors": ["CoA (activation to acyl-CoA)", "O₂ (desaturation reactions)", "NADPH"],
        "sources": ["MetaCyc", "BRENDA EC 6.2.1.3", "PMID:22021740", "HMDB"],
    },
    "mineral": {
        "is_estimated": True,
        "enzymes": [
            {"name": "Superoxide dismutase (SOD1/2)", "action": "Cofactor", "strength": "Strong", "ec": "EC 1.15.1.1",
             "note": "Zinc, copper (SOD1) and manganese (SOD2) are essential catalytic cofactors; SOD1 mutations cause familial ALS.", "source": "MetaCyc; BRENDA"},
            {"name": "Glutathione peroxidase (GPx)", "action": "Cofactor (Selenium)", "strength": "Strong", "ec": "EC 1.11.1.9",
             "note": "Selenium is the active site selenocysteine in GPx; essential for neuronal antioxidant defence.", "source": "MetaCyc; BRENDA"},
            {"name": "Carbonic anhydrase (CA) / Zinc-dependent enzymes", "action": "Cofactor", "strength": "Strong", "ec": "EC 4.2.1.1",
             "note": "Zinc is essential cofactor for ~300 human enzymes including carboxypeptidase, MMP, HDAC, and carbonic anhydrase.", "source": "BRENDA"},
        ],
        "transporters": [
            {"name": "ZIP / ZnT family transporters (SLC39/SLC30)", "role": "Substrate", "location": "Cellular and BBB",
             "note": "ZIP transporters import zinc into cells; ZnT exporters maintain cytoplasmic zinc homeostasis; important in neurons.", "source": "UniProt"},
            {"name": "DMT1 / SLC11A2 (iron)", "role": "Substrate", "location": "BBB and cellular",
             "note": "Divalent metal transporter 1; primary route for non-transferrin-bound iron and other divalent metals into neurons.", "source": "UniProt P49281"},
        ],
        "cofactors": ["Specific ionic form of the mineral (e.g. Zn²⁺, Mg²⁺, Se²⁻)", "ATP (transporter energy)"],
        "sources": ["MetaCyc", "BRENDA", "UniProt", "HMDB"],
    },
    "terpenoid": {
        "is_estimated": True,
        "enzymes": [
            {"name": "NF-κB (IKK)", "action": "Inhibition", "strength": "Moderate", "ec": "EC 2.7.11.10",
             "note": "Anti-inflammatory mechanism common to many terpenoids and adaptogenic herbs; reduces microglial activation. " + _ESTI_NOTE, "source": "MetaCyc; Literature"},
            {"name": "Nrf2", "action": "Activation", "strength": "Moderate", "ec": "Transcription factor",
             "note": "Antioxidant response element induction; class-typical for terpenoids and plant polyphenols. " + _ESTI_NOTE, "source": "PMID:22178138"},
            {"name": "AChE", "action": "Inhibition", "strength": "Weak", "ec": "EC 3.1.1.7",
             "note": "Mild cholinesterase inhibition reported for some terpenoid-rich extracts. " + _ESTI_NOTE, "source": "Literature"},
        ],
        "transporters": [
            {"name": "P-glycoprotein / ABCB1", "role": "Variable (compound-dependent)", "location": "BBB efflux",
             "note": "Lipophilic terpenoids may cross BBB by passive diffusion; P-gp interaction varies. " + _ESTI_NOTE, "source": "Literature"},
            {"name": "Serum albumin", "role": "Binding / transport protein", "location": "Blood plasma",
             "note": "Lipophilic terpenoid aglycones typically bind albumin for systemic transport.", "source": "Literature"},
        ],
        "cofactors": ["NADPH (terpenoid biosynthesis)", "ATP", "Mg²⁺"],
        "sources": ["MetaCyc", "BRENDA", "PMID:22178138", "Literature"],
    },
    "polyphenol": {
        "is_estimated": True,
        "enzymes": [
            {"name": "COX-1 / COX-2", "action": "Inhibition", "strength": "Moderate", "ec": "EC 1.14.99.1",
             "note": "Anti-inflammatory prostaglandin suppression; common polyphenol class activity. " + _ESTI_NOTE, "source": "BRENDA; Literature"},
            {"name": "NF-κB (IKK)", "action": "Inhibition", "strength": "Moderate", "ec": "EC 2.7.11.10",
             "note": "Reduces neuroinflammatory cytokine transcription; well-established for polyphenols. " + _ESTI_NOTE, "source": "MetaCyc; PMID:22178138"},
            {"name": "Nrf2", "action": "Activation", "strength": "Moderate", "ec": "Transcription factor",
             "note": "Antioxidant response element induction; upregulates HO-1, NQO1, GSH synthesis. " + _ESTI_NOTE, "source": "MetaCyc; PMID:22178138"},
        ],
        "transporters": [
            {"name": "P-glycoprotein / ABCB1", "role": "Substrate / Inhibitor", "location": "BBB efflux",
             "note": "Polyphenols are often P-gp substrates and may inhibit P-gp efflux of other compounds. " + _ESTI_NOTE, "source": "PMID:14642381"},
            {"name": "MRP2 / ABCC2", "role": "Inhibitor", "location": "Hepatic and intestinal efflux",
             "note": "Some polyphenols modulate MRP2 affecting systemic bioavailability. " + _ESTI_NOTE, "source": "PMID:16081031"},
        ],
        "cofactors": ["NADPH", "GSH", "Fe²⁺ / Cu²⁺ (redox chelation)"],
        "sources": ["BRENDA", "MetaCyc", "PMID:22178138", "Literature"],
    },
    "organosulfur": {
        "is_estimated": True,
        "enzymes": [
            {"name": "Nrf2 / Keap1", "action": "Activation (via Keap1 cysteine alkylation)", "strength": "Strong", "ec": "Transcription regulator",
             "note": "Organosulfur compounds covalently modify Keap1 cysteines; releases Nrf2 for ARE-driven antioxidant gene expression.", "source": "MetaCyc; PMID:22218793"},
            {"name": "COX-2", "action": "Inhibition", "strength": "Moderate", "ec": "EC 1.14.99.1",
             "note": "Organosulfur compounds from Allium species and cruciferous vegetables suppress COX-2 expression via NF-κB.", "source": "BRENDA; Literature"},
            {"name": "HDACs (Histone deacetylases)", "action": "Inhibition", "strength": "Moderate", "ec": "EC 3.5.1.98",
             "note": "Epigenetic derepression of neuroprotective and tumour-suppressor genes.", "source": "PMID:15033567"},
        ],
        "transporters": [
            {"name": "P-glycoprotein / ABCB1", "role": "Weak Substrate", "location": "BBB efflux",
             "note": "Organosulfur-GSH conjugates may be exported by MRP/ABCC transporters.", "source": "Literature"},
        ],
        "cofactors": ["GSH (conjugation substrate)", "Cysteine thiols (target of electrophilic organosulfur compounds)", "NADPH"],
        "sources": ["MetaCyc", "PMID:22218793", "PMID:15033567", "BRENDA"],
    },
    "phospholipid": {
        "is_estimated": True,
        "enzymes": [
            {"name": "Phospholipase A2 (PLA2)", "action": "Substrate", "strength": "Strong", "ec": "EC 3.1.1.4",
             "note": "PLA2 hydrolyses the sn-2 fatty acyl chain of phospholipids; releases DHA, arachidonate and other bioactive lipids.", "source": "MetaCyc; BRENDA"},
            {"name": "Choline acetyltransferase (ChAT)", "action": "Cofactor (via choline supply)", "strength": "Moderate", "ec": "EC 2.3.1.6",
             "note": "Phosphatidylcholine provides choline for ACh synthesis; PC-DHA supports membrane composition of cholinergic neurons.", "source": "MetaCyc"},
            {"name": "Cytidine diphosphocholine pathway (CDP-choline)", "action": "Substrate", "strength": "Strong", "ec": "EC 2.7.7.15",
             "note": "Phospholipid biosynthesis via Kennedy pathway; required for neuronal membrane maintenance and synaptic function.", "source": "MetaCyc"},
        ],
        "transporters": [
            {"name": "Mfsd2a (Major Facilitator Superfamily Domain 2a)", "role": "Substrate (LPC form)", "location": "BBB",
             "note": "Primary LPC-phospholipid transporter at BBB; critical for brain DHA homeostasis via LPC-DHA transport.", "source": "PMID:24828044"},
            {"name": "ABCA1 (ABC transporter)", "role": "Substrate", "location": "Cellular efflux",
             "note": "Mediates phospholipid efflux; involved in cholesterol and lipid homeostasis in neurons.", "source": "UniProt O95477"},
        ],
        "cofactors": ["CDP-choline (activation)", "CoA (acylation)", "CTP (cytidine triphosphate)"],
        "sources": ["MetaCyc", "BRENDA EC 3.1.1.4", "PMID:24828044", "UniProt O95477"],
    },
    "drug_general": {
        "is_estimated": True,
        "enzymes": [
            {"name": "CYP450 enzymes (CYP2D6, CYP3A4)", "action": "Substrate", "strength": "Moderate", "ec": "EC 1.14.14.1",
             "note": "Most approved CNS drugs are metabolised by hepatic CYP450; polymorphisms affect plasma concentration and response. " + _ESTI_NOTE, "source": "DrugBank; BRENDA"},
        ],
        "transporters": [
            {"name": "P-glycoprotein / ABCB1", "role": "Substrate (variable)", "location": "BBB efflux",
             "note": "CNS drug P-gp interaction is compound-specific; consult DrugBank entry for detailed transporter data.", "source": "DrugBank"},
        ],
        "cofactors": ["Compound-specific — consult DrugBank or BRENDA for precise data"],
        "sources": ["DrugBank", "BRENDA", "FDA Prescribing Information"],
    },
    "general": {
        "is_estimated": True,
        "enzymes": [
            {"name": "COX-2", "action": "Inhibition", "strength": "Weak", "ec": "EC 1.14.99.1",
             "note": "Estimated anti-inflammatory activity; strength unconfirmed for this compound. " + _ESTI_NOTE, "source": "Literature"},
            {"name": "Nrf2", "action": "Activation", "strength": "Weak", "ec": "Transcription factor",
             "note": "Estimated antioxidant response; requires compound-specific validation. " + _ESTI_NOTE, "source": "Literature"},
        ],
        "transporters": [
            {"name": "P-glycoprotein / ABCB1", "role": "Unknown / Compound-specific", "location": "BBB efflux",
             "note": "BBB transporter interaction not confirmed; consult ChEMBL or PubMed for this compound. " + _ESTI_NOTE, "source": "Literature"},
        ],
        "cofactors": ["Not confirmed — consult BRENDA or HMDB for this compound class"],
        "sources": ["Class estimate — consult BRENDA, MetaCyc, ChEMBL for compound-specific data"],
    },
}


def get_enzyme_entry(compound_name, compound_type=""):
    direct = ENZYME_DATA.get(compound_name)
    if direct:
        return direct
    alias_key = ENZYME_ALIASES.get(compound_name)
    if alias_key:
        aliased = ENZYME_DATA.get(alias_key)
        if aliased:
            return aliased
    ct = compound_type.lower()
    if any(x in ct for x in ['flavon', 'flavan', 'anthocyan', 'isoflavon', 'chalcon', 'xanthone']):
        return CLASS_ENZYME_TEMPLATES['flavonoid']
    if 'vitamin b1' in ct or 'thiamine' in ct:
        return CLASS_ENZYME_TEMPLATES['vitamin_b_thiamine']
    if 'vitamin b2' in ct or 'riboflavin' in ct:
        return CLASS_ENZYME_TEMPLATES['vitamin_b_riboflavin']
    if 'vitamin b3' in ct or 'niacin' in ct or ('nad' in ct and 'precursor' in ct):
        return CLASS_ENZYME_TEMPLATES['vitamin_b_niacin']
    if 'vitamin b5' in ct or 'pantothen' in ct or 'coa precursor' in ct:
        return CLASS_ENZYME_TEMPLATES['vitamin_b5']
    if 'vitamin b6' in ct or 'pyridox' in ct:
        return CLASS_ENZYME_TEMPLATES['vitamin_b6']
    if 'vitamin b9' in ct or 'folat' in ct or 'one-carbon' in ct:
        return CLASS_ENZYME_TEMPLATES['vitamin_b9']
    if 'vitamin c' in ct or 'ascorb' in ct:
        return CLASS_ENZYME_TEMPLATES['vitamin_c']
    if 'vitamin e' in ct or 'tocopherol' in ct or 'tocotrienol' in ct:
        return CLASS_ENZYME_TEMPLATES['vitamin_e']
    if 'carotenoid' in ct or 'lutein' in ct or 'lycopene' in ct or 'zeaxanthin' in ct or 'astaxanthin' in ct:
        return CLASS_ENZYME_TEMPLATES['carotenoid']
    if 'amino acid' in ct or 'cysteine' in ct or 'glutathione' in ct or 'gsh precursor' in ct or 'tripeptide' in ct or 'dipeptide' in ct:
        return CLASS_ENZYME_TEMPLATES['amino_acid']
    if 'alkaloid' in ct:
        return CLASS_ENZYME_TEMPLATES['alkaloid']
    if 'phospholipid' in ct or 'phosphatidyl' in ct:
        return CLASS_ENZYME_TEMPLATES['phospholipid']
    if 'fatty acid' in ct or 'omega' in ct or 'dha' in ct or 'epa' in ct or 'lipid' in ct:
        return CLASS_ENZYME_TEMPLATES['fatty_acid']
    if 'mineral' in ct or 'zinc' in ct or 'magnesium' in ct or 'selenium' in ct or 'lithium' in ct or 'trace' in ct:
        return CLASS_ENZYME_TEMPLATES['mineral']
    if 'organosulfur' in ct or 'allicin' in ct or 'isothiocyanate' in ct:
        return CLASS_ENZYME_TEMPLATES['organosulfur']
    if 'terpenoid' in ct or 'terpene' in ct or 'adapto' in ct or 'mushroom' in ct or 'herbal' in ct or 'withanolide' in ct or 'triterpenoid' in ct:
        return CLASS_ENZYME_TEMPLATES['terpenoid']
    if 'polyphenol' in ct or 'phenolic' in ct or 'stilbene' in ct or 'lignin' in ct or 'lignan' in ct or 'neolignan' in ct:
        return CLASS_ENZYME_TEMPLATES['polyphenol']
    if 'drug' in ct or 'ache' in ct or 'mao-b inhibitor drug' in ct or 'hdac inhibitor' in ct or 'approved' in ct:
        return CLASS_ENZYME_TEMPLATES['drug_general']
    return CLASS_ENZYME_TEMPLATES['general']


PATHWAY_METABOLITES = {
    "Nrf2/GSH": ["GSH", "HO-1", "NQO1"],
    "NF-kB": ["TNF-alpha", "IL-6", "IL-1beta"],
    "AMPK": ["ATP", "PGC-1alpha"],
    "SIRT1": ["NAD+", "PGC-1alpha"],
    "PGC-1alpha": ["ATP", "Mitochondrial biogenesis"],
    "PI3K/Akt": ["BDNF", "Cell survival"],
    "mTOR inhibition": ["Autophagy"],
    "Autophagy": ["Protein clearance"],
    "Autophagy/mitophagy": ["Protein clearance", "ATP"],
    "BDNF/TrkB": ["BDNF", "Synaptic plasticity"],
    "AChE inhibition": ["Acetylcholine"],
    "ETC/Mitochondria": ["ATP", "ROS"],
    "NMDA modulation": ["Ca2+", "Synaptic plasticity"],
    "CREB": ["BDNF", "cAMP"],
    "MAPK/ERK": ["BDNF"],
    "GSK-3beta inhibition": ["beta-catenin"],
    "NAD+/SIRT1": ["NAD+", "SIRT1"],
    "Abeta aggregation inhibition": ["Abeta"],
    "alpha-synuclein inhibition": ["alpha-synuclein"],
    "GABA modulation": ["GABA"],
    "NGF induction": ["NGF"],
    "Myelin synthesis": ["Myelin"],
    "Methionine cycle": ["Homocysteine", "GSH"],
    "Ca2+ homeostasis": ["Ca2+"],
    "Senolytic pathway": ["Senescent cells"],
    "VDR signaling": ["Ca2+", "BDNF"],
    "HPA axis modulation": ["Cortisol"],
    "Membrane fluidity": ["Membrane phospholipids"],
    "CDP-choline pathway": ["Acetylcholine", "Choline"],
    "Synaptogenesis": ["Synaptic plasticity"],
    "PCr/ATP buffer": ["ATP", "PCr"],
    "Beta-oxidation": ["ATP"],
    "Wnt/beta-catenin": ["beta-catenin"],
    "Selenoprotein P": ["GPx"],
    "GPx activity": ["GSH"],
    "MT1/MT2 receptor": ["Melatonin signaling"],
    "Circadian rhythm": ["Sleep quality"],
    "Estrogen receptor (ERb)": ["BDNF"],
    "Resolvin/Protectin synthesis": ["Resolvin"],
    "Iron chelation": ["Fe3+"],
    "PDE1 inhibition": ["cAMP"],
    "Cerebral blood flow": ["Dopamine"],
    "FAD/FMN cofactor": ["FAD", "ATP"],
    "Glycolysis / TCA cycle": ["ATP", "Pyruvate"],
    "Lipid peroxidation prevention": ["Membrane phospholipids"],
    "Singlet oxygen quenching": ["ROS"],
    "Radical scavenging": ["ROS"],
    "Myelination": ["Myelin"],
    "Glycine receptor": ["GABA"],
    "eIF5A hypusination": ["Polyamines"],
    "Serotonin reuptake inhibition": ["Serotonin"],
    "Dopamine regulation": ["Dopamine"],
    "NLRP3 inflammasome": ["IL-1beta"],
    "NLRP3 inflammasome inhibition": ["IL-1beta"],
    "Bioavailability enhancement": ["Bioavailability"],
    "Phosphatidylcholine synthesis": ["Membrane phospholipids"],
    "Anti-glycation": ["AGEs"],
    "Metal chelation": ["Fe3+"],
    "Monoamine synthesis": ["Serotonin", "Dopamine"],
    "GSH precursor": ["GSH"],
    "Membrane integrity": ["Membrane phospholipids"],
}

METABOLITE_LABELS = {
    "GSH": "Glutathione (antioxidant)",
    "NAD+": "NAD+ (energy/DNA repair)",
    "ATP": "ATP (neuronal energy)",
    "ROS": "ROS (oxidative load)",
    "BDNF": "BDNF (neuroplasticity)",
    "NGF": "NGF (neuron survival)",
    "Abeta": "Amyloid-beta (AD marker)",
    "alpha-synuclein": "alpha-Syn (PD marker)",
    "Acetylcholine": "ACh (memory/cognition)",
    "Dopamine": "Dopamine (motor/reward)",
    "GABA": "GABA (inhibitory neurotrans.)",
    "Serotonin": "Serotonin (mood/circadian)",
    "IL-6": "IL-6 (neuroinflammation)",
    "TNF-alpha": "TNF-alpha (neuroinflammation)",
    "IL-1beta": "IL-1beta (inflammasome)",
    "Cortisol": "Cortisol (HPA stress)",
    "Ca2+": "Ca2+ (synaptic signaling)",
    "Myelin": "Myelin (axonal integrity)",
    "Homocysteine": "Homocysteine (vascular risk)",
    "Senescent cells": "Senescent cells (clearance)",
    "Mitochondrial biogenesis": "Mitochondrial biogenesis",
    "Synaptic plasticity": "Synaptic plasticity (LTP)",
}


COMPOUND_CLASS_TEMPLATES = {
    "Flavonoid": {
        "compound_type": "Flavonoid / Polyphenol (Class-estimated — not in database)",
        "bbb": "Medium", "antioxidant": 7, "anti_inflammatory": 7, "mitochondrial_support": 6,
        "aggregation_modulation": 7, "cognitive_enhancement": 6, "neurogenesis": 5, "synaptic_plasticity": 6,
        "als": "Low", "alzheimers": "Med", "parkinsons": "Med", "huntingtons": "Low",
        "pathways": ["Nrf2/GSH", "NF-kB", "PI3K/Akt", "Radical scavenging"],
        "metabolites": ["GSH", "IL-6", "BDNF", "ROS"],
        "brain_regions": ["Hippocampus", "Cortex"],
    },
    "Polyphenol": {
        "compound_type": "Polyphenol / Phenolic (Class-estimated — not in database)",
        "bbb": "Medium", "antioxidant": 8, "anti_inflammatory": 7, "mitochondrial_support": 6,
        "aggregation_modulation": 7, "cognitive_enhancement": 6, "neurogenesis": 5, "synaptic_plasticity": 6,
        "als": "Low", "alzheimers": "Med", "parkinsons": "Low", "huntingtons": "Low",
        "pathways": ["Nrf2/GSH", "NF-kB", "Radical scavenging"],
        "metabolites": ["GSH", "ROS", "IL-6"],
        "brain_regions": ["Hippocampus", "Cortex"],
    },
    "Amino Acid": {
        "compound_type": "Amino Acid / Neurotransmitter Precursor (Class-estimated — not in database)",
        "bbb": "High", "antioxidant": 5, "anti_inflammatory": 5, "mitochondrial_support": 6,
        "aggregation_modulation": 4, "cognitive_enhancement": 6, "neurogenesis": 5, "synaptic_plasticity": 6,
        "als": "Low", "alzheimers": "Low", "parkinsons": "Low", "huntingtons": "Low",
        "pathways": ["Neurotransmitter synthesis", "Methionine cycle", "Nrf2/GSH"],
        "metabolites": ["GSH", "Serotonin", "Dopamine"],
        "brain_regions": ["Cortex", "Hippocampus"],
    },
    "Vitamin": {
        "compound_type": "Vitamin / Essential Micronutrient (Class-estimated — not in database)",
        "bbb": "High", "antioxidant": 7, "anti_inflammatory": 6, "mitochondrial_support": 7,
        "aggregation_modulation": 5, "cognitive_enhancement": 6, "neurogenesis": 5, "synaptic_plasticity": 5,
        "als": "Low", "alzheimers": "Low", "parkinsons": "Low", "huntingtons": "Low",
        "pathways": ["Nrf2/GSH", "Methionine cycle", "ETC/Mitochondria"],
        "metabolites": ["GSH", "NAD+", "ATP"],
        "brain_regions": ["Cortex", "Hippocampus"],
    },
    "Terpenoid": {
        "compound_type": "Terpenoid / Terpene (Class-estimated — not in database)",
        "bbb": "High", "antioxidant": 7, "anti_inflammatory": 7, "mitochondrial_support": 6,
        "aggregation_modulation": 6, "cognitive_enhancement": 6, "neurogenesis": 5, "synaptic_plasticity": 6,
        "als": "Low", "alzheimers": "Med", "parkinsons": "Med", "huntingtons": "Low",
        "pathways": ["Nrf2/GSH", "NF-kB", "AMPK"],
        "metabolites": ["GSH", "IL-6", "ATP"],
        "brain_regions": ["Hippocampus", "Cortex"],
    },
    "Alkaloid": {
        "compound_type": "Alkaloid (Class-estimated — not in database)",
        "bbb": "High", "antioxidant": 6, "anti_inflammatory": 6, "mitochondrial_support": 5,
        "aggregation_modulation": 6, "cognitive_enhancement": 7, "neurogenesis": 5, "synaptic_plasticity": 6,
        "als": "Low", "alzheimers": "Med", "parkinsons": "Med", "huntingtons": "Low",
        "pathways": ["Cholinergic signaling", "Dopamine regulation", "NF-kB"],
        "metabolites": ["Acetylcholine", "Dopamine", "BDNF"],
        "brain_regions": ["Cortex", "Hippocampus", "Striatum"],
    },
    "Fatty Acid": {
        "compound_type": "Fatty Acid / Lipid (Class-estimated — not in database)",
        "bbb": "High", "antioxidant": 6, "anti_inflammatory": 7, "mitochondrial_support": 7,
        "aggregation_modulation": 5, "cognitive_enhancement": 6, "neurogenesis": 5, "synaptic_plasticity": 6,
        "als": "Low", "alzheimers": "Med", "parkinsons": "Low", "huntingtons": "Low",
        "pathways": ["Membrane fluidity", "NF-kB", "BDNF/TrkB"],
        "metabolites": ["BDNF", "Membrane phospholipids", "IL-6"],
        "brain_regions": ["Cortex", "Hippocampus"],
    },
    "Drug": {
        "compound_type": "Pharmaceutical Drug / Small Molecule (Class-estimated — not in database)",
        "bbb": "High", "antioxidant": 5, "anti_inflammatory": 7, "mitochondrial_support": 5,
        "aggregation_modulation": 7, "cognitive_enhancement": 7, "neurogenesis": 5, "synaptic_plasticity": 6,
        "als": "Low", "alzheimers": "Med", "parkinsons": "Med", "huntingtons": "Low",
        "pathways": ["Nrf2/GSH", "NF-kB", "Cholinergic signaling"],
        "metabolites": ["Acetylcholine", "BDNF", "GSH"],
        "brain_regions": ["Hippocampus", "Cortex"],
    },
    "General": {
        "compound_type": "Bioactive Compound (Class-estimated — not in database)",
        "bbb": "Medium", "antioxidant": 5, "anti_inflammatory": 5, "mitochondrial_support": 5,
        "aggregation_modulation": 5, "cognitive_enhancement": 5, "neurogenesis": 5, "synaptic_plasticity": 5,
        "als": "Low", "alzheimers": "Low", "parkinsons": "Low", "huntingtons": "Low",
        "pathways": ["Nrf2/GSH", "NF-kB"],
        "metabolites": ["GSH", "ATP"],
        "brain_regions": ["Cortex", "Hippocampus"],
    },
}


def infer_compound_class(name):
    n = name.lower().strip()
    flavonoid_words = ["flavone", "flavonol", "flavanone", "flavonoid", "chalcone", "flavanol",
                       "anthocyanin", "isoflavone", "aurone", "biflavonoid"]
    if any(w in n for w in flavonoid_words):
        return "Flavonoid"
    if any(n.endswith(w) for w in ["flavone", "flavonol", "flavanone", "flavanol"]):
        return "Flavonoid"
    polyphenol_words = ["phenol", "catechol", "stilbene", "tannin", "coumarin", "lignin",
                        "phenolic", "caffeic", "ferulic", "rosmarin", "gallic", "ellagic",
                        "chlorogenic", "protocatechuic", "coumaric", "sinapic"]
    if any(w in n for w in polyphenol_words):
        return "Polyphenol"
    amino_words = ["alanine", "glycine", "serine", "cysteine", "methionine", "tyrosine",
                   "tryptophan", "glutamine", "glutamate", "aspartate", "histidine", "lysine",
                   "arginine", "proline", "valine", "leucine", "isoleucine", "phenylalanine",
                   "threonine", "ornithine", "carnosine", "taurine", "amino acid", "peptide"]
    if any(w in n for w in amino_words):
        return "Amino Acid"
    vitamin_words = ["vitamin", "riboflavin", "thiamin", "niacin", "pantothen", "pyridoxin",
                     "biotin", "folate", "folic", "cobalamin", "ascorbic", "tocopherol"]
    if any(w in n for w in vitamin_words):
        return "Vitamin"
    terpenoid_words = ["terpene", "terpen", "pinene", "limonene", "linalool", "carvacrol",
                       "thymol", "menthol", "camphor", "geraniol", "borneol", "terpineol",
                       "ursolic", "oleanolic", "betulin", "lupeol", "sterol", "steroid"]
    if any(w in n for w in terpenoid_words):
        return "Terpenoid"
    if n.endswith("ene") or n.endswith("terpene"):
        return "Terpenoid"
    fatty_words = ["fatty acid", "omega-", "linoleic", "linolenic", "oleic", "palmitic",
                   "stearic", "arachidonic", "eicosa", "docosa"]
    if any(w in n for w in fatty_words):
        return "Fatty Acid"
    alkaloid_words = ["alkaloid", "caffeine", "morphine", "codeine", "quinine", "berberine",
                      "piperine", "capsaicin", "colchicine", "vincamine", "galantamine",
                      "huperzine", "physostigmine", "rivastigmine", "donepezil"]
    if any(w in n for w in alkaloid_words):
        return "Alkaloid"
    drug_words = ["mab", "tidine", "prazole", "statin", "cillin", "mycin", "cycline",
                  "olol", "pril", "sartan", "zepine", "zolam", "pam", "xam"]
    if any(n.endswith(w) for w in drug_words):
        return "Drug"
    if n.endswith("ol") or n.endswith("in") or n.endswith("one"):
        return "Polyphenol"
    if n.endswith("ine"):
        return "Alkaloid"
    return "General"


def generate_estimated_data(compound_name):
    cls = infer_compound_class(compound_name)
    template = dict(COMPOUND_CLASS_TEMPLATES[cls])
    template["is_estimated"] = True
    template["estimated_class"] = cls
    return template


def fuzzy_match(query, names, n=6, cutoff=0.35):
    q = query.lower().strip()
    matches = difflib.get_close_matches(q, [n.lower() for n in names], n=n, cutoff=cutoff)
    result = []
    for m in matches:
        for name in names:
            if name.lower() == m:
                result.append(name)
                break
    return result


def score_color_hex(score):
    if score >= 70:
        return ("#1B6B45", "#E6F5EE")
    elif score >= 40:
        return ("#9B5C00", "#FFF3E0")
    else:
        return ("#9B2335", "#FDE8EA")


def bbb_color_hex(bbb):
    return {
        "High": ("#1B6B45", "#E6F5EE"),
        "Medium": ("#9B5C00", "#FFF3E0"),
        "Low-Med": ("#9B5C00", "#FFF3E0"),
        "Low": ("#9B2335", "#FDE8EA"),
    }.get(bbb, ("#64748B", "#F1F5F9"))


def disease_style(level):
    return {
        "High": ("#1B6B45", "#E6F5EE", "High", "Strong research evidence"),
        "Med":  ("#9B5C00", "#FFF3E0", "Moderate", "Emerging / supporting evidence"),
        "Low":  ("#9B2335", "#FDE8EA", "Low", "Limited or indirect evidence"),
    }.get(level, ("#64748B", "#F1F5F9", level, ""))


def pill(text, fg, bg):
    return f'<span class="pill" style="color:{fg};background:{bg};">{text}</span>'


def make_radar(data):
    dims = [
        "Antioxidant", "Anti-Inflammatory", "Mitochondrial",
        "Aggregation Mod.", "Cognitive", "Neurogenesis", "Synaptic Plasticity",
    ]
    dims_full = [
        "Antioxidant Activity", "Anti-Inflammatory", "Mitochondrial Support",
        "Aggregation Modulation", "Cognitive Enhancement",
        "Neurogenesis Support", "Synaptic Plasticity",
    ]
    vals = [
        data['antioxidant'], data['anti_inflammatory'], data['mitochondrial_support'],
        data['aggregation_modulation'], data['cognitive_enhancement'],
        data['neurogenesis'], data['synaptic_plasticity'],
    ]
    vals_c = vals + [vals[0]]
    dims_c = dims + [dims[0]]
    dims_full_c = dims_full + [dims_full[0]]

    ref5 = [5] * len(dims)
    ref5_c = ref5 + [ref5[0]]

    ref8 = [8] * len(dims)
    ref8_c = ref8 + [ref8[0]]

    fig = go.Figure()

    fig.add_trace(go.Scatterpolar(
        r=ref8_c, theta=dims_c, mode='lines',
        line=dict(color='#E2E8F2', width=0.8, dash='dot'),
        showlegend=False, hoverinfo='skip'
    ))
    fig.add_trace(go.Scatterpolar(
        r=ref5_c, theta=dims_c, mode='lines',
        line=dict(color='#94A3B8', width=1.2, dash='dash'),
        showlegend=True, hoverinfo='skip', name='Reference (5)'
    ))
    fig.add_trace(go.Scatterpolar(
        r=vals_c, theta=dims_c, fill='toself',
        fillcolor='rgba(15,81,50,0.32)',
        line=dict(color='#0F5132', width=2.8),
        marker=dict(size=8, color='#0F5132', line=dict(width=2.5, color='white')),
        name='Score',
        customdata=dims_full_c,
        hovertemplate='<b>%{customdata}</b><br>Score: <b>%{r}/10</b><extra></extra>'
    ))

    fig.update_layout(
        polar=dict(
            bgcolor='#FAFCFF',
            hole=0.08,
            radialaxis=dict(
                visible=True, range=[0, 10],
                tickvals=[2, 4, 6, 8, 10],
                ticktext=['2', '4', '6', '8', '10'],
                tickfont=dict(size=9, color='#94A3B8', family='Inter, sans-serif'),
                gridcolor='#E2E8F0', linecolor='#CBD5E1',
                tickangle=0,
                tickprefix='  ',
            ),
            angularaxis=dict(
                tickfont=dict(size=11, color='#1E293B', family='Inter, sans-serif'),
                gridcolor='#E2E8F0', linecolor='#CBD5E1',
                direction='clockwise',
                ticklabelstep=1,
            )
        ),
        showlegend=True,
        legend=dict(
            orientation='h', yanchor='bottom', y=-0.12, xanchor='center', x=0.5,
            font=dict(size=10, color='#475569')
        ),
        margin=dict(l=88, r=88, t=50, b=55),
        paper_bgcolor='white',
        height=440,
    )
    return fig


def make_network(compound_name, data):
    pathways = data.get('pathways', [])
    raw_metabolites = data.get('metabolites', [])
    brain_regions = data.get('brain_regions', [])[:4]

    disease_relevance = []
    if data.get('als') in ['High', 'Med']:
        disease_relevance.append("ALS")
    if data.get('alzheimers') in ['High', 'Med']:
        disease_relevance.append("Alzheimer's")
    if data.get('parkinsons') in ['High', 'Med']:
        disease_relevance.append("Parkinson's")
    if data.get('huntingtons') in ['High', 'Med']:
        disease_relevance.append("Huntington's")

    expanded_mets = set(raw_metabolites)
    for pw in pathways:
        for met in PATHWAY_METABOLITES.get(pw, [])[:2]:
            expanded_mets.add(met)
    metabolites = list(expanded_mets)[:10]

    def ring_positions(items, radius, offset_deg=0):
        n = len(items)
        positions = {}
        for i, item in enumerate(items):
            angle = math.radians(90 + offset_deg - i * 360 / n) if n > 0 else 0
            positions[item] = (radius * math.cos(angle), radius * math.sin(angle))
        return positions

    all_nodes = {}
    all_nodes[compound_name] = (0, 0)

    pw_pos = ring_positions(pathways, 1.85, offset_deg=15)
    all_nodes.update(pw_pos)

    met_pos = ring_positions(metabolites, 3.6, offset_deg=-10)
    all_nodes.update(met_pos)

    outer = brain_regions + disease_relevance
    total_outer = len(outer)
    outer_pos = {}
    for i, item in enumerate(outer):
        angle = math.radians(90 - i * 360 / total_outer) if total_outer > 0 else 0
        outer_pos[item] = (5.3 * math.cos(angle), 5.3 * math.sin(angle))
    all_nodes.update(outer_pos)

    EDGE_STYLES = {
        'cp': {'color': '#3B82F6', 'width': 1.8, 'label': 'Compound → Pathway'},
        'pm': {'color': '#10B981', 'width': 1.4, 'label': 'Pathway → Metabolite'},
        'md': {'color': '#EF4444', 'width': 1.4, 'label': 'Metabolite → Disease'},
        'pb': {'color': '#F59E0B', 'width': 1.2, 'label': 'Pathway → Brain Region'},
    }

    typed_edges = []
    for pw in pathways:
        typed_edges.append((compound_name, pw, 'cp'))
    for pw in pathways:
        for met in PATHWAY_METABOLITES.get(pw, [])[:2]:
            if met in all_nodes:
                typed_edges.append((pw, met, 'pm'))
    for met in raw_metabolites:
        if met in all_nodes:
            for d in disease_relevance:
                if d in ["ALS", "Huntington's"] and met in ["ATP", "PCr", "alpha-synuclein", "GSH"]:
                    typed_edges.append((met, d, 'md'))
                elif d == "Alzheimer's" and met in ["Abeta", "Acetylcholine", "BDNF", "GSH"]:
                    typed_edges.append((met, d, 'md'))
                elif d == "Parkinson's" and met in ["alpha-synuclein", "Dopamine", "GSH"]:
                    typed_edges.append((met, d, 'md'))
    for pw in pathways[:2]:
        for br in brain_regions[:2]:
            typed_edges.append((pw, br, 'pb'))

    NODE_STYLES = {
        'compound':     {'color': '#FAFCFF', 'size': 38, 'border': '#F0A500', 'border_width': 3.5, 'symbol': 'circle'},
        'pathway':      {'color': '#0B4D2E', 'size': 20, 'border': '#34D399', 'border_width': 2.2, 'symbol': 'square'},
        'metabolite':   {'color': '#3B1F6E', 'size': 16, 'border': '#A78BFA', 'border_width': 2.0, 'symbol': 'diamond'},
        'disease':      {'color': '#7F1D1D', 'size': 19, 'border': '#FCA5A5', 'border_width': 2.2, 'symbol': 'triangle-up'},
        'brain_region': {'color': '#6B3A10', 'size': 16, 'border': '#FCD34D', 'border_width': 2.0, 'symbol': 'circle'},
    }

    def ntype(name):
        if name == compound_name: return 'compound'
        if name in pathways: return 'pathway'
        if name in metabolites: return 'metabolite'
        if name in disease_relevance: return 'disease'
        return 'brain_region'

    traces = []

    for etype, style in EDGE_STYLES.items():
        ex, ey = [], []
        for u, v, et in typed_edges:
            if et == etype and u in all_nodes and v in all_nodes:
                x0, y0 = all_nodes[u]
                x1, y1 = all_nodes[v]
                ex += [x0, x1, None]
                ey += [y0, y1, None]
        if ex:
            traces.append(go.Scatter(
                x=ex, y=ey, mode='lines',
                line=dict(width=style['width'] * 1.1, color=style['color']),
                opacity=0.72,
                hoverinfo='none', showlegend=False
            ))

    group_legend = {
        'compound':     'Compound (circle)',
        'pathway':      'Pathway (square)',
        'metabolite':   'Metabolite (diamond)',
        'disease':      'Disease (triangle)',
        'brain_region': 'Brain Region (circle)',
    }
    seen_legend = set()

    for name, (x, y) in all_nodes.items():
        nt = ntype(name)
        nd = NODE_STYLES[nt]
        show_legend = nt not in seen_legend
        seen_legend.add(nt)
        label = METABOLITE_LABELS.get(name, name)
        hover = f"<b>{name}</b><br><i style='color:#aaa;font-size:11px;'>{label}</i>"
        traces.append(go.Scatter(
            x=[x], y=[y], mode='markers',
            marker=dict(
                size=nd['size'], color=nd['color'],
                line=dict(width=nd['border_width'], color=nd['border']),
                opacity=0.97,
                symbol=nd['symbol']
            ),
            hovertext=hover,
            hoverinfo='text',
            name=group_legend[nt],
            legendgroup=nt,
            showlegend=show_legend
        ))

    arrow_annotations = []
    for u, v, etype in typed_edges:
        if u in all_nodes and v in all_nodes:
            x0, y0 = all_nodes[u]
            x1, y1 = all_nodes[v]
            ec = EDGE_STYLES[etype]['color']
            angle = math.atan2(y1 - y0, x1 - x0)
            node_r = NODE_STYLES[ntype(v)]['size'] * 0.028
            ax1 = x1 - node_r * math.cos(angle)
            ay1 = y1 - node_r * math.sin(angle)
            arrow_annotations.append(dict(
                x=ax1, y=ay1, ax=x0, ay=y0,
                xref='x', yref='y', axref='x', ayref='y',
                text='', showarrow=True,
                arrowhead=2, arrowsize=1.0,
                arrowwidth=EDGE_STYLES[etype]['width'] * 0.9,
                arrowcolor=ec,
                opacity=0.75
            ))

    label_annotations = []
    for name, (x, y) in all_nodes.items():
        nt = ntype(name)
        nd = NODE_STYLES[nt]
        display_text = name if len(name) <= 16 else name[:14] + '..'
        if nt == 'compound':
            label_annotations.append(dict(
                x=x, y=y,
                text=f'<b>{display_text}</b>',
                font=dict(size=11, color='#0B2040', family='Inter, sans-serif'),
                showarrow=False, xref='x', yref='y',
                bgcolor='rgba(0,0,0,0)', borderpad=2, align='center',
            ))
        else:
            offset = nd['size'] * 0.022 + 0.48
            ty = y + offset if y >= 0 else y - offset
            label_annotations.append(dict(
                x=x, y=ty,
                text=display_text,
                font=dict(size=9.5, color='#0D2137', family='Inter, sans-serif'),
                showarrow=False, xref='x', yref='y',
                bgcolor='rgba(255,255,255,0.88)',
                bordercolor='rgba(148,163,184,0.45)',
                borderpad=2.5, borderwidth=0.8, align='center',
            ))

    for r, dash, lw in [(1.85, 'dot', 1.5), (3.6, 'dot', 1.2), (5.3, 'dot', 1.1)]:
        theta = [math.radians(a) for a in range(361)]
        traces.append(go.Scatter(
            x=[r * math.cos(t) for t in theta],
            y=[r * math.sin(t) for t in theta],
            mode='lines', line=dict(color='#94A3B8', width=lw, dash=dash),
            hoverinfo='none', showlegend=False
        ))

    ring_label_annotations = []
    for r, lbl in [(1.85, 'Pathways'), (3.6, 'Metabolites'), (5.3, 'Brain & Disease')]:
        ring_label_annotations.append(dict(
            x=r * math.cos(math.radians(-28)),
            y=r * math.sin(math.radians(-28)),
            text=f'<i>{lbl}</i>',
            font=dict(size=8.5, color='#64748B', family='Inter, sans-serif'),
            showarrow=False, xref='x', yref='y',
            bgcolor='rgba(255,255,255,0.75)',
            borderpad=2, align='center',
        ))

    annotations = arrow_annotations + label_annotations + ring_label_annotations

    fig = go.Figure(data=traces)
    fig.update_layout(
        annotations=annotations,
        showlegend=True,
        legend=dict(
            orientation='h', yanchor='bottom', y=-0.32, xanchor='center', x=0.5,
            font=dict(size=10.5, color='#475569', family='Inter, sans-serif'),
            bgcolor='rgba(255,255,255,0)', itemsizing='constant',
            traceorder='normal'
        ),
        hovermode='closest',
        dragmode='pan',
        margin=dict(l=20, r=20, t=24, b=72),
        height=620,
        paper_bgcolor='white',
        plot_bgcolor='#FAFCFF',
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, showline=False, range=[-8.5, 8.5]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, showline=False, range=[-8.5, 8.5],
                   scaleanchor='x', scaleratio=1),
    )
    return fig


def generate_assessment(compound_name, data, score, enzyme_entry=None):
    bbb = data.get('bbb', 'Unknown')
    pathways = data.get('pathways', [])
    dimensions = {
        'antioxidant': ('Antioxidant Activity', data.get('antioxidant', 0)),
        'anti_inflammatory': ('Anti-Inflammatory', data.get('anti_inflammatory', 0)),
        'mitochondrial_support': ('Mitochondrial Support', data.get('mitochondrial_support', 0)),
        'aggregation_modulation': ('Aggregation Modulation', data.get('aggregation_modulation', 0)),
        'cognitive_enhancement': ('Cognitive Enhancement', data.get('cognitive_enhancement', 0)),
        'neurogenesis': ('Neurogenesis Support', data.get('neurogenesis', 0)),
        'synaptic_plasticity': ('Synaptic Plasticity', data.get('synaptic_plasticity', 0)),
    }
    top_dims = sorted(dimensions.items(), key=lambda x: x[1][1], reverse=True)[:3]
    top_dim_labels = [v[0] for _, v in top_dims if v[1] >= 5]

    diseases = {
        'ALS': data.get('als', 'Low'),
        "Alzheimer's disease": data.get('alzheimers', 'Low'),
        "Parkinson's disease": data.get('parkinsons', 'Low'),
        "Huntington's disease": data.get('huntingtons', 'Low'),
    }
    high_rel = [d for d, v in diseases.items() if v == 'High']
    med_rel  = [d for d, v in diseases.items() if v == 'Med']

    if score >= 80:
        strength_word = "strong"
        evidence_phrase = "supported by convergent mechanistic evidence across multiple neuroprotective dimensions"
    elif score >= 60:
        strength_word = "moderate"
        evidence_phrase = "with established mechanistic targets across several neuroprotective dimensions"
    elif score >= 40:
        strength_word = "emerging"
        evidence_phrase = "with selective activity in specific neuroprotective pathways, primarily from in vitro evidence"
    else:
        strength_word = "preliminary"
        evidence_phrase = "currently supported by limited in vitro or animal model data"

    bbb_phrase = {
        'High': "readily crosses the blood-brain barrier, enabling direct CNS engagement at physiological doses",
        'Moderate': "shows moderate blood-brain barrier penetration; bioavailability-enhanced formulations may improve CNS delivery",
        'Low': "demonstrates limited blood-brain barrier penetration, which may restrict direct CNS activity",
        'Very Low': "has very limited CNS bioavailability due to restricted blood-brain barrier penetration; targeted delivery strategies are likely required",
    }.get(bbb, "has uncertain blood-brain barrier penetration status based on current data")

    paras = []

    if top_dim_labels:
        dim_str = ', '.join(top_dim_labels[:-1]) + (' and ' + top_dim_labels[-1] if len(top_dim_labels) > 1 else top_dim_labels[0] if top_dim_labels else '')
        paras.append(
            f"{compound_name} demonstrates {strength_word} neuroprotective potential ({evidence_phrase}). "
            f"The compound {bbb_phrase}. "
            f"Its primary mechanistic strengths lie in {dim_str}."
        )
    else:
        paras.append(
            f"{compound_name} demonstrates {strength_word} neuroprotective potential ({evidence_phrase}). "
            f"The compound {bbb_phrase}."
        )

    if high_rel or med_rel:
        rel_parts = []
        if high_rel:
            rel_parts.append(f"strongest evidence for {' and '.join(high_rel)}")
        if med_rel:
            rel_parts.append(f"moderate relevance for {' and '.join(med_rel)}")
        disease_sent = f"Current peer-reviewed research indicates {', and '.join(rel_parts)}."
        if pathways:
            pathway_str = ', '.join(pathways[:3])
            extra = f" and {len(pathways)-3} additional pathways" if len(pathways) > 3 else ""
            disease_sent += f" Key mechanistic pathways include {pathway_str}{extra}, suggesting pleiotropic neuroprotective activity."
        paras.append(disease_sent)
    else:
        if pathways:
            paras.append(
                f"Current evidence does not establish strong disease-specific relevance for the major neurodegenerative conditions modelled, "
                f"though general neuroprotective mechanisms are documented across {', '.join(pathways[:3])} and related pathways."
            )
        else:
            paras.append("Current evidence does not establish strong relevance for specific neurodegenerative conditions at this time.")

    if enzyme_entry and enzyme_entry.get('enzymes'):
        enzymes = enzyme_entry['enzymes']
        inhibited = [e['name'] for e in enzymes if 'inhib' in e.get('action', '').lower()]
        activated = [e['name'] for e in enzymes if any(w in e.get('action', '').lower() for w in ('activ', 'upregul', 'induct'))]
        cofactors  = enzyme_entry.get('cofactors', [])
        transporters = enzyme_entry.get('transporters', [])
        enz_parts = []
        if inhibited:
            enz_parts.append(f"inhibits {', '.join(inhibited[:3])}")
        if activated:
            enz_parts.append(f"activates or upregulates {', '.join(activated[:2])}")
        if enz_parts:
            mol_sent = f"At the molecular level, {compound_name} " + ", while it ".join(enz_parts) + "."
        else:
            mol_sent = ""
        if transporters:
            tnames = [t['name'].split('/')[0].strip() for t in transporters[:2]]
            mol_sent += f" Relevant transporters include {' and '.join(tnames)}."
        if cofactors:
            mol_sent += f" Key cofactors involved: {', '.join(cofactors[:3])}."
        if mol_sent:
            paras.append(mol_sent)

    if score >= 75:
        paras.append(
            f"Overall, {compound_name} represents a scientifically compelling candidate for translational neuroprotection research. "
            "Optimised bioavailability formulations — including nanoparticle encapsulation, liposomal delivery, or structural analogues — "
            "may enhance CNS exposure. Rigorously controlled clinical trials are warranted to establish therapeutic dosing and safety."
        )
    elif score >= 50:
        paras.append(
            f"Overall, {compound_name} warrants continued investigation as a neuroprotective agent. "
            "Further mechanistic characterisation in translational models and bioavailability optimisation are recommended before clinical application."
        )
    else:
        paras.append(
            f"Overall, {compound_name} may contribute modestly as part of a multimodal brain health strategy. "
            "Its neuroprotective profile is currently best characterised in preclinical settings; additional research is needed to confirm human relevance."
        )

    if enzyme_entry and enzyme_entry.get('sources'):
        src_list = enzyme_entry['sources'][:4]
        paras.append(f"Primary data sources: {'; '.join(src_list)}.")
    else:
        paras.append("Data sources: Published neuroscience literature (PubMed/PMC), ChEMBL, DrugBank.")

    return paras


def render_external_data(compound_name):
    """
    Render a panel showing live data fetched from PubChem and ChEMBL.
    Data is fetched automatically (cached in-process via lru_cache).
    Called from render_report.
    """
    st.divider()
    st.subheader("Live External Database Lookup")
    st.caption(
        "Physicochemical properties, BBB penetration prediction, and mechanism-of-action data "
        "fetched in real time from PubChem and ChEMBL. Results are sourced directly from "
        "the primary databases and clearly labelled."
    )

    with st.spinner("Querying PubChem and ChEMBL..."):
        pc = fetch_pubchem(compound_name)
        ch = fetch_chembl(compound_name)

    st.success("Live External Data — PubChem & ChEMBL  |  Fetched directly from primary databases")

    col_pc, col_ch = st.columns(2)

    with col_pc:
        st.markdown("**PubChem Data**")
        if "error" in pc:
            st.warning(f"PubChem: {pc['error']}")
        else:
            cid_link = f'<a href="{pc.get("pubchem_url","#")}" target="_blank" style="color:#1A6FBD;">CID {pc.get("cid","")}</a>'
            st.markdown(
                f'<p style="font-size:0.82rem;color:#1A2B45;margin:0 0 8px 0;">'
                f'<b>PubChem Compound:</b> {cid_link}</p>',
                unsafe_allow_html=True
            )
            m1, m2, m3 = st.columns(3)
            m1.metric("Formula",  pc.get("formula", "N/A"))
            m2.metric("MW (Da)",  pc.get("mw", "N/A"))
            m3.metric("XLogP",    pc.get("xlogp", "N/A"))

            prop_rows = [
                ("Molecular Formula",   pc.get("formula",         "N/A")),
                ("Molecular Weight",    f"{pc.get('mw','N/A')} Da"     if pc.get("mw") != "N/A" else "N/A"),
                ("XLogP (lipophilicity)", str(pc.get("xlogp",    "N/A"))),
                ("TPSA",                f"{pc.get('tpsa','N/A')} Å²"   if pc.get("tpsa") != "N/A" else "N/A"),
                ("H-Bond Donors",       str(pc.get("hbd",         "N/A"))),
                ("H-Bond Acceptors",    str(pc.get("hba",         "N/A"))),
                ("Rotatable Bonds",     str(pc.get("rotatable_bonds","N/A"))),
                ("Heavy Atom Count",    str(pc.get("heavy_atoms", "N/A"))),
                ("InChIKey",            pc.get("inchikey",        "N/A")),
            ]
            rows_html = "".join(
                f'<tr><td style="color:#64748B;font-size:0.78rem;padding:3px 10px 3px 0;">{k}</td>'
                f'<td style="color:#1A2B45;font-size:0.82rem;font-weight:600;padding:3px 0;">{v}</td></tr>'
                for k, v in prop_rows
            )
            st.markdown(
                f'<table style="border-collapse:collapse;width:100%;margin-bottom:10px;">{rows_html}</table>',
                unsafe_allow_html=True
            )

            syns = pc.get("synonyms", [])
            if syns:
                st.markdown(
                    f'<p style="font-size:0.78rem;color:#64748B;margin:0 0 6px 0;">'
                    f'<b>Known synonyms:</b> {", ".join(syns)}</p>',
                    unsafe_allow_html=True
                )

            bbb = pc.get("bbb_prediction", {})
            if bbb and "level" in bbb:
                lvl      = bbb["level"]
                bbb_colours = {
                    "High": ("#14532D", "#F0FFF4"),
                    "Medium": ("#92710C", "#FFFBEA"),
                    "Low-Med": ("#7C3009", "#FFF7ED"),
                    "Low": ("#7F1D1D", "#FFF1F2"),
                }
                fg, bg = bbb_colours.get(lvl, ("#334155", "#F8FAFC"))
                st.markdown(
                    f'<div style="background:{bg};border:1.5px solid {fg};border-radius:6px;'
                    f'padding:10px 14px;margin:8px 0;">'
                    f'<span style="font-size:0.73rem;font-weight:700;color:{fg};letter-spacing:0.8px;">'
                    f'COMPUTED BBB PENETRATION — {lvl.upper()}</span>'
                    f'<p style="font-size:0.8rem;color:#1A2B45;margin:4px 0 4px 0;line-height:1.55;">'
                    f'{bbb["rationale"]}</p>'
                    f'<p style="font-size:0.7rem;color:#64748B;margin:2px 0 0 0;">'
                    f'Method: CNS-MPO-inspired physicochemical rules. '
                    f'Ref: {bbb.get("ref","")}</p>'
                    f'</div>',
                    unsafe_allow_html=True
                )

            ro5 = bbb.get("ro5_violations", "N/A")
            if ro5 != "N/A":
                ro5_label = "Pass" if int(ro5) <= 1 else f"{ro5} violation(s)"
                ro5_colour = "#15803D" if int(ro5) <= 1 else "#B45309"
                st.markdown(
                    f'<p style="font-size:0.78rem;color:#1A2B45;margin:0 0 4px 0;">'
                    f'<b>Lipinski Rule of Five:</b> '
                    f'<span style="color:{ro5_colour};font-weight:700;">{ro5_label}</span></p>',
                    unsafe_allow_html=True
                )

    with col_ch:
        st.markdown("**ChEMBL Data**")
        if "error" in ch:
            st.warning(f"ChEMBL: {ch['error']}")
        else:
            ch_link = (
                f'<a href="{ch.get("chembl_url","#")}" target="_blank" style="color:#1A6FBD;">'
                f'{ch.get("chembl_id","")}</a>'
            )
            st.markdown(
                f'<p style="font-size:0.82rem;color:#1A2B45;margin:0 0 8px 0;">'
                f'<b>ChEMBL ID:</b> {ch_link}</p>',
                unsafe_allow_html=True
            )

            ch_m1, ch_m2, ch_m3 = st.columns(3)
            ch_m1.metric("Max Phase", ch.get("max_phase", "N/A"))
            ch_m2.metric("QED Score", f"{float(ch.get('qed',0)):.2f}" if ch.get("qed") not in ("N/A", None) else "N/A")
            ch_m3.metric("MW (Da)",   ch.get("mw", "N/A"))

            max_phase = ch.get("max_phase")
            phase_map = {
                4: ("Approved Drug", "#14532D", "#F0FFF4"),
                3: ("Phase 3 Clinical Trial", "#1A4E8A", "#EFF6FF"),
                2: ("Phase 2 Clinical Trial", "#1A4E8A", "#EFF6FF"),
                1: ("Phase 1 Clinical Trial", "#64748B", "#F8FAFC"),
                0: ("Preclinical / Research Only", "#64748B", "#F8FAFC"),
            }
            if max_phase is not None:
                try:
                    _phase_key = int(float(max_phase))
                except (TypeError, ValueError):
                    _phase_key = -1
                plabel, pfg, pbg = phase_map.get(_phase_key, (f"Phase {max_phase}", "#64748B", "#F8FAFC"))
                st.markdown(
                    f'<div style="background:{pbg};border:1px solid {pfg};border-radius:5px;'
                    f'padding:4px 10px;display:inline-block;margin-bottom:8px;">'
                    f'<span style="color:{pfg};font-size:0.76rem;font-weight:700;">'
                    f'Drug Status: {plabel}</span></div>',
                    unsafe_allow_html=True
                )

            ch_prop_rows = [
                ("Molecule Type",     ch.get("mol_type",   "N/A")),
                ("ALogP",             str(ch.get("alogp",  "N/A"))),
                ("MW (free base)",    f"{ch.get('mw','N/A')} Da" if ch.get("mw") not in ("N/A", None) else "N/A"),
                ("TPSA",              f"{ch.get('tpsa','N/A')} Å²" if ch.get("tpsa") not in ("N/A", None) else "N/A"),
                ("H-Bond Donors",     str(ch.get("hbd",    "N/A"))),
                ("H-Bond Acceptors",  str(ch.get("hba",    "N/A"))),
                ("Ro5 Violations",    str(ch.get("ro5_violations","N/A"))),
                ("QED Score",         f"{float(ch.get('qed',0)):.3f}" if ch.get("qed") not in ("N/A", None) else "N/A"),
                ("Oral",              "Yes" if ch.get("oral") else ("No" if ch.get("oral") is False else "N/A")),
            ]
            rows_html2 = "".join(
                f'<tr><td style="color:#64748B;font-size:0.78rem;padding:3px 10px 3px 0;">{k}</td>'
                f'<td style="color:#1A2B45;font-size:0.82rem;font-weight:600;padding:3px 0;">{v}</td></tr>'
                for k, v in ch_prop_rows
            )
            st.markdown(
                f'<table style="border-collapse:collapse;width:100%;margin-bottom:10px;">{rows_html2}</table>',
                unsafe_allow_html=True
            )

            qed = ch.get("qed")
            if qed not in ("N/A", None):
                qed_val = float(qed)
                if qed_val >= 0.6:
                    qed_label, qed_colour = "Drug-like (QED ≥ 0.6)", "#15803D"
                elif qed_val >= 0.35:
                    qed_label, qed_colour = "Moderately drug-like (QED 0.35–0.6)", "#B45309"
                else:
                    qed_label, qed_colour = "Low drug-likeness (QED < 0.35)", "#B91C1C"
                st.markdown(
                    f'<p style="font-size:0.78rem;color:#1A2B45;margin:0 0 4px 0;">'
                    f'<b>Drug-likeness (QED):</b> '
                    f'<span style="color:{qed_colour};font-weight:700;">{qed_label}</span></p>',
                    unsafe_allow_html=True
                )

            mechanisms = ch.get("mechanisms", [])
            if mechanisms:
                st.markdown(
                    '<div style="font-size:0.73rem;font-weight:700;color:#64748B;letter-spacing:1.2px;'
                    'text-transform:uppercase;margin:12px 0 6px 0;">Mechanism of Action (ChEMBL)</div>',
                    unsafe_allow_html=True
                )
                for m in mechanisms:
                    st.markdown(
                        f'<div style="background:#F0F7FF;border-left:3px solid #1A6FBD;'
                        f'border-radius:0 5px 5px 0;padding:5px 10px;margin-bottom:5px;">'
                        f'<span style="font-size:0.82rem;font-weight:600;color:#1A2B45;">'
                        f'{m.get("target","")}</span>'
                        f'<span style="font-size:0.78rem;color:#64748B;"> — {m.get("action","N/A")}</span>'
                        f'</div>',
                        unsafe_allow_html=True
                    )

    st.markdown(
        '<p style="font-size:0.72rem;color:#94A3B8;margin:12px 0 0 0;border-top:1px solid #E2E8F0;padding-top:8px;">'
        'Live data fetched via PubChem REST API (pubchem.ncbi.nlm.nih.gov) and ChEMBL REST API (ebi.ac.uk/chembl). '
        'BBB penetration is a computational prediction based on physicochemical rules '
        '(Wager et al. 2010, ACS Chem Neurosci; Lipinski et al. 2001, Adv Drug Deliv Rev). '
        'It is not a direct experimental measurement. Data is provided for educational and research use only.</p>',
        unsafe_allow_html=True
    )



def render_report(compound_name, data):
    score = neuro_score(data)
    bbb = data.get('bbb', 'N/A')
    sc_fg, sc_bg = score_color_hex(score)
    bb_fg, bb_bg = bbb_color_hex(bbb)
    score_label = "Strong" if score >= 70 else ("Moderate" if score >= 40 else "Limited")

    _is_ml   = bool(data.get("ml_predicted"))
    _is_est  = bool(data.get("is_estimated")) and not _is_ml

    if _is_est:
        cls_lbl = data.get('estimated_class', 'General')
        st.markdown(
            f'<div class="estimated-banner">'
            f'<span class="estimated-icon">&#9888;</span>'
            f'<span class="estimated-text">'
            f'<strong>Compound not found in database \u2014 showing class-based estimates</strong><br>'
            f'&ldquo;<b>{compound_name}</b>&rdquo; was identified as a <b>{cls_lbl}</b> compound based on name analysis. '
            f'Scores and pathways below reflect typical values for this compound class from published neuroprotection literature. '
            f'For precise data, consult PubMed, ChEMBL, or DrugBank directly.'
            f'</span></div>',
            unsafe_allow_html=True
        )
    elif _is_ml:
        r2 = data.get("model_cv_r2", 0.0)
        chembl_id = data.get("chembl_id", "")
        _live_ml = bool(data.get("_live_ml"))
        _cid = data.get("pubchem_cid", "")
        if _live_ml:
            _ml_body = (
                f'<strong>Live ML-Estimated Profile:</strong> <b>{compound_name}</b> was not in the '
                f'BrainSafe AI database. Its physicochemical properties were fetched live from '
                f'PubChem (CID {_cid}) and mechanism-of-action data from ChEMBL'
                + (f' ({chembl_id})' if chembl_id else '') +
                f'. A Random Forest model trained on 129 curated compounds '
                f'(5-fold CV R\u00b2 = {r2:.2f}) then predicted the 7-dimension neuroprotective profile. '
                f'This is a <b>live computational estimate</b> — validate in PubMed, ChEMBL, or DrugBank.'
            )
        else:
            _ml_body = (
                f'<strong>ML-Predicted Profile:</strong> Bioactivity scores for <b>{compound_name}</b> '
                f'were pre-computed by a Random Forest model trained on 129 manually curated compounds '
                f'(5-fold CV R\u00b2 = {r2:.2f}). Compound sourced from ChEMBL drug indication data '
                f'for neurological conditions'
                + (f' (ChEMBL: {chembl_id})' if chembl_id else '') +
                f'. This is a <b>computational prediction</b> — use as a starting point and '
                f'cross-reference with primary databases for validation.'
            )
        st.markdown(
            f'<div style="background:#EFF6FF;border:1.5px solid #93C5FD;border-radius:8px;'
            f'padding:10px 16px;margin:0 0 10px 0;display:flex;gap:12px;align-items:flex-start;">'
            f'<span style="font-size:1.1rem;color:#1D4ED8;flex-shrink:0;">&#9881;</span>'
            f'<span style="font-size:0.83rem;color:#1E40AF;line-height:1.6;">{_ml_body}</span></div>',
            unsafe_allow_html=True
        )

    nps_pill = pill(f"Neuroprotective Score: {score}/100 \u2014 {score_label}", sc_fg, sc_bg)
    bbb_pill = pill(f"BBB Penetration: {bbb}", bb_fg, bb_bg)
    if _is_est:
        extra_pill = pill("CLASS ESTIMATE", "#7A5500", "#FFFBEE")
    elif _is_ml:
        extra_pill = pill("ML-PREDICTED", "#1D4ED8", "#DBEAFE")
    else:
        extra_pill = ""
    st.markdown(
        f'<div class="compound-banner">'
        f'<p class="cmpd-name">{compound_name}</p>'
        f'<p class="cmpd-type">{data.get("compound_type", "")}</p>'
        f'<p style="margin:12px 0 0 0;padding:0;">{nps_pill}{bbb_pill}{extra_pill}</p>'
        f'</div>',
        unsafe_allow_html=True
    )

    if _is_est:
        _prov_bg = "#FFFBEB"; _prov_border = "#FDE68A"; _prov_color = "#92710C"
        _prov_icon = "&#9888;&#xFE0E;"
        _prov_label = "Class-Based Inference"
        _cls = data.get('estimated_class', 'General')
        _prov_body = (
            f"Pathway, metabolite, and scoring data for <b>{compound_name}</b> are <b>class-based estimates</b> "
            f"derived from published neuroprotection literature for the <b>{_cls}</b> compound class (MetaCyc, BRENDA, PubMed/PMC). "
            f"Individual compound-level experimental data have not been separately curated for this entry. "
            f"Results should be interpreted as indicative rather than definitive. "
            f"Cross-reference with KEGG, PubChem, or ChEMBL for compound-specific validation."
        )
    elif _is_ml:
        _prov_bg = "#EFF6FF"; _prov_border = "#93C5FD"; _prov_color = "#1E3A8A"
        _prov_icon = "&#9881;"
        _chembl_id = data.get("chembl_id", "")
        _r2 = data.get("model_cv_r2", 0.0)
        _live_ml2 = bool(data.get("_live_ml"))
        if _live_ml2:
            _prov_label = "Live PubChem + ChEMBL + ML Prediction"
            _prov_body = (
                f"<b>{compound_name}</b> was queried live from PubChem and ChEMBL in real time. "
                f"Its physicochemical properties (MW, XLogP, TPSA, HBD) were used to compute the "
                f"CNS-MPO-based BBB level. Mechanism-of-action data from ChEMBL was parsed to infer "
                f"disease relevance features. A <b>MultiOutput Random Forest model</b> trained on 129 "
                f"curated BrainSafe AI compounds (5-fold CV R\u00b2 = {_r2:.2f}) then predicted the "
                f"7-dimension neuroprotective profile. "
                + (f"PubChem CID: <b>{data.get('pubchem_cid','')}</b>. " if data.get("pubchem_cid") else "") +
                (f"ChEMBL ID: <b>{_chembl_id}</b>. " if _chembl_id else "") +
                f"This is a <b>live ML estimate</b> — validate findings in PubMed or ChEMBL."
            )
        else:
            _prov_label = "ChEMBL Drug Indication + ML-Predicted Data"
            _prov_body = (
                f"7-dimension bioactivity scores for <b>{compound_name}</b> were <b>predicted by a MultiOutput "
                f"Random Forest model</b> trained on 129 manually curated BrainSafe AI compounds "
                f"(5-fold CV R\u00b2 = {_r2:.2f}). This compound was identified from ChEMBL drug "
                f"indication data for neurological conditions (ALS, Alzheimer's, Parkinson's, Huntington's). "
                f"Enzyme/receptor data are <b>real IC50/Ki values from ChEMBL bioassays</b> (Homo sapiens). "
                + (f"ChEMBL molecule ID: <b>{_chembl_id}</b>. " if _chembl_id else "") +
                f"Use as a computational hypothesis-generation profile; validate in primary databases."
            )
    else:
        _prov_bg = "#F0FDF4"; _prov_border = "#86EFAC"; _prov_color = "#14532D"
        _prov_icon = "&#10003;"
        _prov_label = "Literature-Curated Data"
        _sources = data.get('sources', ['PubMed/PMC', 'ChEMBL', 'DrugBank'])
        _src_str = ', '.join(_sources[:4]) if isinstance(_sources, list) else str(_sources)
        _prov_body = (
            f"All pathways, metabolites, biomarkers, and dimension scores for <b>{compound_name}</b> are "
            f"<b>manually curated from peer-reviewed literature</b> (PubMed/PMC 2020–2026) and established "
            f"compound databases (ChEMBL, DrugBank, BRENDA, MetaCyc, UniProt). "
            f"Each data point reflects experimentally reported or clinically documented evidence. "
            f"This entry is <b>not a computational prediction</b>."
        )
    st.markdown(
        f'<div style="background:{_prov_bg};border:1px solid {_prov_border};border-radius:8px;'
        f'padding:10px 16px;margin:12px 0 8px 0;display:flex;gap:12px;align-items:flex-start;">'
        f'<span style="font-size:1.1rem;color:{_prov_color};flex-shrink:0;">{_prov_icon}</span>'
        f'<span style="font-size:0.83rem;color:{_prov_color};line-height:1.6;">'
        f'<strong>{_prov_label}:</strong> {_prov_body}</span></div>',
        unsafe_allow_html=True
    )

    render_external_data(compound_name)

    col_left, col_right = st.columns([1, 1], gap="large")

    with col_left:
        with st.container(border=True):
            st.markdown('<div class="section-label">Brain Health Radar — 7 Dimensions</div>', unsafe_allow_html=True)
            st.plotly_chart(make_radar(data), use_container_width=True, config={'displayModeBar': False})

    with col_right:
        with st.container(border=True):
            st.markdown('<div class="section-label">Bioactivity Profile</div>', unsafe_allow_html=True)
            metrics = [
                ("Antioxidant Activity", data['antioxidant']),
                ("Anti-Inflammatory", data['anti_inflammatory']),
                ("Mitochondrial Support", data['mitochondrial_support']),
                ("Aggregation Modulation", data['aggregation_modulation']),
                ("Cognitive Enhancement", data['cognitive_enhancement']),
                ("Neurogenesis Support", data['neurogenesis']),
                ("Synaptic Plasticity", data['synaptic_plasticity']),
            ]
            rows = ""
            for label, val in metrics:
                pct = val * 10
                c = "#1B6B45" if pct >= 70 else ("#9B5C00" if pct >= 40 else "#9B2335")
                rows += f"""<tr>
                    <td style="width:50%;font-weight:500;">{label}</td>
                    <td style="width:35%;">
                        <div style="background:#EEF1F6;border-radius:4px;height:9px;">
                          <div style="width:{pct}%;background:{c};height:9px;border-radius:4px;"></div>
                        </div>
                    </td>
                    <td style="width:15%;text-align:right;font-weight:700;color:{c};">{val}</td>
                </tr>"""
            st.markdown(
                f'<table><thead><tr><th>Dimension</th><th>Score</th><th style="text-align:right;">/10</th></tr>'
                f'</thead><tbody>{rows}</tbody></table>',
                unsafe_allow_html=True
            )

    with st.container(border=True):
        if _is_est:
            _net_badge = (
                '<span style="background:#FEF3C7;color:#92710C;font-size:0.68rem;font-weight:700;'
                'padding:2px 7px;border-radius:8px;vertical-align:middle;margin-left:6px;">Class-Inferred Edges</span>'
            )
        elif _is_ml:
            _net_badge = (
                '<span style="background:#DBEAFE;color:#1D4ED8;font-size:0.68rem;font-weight:700;'
                'padding:2px 7px;border-radius:8px;vertical-align:middle;margin-left:6px;">ML-Predicted Edges</span>'
            )
        else:
            _net_badge = (
                '<span style="background:#DCFCE7;color:#14532D;font-size:0.68rem;font-weight:700;'
                'padding:2px 7px;border-radius:8px;vertical-align:middle;margin-left:6px;">Literature-Curated Edges</span>'
            )
        st.markdown(
            f'<div class="section-label">Pathway and Metabolite Network — Interactive{_net_badge}</div>',
            unsafe_allow_html=True
        )
        st.markdown(
            '<p style="font-size:0.82rem;color:#8896AD;margin:0 0 10px 0;">'
            'Concentric graph showing mechanistic associations: <b>Compound → Pathways</b> (inner ring, curated from MetaCyc/KEGG), '
            '<b>Pathways → Metabolites/Biomarkers</b> (middle ring, from HMDB/MetaCyc), and '
            '<b>Metabolites → Disease relevance</b> (outer ring, from PubMed/PMC). '
            'Edges represent documented biological associations — not algorithmic predictions. '
            'Hover any node for details. Scroll to zoom; drag to pan.</p>',
            unsafe_allow_html=True
        )
        st.plotly_chart(make_network(compound_name, data), use_container_width=True, config={
            'scrollZoom': True,
            'displayModeBar': True,
            'modeBarButtonsToKeep': ['zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d'],
            'displaylogo': False,
        })

    enzyme_entry = get_enzyme_entry(compound_name, data.get('compound_type', ''))

    if enzyme_entry and (enzyme_entry.get('enzymes') or enzyme_entry.get('transporters')):
        is_enz_estimated = enzyme_entry.get('is_estimated', False)
        with st.container(border=True):
            est_badge = (
                ' <span style="background:#FFF3CD;color:#7A5500;font-size:0.7rem;font-weight:700;'
                'padding:2px 8px;border-radius:8px;vertical-align:middle;margin-left:8px;">CLASS ESTIMATE</span>'
                if is_enz_estimated else ''
            )
            st.markdown(
                f'<div class="section-label">Enzyme, Transporter &amp; Cofactor Profile — Molecular Pharmacology{est_badge}</div>',
                unsafe_allow_html=True
            )
            if is_enz_estimated:
                st.markdown(
                    '<p style="font-size:0.82rem;color:#92710C;background:#FFFBEB;border:1px solid #FDE68A;'
                    'border-radius:6px;padding:8px 12px;margin:0 0 14px 0;">'
                    'Compound-specific kinetic data is not yet in the curated database. The enzyme interactions, '
                    'transporters, and cofactors below are <b>class-based estimates</b> derived from published '
                    'literature for this compound type. Consult BRENDA, MetaCyc, or PubMed for compound-level validation.</p>',
                    unsafe_allow_html=True
                )
            elif _is_ml and enzyme_entry.get('sources'):
                st.markdown(
                    '<p style="font-size:0.82rem;color:#1E40AF;background:#EFF6FF;border:1px solid #93C5FD;'
                    'border-radius:6px;padding:8px 12px;margin:0 0 14px 0;">'
                    'Enzyme interactions and kinetic data (IC50/Ki) below are <b>real experimental values from '
                    'ChEMBL bioassays</b> (Homo sapiens, Mendez et al. 2019). '
                    'Target selection reflects bioassay availability for this compound in ChEMBL.</p>',
                    unsafe_allow_html=True
                )
            else:
                st.markdown(
                    '<p style="font-size:0.82rem;color:#8896AD;margin:0 0 14px 0;">'
                    'Curated enzyme interactions, transporter pharmacology, and cofactor involvement sourced from '
                    'BRENDA, MetaCyc, UniProt, DrugBank, ChEMBL, and peer-reviewed PubMed literature.</p>',
                    unsafe_allow_html=True
                )
            ecol, tcol = st.columns([1.1, 1], gap="large")

            with ecol:
                st.markdown('<div style="font-size:0.73rem;font-weight:700;color:#64748B;letter-spacing:1.2px;text-transform:uppercase;margin-bottom:8px;">Enzyme Interactions</div>', unsafe_allow_html=True)
                if enzyme_entry.get('enzymes'):
                    action_colors = {
                        'Inhibition': ('#7F1D1D', '#FEF2F2'),
                        'Activation': ('#14532D', '#F0FDF4'),
                        'Substrate': ('#1E3A5F', '#EFF6FF'),
                        'Cofactor': ('#44337A', '#F5F3FF'),
                        'Allosteric modulation': ('#7C3009', '#FFF7ED'),
                        'Induction': ('#14532D', '#F0FDF4'),
                        'Upregulation': ('#14532D', '#F0FDF4'),
                        'Modulation': ('#374151', '#F9FAFB'),
                    }
                    strength_icons = {'Strong': '■■■', 'Moderate': '■■□', 'Weak': '■□□'}
                    rows_e = ""
                    for enz in enzyme_entry['enzymes']:
                        action = enz.get('action', '').split('(')[0].strip()
                        action_key = next((k for k in action_colors if k.lower() in action.lower()), 'Inhibition')
                        fg, bg = action_colors.get(action_key, ('#1E3A5F', '#EFF6FF'))
                        strength = enz.get('strength', 'Moderate')
                        icon = strength_icons.get(strength, '■■□')
                        ec = f" <span style='color:#94A3B8;font-size:0.75rem;'>({enz['ec']})</span>" if enz.get('ec') else ""
                        src = f"<div style='color:#94A3B8;font-size:0.73rem;margin-top:1px;'>{enz.get('source','')}</div>" if enz.get('source') else ""
                        note = f"<div style='color:#475569;font-size:0.77rem;margin-top:2px;'>{enz.get('note','')}</div>" if enz.get('note') else ""
                        rows_e += (
                            f"<tr>"
                            f"<td style='width:28%;font-weight:600;font-size:0.84rem;padding:7px 6px;vertical-align:top;'>{enz['name']}{ec}</td>"
                            f"<td style='width:18%;padding:7px 6px;vertical-align:top;'>"
                            f"<span style='background:{bg};color:{fg};padding:2px 8px;border-radius:10px;font-size:0.75rem;font-weight:700;'>{action}</span>"
                            f"</td>"
                            f"<td style='width:12%;padding:7px 6px;vertical-align:top;color:#F0A500;font-size:0.75rem;font-weight:700;'>{icon}</td>"
                            f"<td style='width:42%;padding:7px 6px;vertical-align:top;'>{note}{src}</td>"
                            f"</tr>"
                        )
                    st.markdown(
                        f'<table style="width:100%;border-collapse:collapse;">'
                        f'<thead><tr style="border-bottom:2px solid #E2E8F0;">'
                        f'<th style="text-align:left;padding:5px 6px;font-size:0.72rem;color:#94A3B8;">Enzyme / Target</th>'
                        f'<th style="text-align:left;padding:5px 6px;font-size:0.72rem;color:#94A3B8;">Action</th>'
                        f'<th style="text-align:left;padding:5px 6px;font-size:0.72rem;color:#94A3B8;">Evidence</th>'
                        f'<th style="text-align:left;padding:5px 6px;font-size:0.72rem;color:#94A3B8;">Notes</th>'
                        f'</tr></thead><tbody>{rows_e}</tbody></table>',
                        unsafe_allow_html=True
                    )

            with tcol:
                if enzyme_entry.get('transporters'):
                    st.markdown('<div style="font-size:0.73rem;font-weight:700;color:#64748B;letter-spacing:1.2px;text-transform:uppercase;margin-bottom:8px;">Transporter Pharmacology</div>', unsafe_allow_html=True)
                    rows_t = ""
                    for tr in enzyme_entry['transporters']:
                        role = tr.get('role', '')
                        loc = tr.get('location', '')
                        note = tr.get('note', '')
                        rows_t += (
                            f"<tr>"
                            f"<td style='width:32%;font-weight:600;font-size:0.84rem;padding:7px 6px;vertical-align:top;'>{tr['name']}</td>"
                            f"<td style='width:22%;padding:7px 6px;vertical-align:top;color:#1E3A5F;font-size:0.82rem;font-weight:600;'>{role}</td>"
                            f"<td style='width:22%;padding:7px 6px;vertical-align:top;color:#475569;font-size:0.78rem;'>{loc}</td>"
                            f"<td style='width:24%;padding:7px 6px;vertical-align:top;color:#475569;font-size:0.77rem;'>{note}</td>"
                            f"</tr>"
                        )
                    st.markdown(
                        f'<table style="width:100%;border-collapse:collapse;">'
                        f'<thead><tr style="border-bottom:2px solid #E2E8F0;">'
                        f'<th style="text-align:left;padding:5px 6px;font-size:0.72rem;color:#94A3B8;">Transporter</th>'
                        f'<th style="text-align:left;padding:5px 6px;font-size:0.72rem;color:#94A3B8;">Role</th>'
                        f'<th style="text-align:left;padding:5px 6px;font-size:0.72rem;color:#94A3B8;">Location</th>'
                        f'<th style="text-align:left;padding:5px 6px;font-size:0.72rem;color:#94A3B8;">Note</th>'
                        f'</tr></thead><tbody>{rows_t}</tbody></table>',
                        unsafe_allow_html=True
                    )

                if enzyme_entry.get('cofactors'):
                    st.markdown('<div style="font-size:0.73rem;font-weight:700;color:#64748B;letter-spacing:1.2px;text-transform:uppercase;margin:16px 0 8px 0;">Cofactors &amp; Key Molecules</div>', unsafe_allow_html=True)
                    cf_html = "".join(
                        f'<span style="display:inline-block;background:#F0F9FF;color:#0C4A6E;border:1px solid #BAE6FD;'
                        f'border-radius:12px;padding:3px 12px;font-size:0.79rem;font-weight:600;margin:3px 4px 3px 0;">{cf}</span>'
                        for cf in enzyme_entry['cofactors']
                    )
                    st.markdown(cf_html, unsafe_allow_html=True)

    col_a, col_b = st.columns([1, 1], gap="large")

    with col_a:
        with st.container(border=True):
            if _is_est:
                _src_badge = (
                    '<span style="background:#FEF3C7;color:#92710C;font-size:0.68rem;font-weight:700;'
                    'padding:2px 7px;border-radius:8px;vertical-align:middle;margin-left:6px;">Class-Inferred</span>'
                )
                _path_desc = "Class-representative pathways from MetaCyc / KEGG / PubMed for this compound type."
            elif _is_ml:
                _src_badge = (
                    '<span style="background:#DBEAFE;color:#1D4ED8;font-size:0.68rem;font-weight:700;'
                    'padding:2px 7px;border-radius:8px;vertical-align:middle;margin-left:6px;">ML-Inferred</span>'
                )
                _path_desc = "Pathways inferred from ChEMBL bioassay target associations and disease indication data."
            else:
                _src_badge = (
                    '<span style="background:#DCFCE7;color:#14532D;font-size:0.68rem;font-weight:700;'
                    'padding:2px 7px;border-radius:8px;vertical-align:middle;margin-left:6px;">Literature-Curated</span>'
                )
                _path_desc = "Experimentally documented mechanistic pathways from MetaCyc, KEGG, and PubMed/PMC."
            st.markdown(
                f'<div class="section-label">Molecular Pathways{_src_badge}</div>'
                f'<p style="font-size:0.78rem;color:#8896AD;margin:2px 0 8px 0;">{_path_desc}</p>',
                unsafe_allow_html=True
            )
            pw_html = "".join(f'<span class="pathway-chip">{pw}</span>' for pw in data.get('pathways', []))
            st.markdown(pw_html, unsafe_allow_html=True)

            with st.spinner("Fetching live KEGG pathway links..."):
                kegg_data = fetch_kegg_pathways(compound_name)

            if "error" not in kegg_data and kegg_data.get("pathways"):
                kid      = kegg_data["kegg_id"]
                kname    = kegg_data.get("kegg_name", compound_name)
                kurl     = kegg_data["kegg_url"]
                kpaths   = kegg_data["pathways"]
                hsa_pws  = [p for p in kpaths if p["category"] == "Human (HSA)"]
                map_pws  = [p for p in kpaths if p["category"] == "Reference Map"]

                st.markdown(
                    f'<div style="margin:10px 0 4px 0;font-size:0.73rem;font-weight:700;color:#64748B;'
                    f'letter-spacing:1.1px;text-transform:uppercase;">Live KEGG Pathway Associations</div>'
                    f'<p style="font-size:0.75rem;color:#94A3B8;margin:0 0 6px 0;">'
                    f'KEGG compound: <a href="{kurl}" target="_blank" style="color:#1A6FBD;">{kid} — {kname}</a> '
                    f'&nbsp;&bull;&nbsp; '
                    f'<a href="https://www.kegg.jp/kegg-bin/show_compound?map=compound&query={kid}" '
                    f'target="_blank" style="color:#1A6FBD;">View in KEGG</a>'
                    f'</p>',
                    unsafe_allow_html=True
                )
                if hsa_pws:
                    links = " &nbsp; ".join(
                        f'<a href="{p["url"]}" target="_blank" '
                        f'style="display:inline-block;background:#EFF6FF;color:#1D4ED8;font-size:0.73rem;'
                        f'font-weight:600;padding:3px 9px;border-radius:20px;border:1px solid #BFDBFE;'
                        f'margin:2px 2px;text-decoration:none;">'
                        f'{p["name"]} <span style="opacity:0.6;font-size:0.65rem;">({p["id"]})</span></a>'
                        for p in hsa_pws
                    )
                    st.markdown(
                        f'<div style="margin-bottom:4px;"><span style="font-size:0.7rem;color:#1D4ED8;'
                        f'font-weight:700;margin-right:6px;">Human pathways:</span>{links}</div>',
                        unsafe_allow_html=True
                    )
                if map_pws:
                    links2 = " &nbsp; ".join(
                        f'<a href="{p["url"]}" target="_blank" '
                        f'style="display:inline-block;background:#F0FDF4;color:#15803D;font-size:0.73rem;'
                        f'font-weight:600;padding:3px 9px;border-radius:20px;border:1px solid #BBF7D0;'
                        f'margin:2px 2px;text-decoration:none;">'
                        f'{p["name"]} <span style="opacity:0.6;font-size:0.65rem;">({p["id"]})</span></a>'
                        for p in map_pws
                    )
                    st.markdown(
                        f'<div style="margin-bottom:6px;"><span style="font-size:0.7rem;color:#15803D;'
                        f'font-weight:700;margin-right:6px;">Reference maps:</span>{links2}</div>',
                        unsafe_allow_html=True
                    )
                st.markdown(
                    '<p style="font-size:0.68rem;color:#94A3B8;margin:2px 0 0 0;">'
                    'Source: KEGG REST API (rest.kegg.jp). Live data fetched this session. '
                    'Ref: Kanehisa M &amp; Goto S (2000). Nucleic Acids Res 28(1):27–30.</p>',
                    unsafe_allow_html=True
                )
            elif "error" in kegg_data:
                st.caption(f"KEGG: {kegg_data['error']}")

            st.markdown(
                f'<div class="section-label" style="margin-top:18px;">Key Metabolites and Biomarkers{_src_badge}</div>'
                f'<p style="font-size:0.78rem;color:#8896AD;margin:2px 0 8px 0;">'
                f'{"Typical metabolic markers for this compound class (HMDB, MetaCyc)." if _is_est else "Experimentally linked metabolites and biomarkers from HMDB, MetaCyc, and peer-reviewed studies."}'
                f'</p>',
                unsafe_allow_html=True
            )
            met_html = ""
            for met in data.get('metabolites', []):
                lbl = METABOLITE_LABELS.get(met, met)
                met_html += f'<span class="met-chip" title="{lbl}">{met}</span>'
            st.markdown(met_html, unsafe_allow_html=True)

            st.markdown('<div class="section-label" style="margin-top:18px;">Brain Regions Affected</div>', unsafe_allow_html=True)
            br_html = "".join(f'<span class="br-chip">{br}</span>' for br in data.get('brain_regions', []))
            st.markdown(br_html, unsafe_allow_html=True)

    with col_b:
        with st.container(border=True):
            st.markdown('<div class="section-label">Neurodegenerative Disease Relevance</div>', unsafe_allow_html=True)
            disease_data = [
                ("ALS", data.get('als', 'Low')),
                ("Alzheimer's Disease", data.get('alzheimers', 'Low')),
                ("Parkinson's Disease", data.get('parkinsons', 'Low')),
                ("Huntington's Disease", data.get('huntingtons', 'Low')),
            ]
            rows_d = ""
            for dname, dlevel in disease_data:
                fg, bg, lbl, desc = disease_style(dlevel)
                rows_d += (
                    f"<tr><td style='font-weight:500;width:38%;'>{dname}</td>"
                    f"<td style='width:20%;'><span class='pill' style='color:{fg};background:{bg};'>{lbl}</span></td>"
                    f"<td style='color:#64748B;font-size:0.83rem;'>{desc}</td></tr>"
                )
            st.markdown(
                f'<table><thead><tr><th>Condition</th><th>Evidence</th><th>Assessment</th></tr></thead>'
                f'<tbody>{rows_d}</tbody></table>',
                unsafe_allow_html=True
            )

    with st.container(border=True):
        st.markdown('<div class="section-label">Overall Brain Health Assessment</div>', unsafe_allow_html=True)
        assessment_paras = generate_assessment(compound_name, data, score, enzyme_entry)
        for para in assessment_paras:
            if para.startswith("Primary data sources") or para.startswith("Data sources"):
                st.markdown(
                    f'<p style="font-size:0.77rem;color:#94A3B8;margin:6px 0 0 0;border-top:1px solid #F0F4FA;padding-top:8px;">{para}</p>',
                    unsafe_allow_html=True
                )
            else:
                st.markdown(
                    f'<p style="font-size:0.89rem;color:#1A2B45;line-height:1.75;margin:0 0 10px 0;">{para}</p>',
                    unsafe_allow_html=True
                )

    st.markdown(
        '<div class="disclaimer"><b>Educational Use Only.</b> '
        'This report integrates data from published neuroscience literature (PubMed/PMC), BRENDA enzyme database, '
        'MetaCyc metabolic pathway database, UniProt protein database, ChEMBL, DrugBank, and HMDB '
        'for educational and research purposes only. '
        'It does not constitute medical advice. Consult a qualified healthcare professional before '
        'using any supplement, nutraceutical, or pharmaceutical compound. '
        'This tool is part of the BrainSafe AI module of SAI-Net — Sri Sathya Sai Institute of Higher Learning.</div>',
        unsafe_allow_html=True
    )


_logo_html = (
    f'<img src="data:image/png;base64,{LOGO_B64}" class="header-logo-img">'
    if LOGO_B64 else
    '<div class="header-logo-img" style="display:flex;align-items:center;justify-content:center;background:#1A3A5C;font-size:1.6rem;color:#F0A500;font-weight:800;">B</div>'
)

_sssihl_logo_html = (
    f'<img src="data:image/png;base64,{SSSIHL_LOGO_B64}" class="header-inst-logo">'
    if SSSIHL_LOGO_B64 else ""
)

st.markdown(f"""
<div class="site-header">
    <div class="header-logo-wrap">
        {_logo_html}
        <div class="header-divider"></div>
        <div class="header-text-block">
            <div class="header-title">Brain<span>Safe</span> AI</div>
            <div class="header-sub">Brain Health Compound Explorer: a SAI-Net Translational Module</div>
            <div class="header-tags">
                <span class="header-tag-pill">Flavonoids</span>
                <span class="header-tag-pill">Vitamins</span>
                <span class="header-tag-pill">Supplements</span>
                <span class="header-tag-pill">Nutraceuticals</span>
                <span class="header-tag-pill">Approved Drugs</span>
                <span class="header-tag-pill">BBB Penetration</span>
                <span class="header-tag-pill">GSH / NAD+ / ATP Pathways</span>
                <span class="header-tag-pill">ALS</span>
                <span class="header-tag-pill">Alzheimer's</span>
                <span class="header-tag-pill">Parkinson's</span>
                <span class="header-tag-pill">Huntington's</span>
            </div>
        </div>
        <div class="header-divider"></div>
        {_sssihl_logo_html}
    </div>
</div>
""", unsafe_allow_html=True)

tab_search, tab_about = st.tabs(["Compound Search", "About SAI-Net"])

with tab_search:
    st.markdown('<div style="height:10px;"></div>', unsafe_allow_html=True)
    with st.container(border=True):
        st.markdown(
            '<p style="font-size:1.05rem;font-weight:700;color:#0D2137;margin:0 0 3px 0;">Search a Compound</p>'
            '<p style="font-size:0.84rem;color:#4A5568;margin:0 0 14px 0;">'
            'Select from the dropdown (325 compounds) or type any name below. '
            'Compounds not in the database are first looked up live in PubChem + ChEMBL '
            'and profiled by the BrainSafe AI Random Forest model. '
            'Supports abbreviations: EGCG, NAC, CoQ10, BHB, ALA, etc.</p>',
            unsafe_allow_html=True
        )
        sel_col, btn_col = st.columns([4, 1])
        with sel_col:
            selected = st.selectbox("Select compound", [""] + COMPOUND_NAMES, label_visibility="collapsed")
        with btn_col:
            search_btn = st.button("View Report", use_container_width=True, type="primary")

        st.markdown('<div style="height:6px;"></div>', unsafe_allow_html=True)

        free_col, free_btn_col = st.columns([4, 1])
        with free_col:
            free_text = st.text_input(
                "Free text",
                placeholder="Type any compound — e.g. curcumin, apigenin, leucine, diosmin, myristicin...",
                label_visibility="collapsed"
            )
        with free_btn_col:
            free_btn = st.button("Search / Estimate", use_container_width=True)

    if search_btn and selected:
        if selected in COMPOUNDS:
            render_report(selected, COMPOUNDS[selected])

    elif free_btn and free_text:
        matches = fuzzy_match(free_text, COMPOUND_NAMES, cutoff=0.80)
        if matches:
            best = matches[0]
            others = [m for m in matches[1:] if m != best]
            if others:
                st.info(f"Showing best match: **{best}**. Other close matches: {', '.join(others)}")
            render_report(best, COMPOUNDS[best])
        else:
            close = difflib.get_close_matches(free_text.lower(), [n.lower() for n in COMPOUND_NAMES], n=3, cutoff=0.2)
            suggestions = [n for c in close for n in COMPOUND_NAMES if n.lower() == c]
            if suggestions:
                st.info(
                    f"No exact match for **'{free_text}'** in the database. "
                    f"Closest entries: {', '.join(suggestions)}. "
                    f"Querying PubChem + ChEMBL for a live ML-estimated profile..."
                )
            with st.spinner("Querying PubChem and ChEMBL for ML prediction..."):
                ml_entry = predict_unknown_via_ml(free_text.strip())
            if ml_entry:
                st.success(
                    f"**'{free_text}'** identified in PubChem (CID {ml_entry.get('pubchem_cid', '?')}). "
                    f"7-dimension profile predicted by the BrainSafe AI Random Forest model "
                    f"using live physicochemical data and ChEMBL mechanism-of-action."
                )
                render_report(free_text.strip().title(), ml_entry)
            else:
                st.info(
                    f"**'{free_text}'** was not found in PubChem or ChEMBL. "
                    f"Generating a class-based scientific estimate from compound name analysis."
                )
                estimated = generate_estimated_data(free_text.strip())
                render_report(free_text.strip().title(), estimated)


# with tab_browse:
#     st.markdown('<div style="height:10px;"></div>', unsafe_allow_html=True)
# 
#     n_cur = DB_STATS["curated"]
#     n_ml  = DB_STATS["ml"]
#     n_tot = DB_STATS["total"]
#     bcol1, bcol2, bcol3 = st.columns(3)
#     bcol1.metric("Literature-Curated", n_cur, help="Manually curated from PubMed/PMC and primary databases")
#     bcol2.metric("ML-Predicted (ChEMBL)", n_ml,
#                  help="Profiles generated by Random Forest trained on curated compounds; ERC from ChEMBL bioassays" if n_ml else "Run python3 ml_expander.py to generate")
#     bcol3.metric("Total Compounds", n_tot)
# 
#     if n_ml == 0:
#         st.info(
#             "The ML-expanded database has not been generated yet. "
#             "The Random Forest expansion pipeline generates profiles for 400–700 additional "
#             "ChEMBL-indicated neuroprotection compounds. "
#             "Run `python3 ml_expander.py` in the terminal to generate (takes ~3 minutes)."
#         )
# 
#     source_filter = st.radio(
#         "Show compounds:",
#         ["All", "Literature-Curated only", "ML-Predicted only"],
#         horizontal=True,
#         label_visibility="visible",
#     )
# 
#     type_groups_browse: dict[str, list[str]] = {}
#     for name in COMPOUND_NAMES:
#         d = COMPOUNDS[name]
#         is_ml_entry = bool(d.get("ml_predicted"))
#         if source_filter == "Literature-Curated only" and is_ml_entry:
#             continue
#         if source_filter == "ML-Predicted only" and not is_ml_entry:
#             continue
#         ctype   = d.get('compound_type', 'Other')
#         primary = ctype.split('/')[0].strip()
#         type_groups_browse.setdefault(primary, []).append(name)
# 
#     for group_name in sorted(type_groups_browse.keys()):
#         items = sorted(type_groups_browse[group_name])
#         with st.expander(f"{group_name}  ({len(items)} compounds)"):
#             cols = st.columns(3)
#             for i, name in enumerate(items):
#                 d = COMPOUNDS[name]
#                 score     = neuro_score(d)
#                 fg, bg    = score_color_hex(score)
#                 bbb       = d.get('bbb', '')
#                 bb_fg, bb_bg = bbb_color_hex(bbb)
#                 is_ml_c   = bool(d.get("ml_predicted"))
#                 src_badge = (
#                     '<span style="background:#DBEAFE;color:#1D4ED8;font-size:0.65rem;font-weight:700;'
#                     'padding:1px 6px;border-radius:8px;margin-left:4px;">ML</span>'
#                     if is_ml_c else
#                     '<span style="background:#DCFCE7;color:#14532D;font-size:0.65rem;font-weight:700;'
#                     'padding:1px 6px;border-radius:8px;margin-left:4px;">Curated</span>'
#                 )
#                 cols[i % 3].markdown(
#                     f'<div style="padding:6px 2px;border-bottom:1px solid #F0F3F8;">'
#                     f'<span style="font-weight:700;color:#0D2137;font-size:0.9rem;">{name}</span>{src_badge}<br>'
#                     f'<span class="pill" style="color:{fg};background:{bg};font-size:0.72rem;padding:2px 8px;">Score: {score}</span>'
#                     f'<span class="pill" style="color:{bb_fg};background:{bb_bg};font-size:0.72rem;padding:2px 8px;">BBB: {bbb}</span>'
#                     f'</div>',
#                     unsafe_allow_html=True
#                 )
# 
# 
with tab_about:
    st.markdown('<div style="height:10px;"></div>', unsafe_allow_html=True)
    col_logo, col_text = st.columns([1, 3], gap="large")
    with col_logo:
        if LOGO_B64:
            st.markdown(
                f'<img src="data:image/png;base64,{LOGO_B64}" '
                f'style="width:220px;border-radius:12px;box-shadow:0 4px 16px rgba(13,33,55,0.15);display:block;margin:0 auto;">',
                unsafe_allow_html=True
            )
        else:
            st.markdown(
                '<div style="width:220px;height:220px;background:#0D2137;border-radius:12px;'
                'display:flex;align-items:center;justify-content:center;'
                'color:#F0A500;font-size:2rem;font-weight:800;">B</div>',
                unsafe_allow_html=True
            )
    with col_text:
        st.markdown("""
        <div style="padding: 6px 0 20px 0;">
          <h2 style="font-family: 'Inter', sans-serif; font-size: 1.75rem; font-weight: 700;
                     color: #0D2137; margin: 0 0 16px 0; letter-spacing: -0.015em;">
            About Brain<span style="color: #F0A500;">Safe</span> AI
          </h2>
          <div style="border-left: 3px solid #F0A500; padding-left: 14px;">
            <p style="font-family: 'Inter', sans-serif; font-size: 0.98rem; font-weight: 700;
                      color: #B45309; letter-spacing: 0.06em; text-transform: uppercase;
                      margin: 0 0 6px 0;">
              Science for Society
            </p>
            <p style="font-family: 'Inter', sans-serif; font-size: 0.95rem; color: #475569;
                      margin: 0; line-height: 1.65; font-style: italic;">
              Bridging neuroscience with public health education.
            </p>
          </div>
        </div>
        """, unsafe_allow_html=True)

    with st.container(border=True):
        st.markdown("""
        ### SAI-Net Connection

        BrainSafe AI is a **translational educational module** derived from **SAI-Net (Structure-Activity Intelligence Network)** —
        a computational neuropharmacology framework for multiomics drug discovery in neurodegenerative diseases.

        SAI-Net integrates structure-activity relationships, network pharmacology, proteomics, transcriptomics,
        and metabolomics to identify candidate therapeutic compounds and mechanisms across ALS, Alzheimer's,
        Parkinson's, and Huntington's disease spectra.

        **Future Integration:** SAI-Net analytical APIs for real-time multiomics predictions, protein docking,
        and network pharmacology will be incorporated in subsequent versions of this platform.
        """)

    with st.container(border=True):
        st.markdown("""
        ### 100 Years of Selfless Service

        This offering celebrates **100 years of selfless service and unconditional love to our society**
        by **Bhagawan Sri Sathya Sai Baba** (1926 - 2011).
        """)
        st.markdown("""
        <div style="display:flex;flex-direction:column;gap:14px;margin-top:8px;">
          <div style="border-left:4px solid #F0A500;background:#FFFBF0;
                      border-radius:0 10px 10px 0;padding:16px 20px;
                      box-shadow:0 1px 6px rgba(240,165,0,0.10);">
            <p style="font-size:1.05rem;font-style:italic;color:#5A4000;
                      margin:0 0 8px 0;line-height:1.65;">
              &#8220;Let the world achieve the glory of becoming a family &#8212; through Love.&#8221;
            </p>
            <p style="font-size:0.78rem;font-weight:700;color:#B45309;
                      letter-spacing:0.06em;text-transform:uppercase;margin:0;">
              &mdash; Bhagawan Sri Sathya Sai Baba
            </p>
          </div>
          <div style="border-left:4px solid #1A3A5C;background:#F0F5FB;
                      border-radius:0 10px 10px 0;padding:16px 20px;
                      box-shadow:0 1px 6px rgba(13,33,55,0.07);">
            <p style="font-size:1.05rem;font-style:italic;color:#1A3A5C;
                      margin:0 0 8px 0;line-height:1.65;">
              &#8220;True knowledge is that which makes man work for the welfare of humanity.&#8221;
            </p>
            <p style="font-size:0.78rem;font-weight:700;color:#1A3A5C;
                      letter-spacing:0.06em;text-transform:uppercase;margin:0;opacity:0.7;">
              &mdash; Bhagawan Sri Sathya Sai Baba
            </p>
          </div>
        </div>
        """, unsafe_allow_html=True)

    with st.container(border=True):
        st.markdown("""
        <div style="margin-bottom:16px;">
          <p style="font-size:1.05rem;font-weight:700;color:#0D2137;margin:0 0 14px 0;
                    letter-spacing:-0.01em;">Research Team</p>
          <div style="display:flex;flex-direction:column;gap:10px;">
            <div style="display:flex;align-items:center;gap:14px;
                        background:#F7F9FD;border:1px solid #E4EAF4;
                        border-radius:10px;padding:14px 18px;">
              <div style="background:#0D2137;color:#F0A500;font-size:0.7rem;font-weight:800;
                          padding:4px 10px;border-radius:20px;letter-spacing:0.06em;
                          text-transform:uppercase;white-space:nowrap;flex-shrink:0;">
                Principal Investigator
              </div>
              <div>
                <p style="margin:0;font-size:0.95rem;font-weight:700;color:#0D2137;">
                  Prof. Venketesh Sivaramakrishnan
                </p>
                <p style="margin:0;font-size:0.80rem;color:#64748B;">
                  Sri Sathya Sai Institute of Higher Learning
                </p>
              </div>
            </div>
            <div style="display:flex;align-items:center;gap:14px;
                        background:#F7F9FD;border:1px solid #E4EAF4;
                        border-radius:10px;padding:14px 18px;">
              <div style="background:#1A3A5C;color:#F0A500;font-size:0.7rem;font-weight:800;
                          padding:4px 10px;border-radius:20px;letter-spacing:0.06em;
                          text-transform:uppercase;white-space:nowrap;flex-shrink:0;">
                Developer
              </div>
              <div>
                <p style="margin:0;font-size:0.95rem;font-weight:700;color:#0D2137;">
                  Krishnasalini Gunanathan
                </p>
                <p style="margin:0;font-size:0.80rem;color:#64748B;">
                  Sri Sathya Sai Institute of Higher Learning
                </p>
              </div>
            </div>
          </div>
        </div>
        <hr style="border:none;border-top:1px solid #E4EAF4;margin:16px 0;">
        <p style="font-size:0.87rem;color:#475569;margin:0 0 10px 0;line-height:1.7;">
          <strong>Educational Purpose Only</strong> — This tool is not intended as medical advice.
          Consult a qualified healthcare professional before using any compound discussed here.
        </p>
        <p style="font-size:0.87rem;font-weight:700;color:#0D2137;margin:0 0 6px 0;">
          Data Sources &amp; Acknowledgements
        </p>
        <p style="font-size:0.85rem;color:#475569;line-height:1.75;margin:0;">
          Compound bioactivity and pharmacology data were curated from <strong>ChEMBL</strong>
          (European Molecular Biology Laboratory) and <strong>DrugBank</strong>. Neuropharmacological
          literature was sourced from <strong>PubMed / PMC</strong> (2020–2026). Enzyme kinetics,
          reaction mechanisms, and EC classifications were drawn from <strong>BRENDA</strong>
          (Braunschweig Enzyme Database) and <strong>MetaCyc</strong> (Pathway/Genome Database).
          Protein sequences and transporter annotations were obtained from <strong>UniProt</strong>
          and the <strong>Human Metabolome Database (HMDB)</strong>. Blood-brain barrier and CNS
          penetration data were referenced from <strong>FDA prescribing information</strong> and
          published clinical literature. Neurodegenerative disease relevance data were
          cross-referenced with the <strong>ALS Association</strong>, <strong>Alzheimer's
          Association</strong>, <strong>Parkinson's Foundation</strong>, and <strong>HDSA</strong>
          (Huntington's Disease Society of America) research resources.<br><br>
          All curated data represent the state of scientific knowledge available at the time of tool
          development and should be independently verified before any research or clinical use.
        </p>
        """, unsafe_allow_html=True)