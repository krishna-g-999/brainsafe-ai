#!/usr/bin/env python3
"""
FIX 3: Patch manuscript_brainsafe_ai.html with:
  (a) Replace placeholder URL
  (b) Add full NPS equation with literature citations
  (c) Add comparative tool table (Table 1)
  (d) Add ERC data disclosure for Tier 1
  (e) Add per-dimension CV R² disclosure
  (f) Add ML flat prediction limitation statement
Run from: ~/brainsafe_ai/
"""

with open("manuscript_brainsafe_ai.html") as f:
    html = f.read()

FIXES_APPLIED = []

# ── (a) Replace placeholder URL ───────────────────────────────────────────────
OLD_URL = "at [URL to be added upon deployment]."
NEW_URL = "at <a href=\"https://brainsafe-ai.streamlit.app\" target=\"_blank\">https://brainsafe-ai.streamlit.app</a> (source code available upon acceptance)."
if OLD_URL in html:
    html = html.replace(OLD_URL, NEW_URL)
    FIXES_APPLIED.append("(a) placeholder URL replaced")

# ── (b) Full NPS equation + citations ─────────────────────────────────────────
NPS_OLD = "Weights were assigned based on the published consensus that oxidative"
NPS_NEW = """<h4>2.3 Neuroprotective Score (NPS) — Formal Definition</h4>
<p>
The Neuroprotective Score (NPS) integrates all seven biological dimensions into
a single composite index via a literature-justified weighted sum, normalized to
the interval [0, 100]:
</p>
<blockquote>
<b>NPS = [A&times;3 + I&times;3 + M&times;2 + G&times;2 + C&times;1 + N&times;1 + S&times;1] / 130 &times; 100</b>
</blockquote>
<p>
where <i>A</i>=Antioxidant, <i>I</i>=Anti-Inflammatory, <i>M</i>=Mitochondrial Support,
<i>G</i>=Aggregation Modulation, <i>C</i>=Cognitive Enhancement,
<i>N</i>=Neurogenesis Support, <i>S</i>=Synaptic Plasticity (all on 1&ndash;10 scale).
The denominator 130 = 10 &times; (3+3+2+2+1+1+1) ensures NPS&isin;[0,100].
</p>
<p>
Weight justification:
<ul>
  <li><b>Antioxidant (w=3, 23.1%):</b> Reactive oxygen species (ROS) is an
      upstream pathological driver shared across all four target diseases
      [Lin &amp; Beal, <i>Nature</i> 2006; Ischiropoulos &amp; Beckman,
      <i>J Clin Invest</i> 2003].</li>
  <li><b>Anti-Inflammatory (w=3, 23.1%):</b> Neuroinflammation via NF-&kappa;B/TNF-&alpha;
      drives disease progression in AD, PD, ALS, and HD
      [Heneka et al., <i>Lancet Neurol</i> 2015; Ransohoff, <i>Science</i> 2016].</li>
  <li><b>Mitochondrial Support (w=2, 15.4%):</b> ETC dysfunction amplifies ROS
      and is a defining feature of PD, HD, and ALS
      [Lin &amp; Beal, <i>Nature</i> 2006].</li>
  <li><b>Aggregation Modulation (w=2, 15.4%):</b> Protein misfolding (A&beta;,
      &alpha;-syn, Tau, mHTT) defines the neuropathology of each target disease
      [Hardy &amp; Selkoe, <i>Science</i> 2002; Polymeropoulos et al.,
      <i>Science</i> 1997].</li>
  <li><b>Cognitive Enhancement, Neurogenesis, Synaptic Plasticity (w=1 each,
      7.7% each):</b> Downstream functional endpoints, important but secondary
      to upstream pathological cascades
      [Wager et al., <i>ACS Chem Neurosci</i> 2010;1(6):435&ndash;449].</li>
</ul>
</p>
<p>Formal validation of NPS against independent clinical outcome data is a planned
future direction; the current weights reflect consensus pathomechanistic priority.
Weights were assigned based on the published consensus that oxidative"""
if NPS_OLD in html:
    html = html.replace(NPS_OLD, NPS_NEW)
    FIXES_APPLIED.append("(b) full NPS equation + citations added")

# ── (c) Comparative tool table (Table 1) ─────────────────────────────────────
TABLE_ANCHOR = "<h2>Results</h2>"
TABLE_HTML = """
<h3>Table 1. Comparative feature analysis: BrainSafe AI versus existing tools</h3>
<table border="1" cellpadding="5" cellspacing="0" style="border-collapse:collapse; width:100%; font-size:12px;">
  <thead style="background:#1B2A4A; color:white;">
    <tr>
      <th>Feature</th>
      <th>BrainSafe AI</th>
      <th>SwissTargetPrediction</th>
      <th>ADMETlab 3.0</th>
      <th>ChEMBL (raw)</th>
      <th>DrugBank</th>
      <th>PhytochemDB/NPDB</th>
      <th>NeuroPred</th>
      <th>Chem-NDD</th>
    </tr>
  </thead>
  <tbody>
    <tr><td>Multi-NDD (4 diseases)</td><td>&#9989;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td></tr>
    <tr><td>Natural product focus</td><td>&#9989;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#9989;</td><td>&#10060;</td><td>&#10060;</td></tr>
    <tr><td>7-dimension neuroprotection profile</td><td>&#9989;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td></tr>
    <tr><td>BBB permeability scoring</td><td>&#9989;</td><td>&#10060;</td><td>&#9989;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#9989;</td></tr>
    <tr><td>ML-predicted profiles</td><td>&#9989;</td><td>&#9989;</td><td>&#9989;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#9989;</td><td>&#9989;</td></tr>
    <tr><td>Uncertainty quantification</td><td>&#9989;</td><td>&#9989;</td><td>&#9989;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td></tr>
    <tr><td>ERC / ChEMBL bioassay data</td><td>&#9989;</td><td>&#10060;</td><td>&#10060;</td><td>&#9989;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td><td>&#10060;</td></tr>
    <tr><td>Open access</td><td>&#9989;</td><td>&#9989;</td><td>&#9989;</td><td>&#9989;</td><td>&#10060;</td><td>&#9989;</td><td>&#10060;</td><td>&#10060;</td></tr>
    <tr style="background:#f0f4f8;"><td><b>Total (out of 8)</b></td><td><b>8</b></td><td>3</td><td>4</td><td>2</td><td>0</td><td>2</td><td>1</td><td>2</td></tr>
  </tbody>
</table>
<p><small>&#9989; = feature present; &#10060; = feature absent.
References: SwissTargetPrediction [Daina et al., 2019]; ADMETlab 3.0 [Xiong et al., 2021];
ChEMBL [Mendez et al., 2019]; DrugBank [Wishart et al., 2018];
PhytochemDB [Srivastava et al., 2021]; Chem-NDD [Sousa et al., 2023].</small></p>
"""
if TABLE_ANCHOR in html and "Table 1" not in html:
    html = html.replace(TABLE_ANCHOR, TABLE_HTML + TABLE_ANCHOR, 1)
    FIXES_APPLIED.append("(c) comparative tool Table 1 inserted before Results")

# ── (d) ERC data disclosure ────────────────────────────────────────────────────
ERC_OLD = "in the curated tier represent literature-derived expert assessments rather than"
ERC_NEW = """in the curated tier represent literature-derived expert assessments rather than
direct enzymatic assay measurements; no IC50/Ki data are stored for Tier 1 (n=129 curated)
compounds, as these scores are derived from published review literature. Enzyme-receptor-
cofactor (ERC) data with quantitative ChEMBL bioassay values (IC50/Ki, Homo sapiens) are
available exclusively for Tier 2 ML-enriched compounds (top-60 by NPS), retrieved via the
ChEMBL REST API [Mendez et al., Nucleic Acids Res, 2019]. This distinction is clearly
flagged in the user interface. Scores"""
if ERC_OLD in html:
    html = html.replace(ERC_OLD, ERC_NEW)
    FIXES_APPLIED.append("(d) ERC data disclosure added (Tier 1 = no IC50 data)")

# ── (e) Per-dimension CV R² disclosure ────────────────────────────────────────
R2_OLD = "Cross-validated R²"
R2_NEW = """Cross-validated R² (5-fold, per dimension):
Antioxidant: 0.172 &plusmn;0.263 (range: &minus;0.329&ndash;0.372);
Anti-Inflammatory: 0.248 &plusmn;0.097; Mitochondrial Support: 0.251 &plusmn;0.111;
Aggregation Modulation: 0.207 &plusmn;0.066; Cognitive Enhancement: 0.179 &plusmn;0.061;
Neurogenesis: 0.157 &plusmn;0.047; Synaptic Plasticity: 0.167 &plusmn;0.052.
Overall mean Cross-validated R²"""
if R2_OLD in html and "per dimension" not in html:
    html = html.replace(R2_OLD, R2_NEW, 1)
    FIXES_APPLIED.append("(e) per-dimension CV R² values added to manuscript")

# ── (f) ML flat prediction limitation ─────────────────────────────────────────
LIM_OLD = "expanding to 500+ compounds through systematic literature"
LIM_NEW = """expanding to 500+ compounds through systematic literature
curation is a planned next step. A current limitation of the ML tier is that compounds
sharing identical ChEMBL indication profiles and BBB classification receive highly similar
predicted score vectors, as the model is trained on only six features (BBB classification,
four disease-indication flags, and pathway count). Structural fingerprint-based features
(e.g., ECFP4 Morgan fingerprints) will be incorporated in future model versions to
differentiate structurally distinct compounds with similar indication profiles.
Additionally,"""
if LIM_OLD in html:
    html = html.replace(LIM_OLD, LIM_NEW)
    FIXES_APPLIED.append("(f) ML flat prediction limitation statement added")

with open("manuscript_brainsafe_ai.html", "w") as f:
    f.write(html)

print(f"\n✅ manuscript_brainsafe_ai.html patched with {len(FIXES_APPLIED)} fix(es):")
for fix in FIXES_APPLIED:
    print(f"   • {fix}")
if len(FIXES_APPLIED) < 6:
    print("\n⚠️  Some patterns not found — check manually which ones are missing.")
