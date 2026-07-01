"""
model_config.py  –  BrainSafe AI v6
Central configuration: all constants, mappings, and lookup tables.

FIXES vs v5:
  1. Renamed _BBB_MAP / _DIS_MAP to BBB_MAP / DIS_MAP (removes NameError in
     build_ml_predictor where private-name versions were used but the public
     names were called).
  2. Defines TARGETTOPATHWAY (was referenced in predict_unknown_via_ml but
     never defined anywhere — caused a silent NameError on every live prediction).
  3. Defines DIMENSION_COLS as the canonical ordered list of the 7 dimensions.
  4. Adds SCORE_MAX = 100 so scale is enforced from one place.
"""

# ---------------------------------------------------------------------------
# Canonical dimension column names (7 mechanistic dimensions)
# ---------------------------------------------------------------------------
DIMENSION_COLS: list[str] = [
    "antioxidant",
    "anti_inflammatory",
    "mitochondrial_support",
    "aggregation_modulation",
    "cognitive_enhancement",
    "neurogenesis",
    "synaptic_plasticity",
]

SCORE_MAX: float = 100.0   # all dimension scores are on 0–100 scale
SCORE_MIN: float = 0.0

DIM_LABELS: dict[str, str] = {
    "antioxidant":            "Antioxidant",
    "anti_inflammatory":      "Anti-inflammatory",
    "mitochondrial_support":  "Mitochondrial Support",
    "aggregation_modulation": "Aggregation Modulation",
    "cognitive_enhancement":  "Cognitive Enhancement",
    "neurogenesis":           "Neurogenesis",
    "synaptic_plasticity":    "Synaptic Plasticity",
}

N_FEATURE_TOTAL: int = 93  # 50 (ECFP PCA) + 32 (ChemBERTa PCA) + 4 (physicochemical) + 1 (BBB)


# ---------------------------------------------------------------------------
# BBB permeability encoding (ordinal, for feature vector)
# ---------------------------------------------------------------------------
BBB_MAP: dict[str, int] = {
    "Low":     0,
    "Low-Med": 1,
    "Medium":  2,
    "High":    3,
}

BBB_LABELS: dict[int, str] = {v: k for k, v in BBB_MAP.items()}


# ---------------------------------------------------------------------------
# Disease relevance encoding (ordinal, for feature vector)
# ---------------------------------------------------------------------------
DIS_MAP: dict[str, int] = {
    "Low":  0,
    "Med":  1,
    "High": 2,
}

DIS_LABELS: dict[int, str] = {v: k for k, v in DIS_MAP.items()}


# ---------------------------------------------------------------------------
# Disease keyword matching (for live compound disease relevance prediction)
# ---------------------------------------------------------------------------
DIS_KW: dict[str, list[str]] = {
    "alzheimers": [
        "alzheimer", "cholinesterase", "acetylcholinesterase", "bace", "bace1",
        "amyloid", "tau", "gsk3", "app", "presenilin", "donepezil", "memantine",
        "galantamine", "rivastigmine", "aricept", "nmda antagonist",
    ],
    "parkinsons": [
        "parkinson", "dopaminerg", "levodopa", "alpha-synuclein", "synuclein",
        "mao-b", "selegiline", "rasagiline", "pramipexole", "ropinirole",
        "lrrk2", "pink1", "parkin", "dj-1", "substantia nigra",
    ],
    "als": [
        "amyotrophic", "als", "motor neuron", "sod1", "tdp-43", "fus",
        "c9orf72", "riluzole", "edaravone", "excitotoxic", "superoxide dismutase",
    ],
    "huntingtons": [
        "huntington", "striatum", "hdac", "tetrabenazine", "mutant htt",
        "polyglutamine", "pde10a", "sigma-1", "mtor inhibit",
    ],
}


# ---------------------------------------------------------------------------
# TARGET → PATHWAY mapping
# FIX: this was referenced in predict_unknown_via_ml() but never defined,
# causing a NameError on every live compound prediction.
# ---------------------------------------------------------------------------
TARGETTOPATHWAY: dict[str, str] = {
    # AD targets
    "acetylcholinesterase":    "cholinergic",
    "acetylcholine esterase":  "cholinergic",
    "ache":                    "cholinergic",
    "bace1":                   "amyloid",
    "bace-1":                  "amyloid",
    "beta secretase":          "amyloid",
    "gsk3":                    "tau",
    "gsk-3":                   "tau",
    "glycogen synthase kinase": "tau",
    "tau":                     "tau",
    "amyloid":                 "amyloid",
    "app":                     "amyloid",
    # PD targets
    "mao-b":                   "dopaminergic",
    "maob":                    "dopaminergic",
    "monoamine oxidase b":     "dopaminergic",
    "alpha-synuclein":         "aggregation",
    "lrrk2":                   "dopaminergic",
    "tyrosine hydroxylase":    "dopaminergic",
    "dopamine":                "dopaminergic",
    # ALS / general
    "sod1":                    "oxidative",
    "superoxide dismutase":    "oxidative",
    "tdp-43":                  "aggregation",
    "fus":                     "aggregation",
    "glutamate receptor":      "excitotoxicity",
    "nmda":                    "excitotoxicity",
    "ampa":                    "synaptic",
    # HD
    "hdac":                    "epigenetic",
    "histone deacetylase":     "epigenetic",
    "pde10a":                  "cAMP",
    # Shared neuroprotective pathways
    "nf-kb":                   "inflammation",
    "nf-κb":                   "inflammation",
    "nfkb":                    "inflammation",
    "tnf-alpha":               "inflammation",
    "interleukin":             "inflammation",
    "il-6":                    "inflammation",
    "il-1b":                   "inflammation",
    "cox-2":                   "inflammation",
    "cyclooxygenase":          "inflammation",
    "nrf2":                    "antioxidant",
    "hmox1":                   "antioxidant",
    "heme oxygenase":          "antioxidant",
    "glutathione":             "antioxidant",
    "complex i":               "mitochondrial",
    "mitochondrial":           "mitochondrial",
    "bdnf":                    "neurogenesis",
    "trkb":                    "neurogenesis",
    "ngf":                     "neurogenesis",
    "wnt":                     "neurogenesis",
    "creb":                    "neurogenesis",
    "mtor":                    "autophagy",
    "autophagy":               "autophagy",
    "p62":                     "autophagy",
    "beclin":                  "autophagy",
}


# ---------------------------------------------------------------------------
# Polyphenol / compound class type keywords (for feature engineering)
# ---------------------------------------------------------------------------
POLYPHENOL_TYPES: frozenset[str] = frozenset({
    "flavonoid", "polyphenol", "catechin", "stilbene", "terpene",
    "carotenoid", "vitamin", "phenolic", "alkaloid", "curcuminoid",
    "anthocyanin", "isoflavone", "chalcone", "xanthone", "lignan",
    "coumarin", "resveratrol", "quercetin",
})

NEURO_KWS: frozenset[str] = frozenset({
    "bdnf", "trkb", "wnt", "ngf", "neurogenesis", "hippocampus",
    "creb", "notch", "shh", "vegf", "fgf", "sox2", "nestin",
    "neurite", "synaptogenesis", "ltp", "synaptic", "plasticity",
})


# ---------------------------------------------------------------------------
# BBB estimation rules from physicochemical properties
# (used in predict_unknown_via_ml for live compound lookup)
# ---------------------------------------------------------------------------
def estimate_bbb_class(mw: float, logp: float, tpsa: float, hbd: float) -> tuple[str, int]:
    """
    Estimate BBB permeability class from Lipinski-like physicochemical parameters.
    Based on CNS-MPO heuristics (Wager et al. 2010).

    Returns
    -------
    (bbb_str, bbb_int)  e.g. ("High", 3)
    """
    if mw <= 360 and 1.0 <= logp <= 3.0 and tpsa <= 60 and hbd <= 1:
        return "High",    3
    if mw <= 450 and 0.0 <= logp <= 4.0 and tpsa <= 90 and hbd <= 3:
        return "Medium",  2
    if mw <= 500 and tpsa <= 120:
        return "Low-Med", 1
    return "Low", 0


# ---------------------------------------------------------------------------
# Model file count (used in Methods section and app loading)
# FIX: v5 stated 36 but 4 models × 7 dims + 4 models for NPS = 32.
# We now also count 4 PCA/scaler objects = 36 total serialised objects.
# ---------------------------------------------------------------------------
N_MODEL_FILES: int = 32    # 4 models × (7 dims + 1 NPS composite)
N_TRANSFORM_FILES: int = 4 # ECFP PCA, ChemBERTa PCA, scaler, BBB encoder
N_TOTAL_SERIALISED: int = N_MODEL_FILES + N_TRANSFORM_FILES  # = 36
