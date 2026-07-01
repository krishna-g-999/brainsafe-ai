"""
brain_region_mapper.py  –  BrainSafe AI v6
Maps a compound's 7-axis neuroprotective profile to brain region specificity.

This NEW module implements the "Brain Region Map" vision feature.
Each of the 7 neuroprotective dimensions has differential importance
in specific brain regions based on known neuroanatomy and disease pathology.

Data sources:
  - Allen Brain Atlas human gene expression data (allen.brain-map.org)
  - Disease-region associations from key review papers:
    AD:  Masters et al. 2015 (Nat Rev Dis Primers); Mufson et al. 2008
    PD:  Kalia & Lang 2015 (Lancet); Bose & Beal 2016
    ALS: Taylor et al. 2016 (Nature)
    HD:  Ross & Tabrizi 2011 (Lancet Neurol)
"""

from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Brain regions
# ---------------------------------------------------------------------------
BRAIN_REGIONS: dict[str, dict] = {
    "hippocampus": {
        "display":     "Hippocampus",
        "role":        "Memory formation, spatial navigation, neurogenesis hub",
        "ndd_primary": ["alzheimers"],
        "coordinates": {"x": 310, "y": 185},   # SVG coordinates for brain map
    },
    "entorhinal_cortex": {
        "display":     "Entorhinal cortex",
        "role":        "Memory consolidation, early site of tau pathology in AD",
        "ndd_primary": ["alzheimers"],
        "coordinates": {"x": 275, "y": 210},
    },
    "prefrontal_cortex": {
        "display":     "Prefrontal cortex",
        "role":        "Executive function, working memory, decision-making",
        "ndd_primary": ["alzheimers", "huntingtons"],
        "coordinates": {"x": 180, "y": 130},
    },
    "substantia_nigra": {
        "display":     "Substantia nigra",
        "role":        "Dopamine production, motor control; primary PD target",
        "ndd_primary": ["parkinsons"],
        "coordinates": {"x": 290, "y": 280},
    },
    "striatum": {
        "display":     "Striatum (caudate/putamen)",
        "role":        "Motor control, reward; primary HD target",
        "ndd_primary": ["huntingtons", "parkinsons"],
        "coordinates": {"x": 255, "y": 195},
    },
    "motor_cortex": {
        "display":     "Motor cortex",
        "role":        "Voluntary movement; upper motor neuron origin in ALS",
        "ndd_primary": ["als"],
        "coordinates": {"x": 210, "y": 115},
    },
    "spinal_cord": {
        "display":     "Spinal cord (anterior horn)",
        "role":        "Lower motor neurons; degeneration site in ALS",
        "ndd_primary": ["als"],
        "coordinates": {"x": 320, "y": 370},
    },
    "cerebellum": {
        "display":     "Cerebellum",
        "role":        "Motor coordination, balance; involved in ataxic forms of NDD",
        "ndd_primary": [],
        "coordinates": {"x": 375, "y": 310},
    },
    "amygdala": {
        "display":     "Amygdala",
        "role":        "Emotional processing, fear memory; affected in AD",
        "ndd_primary": ["alzheimers"],
        "coordinates": {"x": 295, "y": 225},
    },
    "brainstem": {
        "display":     "Brainstem",
        "role":        "Autonomic function; Lewy body pathology in early PD",
        "ndd_primary": ["parkinsons"],
        "coordinates": {"x": 320, "y": 315},
    },
}


# ---------------------------------------------------------------------------
# Dimension → brain region importance weights
# Each entry: (region_key, weight_0_to_1, rationale)
# ---------------------------------------------------------------------------
DIM_TO_REGION: dict[str, list[tuple[str, float, str]]] = {

    "antioxidant": [
        ("hippocampus",        0.9, "High metabolic rate → high ROS vulnerability; "
                                     "oxidative stress drives hippocampal neurodegeneration in AD"),
        ("substantia_nigra",   0.9, "Dopamine metabolism generates H2O2; "
                                     "complex I ROS → SN neurodegeneration in PD"),
        ("motor_cortex",       0.7, "Upper motor neurons vulnerable to SOD1 ROS in ALS"),
        ("striatum",           0.6, "Mitochondrial ROS drives striatal neurodegeneration in HD"),
        ("cerebellum",         0.5, "Moderate antioxidant relevance"),
        ("prefrontal_cortex",  0.6, "AD-related oxidative stress in PFC"),
    ],

    "anti_inflammatory": [
        ("hippocampus",        0.9, "Microglial activation and IL-6/TNF-α "
                                     "neuroinflammation in AD hippocampus"),
        ("substantia_nigra",   0.9, "NLRP3 inflammasome activation drives "
                                     "dopaminergic neurodegeneration in PD"),
        ("striatum",           0.8, "Neuroinflammation central to HD pathology"),
        ("prefrontal_cortex",  0.7, "NF-κB-driven inflammation in cortex"),
        ("motor_cortex",       0.7, "Neuroinflammation in ALS upper motor neurons"),
        ("spinal_cord",        0.8, "Microglial/astrocyte activation in ALS spinal cord"),
    ],

    "mitochondrial_support": [
        ("substantia_nigra",   0.95, "Complex I deficiency is THE primary mitochondrial "
                                      "dysfunction in PD; Beal 2016"),
        ("hippocampus",        0.8,  "Mitochondrial dysfunction in AD hippocampal neurons"),
        ("motor_cortex",       0.85, "Mitochondrial stress in ALS upper motor neurons"),
        ("spinal_cord",        0.85, "Mitochondrial dysfunction in ALS motor neurons"),
        ("striatum",           0.75, "Impaired mitochondrial function in HD striatum"),
        ("cerebellum",         0.5,  "Moderate mitochondrial relevance"),
    ],

    "aggregation_modulation": [
        ("hippocampus",        0.95, "Amyloid-β plaques and tau neurofibrillary "
                                      "tangles begin in hippocampus/entorhinal cortex in AD"),
        ("entorhinal_cortex",  0.95, "Earliest site of tau pathology (Braak stages I–II)"),
        ("substantia_nigra",   0.85, "α-synuclein Lewy bodies in dopaminergic neurons"),
        ("striatum",           0.9,  "Mutant HTT aggregates selectively toxic to striatal neurons"),
        ("motor_cortex",       0.8,  "TDP-43 / FUS aggregation in ALS upper motor neurons"),
        ("spinal_cord",        0.8,  "SOD1/TDP-43 aggregation in ALS lower motor neurons"),
        ("amygdala",           0.75, "Amyloid pathology extends to amygdala in advanced AD"),
    ],

    "cognitive_enhancement": [
        ("hippocampus",        0.95, "Cholinergic, BDNF, and synaptic plasticity mechanisms "
                                      "for memory; hippocampus is primary cognitive target"),
        ("entorhinal_cortex",  0.85, "Memory consolidation gateway to hippocampus"),
        ("prefrontal_cortex",  0.85, "Working memory, attention, executive function"),
        ("amygdala",           0.6,  "Emotional learning, contextual memory"),
        ("striatum",           0.55, "Procedural learning and habit formation"),
        ("cerebellum",         0.45, "Cognitive cerebellar syndrome in some NDDs"),
    ],

    "neurogenesis": [
        ("hippocampus",        0.95, "Adult hippocampal neurogenesis (dentate gyrus) "
                                      "is the primary neurogenesis site; BDNF/Wnt/Notch dependent"),
        ("prefrontal_cortex",  0.5,  "Limited adult neurogenesis; BDNF signalling important"),
        ("striatum",           0.4,  "Limited striatal neurogenesis; relevant in HD models"),
        ("cerebellum",         0.3,  "Cerebellar granule cell renewal in some models"),
    ],

    "synaptic_plasticity": [
        ("hippocampus",        0.95, "LTP and LTD at CA1-CA3 synapses; "
                                      "primary site for memory synaptic plasticity"),
        ("prefrontal_cortex",  0.8,  "Cortical synaptic plasticity for learning"),
        ("striatum",           0.75, "Corticostriatal synaptic plasticity in HD/PD"),
        ("cerebellum",         0.7,  "Cerebellar LTD for motor learning"),
        ("amygdala",           0.6,  "Fear conditioning synaptic plasticity"),
        ("motor_cortex",       0.6,  "Motor learning plasticity"),
    ],
}


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------
@dataclass
class RegionScore:
    region_key:  str
    display:     str
    score:       float    # 0–100 activation score
    role:        str
    ndd_primary: list[str]
    coordinates: dict
    top_mechanism: str


def compute_region_scores(dimension_scores: dict[str, float]) -> list[RegionScore]:
    """
    Compute brain region activation scores from 7-axis dimension scores.

    Parameters
    ----------
    dimension_scores : dict mapping dimension name → score (0–100)

    Returns
    -------
    List of RegionScore objects sorted by score descending.
    """
    region_totals: dict[str, dict] = {k: {"weighted": 0.0, "weight_sum": 0.0,
                                          "best_mech": ""}
                                       for k in BRAIN_REGIONS}

    for dim, dim_score in dimension_scores.items():
        if dim not in DIM_TO_REGION:
            continue
        for region_key, importance, rationale in DIM_TO_REGION[dim]:
            contribution = dim_score * importance
            region_totals[region_key]["weighted"]    += contribution
            region_totals[region_key]["weight_sum"]  += importance
            if importance >= 0.8 and dim_score >= 50.0:
                if not region_totals[region_key]["best_mech"]:
                    region_totals[region_key]["best_mech"] = rationale

    # Normalise by the global maximum weighted sum so that the most-activated
    # region scores 100. This rewards regions (like hippocampus) that are
    # relevant across ALL 7 dimensions rather than just a few high-scoring ones.
    max_weighted = max((d["weighted"] for d in region_totals.values()), default=1.0)
    if max_weighted == 0:
        max_weighted = 1.0

    results = []
    for region_key, info in BRAIN_REGIONS.items():
        d = region_totals[region_key]
        score = round(min(100.0, d["weighted"] / max_weighted * 100.0), 1)
        results.append(RegionScore(
            region_key    = region_key,
            display       = info["display"],
            score         = score,
            role          = info["role"],
            ndd_primary   = info["ndd_primary"],
            coordinates   = info["coordinates"],
            top_mechanism = d["best_mech"] or f"Activity in {info['display']}",
        ))

    results.sort(key=lambda x: -x.score)
    return results


def get_primary_regions(dimension_scores: dict[str, float],
                        top_n: int = 4) -> list[RegionScore]:
    """Return top N most activated brain regions for the compound."""
    return compute_region_scores(dimension_scores)[:top_n]


def format_region_summary(regions: list[RegionScore]) -> list[dict]:
    """Format region scores for Streamlit display."""
    return [
        {
            "region":    r.display,
            "score":     r.score,
            "role":      r.role,
            "diseases":  r.ndd_primary,
            "mechanism": r.top_mechanism,
            "color":     "green" if r.score >= 65 else "orange" if r.score >= 35 else "red",
        }
        for r in regions
    ]
