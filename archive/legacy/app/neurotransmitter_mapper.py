"""
neurotransmitter_mapper.py  –  BrainSafe AI v6
Maps enzyme/transporter interactions to net neurotransmitter effects.

This is a NEW module implementing the "Neurotransmitter Effect Panel" vision
feature. It bridges the existing enzyme/transporter data (in
CLASS_ENZYME_TEMPLATES) to the 5 key neurotransmitters relevant to NDDs:
    Dopamine, Serotonin, Acetylcholine, GABA, Glutamate.

For each compound, it computes a direction (↑ / ↓ / →) and confidence
(Strong / Moderate / Weak) per neurotransmitter based on all known
enzyme interactions.
"""

from dataclasses import dataclass, field
from typing import Literal

# ---------------------------------------------------------------------------
# Type definitions
# ---------------------------------------------------------------------------
Direction = Literal["increases", "decreases", "modulates", "unknown"]
Strength  = Literal["Strong", "Moderate", "Weak"]


@dataclass
class NTEffect:
    """Effect on one neurotransmitter."""
    neurotransmitter: str
    direction:        Direction
    strength:         Strength
    mechanism:        str          # brief explanation
    relevant_to:      list[str]    # NDDs where this matters most
    display_arrow:    str = field(init=False)

    def __post_init__(self):
        self.display_arrow = {
            "increases":  "↑",
            "decreases":  "↓",
            "modulates":  "⟷",
            "unknown":    "→",
        }.get(self.direction, "→")


# ---------------------------------------------------------------------------
# Enzyme action → neurotransmitter effect rules
# ---------------------------------------------------------------------------
# Format: "enzyme_keyword": {neurotransmitter: (direction, strength, mechanism)}
# These rules are grounded in established neuropharmacology.
ENZYME_TO_NT_RULES: dict[str, dict[str, tuple]] = {

    # Dopamine
    "mao-b":     {"dopamine": ("increases", "Moderate",
                               "MAO-B inhibition reduces dopamine catabolism, "
                               "increasing synaptic dopamine availability.")},
    "maob":      {"dopamine": ("increases", "Moderate",
                               "MAO-B inhibition → ↑ dopamine.")},
    "monoamine oxidase b": {"dopamine": ("increases", "Strong",
                               "Selective MAO-B inhibition is the mechanism of "
                               "selegiline; raises dopamine in striatum.")},
    "tyrosine hydroxylase": {"dopamine": ("increases", "Moderate",
                               "Activation of TH (rate-limiting enzyme) → ↑ dopamine synthesis.")},
    "dopa decarboxylase":   {"dopamine": ("increases", "Strong",
                               "DOPA decarboxylase converts L-DOPA → dopamine (PD mechanism).")},
    "aadc":               {"dopamine": ("increases", "Strong",
                               "AADC (DOPA decarboxylase) requires PLP cofactor; "
                               "converts L-DOPA and 5-HTP to dopamine and serotonin.")},
    "dopamine transporter": {"dopamine": ("increases", "Moderate",
                               "DAT inhibition raises synaptic dopamine (stimulant mechanism).")},
    "comt":               {"dopamine": ("increases", "Moderate",
                               "COMT inhibition reduces dopamine breakdown; "
                               "used adjunctively in PD therapy.")},

    # Serotonin
    "mao-a":     {"serotonin": ("increases", "Moderate",
                                "MAO-A inhibition reduces serotonin catabolism.")},
    "maoa":      {"serotonin": ("increases", "Moderate",
                                "MAO-A preferentially deaminates serotonin and norepinephrine.")},
    "serotonin transporter": {"serotonin": ("increases", "Strong",
                                "SERT inhibition (SSRI mechanism) raises synaptic serotonin.")},
    "tryptophan hydroxylase": {"serotonin": ("increases", "Moderate",
                                "TPH is the rate-limiting enzyme for serotonin synthesis.")},
    "5-ht3":     {"serotonin": ("modulates", "Moderate",
                                "5-HT3 antagonism modulates serotonergic signalling.")},

    # Acetylcholine
    "acetylcholinesterase": {"acetylcholine": ("increases", "Strong",
                              "AChE inhibition prevents ACh breakdown, raising synaptic "
                              "acetylcholine — mechanism of donepezil, galantamine.")},
    "ache":     {"acetylcholine": ("increases", "Strong",
                              "AChE inhibition → ↑ acetylcholine (cholinergic enhancement).")},
    "choline acetyltransferase": {"acetylcholine": ("increases", "Moderate",
                              "ChAT activation → ↑ ACh synthesis.")},
    "muscarinic": {"acetylcholine": ("modulates", "Weak",
                              "Muscarinic receptor modulation affects cholinergic signalling.")},
    "chat":     {"acetylcholine": ("increases", "Moderate",
                              "ChAT converts choline + acetyl-CoA → ACh.")},

    # GABA
    "gaba transaminase": {"gaba": ("increases", "Moderate",
                          "GABA-T inhibition → ↑ GABA levels (mechanism of vigabatrin).")},
    "gaba-t":   {"gaba": ("increases", "Moderate",
                          "GABA transaminase inhibition elevates brain GABA.")},
    "glutamate decarboxylase": {"gaba": ("increases", "Moderate",
                          "GAD converts glutamate → GABA; activation raises GABA.")},
    "gad":      {"gaba": ("increases", "Moderate",
                          "GAD (glutamate decarboxylase) is the GABA synthesis enzyme.")},
    "gaba-a":   {"gaba": ("increases", "Weak",
                          "GABA-A receptor modulation (benzodiazepine-like mechanism).")},
    "benzodiazepine": {"gaba": ("increases", "Strong",
                          "Positive allosteric modulation of GABA-A receptor.")},

    # Glutamate (excitatory — most effects are reductive for neuroprotection)
    "nmda":     {"glutamate": ("decreases", "Moderate",
                               "NMDA receptor antagonism reduces excitotoxic glutamate "
                               "activity (mechanism of memantine in AD).")},
    "glutamate receptor": {"glutamate": ("decreases", "Moderate",
                               "Glutamate receptor modulation; antagonism is neuroprotective.")},
    "ampa":     {"glutamate": ("modulates", "Weak",
                               "AMPA receptor modulation affects fast excitatory transmission.")},
    "glast":    {"glutamate": ("decreases", "Moderate",
                               "Glutamate transporter activation reduces synaptic glutamate.")},
    "glt-1":    {"glutamate": ("decreases", "Moderate",
                               "GLT-1 (EAAT2) upregulation clears synaptic glutamate; "
                               "target of riluzole in ALS.")},
    "xct":      {"glutamate": ("modulates", "Weak",
                               "xCT cystine/glutamate antiporter modulates glutamate release.")},
}

# NDD relevance of each neurotransmitter
NT_DISEASE_RELEVANCE: dict[str, list[str]] = {
    "dopamine":      ["parkinsons", "huntingtons"],
    "serotonin":     ["alzheimers", "parkinsons"],
    "acetylcholine": ["alzheimers"],
    "gaba":          ["huntingtons", "als"],
    "glutamate":     ["als", "alzheimers", "huntingtons"],
}

# Plain-language NT names
NT_DISPLAY_NAMES: dict[str, str] = {
    "dopamine":      "Dopamine",
    "serotonin":     "Serotonin",
    "acetylcholine": "Acetylcholine",
    "gaba":          "GABA (calming signal)",
    "glutamate":     "Glutamate (excitatory signal)",
}


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------
def compute_nt_effects(enzyme_entries: list[dict]) -> list[NTEffect]:
    """
    Compute neurotransmitter effects from a list of enzyme interaction records.

    Parameters
    ----------
    enzyme_entries : list of dicts with at minimum:
        {"name": str, "action": str, "strength": str}
        (from CLASS_ENZYME_TEMPLATES or compound-specific ENZYME_DATA)

    Returns
    -------
    list of NTEffect objects, one per affected neurotransmitter.
    """
    # Accumulate effects per NT
    nt_effects: dict[str, list[tuple[str, str, str]]] = {}   # NT → [(dir, strength, mech)]

    for entry in enzyme_entries:
        enzyme_name   = (entry.get("name", "") or "").lower()
        enzyme_action = (entry.get("action", "") or "").lower()
        entry_strength = entry.get("strength", "Weak")

        for keyword, nt_map in ENZYME_TO_NT_RULES.items():
            if keyword in enzyme_name:
                for nt, (direction, rule_strength, mechanism) in nt_map.items():
                    # Allow the compound entry's strength to override the rule strength
                    # if it provides more specific data
                    actual_strength = entry_strength if entry_strength in ("Strong", "Moderate", "Weak") else rule_strength
                    if nt not in nt_effects:
                        nt_effects[nt] = []
                    nt_effects[nt].append((direction, actual_strength, mechanism))

    # Consolidate: if multiple enzymes affect the same NT, take dominant direction
    results: list[NTEffect] = []
    strength_order = {"Strong": 3, "Moderate": 2, "Weak": 1}

    for nt, effect_list in nt_effects.items():
        # Count direction votes, weighted by strength
        direction_votes: dict[str, float] = {}
        best_mechanism = ""
        best_strength  = "Weak"

        for direction, strength, mechanism in effect_list:
            weight = strength_order.get(strength, 1)
            direction_votes[direction] = direction_votes.get(direction, 0) + weight
            if strength_order.get(strength, 1) >= strength_order.get(best_strength, 1):
                best_strength  = strength
                best_mechanism = mechanism

        dominant_dir = max(direction_votes, key=direction_votes.__getitem__)

        results.append(NTEffect(
            neurotransmitter = nt,
            direction        = dominant_dir,
            strength         = best_strength,
            mechanism        = best_mechanism,
            relevant_to      = NT_DISEASE_RELEVANCE.get(nt, []),
        ))

    # Sort by NDD relevance count (most broadly relevant first)
    results.sort(key=lambda x: -len(x.relevant_to))
    return results


def format_nt_summary(effects: list[NTEffect],
                      patient_mode: bool = False) -> list[dict]:
    """
    Format NT effects for display.

    Parameters
    ----------
    effects      : list of NTEffect
    patient_mode : bool – if True, return plain-language descriptions

    Returns
    -------
    list of dicts ready for the Streamlit UI
    """
    rows = []
    for e in effects:
        display_name = NT_DISPLAY_NAMES.get(e.neurotransmitter, e.neurotransmitter.title())
        if patient_mode:
            direction_text = {
                "increases":  f"may increase {display_name}",
                "decreases":  f"may decrease {display_name}",
                "modulates":  f"may modulate {display_name}",
                "unknown":    f"has uncertain effects on {display_name}",
            }.get(e.direction, "")
            disease_note = ""
            if e.relevant_to:
                d_names = {
                    "alzheimers": "Alzheimer's", "parkinsons": "Parkinson's",
                    "als": "ALS", "huntingtons": "Huntington's",
                }
                disease_note = "Relevant to: " + ", ".join(
                    d_names.get(d, d) for d in e.relevant_to
                )
            rows.append({
                "neurotransmitter": display_name,
                "arrow":            e.display_arrow,
                "direction":        direction_text,
                "strength":         e.strength,
                "note":             disease_note,
                "mechanism":        e.mechanism,
            })
        else:
            rows.append({
                "neurotransmitter": display_name,
                "arrow":            e.display_arrow,
                "direction":        e.direction,
                "strength":         e.strength,
                "mechanism":        e.mechanism,
                "relevant_to":      e.relevant_to,
            })
    return rows
