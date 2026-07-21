"""
BS_druggability.py — BrainSafe AI (BS)
CNS-weighted druggability scoring from SMILES, grounded in published
medicinal-chemistry scales (no training required — these are validated,
peer-reviewed desirability metrics, so combining them analytically IS the
correct scientific approach).

Components (all computed offline from SMILES via RDKit):
  - QED            Quantitative Estimate of Drug-likeness.
                   Bickerton GR et al. (2012) Nat Chem 4:90-98. doi:10.1038/nchem.1243
  - Lipinski Ro5   MW<=500, cLogP<=5, HBD<=5, HBA<=10.
                   Lipinski CA et al. (2001) Adv Drug Deliv Rev 46:3-26.
  - Veber          rotatable bonds<=10 AND TPSA<=140 (oral bioavailability).
                   Veber DF et al. (2002) J Med Chem 45:2615-2623.
  - CNS-MPO        6-property CNS desirability (0-6). Weighted heavily here
                   because a brain-health compound must reach the CNS to matter.
                   Wager TT et al. (2010) ACS Chem Neurosci 1:435-449. doi:10.1021/cn100008c
  - PAINS          pan-assay interference structural alerts (flag only).
                   Baell JB & Holloway GA (2010) J Med Chem 53:2719-2740.

Public API:
    compute_druggability(smiles: str) -> dict
        Returns overall 0-100 druggability (CNS-weighted), a general
        (non-CNS) drug-likeness score, every sub-score, per-rule pass/fail,
        a category label, a plain-language rationale, and citations.

Run as a script to batch-annotate BS_compounds_full.json:
    python BS_druggability.py            # writes druggability fields in place (with backup)
    python BS_druggability.py --selftest # sanity-check on reference compounds
"""

from __future__ import annotations
import json
import os
import sys
from pathlib import Path

# CNS-weighted composite weights (sum = 1.0). CNS-MPO is weighted highest
# because brain penetration is the gating requirement for a neuro therapeutic.
_W_QED, _W_RO5, _W_VEBER, _W_CNSMPO = 0.30, 0.20, 0.15, 0.35
# General (non-CNS) oral drug-likeness weights (sum = 1.0).
_GW_QED, _GW_RO5, _GW_VEBER = 0.55, 0.25, 0.20

_DIM_DECIMALS = 1

# ---------------------------------------------------------------------------
# RDKit availability (graceful: never crash the app if RDKit is unavailable)
# ---------------------------------------------------------------------------
try:
    from rdkit import Chem
    from rdkit.Chem import QED, Crippen, Descriptors, rdMolDescriptors
    from rdkit.Chem.Lipinski import NumHDonors, NumHAcceptors
    _RDKIT_OK = True
    _RDKIT_ERR = ""
except Exception as _e:  # pragma: no cover
    _RDKIT_OK = False
    _RDKIT_ERR = str(_e)

# PAINS catalog is optional (older/newer RDKit layouts differ)
_PAINS_CATALOG = None
if _RDKIT_OK:
    try:
        from rdkit.Chem import FilterCatalog
        _params = FilterCatalog.FilterCatalogParams()
        _params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        _PAINS_CATALOG = FilterCatalog.FilterCatalog(_params)
    except Exception:
        _PAINS_CATALOG = None


# ---------------------------------------------------------------------------
# CNS-MPO desirability functions (Wager et al. 2010)
# Each returns a 0..1 desirability ("T") score. Piecewise-linear ramps as
# published. pKa is omitted (no reliable offline predictor) — CNS-MPO is
# therefore computed over the 5 structure-derived properties and reported on
# a 0-6 scale by rescaling (documented in the returned `cns_mpo_note`).
# ---------------------------------------------------------------------------
def _ramp_down(x: float, hi_good: float, lo_zero: float) -> float:
    """1.0 at/below hi_good, 0.0 at/above lo_zero, linear between."""
    if x <= hi_good:
        return 1.0
    if x >= lo_zero:
        return 0.0
    return (lo_zero - x) / (lo_zero - hi_good)


def _hump(x: float, rise_lo: float, plateau_lo: float,
          plateau_hi: float, fall_hi: float) -> float:
    """0 below rise_lo, ramps to 1 by plateau_lo, 1 across plateau, ramps to
    0 by fall_hi. Used for TPSA (most desirable in a mid window)."""
    if x <= rise_lo or x >= fall_hi:
        return 0.0
    if plateau_lo <= x <= plateau_hi:
        return 1.0
    if x < plateau_lo:
        return (x - rise_lo) / (plateau_lo - rise_lo)
    return (fall_hi - x) / (fall_hi - plateau_hi)


def _cns_mpo(mw: float, clogp: float, tpsa: float, hbd: int) -> tuple[float, dict]:
    """
    CNS-MPO over 5 structure-derived properties (Wager 2010), rescaled to 0-6.
    cLogD is approximated by cLogP (no offline pKa); pKa term omitted.
    Returns (mpo_0_to_6, per_property_desirabilities).
    """
    t_mw    = _ramp_down(mw,    360.0, 500.0)
    t_clogp = _ramp_down(clogp, 3.0,   5.0)
    t_clogd = _ramp_down(clogp, 2.0,   4.0)          # cLogD ~ cLogP approximation
    t_hbd   = _ramp_down(hbd,   0.5,   3.5)
    t_tpsa  = _hump(tpsa, 20.0, 40.0, 90.0, 120.0)
    parts = {
        "MW": round(t_mw, 3), "cLogP": round(t_clogp, 3),
        "cLogD~cLogP": round(t_clogd, 3), "HBD": round(t_hbd, 3),
        "TPSA": round(t_tpsa, 3),
    }
    mpo5 = t_mw + t_clogp + t_clogd + t_hbd + t_tpsa     # 0..5
    mpo6 = mpo5 * (6.0 / 5.0)                             # rescale to canonical 0..6
    return round(mpo6, 2), parts


# ---------------------------------------------------------------------------
# Main scoring
# ---------------------------------------------------------------------------
def _empty(reason: str) -> dict:
    return {
        "druggability": None, "druggability_general": None,
        "category": "Unknown", "qed": None, "cns_mpo": None,
        "ro5_violations": None, "veber_pass": None, "pains_alerts": None,
        "components": {}, "properties": {}, "rationale": reason,
        "citations": _CITATIONS, "ok": False,
    }


_CITATIONS = [
    "QED: Bickerton GR et al. (2012) Nat Chem 4:90-98.",
    "Lipinski Ro5: Lipinski CA et al. (2001) Adv Drug Deliv Rev 46:3-26.",
    "Veber: Veber DF et al. (2002) J Med Chem 45:2615-2623.",
    "CNS-MPO: Wager TT et al. (2010) ACS Chem Neurosci 1:435-449.",
    "PAINS: Baell JB & Holloway GA (2010) J Med Chem 53:2719-2740.",
]


def compute_druggability(smiles: str | None) -> dict:
    """
    Compute the CNS-weighted druggability score (0-100) for a SMILES string.

    Returns a dict with: druggability (CNS-weighted 0-100), druggability_general
    (non-CNS oral 0-100), category, qed (0-1), cns_mpo (0-6), ro5_violations,
    veber_pass, pains_alerts, a per-component sub-score breakdown, the raw
    physicochemical properties, a plain-language rationale, and citations.
    """
    if not _RDKIT_OK:
        return _empty(f"RDKit unavailable ({_RDKIT_ERR}).")
    if not smiles or str(smiles).strip().lower() in ("n/a", "na", "none", "null", "", "nan"):
        return _empty("No SMILES available for this compound.")

    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return _empty("SMILES could not be parsed by RDKit.")

    # --- raw physicochemical descriptors ---
    mw    = float(Descriptors.MolWt(mol))
    clogp = float(Crippen.MolLogP(mol))
    tpsa  = float(Descriptors.TPSA(mol))
    hbd   = int(NumHDonors(mol))
    hba   = int(NumHAcceptors(mol))
    rotb  = int(rdMolDescriptors.CalcNumRotatableBonds(mol))
    qed   = float(QED.qed(mol))

    # --- Lipinski Ro5 violations ---
    ro5_flags = {
        "MW<=500":   mw <= 500.0,
        "cLogP<=5":  clogp <= 5.0,
        "HBD<=5":    hbd <= 5,
        "HBA<=10":   hba <= 10,
    }
    ro5_violations = sum(1 for ok in ro5_flags.values() if not ok)

    # --- Veber (oral bioavailability) ---
    veber_rotb_ok = rotb <= 10
    veber_tpsa_ok = tpsa <= 140.0
    veber_pass = veber_rotb_ok and veber_tpsa_ok

    # --- CNS-MPO ---
    cns_mpo, cns_parts = _cns_mpo(mw, clogp, tpsa, hbd)

    # --- PAINS structural alerts (flag only; not in the score) ---
    pains_alerts = 0
    pains_names = []
    if _PAINS_CATALOG is not None:
        try:
            matches = _PAINS_CATALOG.GetMatches(mol)
            pains_alerts = len(matches)
            pains_names = [m.GetDescription() for m in matches][:5]
        except Exception:
            pains_alerts = 0

    # --- normalise each component to 0-100 ---
    s_qed   = qed * 100.0
    s_ro5   = (4 - ro5_violations) / 4.0 * 100.0
    s_veber = 100.0 if veber_pass else (50.0 if (veber_rotb_ok or veber_tpsa_ok) else 0.0)
    s_cns   = cns_mpo / 6.0 * 100.0

    # --- composite scores ---
    drug_cns = (_W_QED * s_qed + _W_RO5 * s_ro5 +
                _W_VEBER * s_veber + _W_CNSMPO * s_cns)
    drug_gen = (_GW_QED * s_qed + _GW_RO5 * s_ro5 + _GW_VEBER * s_veber)

    # Mild penalty for PAINS liabilities (developability risk), capped.
    if pains_alerts:
        drug_cns = max(0.0, drug_cns - min(15.0, 7.5 * pains_alerts))
        drug_gen = max(0.0, drug_gen - min(15.0, 7.5 * pains_alerts))

    drug_cns = round(min(100.0, max(0.0, drug_cns)), _DIM_DECIMALS)
    drug_gen = round(min(100.0, max(0.0, drug_gen)), _DIM_DECIMALS)

    if drug_cns >= 70:
        category = "High (drug-like, CNS-favourable)"
    elif drug_cns >= 50:
        category = "Moderate"
    else:
        category = "Low"

    # --- rationale ---
    cns_phrase = ("favourable CNS-MPO" if cns_mpo >= 4 else
                  "borderline CNS-MPO" if cns_mpo >= 3 else "poor CNS-MPO")
    ro5_phrase = ("passes Lipinski Ro5" if ro5_violations <= 1
                  else f"{ro5_violations} Ro5 violations")
    veber_phrase = "passes Veber" if veber_pass else "fails Veber"
    pains_phrase = (f"; {pains_alerts} PAINS alert(s)" if pains_alerts else "")
    rationale = (
        f"QED {qed:.2f}, CNS-MPO {cns_mpo:.1f}/6 ({cns_phrase}), "
        f"{ro5_phrase}, {veber_phrase}{pains_phrase}. "
        f"MW {mw:.0f}, cLogP {clogp:.1f}, TPSA {tpsa:.0f} A^2, HBD {hbd}, "
        f"HBA {hba}, rot-bonds {rotb}."
    )

    return {
        "druggability":         drug_cns,         # CNS-weighted, 0-100 (primary)
        "druggability_general": drug_gen,         # generic oral drug-likeness, 0-100
        "category":             category,
        "qed":                  round(qed, 3),
        "cns_mpo":              cns_mpo,           # 0-6
        "ro5_violations":       ro5_violations,
        "ro5_pass":             ro5_violations <= 1,
        "veber_pass":           bool(veber_pass),
        "pains_alerts":         pains_alerts,
        "pains_names":          pains_names,
        "components": {
            "qed_score":   round(s_qed, 1),
            "ro5_score":   round(s_ro5, 1),
            "veber_score": round(s_veber, 1),
            "cns_score":   round(s_cns, 1),
            "weights_cns": {"qed": _W_QED, "ro5": _W_RO5,
                            "veber": _W_VEBER, "cns_mpo": _W_CNSMPO},
            "cns_mpo_desirabilities": cns_parts,
        },
        "properties": {
            "mw": round(mw, 1), "clogp": round(clogp, 2), "tpsa": round(tpsa, 1),
            "hbd": hbd, "hba": hba, "rotatable_bonds": rotb,
            "ro5_flags": ro5_flags,
        },
        "cns_mpo_note": "CNS-MPO computed over 5 structure-derived properties "
                        "(pKa term omitted: no offline predictor); rescaled to 0-6.",
        "rationale":   rationale,
        "citations":   _CITATIONS,
        "ok":          True,
    }


# ---------------------------------------------------------------------------
# Batch annotation of the app database
# ---------------------------------------------------------------------------
_DB_FILE = "BS_compounds_full.json"


def annotate_database(path: str = _DB_FILE) -> dict:
    """Add druggability fields to every compound in BS_compounds_full.json
    that has a SMILES. Writes in place after backing up. Returns a summary."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{path} not found (run from project root).")
    db = json.loads(p.read_text(encoding="utf-8"))

    backup = p.with_suffix(p.suffix + ".pre_druggability_bak")
    if not backup.exists():
        backup.write_text(json.dumps(db, ensure_ascii=False, indent=2), encoding="utf-8")

    scored, no_smiles, failed = 0, 0, 0
    vals = []
    for name, entry in db.items():
        if str(name).startswith("_") or not isinstance(entry, dict):
            continue
        smi = entry.get("smiles")
        d = compute_druggability(smi)
        if d.get("ok"):
            entry["druggability"]         = d["druggability"]
            entry["druggability_general"] = d["druggability_general"]
            entry["druggability_category"] = d["category"]
            entry["qed"]                  = d["qed"]
            entry["cns_mpo"]              = d["cns_mpo"]
            entry["ro5_violations"]       = d["ro5_violations"]
            entry["veber_pass"]           = d["veber_pass"]
            entry["pains_alerts"]         = d["pains_alerts"]
            entry["druggability_rationale"] = d["rationale"]
            scored += 1
            vals.append(d["druggability"])
        elif not smi:
            no_smiles += 1
        else:
            failed += 1

    p.write_text(json.dumps(db, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = {
        "scored": scored, "no_smiles": no_smiles, "failed_parse": failed,
        "total": len([k for k in db if not str(k).startswith("_")]),
        "mean": round(sum(vals) / len(vals), 1) if vals else None,
        "min": min(vals) if vals else None, "max": max(vals) if vals else None,
        "backup": str(backup),
    }
    return summary


def _selftest():
    """Sanity-check against reference compounds with known profiles."""
    refs = {
        "Donepezil (CNS drug)":      "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2",
        "Caffeine (CNS, small)":     "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "Curcumin (polyphenol)":     "O=C(/C=C/c1ccc(O)c(OC)c1)CC(=O)/C=C/c1ccc(O)c(OC)c1",
        "Sucrose (non-drug, polar)": "OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O",
        "Atorvastatin (big, oral)":  "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
    }
    print(f"{'compound':28} drug(CNS)  gen   QED   MPO  Ro5v  Veber  cat")
    for name, smi in refs.items():
        d = compute_druggability(smi)
        print(f"{name:28} {d['druggability']:6}  {d['druggability_general']:5}  "
              f"{d['qed']:.2f}  {d['cns_mpo']:.1f}   {d['ro5_violations']}    "
              f"{str(d['veber_pass']):5}  {d['category']}")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        _selftest()
    else:
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        print("Annotating", _DB_FILE, "...")
        print(json.dumps(annotate_database(), indent=2))
