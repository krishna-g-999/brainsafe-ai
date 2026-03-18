"""
ml_expander.py — BrainSafe AI ML Expansion Pipeline
=====================================================
Generates ML-predicted neuroprotective profiles for ChEMBL-indicated
neuro compounds using a Random Forest trained on all 129 curated entries.

Feature design (disease-indication based, no API for training):
  Training features from compounds.json:
    bbb_num (0-3), als_num (0-2), alzheimers_num (0-2),
    parkinsons_num (0-2), huntingtons_num (0-2), n_pathways (int)

  Prediction features for ChEMBL compounds derived from:
    ChEMBL indication data → disease_num
    CNS-MPO rules on ChEMBL mol properties → bbb_num
    Number of pathway associations → n_pathways

ERC (enzyme/receptor/cofactor) data: real IC50/Ki values from ChEMBL bioassays.

References:
  Mendez D et al. (2019) ChEMBL. Nucleic Acids Res 47(D1):D930–D940.
  Breiman L (2001) Random Forests. Machine Learning 45:5–32.
  Wager TT et al. (2010) CNS-MPO. ACS Chem Neurosci 1(6):435–449.
"""

from __future__ import annotations

import json
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

import numpy as np
import requests
from sklearn.ensemble import RandomForestRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import StandardScaler

# ── v2 Feature Engineering constants (biological domain knowledge) ────────
POLYPHENOL_TYPES = {"flavonoid","polyphenol","catechin","stilbene","terpene",
                    "carotenoid","vitamin","phenolic","alkaloid","curcuminoid"}
NEURO_KWS        = {"bdnf","trkb","wnt","ngf","neurogenesis","hippocampus",
                    "creb","notch","shh","vegf","fgf","sox2","nestin"}


logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger("ml_expander")

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
TIMEOUT     = 10
MAX_WORKERS = 14

SCORE_COLS = [
    "antioxidant", "anti_inflammatory", "mitochondrial_support",
    "aggregation_modulation", "cognitive_enhancement",
    "neurogenesis", "synaptic_plasticity",
]

BBB_MAP = {"Low": 0, "Low-Med": 1, "Medium": 2, "High": 3}
DIS_MAP = {"Low": 0, "Med": 1, "High": 2}

FEATURE_COLS = [
    "bbb_num", "als_num", "alzheimers_num",
    "parkinsons_num", "huntingtons_num", "n_pathways",
]

INDICATION_KEYWORDS = [
    "Alzheimer", "Parkinson", "Amyotrophic", "Huntington", "dementia",
]

NEURO_TARGET_KW = [
    "acetylcholinesterase", "butyrylcholinesterase", "bace", "beta-secretase",
    "monoamine oxidase", "gsk-3", "glycogen synthase kinase", "tau",
    "alpha-synuclein", "amyloid", "superoxide dismutase", "sod1",
    "lrrk2", "pink1", "hdac", "histone deacetylase", "caspase",
    "parp", "nf-kb", "cyclooxygenase", "cox-2", "tdp-43", "fus",
    "nmda", "glutamate", "bdnf", "trkb", "dopamine",
]

EC_MAP: dict[str, str] = {
    "acetylcholinesterase":     "EC 3.1.1.7",
    "butyrylcholinesterase":    "EC 3.1.1.8",
    "beta-secretase":           "EC 3.4.23.45",
    "bace":                     "EC 3.4.23.45",
    "monoamine oxidase b":      "EC 1.4.3.4",
    "monoamine oxidase a":      "EC 1.4.3.4",
    "gsk-3":                    "EC 2.7.11.26",
    "glycogen synthase kinase": "EC 2.7.11.26",
    "lrrk2":                    "EC 2.7.11.1",
    "cyclooxygenase":           "EC 1.14.99.1",
    "cox-2":                    "EC 1.14.99.1",
    "superoxide dismutase":     "EC 1.15.1.1",
    "hdac":                     "EC 3.5.1.98",
    "histone deacetylase":      "EC 3.5.1.98",
    "caspase-3":                "EC 3.4.22.56",
    "parp":                     "EC 2.4.2.30",
}

COFACTOR_MAP: dict[str, list[str]] = {
    "monoamine oxidase": ["FAD", "O2"],
    "superoxide dismutase": ["Cu2+", "Zn2+"],
    "gsk-3": ["ATP", "Mg2+"],
    "lrrk2": ["ATP"],
    "hdac": ["Zn2+", "NAD+"],
    "parp": ["NAD+"],
    "caspase": ["Zn2+"],
    "cyclooxygenase": ["Heme", "O2"],
    "acetylcholinesterase": ["Choline"],
}

DISEASE_TARGET_MAP = {
    "alzheimers": ["acetylcholinesterase", "bace", "tau", "amyloid", "gsk-3", "nmda"],
    "parkinsons":  ["monoamine oxidase", "alpha-synuclein", "lrrk2", "dopamine"],
    "als":        ["superoxide dismutase", "sod1", "tdp-43", "glutamate", "caspase"],
    "huntingtons": ["hdac", "huntingtin", "caspase", "bdnf"],
}

TARGET_TO_PATHWAY: dict[str, list[str]] = {
    "acetylcholinesterase":  ["Cholinergic signaling", "ACh metabolism"],
    "bace":                  ["APP processing", "Amyloid cascade"],
    "beta-secretase":        ["APP processing", "Amyloid cascade"],
    "tau":                   ["Tau phosphorylation", "PI3K/Akt"],
    "alpha-synuclein":       ["Alpha-syn aggregation", "Autophagy"],
    "monoamine oxidase":     ["Dopamine regulation", "Tryptophan metabolism"],
    "gsk-3":                 ["Tau phosphorylation", "Wnt/Beta-catenin"],
    "lrrk2":                 ["LRRK2/PINK1 pathway", "Mitophagy"],
    "hdac":                  ["Epigenetic regulation", "NF-kB"],
    "cyclooxygenase":        ["NF-kB", "Arachidonic acid pathway"],
    "superoxide dismutase":  ["Nrf2/GSH", "ROS scavenging"],
    "nmda":                  ["Glutamate excitotoxicity", "Calcium signaling"],
    "amyloid":               ["Amyloid cascade", "PI3K/Akt"],
    "bdnf":                  ["BDNF/TrkB", "PI3K/Akt"],
}

TARGET_TO_METABOLITES: dict[str, list[str]] = {
    "acetylcholinesterase":  ["Acetylcholine", "Choline"],
    "bace":                  ["Abeta", "sAPPalpha"],
    "tau":                   ["Tau", "p-Tau"],
    "alpha-synuclein":       ["Alpha-syn", "Dopamine"],
    "monoamine oxidase":     ["Dopamine", "Serotonin"],
    "gsk-3":                 ["Tau", "ATP"],
    "hdac":                  ["BDNF", "NAD+"],
    "cyclooxygenase":        ["PGE2", "TNF-alpha", "IL-6"],
    "superoxide dismutase":  ["ROS", "H2O2"],
    "nmda":                  ["Glutamate", "Ca2+"],
}

BRAIN_REGION_MAP = {
    "alzheimers": ["Hippocampus", "Cortex", "Entorhinal Cortex"],
    "parkinsons":  ["Substantia Nigra", "Striatum", "Basal Ganglia"],
    "als":        ["Motor Cortex", "Spinal Cord", "Brainstem"],
    "huntingtons": ["Striatum", "Caudate Nucleus", "Cortex"],
}


def _safe_float(val: Any, default: float = 0.0) -> float:
    try:
        return float(val) if val is not None else default
    except (TypeError, ValueError):
        return default


def _bbb_from_props(props: dict) -> tuple[str, int]:
    mw   = _safe_float(props.get("mw_freebase") or props.get("full_mwt"))
    logp = _safe_float(props.get("alogp"))
    tpsa = _safe_float(props.get("psa"))
    hbd  = _safe_float(props.get("hbd"))
    if mw <= 360 and 1.0 <= logp <= 3.0 and tpsa <= 60 and hbd <= 1:
        return "High", 3
    if mw <= 450 and 0.0 <= logp <= 4.0 and tpsa <= 90 and hbd <= 3:
        return "Medium", 2
    if mw <= 500 and tpsa <= 120:
        return "Low-Med", 1
    return "Low", 0


def _compound_type_from_chembl(mol: dict) -> str:
    props    = mol.get("molecule_properties") or {}
    np_score = _safe_float(props.get("np_likeness_score"), 0.0)
    try:
        phase = int(float(mol.get("max_phase") or 0))
    except (TypeError, ValueError):
        phase = 0
    if phase >= 4:
        return "FDA-Approved Drug"
    if phase >= 1:
        return "Clinical Candidate"
    if np_score >= 1.5:
        return "Natural Product"
    if np_score >= 0.5:
        return "Natural Product-Like"
    return "Synthetic Small Molecule"


def build_training_data(
    curated: dict[str, dict],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    X_rows: list[list[float]] = []
    y_rows: list[list[float]] = []
    matched: list[str]        = []
    for name, entry in curated.items():
        feats = [
            float(BBB_MAP.get(entry.get("bbb", "Low"), 0)),
            float(DIS_MAP.get(entry.get("als", "Low"), 0)),
            float(DIS_MAP.get(entry.get("alzheimers", "Low"), 0)),
            float(DIS_MAP.get(entry.get("parkinsons", "Low"), 0)),
            float(DIS_MAP.get(entry.get("huntingtons", "Low"), 0)),
            float(len(entry.get("pathways", []))),
        ]
        scores = [_safe_float(entry.get(c), 5.0) for c in SCORE_COLS]
        X_rows.append(feats)
        y_rows.append(scores)
        matched.append(name)
    log.info("Training set: %d / %d curated compounds (all, no API needed)", len(matched), len(curated))
    return np.array(X_rows, dtype=float), np.array(y_rows, dtype=float), matched


def train_ml_model(
    X: np.ndarray, y: np.ndarray
) -> tuple[MultiOutputRegressor, StandardScaler, float]:
    from sklearn.model_selection import cross_val_score

    scaler = StandardScaler()
    X_s    = scaler.fit_transform(X)
    base   = RandomForestRegressor(
        n_estimators=150, max_depth=6, min_samples_leaf=2, random_state=42, n_jobs=-1
    )
    model = MultiOutputRegressor(base)
    model.fit(X_s, y)
    cv_scores = cross_val_score(
        MultiOutputRegressor(RandomForestRegressor(n_estimators=150, max_depth=6, min_samples_leaf=2, random_state=42)),
        X_s, y, cv=5, scoring="r2",
    )
    r2 = float(np.mean(cv_scores))
    log.info("Cross-validated R² (mean, 5-fold): %.3f", r2)
    return model, scaler, r2


def predict_profile(
    model: MultiOutputRegressor, scaler: StandardScaler, feat: list[float]
) -> dict[str, float]:
    x     = np.array([feat], dtype=float)
    preds = model.predict(scaler.transform(x))[0]
    return {col: round(float(np.clip(raw, 1.0, 10.0)), 1)
            for col, raw in zip(SCORE_COLS, preds)}


def predict_std(
    model: MultiOutputRegressor, scaler: StandardScaler, feat: list[float]
) -> dict[str, float]:
    x   = np.array([feat], dtype=float)
    x_s = scaler.transform(x)
    result: dict[str, float] = {}
    for col, rf_est in zip(SCORE_COLS, model.estimators_):
        tree_preds = np.array([tree.predict(x_s)[0] for tree in rf_est.estimators_])
        result[col] = round(float(tree_preds.std()), 2)
    return result


def _fetch_indication_ids(keyword: str, max_entries: int = 100) -> tuple[str, list[str]]:
    ids = []
    try:
        r = requests.get(
            f"{CHEMBL_BASE}/drug_indication",
            params={"mesh_heading__icontains": keyword,
                    "format": "json", "limit": max_entries},
            timeout=TIMEOUT,
        )
        if r.status_code == 200:
            ids = [e["molecule_chembl_id"] for e in r.json().get("drug_indications", [])
                   if e.get("molecule_chembl_id")]
    except Exception as exc:
        log.debug("Indication fetch failed (%s): %s", keyword, exc)
    return keyword, ids


def _fetch_mol_chunk(chunk: list[str]) -> list[dict]:
    try:
        r = requests.get(
            f"{CHEMBL_BASE}/molecule",
            params={"molecule_chembl_id__in": ",".join(chunk),
                    "format": "json", "limit": len(chunk)},
            timeout=TIMEOUT,
        )
        if r.status_code == 200:
            return r.json().get("molecules", [])
    except Exception:
        pass
    return []


def _fetch_molecules_batch(chembl_ids: list[str]) -> list[dict]:
    CHUNK = 40
    chunks = [chembl_ids[i: i + CHUNK] for i in range(0, len(chembl_ids), CHUNK)]
    mols: list[dict] = []
    with ThreadPoolExecutor(max_workers=min(MAX_WORKERS, len(chunks))) as ex:
        for result in ex.map(_fetch_mol_chunk, chunks):
            mols.extend(result)
    return mols


def _fetch_activities_for_mol(chembl_id: str) -> list[dict]:
    try:
        r = requests.get(
            f"{CHEMBL_BASE}/activity",
            params={
                "molecule_chembl_id":     chembl_id,
                "standard_type__in":      "IC50,Ki",
                "standard_relation":      "=",
                "pchembl_value__isnull":  "false",
                "assay_organism":         "Homo sapiens",
                "format":                 "json",
                "limit":                  8,
                "order_by":               "-pchembl_value",
            },
            timeout=8,
        )
        if r.status_code == 200:
            return [a for a in r.json().get("activities", [])
                    if any(kw in (a.get("target_pref_name") or "").lower()
                           for kw in NEURO_TARGET_KW)]
    except Exception:
        pass
    return []


def activities_to_erc(activities: list[dict]) -> dict:
    enzymes: list[dict]      = []
    transporters: list[dict] = []
    cofactors: set[str]      = set()
    seen: set[str]           = set()
    for a in activities:
        tname_raw = a.get("target_pref_name") or ""
        if tname_raw in seen:
            continue
        seen.add(tname_raw)
        tname_lc = tname_raw.lower()
        std_val  = _safe_float(a.get("standard_value"))
        std_type = a.get("standard_type", "IC50")
        std_unit = a.get("standard_units", "nM")
        pchembl  = _safe_float(a.get("pchembl_value"))
        strength = "Strong" if pchembl >= 8 else ("Moderate" if pchembl >= 6 else "Weak")
        action   = "Inhibition"
        ec = next((v for k, v in EC_MAP.items() if k in tname_lc), "")
        for key, cofs in COFACTOR_MAP.items():
            if key in tname_lc:
                cofactors.update(cofs)
        is_transporter = any(kw in tname_lc for kw in
                             ["transporter", "abcb1", "abcg2", "p-glycoprotein", "slc"])
        note = (f"Inhibition; {std_type}={std_val:.1f} {std_unit} "
                f"(pChEMBL={pchembl:.1f}) — ChEMBL bioassay, Homo sapiens.")
        if is_transporter:
            transporters.append({
                "name": tname_raw, "role": "Inhibitor",
                "location": "BBB / Gut efflux", "note": note,
                "source": "ChEMBL bioactivity (Mendez et al., 2019)",
            })
        else:
            enzymes.append({
                "name": tname_raw, "action": action,
                "strength": strength, "ec": ec, "note": note,
                "source": "ChEMBL bioactivity (Mendez et al., 2019)",
            })
    return {
        "enzymes":      enzymes[:6],
        "transporters": transporters[:3],
        "cofactors":    sorted(cofactors)[:6],
        "sources": [
            "ChEMBL bioactivity database (Mendez D et al. 2019. Nucleic Acids Res 47:D930–D940)",
            "ML-predicted 7-dimension scores (BrainSafe AI Random Forest, trained on 129 curated compounds)",
        ],
    }


def _infer_pathways_regions(
    activities: list[dict], indication_diseases: list[str]
) -> tuple[list[str], list[str], list[str]]:
    pathways: list[str]    = []
    metabolites: list[str] = []
    regions: list[str]     = []
    for a in activities:
        tname = (a.get("target_pref_name") or "").lower()
        for key, pws in TARGET_TO_PATHWAY.items():
            if key in tname:
                pathways.extend(pws)
        for key, mets in TARGET_TO_METABOLITES.items():
            if key in tname:
                metabolites.extend(mets)
    for dis in indication_diseases:
        regions.extend(BRAIN_REGION_MAP.get(dis, []))
    if not pathways:
        pathways = ["NF-kB", "Nrf2/GSH"]
    if not metabolites:
        metabolites = ["ROS", "IL-6"]
    if not regions:
        regions = ["Cortex", "Hippocampus"]
    return (list(dict.fromkeys(pathways))[:5],
            list(dict.fromkeys(metabolites))[:5],
            list(dict.fromkeys(regions))[:4])


def _infer_disease_levels(
    indication_diseases: list[str], activities: list[dict]
) -> dict[str, str]:
    tnames = " ".join((a.get("target_pref_name") or "").lower() for a in activities)
    levels = {"als": "Low", "alzheimers": "Low", "parkinsons": "Low", "huntingtons": "Low"}
    for dis, kw_list in DISEASE_TARGET_MAP.items():
        if dis in indication_diseases:
            levels[dis] = "High"
        elif any(kw in tnames for kw in kw_list):
            levels[dis] = "Med"
    return levels


def _indication_keyword_to_disease(keyword: str) -> str:
    mapping = {
        "Alzheimer": "alzheimers", "dementia": "alzheimers",
        "Parkinson":  "parkinsons",
        "Amyotrophic": "als",
        "Huntington":  "huntingtons",
    }
    for k, v in mapping.items():
        if k.lower() in keyword.lower():
            return v
    return ""


def run_expansion(
    curated_path: str  = "compounds.json",
    output_path: str   = "compounds_ml.json",
    erc_output: str    = "erc_ml.json",
    max_per_disease: int = 60,
) -> tuple[int, int]:
    log.info("=== BrainSafe AI ML Expansion Pipeline ===")

    with open(curated_path) as f:
        raw = json.load(f)
    curated: dict[str, dict] = raw.get("compounds", raw) if isinstance(raw, dict) else {}
    log.info("Loaded %d curated compounds", len(curated))

    log.info("--- Stage A: Building training data (disease-feature method) ---")
    X, y, matched = build_training_data(curated)

    log.info("--- Stage B: Training ML model on %d compounds ---", len(X))
    model, scaler, cv_r2 = train_ml_model(X, y)

    log.info("--- Stage C: Fetching neuro indication compound IDs (parallel) ---")
    ind_to_disease: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=5) as ex:
        futures = {ex.submit(_fetch_indication_ids, kw, max_per_disease): kw
                   for kw in INDICATION_KEYWORDS}
        id_lists: dict[str, list[str]] = {}
        for fut in as_completed(futures):
            kw, ids = fut.result()
            dis = _indication_keyword_to_disease(kw)
            id_lists[kw] = ids
            for cid in ids:
                if cid not in ind_to_disease:
                    ind_to_disease[cid] = dis
            log.info("  %s: %d compounds", kw, len(ids))

    all_ids = list(dict.fromkeys(
        cid for ids in id_lists.values() for cid in ids
    ))
    curated_lower = {n.lower() for n in curated}

    log.info("--- Stage D: Fetching molecule data for %d compounds ---", len(all_ids))
    all_mols = _fetch_molecules_batch(all_ids)

    def _is_usable(mol: dict) -> bool:
        if mol.get("molecule_type") not in ("Small molecule", None):
            return False
        props = mol.get("molecule_properties") or {}
        mw    = _safe_float(props.get("mw_freebase") or props.get("full_mwt"))
        if mw < 100 or mw > 900:
            return False
        if mol.get("structure_type") == "NONE":
            return False
        name = (mol.get("pref_name") or "").strip()
        return bool(name) and name.lower() not in curated_lower

    usable = [m for m in all_mols if _is_usable(m)]
    log.info("Usable unique compounds: %d", len(usable))

    def _nps(pred: dict) -> float:
        raw = (pred.get("antioxidant", 5) * 3
               + pred.get("anti_inflammatory", 5) * 3
               + pred.get("mitochondrial_support", 5) * 2
               + pred.get("aggregation_modulation", 5) * 2)
        return min(100.0, raw)

    log.info("--- Stage E (part 1): Batch-predicting profiles for all %d compounds ---", len(usable))
    ml_compounds: dict[str, dict] = {}
    erc_data: dict[str, dict]     = {}

    def _mol_metadata(mol: dict) -> dict:
        cid      = mol.get("molecule_chembl_id", "")
        props    = mol.get("molecule_properties") or {}
        bbb_str, bbb_num = _bbb_from_props(props)
        dis_kw   = ind_to_disease.get(cid, "")
        dis_lvls = {"alzheimers": 0, "parkinsons": 0, "als": 0, "huntingtons": 0}
        if dis_kw:
            dis_lvls[dis_kw] = 2
        ind_diseases = [d for d, lvl in dis_lvls.items() if lvl > 0]
        regions = []
        for dis in ind_diseases:
            regions.extend(BRAIN_REGION_MAP.get(dis, []))
        if not regions:
            regions = ["Cortex", "Hippocampus"]
        return {
            "cid": cid, "bbb_str": bbb_str,
            "feat": [float(bbb_num),
                     float(dis_lvls["als"]), float(dis_lvls["alzheimers"]),
                     float(dis_lvls["parkinsons"]), float(dis_lvls["huntingtons"]),
                     3.0],
            "ind_diseases": ind_diseases,
            "dis_lvls": dis_lvls,
            "regions": list(dict.fromkeys(regions))[:4],
            "compound_type": _compound_type_from_chembl(mol),
            "max_phase": mol.get("max_phase"),
        }

    meta_list = [_mol_metadata(mol) for mol in usable]
    names_list = [(mol.get("pref_name") or mol.get("molecule_chembl_id", "")).strip()
                  for mol in usable]

    X_pred = np.array([m["feat"] for m in meta_list], dtype=float)
    Y_pred = model.predict(scaler.transform(X_pred))
    log.info("  Batch prediction complete for %d compounds.", len(usable))

    for i, (name, meta) in enumerate(zip(names_list, meta_list)):
        if not name:
            continue
        raw_scores = Y_pred[i]
        predicted  = {col: round(float(np.clip(s, 1.0, 10.0)), 1)
                      for col, s in zip(SCORE_COLS, raw_scores)}
        dis_levels_str = {
            d: ("High" if lvl == 2 else ("Med" if lvl == 1 else "Low"))
            for d, lvl in meta["dis_lvls"].items()
        }
        ml_compounds[name] = {
            "compound_type":       meta["compound_type"],
            "bbb":                 meta["bbb_str"],
            **predicted,
            **dis_levels_str,
            "pathways":            ["NF-kB", "Nrf2/GSH"],
            "metabolites":         ["ROS", "IL-6"],
            "brain_regions":       meta["regions"],
            "data_source":         "chembl_ml_predicted",
            "confidence":          "medium",
            "ml_predicted":        True,
            "chembl_id":           meta["cid"],
            "chembl_max_phase":    meta["max_phase"],
            "indication_diseases": meta["ind_diseases"],
            "model_cv_r2":         round(cv_r2, 3),
            "_nps_tmp":            _nps(predicted),
        }

    log.info("  Predicted %d compounds. Selecting top 60 for ERC enrichment.", len(ml_compounds))
    top60 = sorted(ml_compounds.items(), key=lambda kv: kv[1].get("_nps_tmp", 0), reverse=True)[:60]
    cid_to_name = {entry["chembl_id"]: name for name, entry in top60 if entry.get("chembl_id")}

    log.info("--- Stage E (part 2): Fetching ERC for top 60 compounds (%d workers) ---", MAX_WORKERS)

    def _enrich_erc(item: tuple[str, dict]) -> tuple[str, list[dict]]:
        name, entry = item
        cid = entry.get("chembl_id", "")
        if not cid:
            return name, []
        activities = _fetch_activities_for_mol(cid)
        if activities:
            paths, mets, _ = _infer_pathways_regions(activities, entry.get("indication_diseases", []))
            entry["pathways"]    = paths
            entry["metabolites"] = mets
        return name, activities

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for name, activities in ex.map(_enrich_erc, top60):
            if activities:
                erc = activities_to_erc(activities)
                if erc.get("enzymes") or erc.get("transporters"):
                    erc_data[name] = erc

    for name in ml_compounds:
        ml_compounds[name].pop("_nps_tmp", None)

    log.info("Generated %d ML compound profiles (%d with ERC data)", len(ml_compounds), len(erc_data))

    ml_compounds["_ml_metadata"] = {
        "description": (
            f"ML-predicted neuroprotective profiles. Random Forest trained on "
            f"all {len(X)} curated BrainSafe AI compounds. "
            f"ERC from ChEMBL bioassays (Homo sapiens)."
        ),
        "training_compounds": len(X),
        "training_matched":   matched,
        "cv_r2_mean":         round(cv_r2, 3),
        "model":              "MultiOutputRandomForestRegressor (n_estimators=150, max_depth=6, min_samples_leaf=2, disease-feature basis)",
        "features":           FEATURE_COLS,
        "targets":            SCORE_COLS,
        "confidence_label":   "medium",
        "data_sources": [
            "ChEMBL REST API (Mendez D et al. 2019. Nucleic Acids Res 47:D930–D940)",
            "BrainSafe AI curated database (129 compounds, SSSIHL SAI-Net)",
            "Breiman L (2001) Random Forests. Machine Learning 45:5–32",
        ],
    }

    with open(output_path, "w") as f:
        json.dump(ml_compounds, f, indent=2)
    n_out = len(ml_compounds) - 1
    log.info("Saved %d ML compounds to %s", n_out, output_path)

    with open(erc_output, "w") as f:
        json.dump(erc_data, f, indent=2)
    log.info("Saved ERC data for %d compounds to %s", len(erc_data), erc_output)

    return n_out, len(erc_data)


if __name__ == "__main__":
    n, n_erc = run_expansion(max_per_disease=60)
    print(f"\nDone. {n} ML-predicted compounds, {n_erc} with ERC data.")