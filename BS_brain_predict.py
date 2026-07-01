"""
BS_brain_predict.py — GENUINE, evidence-grounded, BBB-gated multi-endpoint
brain-effect engine (the BrainSafe innovation layer).

For a SMILES it returns, all from models trained on MEASURED data:
  * per-endpoint CALIBRATED probability (isotonic on scaffold-CV OOF), class call,
    honest AUROC, applicability-domain flag, and the NEAREST REAL MEASURED ANALOGS
    (so every prediction is backed by actual data, not a black box);
  * BBB-GATED disease scores (a CNS target hit only matters if the compound reaches
    the brain): effective engagement = P(target) x P(BBB);
  * a SAFETY axis (hERG cardiotox liability);
  * the genuine antioxidant model + deterministic druggability/CNS-MPO;
  * an overall, transparent benefit/risk verdict.
"""
from __future__ import annotations
import os, glob, json
import numpy as np

_DIR = os.path.dirname(os.path.abspath(__file__))
_BRAIN = os.path.join(_DIR, "models_brain")
_BUNDLES, _META, _CONF = {}, {}, {}
AD_MIN = 0.30
# Fact-based quality gate: only deploy endpoints whose scaffold-CV MCC clears this
# bar. Excludes near-trivial, severely class-imbalanced receptor models (D2/A2A/
# HT2A/SERT, MCC 0.21-0.44) where binary active/inactive QSAR is ill-posed because
# ChEMBL reports almost only actives. BBB (gate) and hERG (safety) are always kept.
MIN_MCC = 0.45
_ALWAYS = {"BBB", "hERG"}

# Brain-effect semantics
THERAPEUTIC = {"AChE", "BChE", "BACE1", "MAO_B", "MAO_A", "GSK3B", "D2", "A2A", "HT2A", "SERT"}
SAFETY = {"hERG"}
DISEASE_MAP = {
    "Alzheimer's disease":     ["AChE", "BChE", "BACE1"],
    "Parkinson's disease":     ["MAO_B", "D2", "A2A"],
    "Depression / mood":       ["MAO_A", "SERT", "HT2A"],
    "Neuroprotection (tau)":   ["GSK3B"],
}


def _load():
    global _BUNDLES, _META
    if _BUNDLES:
        return
    import joblib
    for f in glob.glob(os.path.join(_BRAIN, "*.joblib")):
        name = os.path.basename(f).replace(".joblib", "")
        try:
            mp = os.path.join(_BRAIN, f"{name}_meta.json")
            meta = json.load(open(mp)) if os.path.exists(mp) else {}
            # quality gate (fact-based): drop weak/ill-posed endpoints
            if name not in _ALWAYS and float(meta.get("mcc", 0.0)) < MIN_MCC:
                continue
            _BUNDLES[name] = joblib.load(f)
            _META[name] = meta
            cp = os.path.join(_BRAIN, f"{name}_conformal.json")
            _CONF[name] = json.load(open(cp)) if os.path.exists(cp) else None
        except Exception:
            pass


def _conformal_set(raw, conf):
    """Mondrian conformal prediction set + class p-values from the raw ensemble score."""
    import bisect
    if not conf:
        return None
    eps = conf.get("eps", 0.10)
    A1, A0 = conf.get("calib_active_nonconf", []), conf.get("calib_inactive_nonconf", [])
    pv_act = ((len(A1) - bisect.bisect_left(A1, 1 - raw)) + 1) / (len(A1) + 1) if A1 else 1.0
    pv_in = ((len(A0) - bisect.bisect_left(A0, raw)) + 1) / (len(A0) + 1) if A0 else 1.0
    pset = [c for c, pv in (("active", pv_act), ("inactive", pv_in)) if pv > eps]
    if pset == ["active"]:
        label = f"Confident active ({int((1-eps)*100)}% conformal)"
    elif pset == ["inactive"]:
        label = f"Confident inactive ({int((1-eps)*100)}%)"
    elif set(pset) == {"active", "inactive"}:
        label = "Uncertain (both plausible)"
    else:
        label = "Atypical (out-of-distribution)"
    return {"set": pset, "p_active": round(pv_act, 3), "p_inactive": round(pv_in, 3), "label": label}


def _features(smiles):
    from BS_predictive_model import morgan, descriptors
    return np.hstack([morgan([smiles]), descriptors([smiles])])


def _query_fp(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024) if m else None


def _endpoint_predict(name, bundle, X, qfp):
    from rdkit import DataStructs
    raw = float(np.mean([m.predict_proba(X)[:, 1][0] for m in bundle["models"].values()]))
    cal = bundle.get("calibrator")
    p = float(cal.predict([raw])[0]) if cal is not None else raw
    thr = bundle.get("threshold", 0.5)
    conformal = _conformal_set(raw, _CONF.get(name))
    ev = bundle.get("evidence", {})
    sims, analogs = [], []
    if qfp is not None and ev.get("fps"):
        sims = np.array(DataStructs.BulkTanimotoSimilarity(qfp, ev["fps"]))
        top = np.argsort(sims)[::-1][:3]
        for i in top:
            analogs.append({"smiles": ev["smiles"][i], "similarity": round(float(sims[i]), 2),
                            "measured": ("active" if ev["label"][i] == 1 else "inactive"),
                            "pchembl": ev["value"][i]})
    ad = round(float(max(sims)), 3) if len(sims) else 0.0
    return {"probability": round(p, 3), "threshold": round(thr, 3),
            "applicability_tanimoto": ad, "in_domain": ad >= AD_MIN, "analogs": analogs,
            "conformal": conformal}


def predict_brain_profile(smiles: str | None) -> dict:
    _load()
    out = {"ok": False, "endpoints": [], "diseases": [], "safety": None,
           "antioxidant": None, "druggability": None, "verdict": "", "bbb_p": None}
    if not smiles or str(smiles).strip().lower() in ("", "nan", "n/a", "none"):
        return out
    from rdkit import Chem
    if Chem.MolFromSmiles(str(smiles)) is None:
        return out
    out["ok"] = True
    X = _features(str(smiles)); qfp = _query_fp(str(smiles))

    P = {}   # endpoint -> calibrated probability
    for name, bundle in _BUNDLES.items():
        try:
            r = _endpoint_predict(name, bundle, X, qfp)
        except Exception:
            continue
        meta = _META.get(name, {})
        P[name] = r["probability"]
        pos = "Penetrant" if name == "BBB" else ("Blocker" if name in SAFETY else "Active")
        neg = "Non-penetrant" if name == "BBB" else ("Low-risk" if name in SAFETY else "Inactive")
        out["endpoints"].append({
            "endpoint": name, "meaning": meta.get("meaning", name),
            "probability": r["probability"], "call": pos if r["probability"] >= r["threshold"] else neg,
            "auroc": meta.get("auroc"), "brier": meta.get("brier"), "n_train": meta.get("n"),
            "in_domain": r["in_domain"], "applicability_tanimoto": r["applicability_tanimoto"],
            "kind": "safety" if name in SAFETY else ("gate" if name == "BBB" else "target"),
            "evidence": r["analogs"], "source": meta.get("source", ""),
            "conformal": r.get("conformal"),
        })
    order = {"BBB": 0, "AChE": 1, "BChE": 2, "BACE1": 3, "GSK3B": 4, "MAO_B": 5,
             "MAO_A": 6, "D2": 7, "A2A": 8, "HT2A": 9, "SERT": 10, "hERG": 20}
    out["endpoints"].sort(key=lambda e: order.get(e["endpoint"], 8))

    bbb_p = P.get("BBB", None)
    out["bbb_p"] = bbb_p

    # --- engagement map: calibrated classifier probs + regression-derived engagement
    # for receptor targets (pKi 5->0, 8->1). Regression contributions are ranking-grade. ---
    ENG = dict(P)
    reg = _receptor_potency(X)
    reg_targets = set()
    for r in reg:
        ENG[r["target"]] = float(np.clip((r["predicted_pKi"] - 5.0) / 3.0, 0.0, 1.0))
        reg_targets.add(r["target"])

    # --- BBB-gated per-disease synthesis ---
    for disease, targets in DISEASE_MAP.items():
        present = [t for t in targets if t in ENG]
        if not present:
            continue
        best_t = max(present, key=lambda t: ENG[t])
        gated = ENG[best_t] * (bbb_p if bbb_p is not None else 1.0)
        out["diseases"].append({
            "disease": disease, "score": round(gated, 3),
            "driver": best_t, "driver_p": round(ENG[best_t], 3),
            "driver_kind": "regression(pKi)" if best_t in reg_targets else "calibrated-classifier",
            "bbb_gated": bbb_p is not None,
            "basis": {t: round(ENG[t], 2) for t in present},
        })
    out["diseases"].sort(key=lambda d: d["score"], reverse=True)

    # --- safety axis (hERG: active == liability) ---
    if "hERG" in P:
        hp = P["hERG"]
        out["safety"] = {"herg_p": round(hp, 3),
                         "risk": "High" if hp >= 0.5 else ("Moderate" if hp >= 0.3 else "Low")}

    # --- antioxidant + druggability ---
    try:
        from BS_predict import predict_antioxidant
        out["antioxidant"] = predict_antioxidant(str(smiles))
    except Exception:
        pass
    try:
        from BS_druggability import compute_druggability
        d = compute_druggability(str(smiles))
        out["druggability"] = d if d.get("ok") else None
    except Exception:
        pass

    # --- receptor binding potency (regression; ranking-grade) ---
    out["receptor_potency"] = _receptor_potency(X)

    # --- clinical/translational precedent (real ChEMBL clinical-phase data) ---
    try:
        from BS_clinical_evidence import clinical_analogs
        out["clinical_evidence"] = clinical_analogs(str(smiles))
    except Exception:
        out["clinical_evidence"] = []

    # --- transparent overall verdict ---
    out["verdict"] = _verdict(out)
    return out


_REG, _REGMETA = {}, {}
def _receptor_potency(X):
    """Predicted pKi/pIC50 for receptor targets via measured-data regression ensembles
    (D2/A2A/5-HT2A/SERT). Ranking-grade; temporal generalisation is weak (reported)."""
    import glob as _glob, joblib as _jl
    rdir = os.path.join(_DIR, "models_brain_reg")
    if not _REG and os.path.isdir(rdir):
        for f in _glob.glob(os.path.join(rdir, "*.joblib")):
            nm = os.path.basename(f).replace(".joblib", "")
            try:
                _REG[nm] = _jl.load(f)
                mp = os.path.join(rdir, f"{nm}_meta.json")
                _REGMETA[nm] = json.load(open(mp)) if os.path.exists(mp) else {}
            except Exception:
                pass
    res = []
    for nm, b in _REG.items():
        try:
            pki = float(np.mean([m.predict(X)[0] for m in b["models"]]))
            mt = _REGMETA.get(nm, {})
            res.append({"target": nm, "meaning": mt.get("meaning", nm),
                        "predicted_pKi": round(pki, 2), "scaffold_cv_r2": mt.get("scaffold_cv_r2"),
                        "spearman": mt.get("spearman"),
                        "temporal_r2": (mt.get("temporal") or {}).get("r2")})
        except Exception:
            continue
    return sorted(res, key=lambda r: -r["predicted_pKi"])


def _verdict(out):
    if not out["diseases"]:
        return "No disease-relevant CNS target activity predicted."
    top = out["diseases"][0]
    drug = out.get("druggability") or {}
    saf = out.get("safety") or {}
    bbb_ok = (out["bbb_p"] or 0) >= 0.5
    parts = []
    if top["score"] >= 0.5:
        parts.append(f"Most promising for **{top['disease']}** "
                     f"(BBB-gated engagement {top['score']*100:.0f}%, via {top['driver'].replace('_','-')}).")
    elif top["score"] >= 0.3:
        parts.append(f"Possible {top['disease']} relevance "
                     f"({top['score']*100:.0f}% BBB-gated engagement) — weak/uncertain.")
    else:
        parts.append("Low predicted CNS-target engagement across the modelled diseases.")
    if not bbb_ok:
        parts.append("Limited brain penetration may blunt central effects.")
    if drug.get("druggability") is not None:
        parts.append(f"Druggability {drug['druggability']:.0f}/100.")
    if saf.get("risk") in ("High", "Moderate"):
        parts.append(f"⚠ {saf['risk']} hERG (cardiotoxicity) risk.")
    ce = out.get("clinical_evidence") or []
    if ce and ce[0]["similarity"] >= 0.5:
        c = ce[0]
        parts.append(f"Clinical precedent: structurally close (T={c['similarity']}) to {c['name']} "
                     f"(phase {c['max_phase']}, {c['disease']}).")
    return " ".join(parts)


def available_endpoints():
    _load(); return sorted(_BUNDLES.keys())


if __name__ == "__main__":
    import sys
    smi = sys.argv[1] if len(sys.argv) > 1 else "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2"  # donepezil
    print(json.dumps(predict_brain_profile(smi), indent=2, default=str)[:2500])
