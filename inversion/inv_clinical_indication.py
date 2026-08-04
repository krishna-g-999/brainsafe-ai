"""H6: do the disease scores predict real clinical indications?

H1 showed the disease layer recovers the disease that a compound's TARGET maps to. That mapping is
this project's own, so H1 established internal consistency rather than external truth. The claim a
user actually relies on is stronger: that a high score for Parkinson's disease says something about
Parkinson's disease as clinicians use the term. Nothing has tested it.

Ground truth here is ChEMBL's drug_indication table, 60,055 records linking drugs to MeSH disease
headings with a maximum clinical phase. It is external to this project and was not used in any
training. Only indications that reached phase 4, meaning approved use, are counted, so the comparison
is against clinical reality rather than aspiration.

The MeSH-to-panel mapping is written out in full below and fixed before any prediction is computed.
It is deliberately narrow: a heading that does not clearly correspond to one of the modelled
conditions is discarded rather than forced into the nearest bucket, because a generous mapping would
manufacture agreement.

Three things could produce apparent success without any real ability, and each is measured:

  permutation null   shuffles which drug carries which indication, destroying the drug-to-disease
                     link while preserving both marginal distributions
  frequency null     always answers with the commonest indications in the set
  memorisation       approved drugs are heavily represented in ChEMBL, so many of these structures
                     were in the training actives. Results are stratified by whether the exact
                     structure appears in the training chemistry, and the stratum that was never
                     seen is the one that carries evidential weight

Read-only with respect to models and data. Writes inversion/results/H6_clinical_indication*.csv
"""
from __future__ import annotations

import json
import pickle
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import urllib3
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
urllib3.disable_warnings()
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))

OUT = ROOT / "inversion" / "results"
CACHE = ROOT / "inversion" / "cache"
OUT.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)
BASE = "https://www.ebi.ac.uk/chembl/api/data"
TOPK = 3
MIN_PHASE = 4.0
N_PERM = 2000
rng = np.random.default_rng(2718)

# Headings that a substring rule would catch wrongly. Wolff-Parkinson-White syndrome is a cardiac
# pre-excitation disorder and has nothing to do with Parkinson's disease; it was found by
# audit_mesh_map.py, which lists what the mapping actually matched instead of assuming it is sound.
MESH_EXCLUDE = ["wolff-parkinson-white"]

# MeSH heading fragment -> panel condition. Fixed before any prediction is computed.
# Narrow by design: an unmatched heading is discarded, never coerced. Generic "dementia" is
# deliberately absent because vascular and frontotemporal dementia are not Alzheimer's disease.
MESH_MAP = {
    "alzheimer": "Alzheimer's disease",
    "parkinson": "Parkinson's disease",
    "amyotrophic lateral sclerosis": "Amyotrophic lateral sclerosis",
    "huntington": "Huntington's disease",
    "depressive": "Depression / anxiety",
    "depression": "Depression / anxiety",
    "anxiety": "Depression / anxiety",
    "obsessive-compulsive": "Depression / anxiety",
    "panic disorder": "Depression / anxiety",
    "schizophren": "Psychosis / schizophrenia",
    "psychotic": "Psychosis / schizophrenia",
    "bipolar": "Psychosis / schizophrenia",
    "epilep": "Epilepsy",
    "seizure": "Epilepsy",
    "status epilepticus": "Epilepsy",
    "pain": "Chronic pain",
    "migraine": "Chronic pain",
    "neuralgia": "Chronic pain",
    "alcoholism": "Addiction",
    "opioid-related disorders": "Addiction",
    "substance withdrawal": "Addiction",
    "tobacco use disorder": "Addiction",
    "smoking cessation": "Addiction",
    "attention deficit": "ADHD",
    "narcolep": "Sleep / wakefulness",
    "sleep initiation and maintenance": "Sleep / wakefulness",
    "insomnia": "Sleep / wakefulness",
}


def get_json(url, params=None, tries=4):
    for _ in range(tries):
        try:
            r = requests.get(url, params=params, timeout=90, verify=False)
            if r.status_code == 200:
                return r.json()
        except Exception:
            pass
        time.sleep(3)
    return None


def fetch_indications():
    f = CACHE / "chembl_indications.json"
    if f.exists():
        return json.loads(f.read_text())
    out, url, page = {}, f"{BASE}/drug_indication.json?limit=1000", 0
    while url:
        j = get_json(url)
        if not j:
            break
        for d in j.get("drug_indications", []):
            try:
                if float(d.get("max_phase_for_ind") or 0) < MIN_PHASE:
                    continue
            except (TypeError, ValueError):
                continue
            mesh = (d.get("mesh_heading") or "").lower()
            if any(x in mesh for x in MESH_EXCLUDE):
                continue
            hit = next((v for k, v in MESH_MAP.items() if k in mesh), None)
            if hit:
                out.setdefault(d["molecule_chembl_id"], set()).add(hit)
        page += 1
        if page % 10 == 0:
            print(f"  page {page}: {len(out)} drugs matched so far", flush=True)
        nxt = j["page_meta"]["next"]
        url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
    out = {k: sorted(v) for k, v in out.items()}
    f.write_text(json.dumps(out))
    return out


def fetch_smiles(ids):
    f = CACHE / "chembl_drug_smiles.json"
    known = json.loads(f.read_text()) if f.exists() else {}
    todo = [i for i in ids if i not in known]
    for k in range(0, len(todo), 40):
        chunk = ",".join(todo[k:k + 40])
        j = get_json(f"{BASE}/molecule.json", {"molecule_chembl_id__in": chunk, "limit": 40})
        if j:
            for m in j.get("molecules", []):
                st = m.get("molecule_structures") or {}
                known[m["molecule_chembl_id"]] = st.get("canonical_smiles")
        if (k // 40) % 5 == 0:
            print(f"  structures resolved: {len(known)}/{len(ids)}", flush=True)
        time.sleep(0.15)
    f.write_text(json.dumps(known))
    return known


def training_structures():
    """Canonical SMILES of every structure the deployed models or the read-across index have seen."""
    seen = set()
    p = ROOT / "models_rf" / "ad_reference.pkl"
    if p.exists():
        with p.open("rb") as fh:
            ref_smiles, _ = pickle.load(fh)
        for s in ref_smiles:
            m = Chem.MolFromSmiles(s)
            if m is not None:
                seen.add(Chem.MolToSmiles(m))
    q = ROOT / "models_rf" / "readacross_index.pkl"
    if q.exists():
        idx = pickle.loads(q.read_bytes())
        for s in idx["smiles"]:
            m = Chem.MolFromSmiles(s)
            if m is not None:
                seen.add(Chem.MolToSmiles(m))
    return seen


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return max(0.0, c - h), min(1.0, c + h)


def accuracy_block(label, sub, col="predicted_top3", hitcol="hit"):
    """Observed accuracy plus both nulls, computed within one stratum."""
    n = len(sub)
    if n == 0:
        return None
    k = int(sub[hitcol].sum())
    lo, hi = wilson(k, n)
    truth = [set(t.split(" | ")) for t in sub.true_indications]
    preds = [set(p.split(" | ")) if p != "NONE" else set() for p in sub[col]]
    perm = np.array([np.mean([bool(preds[i] & t) for i, t in enumerate(rng.permutation(truth))])
                     for _ in range(N_PERM)])
    p_perm = float(((perm >= k / n).sum() + 1) / (N_PERM + 1))
    cnt = Counter(d for t in truth for d in t)
    common = {d for d, _ in cnt.most_common(TOPK)}
    freq = float(np.mean([bool(common & t) for t in truth]))
    return {"stratum": label, "n": n, "top3_accuracy": round(k / n, 4),
            "ci95_low": round(lo, 4), "ci95_high": round(hi, 4),
            "permutation_null": round(float(perm.mean()), 4),
            "permutation_p": round(p_perm, 4),
            "frequency_null": round(freq, 4)}


def main():
    import app  # read-only use of the deployed models and knowledge graph

    print("fetching approved indications from ChEMBL ...", flush=True)
    ind = fetch_indications()
    print(f"  {len(ind)} drugs carry a phase-4 indication that maps to the panel", flush=True)
    if not ind:
        print("nothing to test; aborting")
        return

    print("resolving structures ...", flush=True)
    smi = fetch_smiles(sorted(ind))

    print("indexing training chemistry ...", flush=True)
    seen = training_structures()
    print(f"  {len(seen):,} distinct structures already seen by the models", flush=True)

    # Scoring these drugs across 63 models with a full applicability-domain search takes about an
    # hour, which makes reruns expensive enough to discourage checking. What is cached is the full
    # ranked disease vector per compound, which depends only on the structure, so a correction to
    # the indication map costs seconds rather than another hour. Delete the file to rescore.
    scored = json.loads((CACHE / "H6_scored.json").read_text()) if (CACHE / "H6_scored.json").exists() else {}
    todo = [c for c in ind if c not in scored and smi.get(c)]
    if todo:
        models = app.load_models()
        for i, cid in enumerate(sorted(todo), 1):
            m = Chem.MolFromSmiles(smi[cid])
            if m is None:
                continue
            cs = Chem.MolToSmiles(m)
            r = app.predict_all(cs, models)
            if r is None:
                continue
            _bbb, _neuro, dz = app.disease_scores(r)
            ad = app.assess_domain(cs)
            scored[cid] = {"smiles": cs,
                           "ranking": [[d["disease"], round(float(d["gated"]), 6)] for d in dz],
                           "in_training": int(cs in seen),
                           "ad_max_tanimoto": round(ad["max_sim"], 3) if ad else None}
            if i % 100 == 0:
                print(f"  scored {i}/{len(todo)}", flush=True)
        (CACHE / "H6_scored.json").write_text(json.dumps(scored))
    else:
        print(f"reusing {len(scored)} cached rankings", flush=True)

    rows = []
    for cid, diseases in sorted(ind.items()):
        sc = scored.get(cid)
        if not sc:
            continue
        dz = sc["ranking"]
        top = [d for d, v in dz[:TOPK] if v >= app.MIN_ACTIONABLE_SCORE]
        # the tool reports nothing below MIN_ACTIONABLE_SCORE, which is correct behaviour but makes
        # it silent where a baseline that always answers still scores. The unthresholded ranking is
        # kept as well, so ranking ability can be judged separately from the reporting rule.
        raw = [d for d, _ in dz[:TOPK]]
        rank = next((i + 1 for i, (d, _) in enumerate(dz) if d in diseases), None)
        rows.append({
            "chembl_id": cid, "smiles": sc["smiles"],
            "true_indications": " | ".join(diseases),
            "predicted_top3": " | ".join(top) if top else "NONE",
            "predicted_top3_unthresholded": " | ".join(raw),
            "hit": int(bool(set(top) & set(diseases))),
            "hit_unthresholded": int(bool(set(raw) & set(diseases))),
            "rank_of_true": rank,
            "any_signal": int(bool(top)),
            "in_training": sc["in_training"],
            "top_score": round(dz[0][1], 3),
            "ad_max_tanimoto": sc["ad_max_tanimoto"],
        })

    if not rows:
        print("no scorable drugs; aborting")
        return
    return analyse(pd.DataFrame(rows), app)


def analyse(df, app):
    df.to_csv(OUT / "H6_clinical_indication_predictions.csv", index=False)

    blocks = [accuracy_block("ALL", df),
              accuracy_block("never seen in training", df[df.in_training == 0]),
              accuracy_block("present in training chemistry", df[df.in_training == 1]),
              accuracy_block("tool reported a signal", df[df.any_signal == 1]),
              accuracy_block("never seen AND signal reported",
                             df[(df.in_training == 0) & (df.any_signal == 1)]),
              accuracy_block("ALL, ranking only, reporting threshold removed", df,
                             "predicted_top3_unthresholded", "hit_unthresholded"),
              accuracy_block("never seen, ranking only", df[df.in_training == 0],
                             "predicted_top3_unthresholded", "hit_unthresholded")]
    res = pd.DataFrame([b for b in blocks if b])
    res.to_csv(OUT / "H6_clinical_indication.csv", index=False)

    # where in the ranking the true indication actually falls, for drugs that were ranked at all
    rk = df.rank_of_true.dropna().astype(int)
    if len(rk):
        mrr = float((1.0 / rk).sum() / len(df))
        print(f"\nmedian rank of the true indication: {int(rk.median())} of {len(app.DISEASE_ORDER)}; "
              f"mean reciprocal rank {mrr:.3f}", flush=True)

    # per-disease breakdown, so a single dominant class cannot carry the headline number
    per = []
    for d in app.DISEASE_ORDER:
        sub = df[df.true_indications.str.contains(d, regex=False)]
        if len(sub) >= 5:
            per.append({"indication": d, "n": len(sub),
                        "recovered_in_top3": round(float(
                            sub.predicted_top3.str.contains(d, regex=False).mean()), 3),
                        "recovered_ranking_only": round(float(
                            sub.predicted_top3_unthresholded.str.contains(d, regex=False).mean()), 3),
                        "silent": int((sub.any_signal == 0).sum())})
    if per:
        pd.DataFrame(per).to_csv(OUT / "H6_clinical_per_indication.csv", index=False)

    pd.set_option("display.width", 200)
    print()
    print(res.to_string(index=False))
    if per:
        print()
        print(pd.DataFrame(per).to_string(index=False))

    key = next((b for b in blocks if b and b["stratum"] == "never seen in training"), blocks[0])
    verdict = ("SUPPORTED" if (key["top3_accuracy"] > key["permutation_null"] + 0.05
                               and key["permutation_p"] < 0.05
                               and key["top3_accuracy"] > key["frequency_null"])
               else "WEAKENED" if key["top3_accuracy"] > key["permutation_null"] else "REFUTED")
    print(f"\nVERDICT H6 (judged on the never-seen stratum, n={key['n']}): {verdict}")
    print("wrote", OUT / "H6_clinical_indication.csv")


if __name__ == "__main__":
    main()
