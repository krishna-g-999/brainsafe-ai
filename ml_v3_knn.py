"""
BrainSafe AI v3 — Correct Architecture
NPS: Similarity-weighted KNN (Tanimoto)
BBB: Physicochemical rules
Disease/Pathway: Target mapping + keyword system
"""
import json, numpy as np, joblib, os, time
import requests, urllib.parse
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

# ── Build similarity index from curated compounds ─────────────────────────────
def build_index(compounds_json='compounds.json', ml_json='compounds_ml.json',
                out='models_v3/knn_index.pkl'):
    NPS_DIMS = ['antioxidant','anti_inflammatory','mitochondrial_support',
                'aggregation_modulation','cognitive_enhancement',
                'neurogenesis','synaptic_plasticity']
    fps, names, nps_vals = [], [], []
    for path in [compounds_json, ml_json]:
        with open(path) as f: data = json.load(f)
        items = list(data.items()) if isinstance(data, dict) else list(enumerate(data))
        for name, item in items:
            cname = item.get('name', str(name))
            smiles = item.get('smiles') or item.get('canonical_smiles','')
            # Use stored NPS dims (no API needed — already in DB)
            dim_vals = [float(item.get(d, 0) or 0) for d in NPS_DIMS]
            nps = round(np.mean(dim_vals) * 10, 1)
            if smiles:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        gen = AllChem.GetMorganGenerator(radius=2, fpSize=1024)
                        fp = np.array(gen.GetFingerprint(mol), dtype=np.float32)
                        fps.append(fp); names.append(cname); nps_vals.append(nps)
                except: pass

    os.makedirs('models_v3', exist_ok=True)
    joblib.dump({'fps': fps, 'names': names, 'nps': nps_vals}, out)
    print(f"✅ KNN index built: {len(fps)} compounds with SMILES in DB")
    return fps, names, nps_vals

def tanimoto(a, b):
    """Tanimoto similarity for binary fingerprints."""
    intersection = np.dot(a, b)
    return intersection / (np.sum(a) + np.sum(b) - intersection + 1e-9)

def predict_nps_knn(smiles, k=3, index_path='models_v3/knn_index.pkl'):
    """Similarity-weighted KNN NPS prediction."""
    if not os.path.exists(index_path):
        return None, [], 0.0
    idx = joblib.load(index_path)
    try:
        mol = Chem.MolFromSmiles(smiles)
        gen = AllChem.GetMorganGenerator(radius=2, fpSize=1024)
        query_fp = np.array(gen.GetFingerprint(mol), dtype=np.float32)
    except:
        return None, [], 0.0

    sims = [tanimoto(query_fp, fp) for fp in idx['fps']]
    top_k_idx = np.argsort(sims)[-k:][::-1]

    top_sims  = [sims[i] for i in top_k_idx]
    top_names = [idx['names'][i] for i in top_k_idx]
    top_nps   = [idx['nps'][i] for i in top_k_idx]

    # Weighted mean
    total_sim = sum(top_sims)
    if total_sim < 0.01:
        nps_pred = round(np.mean(top_nps), 1)
    else:
        nps_pred = round(sum(n*s for n,s in zip(top_nps, top_sims)) / total_sim, 1)

    confidence = "High" if top_sims[0] > 0.6 else "Med" if top_sims[0] > 0.4 else "Low"
    neighbours = [f"{n} (sim={s:.2f}, NPS={p})"
                  for n,s,p in zip(top_names, top_sims, top_nps)]
    return nps_pred, neighbours, top_sims[0], confidence

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'build':
        build_index()
    else:
        # Quick test
        idx = build_index()
        test_smiles = {
            'Quercetin':      'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O',
            'Pterostilbene':  'COC1=CC(=CC(=C1)/C=C/C2=CC=C(C=C2)O)OC',
            'Urolithin A':    'C1=CC2=C(C=C1O)C(=O)OC3=C2C=CC(=C3)O',
        }
        for name, smi in test_smiles.items():
            nps, neighbours, max_sim, conf = predict_nps_knn(smi)
            print(f"\n{name}:")
            print(f"  NPS estimate: {nps}/100  (confidence: {conf})")
            for n in neighbours:
                print(f"  Neighbour: {n}")

# ── Final prediction function (used by app.py) ────────────────────────────────
def predict_compound_full(compound_name, smiles=None):
    """
    Full v3 prediction: NPS (KNN) + BBB (rules) + Disease + Pathway + Similarity
    """
    from ml_v3_engine import (get_smiles_pubchem, predict_disease_relevance,
                               predict_pathways, get_chembl_targets)
    import time

    # 1. Get SMILES
    if not smiles:
        smiles, props = get_smiles_pubchem(compound_name)
    else:
        props = {}
    if not smiles:
        return {"error": f"SMILES not found for {compound_name}"}

    # 2. NPS via KNN
    nps, neighbours, max_sim, confidence = predict_nps_knn(smiles)

    # 3. BBB via physicochemical rules
    mw   = float(props.get('MolecularWeight', 500))
    logp = float(props.get('XLogP', 3))
    tpsa = float(props.get('TPSA', 90))
    if mw < 400 and 1 <= logp <= 4 and tpsa < 75:   bbb = 'High'
    elif mw < 500 and logp <= 5 and tpsa < 90:       bbb = 'Medium'
    elif mw < 600 and tpsa < 120:                     bbb = 'Low-Med'
    else:                                              bbb = 'Low'

    # 4. Disease + Pathway
    targets  = get_chembl_targets(compound_name)
    diseases = predict_disease_relevance(compound_name, smiles, targets)
    pathways = predict_pathways(compound_name, targets)

    return {
        'compound':    compound_name,
        'nps':         nps,
        'confidence':  confidence,
        'bbb':         bbb,
        'mw':          mw, 'logp': logp, 'tpsa': tpsa,
        'diseases':    diseases,
        'pathways':    pathways,
        'neighbours':  neighbours,
        'max_sim':     round(max_sim, 3),
    }
