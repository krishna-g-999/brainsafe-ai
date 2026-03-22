"""
BrainSafe AI v3 — SMILES-based NPS + Disease + Pathway Predictor
Features: Morgan fingerprints (2048-bit) + physicochemical descriptors
"""

import requests, json, time, numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score
import joblib, os, pandas as pd

# ── NDD Target → Disease Map ──────────────────────────────────────────────────
NDD_TARGET_DISEASE = {
    # Alzheimer's
    'APP':     'alzheimers', 'BACE1': 'alzheimers', 'BACE2': 'alzheimers',
    'ACHE':    'alzheimers', 'BCHE':  'alzheimers', 'MAPT':  'alzheimers',
    'PSEN1':   'alzheimers', 'PSEN2': 'alzheimers', 'APOE':  'alzheimers',
    'NMDA':    'alzheimers', 'CDK5':  'alzheimers',
    # Parkinson's
    'LRRK2':   'parkinsons', 'SNCA':  'parkinsons', 'PINK1': 'parkinsons',
    'PARK2':   'parkinsons', 'MAOB':  'parkinsons', 'MAOA':  'parkinsons',
    'DRD2':    'parkinsons', 'DRD3':  'parkinsons', 'TH':    'parkinsons',
    # ALS
    'SOD1':    'als', 'TARDBP': 'als', 'FUS': 'als', 'C9ORF72': 'als',
    'SQSTM1':  'als', 'OPTN': 'als',
    # Huntington's
    'HTT':     'huntingtons', 'MTOR': 'huntingtons', 'HDAC':  'huntingtons',
    'PGC1A':   'huntingtons', 'CASP3':'huntingtons',
    # Multi-disease
    'NRF2':    'multi', 'NFKB1': 'multi', 'TNF':   'multi',
    'IL6':     'multi', 'CASP9': 'multi', 'BCL2':  'multi',
    'TP53':    'multi', 'SIRT1': 'multi', 'AMPK':  'multi',
    'MTOR':    'multi', 'BECN1': 'multi',
}

# ── NDD Pathway Map ───────────────────────────────────────────────────────────
NDD_PATHWAYS = {
    'Nrf2/GSH':        ['nrf2','nfe2l2','keap1','gsh','glutathione','hmox1','nqo1'],
    'NF-kB':           ['nfkb','ikk','tnf','il6','il1b','inflammation','neuroinflam'],
    'PI3K/Akt/mTOR':   ['pi3k','akt','mtor','pten','s6k','4ebp1','rapamycin'],
    'Autophagy':       ['autophagy','beclin','atg','lc3','p62','sqstm','ubiquitin'],
    'Mitophagy/PINK1': ['pink1','parkin','mitophagy','mitochondria','drp1','mfn'],
    'Nrf2/HO-1':       ['ho-1','hmox','heme oxygenase','nrf2'],
    'Wnt/β-catenin':   ['wnt','beta-catenin','gsk3','tcf','dishevelled'],
    'BDNF/TrkB':       ['bdnf','trkb','ntrk2','ngf','neurotrophin','creb'],
    'Sirtuin/AMPK':    ['sirt1','sirt3','ampk','nad','nampt','pgc1'],
    'Caspase/Apoptosis':['caspase','casp3','casp9','bcl2','bax','cytochrome c'],
    'Tau/Aggregation': ['tau','mapt','tauopathy','aggregation','fibrils','amyloid'],
    'Alpha-syn':       ['alpha-synuclein','snca','lewy','synuclein'],
    'AChE/Cholinergic':['ache','acetylcholinesterase','choline','cholinergic','bche'],
    'Dopamine':        ['dopamine','drd2','drd3','tyrosine hydroxylase','maob'],
    'Glutamate/NMDA':  ['nmda','glutamate','excitotoxicity','ampa','glu'],
    'GSH/NAD+/ATP':    ['gsh','nad','atp','coenzyme','redox','ros','oxidative'],
}

# ── Feature Extraction ────────────────────────────────────────────────────────
def smiles_to_features(smiles, radius=2, nbits=2048):
    """Morgan fingerprint + physicochemical descriptors from SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        # Morgan fingerprint
        fp = AllChem.GetMorganGenerator(radius=radius, fpSize=nbits).GetFingerprint(mol)
        fp_arr = np.array(fp)
        # Physicochemical
        physchem = np.array([
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.TPSA(mol),
            rdMolDescriptors.CalcNumHBD(mol),
            rdMolDescriptors.CalcNumHBA(mol),
            rdMolDescriptors.CalcNumRotatableBonds(mol),
            rdMolDescriptors.CalcNumAromaticRings(mol),
            Descriptors.NumHeteroatoms(mol),
        ])
        return np.concatenate([fp_arr, physchem])
    except:
        return None

def get_smiles_pubchem(name):
    """Fetch canonical SMILES from PubChem by compound name."""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/property/IsomericSMILES,MolecularWeight,XLogP,TPSA/JSON"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            props = r.json()['PropertyTable']['Properties'][0]
            return props.get('IsomericSMILES') or props.get('SMILES'), props
    except: pass
    return None, {}

def get_chembl_targets(smiles_or_name):
    """Fetch ChEMBL target annotations for a compound."""
    targets = []
    try:
        # Search by name first
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?pref_name__iexact={requests.utils.quote(smiles_or_name)}&format=json"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            if data.get('molecules'):
                chembl_id = data['molecules'][0]['molecule_chembl_id']
                # Get activities
                act_url = f"https://www.ebi.ac.uk/chembl/api/data/activity?molecule_chembl_id={chembl_id}&format=json&limit=50"
                r2 = requests.get(act_url, timeout=10)
                if r2.status_code == 200:
                    for act in r2.json().get('activities', []):
                        tgt = act.get('target_pref_name', '')
                        gene = act.get('target_organism', '')
                        if tgt: targets.append(tgt.upper())
    except: pass
    return targets

def predict_disease_relevance(compound_name, smiles, chembl_targets=[]):
    """
    Predict disease relevance from:
    1. ChEMBL full target name keyword matching (gene symbols + full names)
    2. Compound name/class keyword matching
    3. Compound-class broad rules (flavonoids, polyphenols → multi-disease)
    Returns: dict {disease: 'High'/'Med'/'Low'}
    """
    disease_scores = {'alzheimers': 0, 'parkinsons': 0, 'als': 0, 'huntingtons': 0}

    # ── ChEMBL target full-name keyword map ──
    TARGET_KW_DISEASE = {
        'alzheimers': [
            'acetylcholinesterase','butyrylcholinesterase','beta-secretase','bace',
            'amyloid','presenilin','tau','mapt','cdk5','nmda','glutamate receptor',
            'apolipoprotein','app','abeta','cholinesterase','ache','bche',
        ],
        'parkinsons': [
            'tyrosine hydroxylase','monoamine oxidase','dopamine','lrrk2','alpha-synuclein',
            'synuclein','pink1','parkin','drd2','drd3','maob','maoa','nurr1',
            'dopaminergic','catechol','comt',
        ],
        'als': [
            'superoxide dismutase','sod1','tdp-43','tardbp','fus','c9orf72',
            'optineurin','sqstm1','p62','motor neuron','neurofilament','ubiquilin',
        ],
        'huntingtons': [
            'huntingtin','htt','striatum','hdac','histone deacetylase','mtor',
            'pgc-1','caspase','mutant','polyglutamine','caudate',
        ],
        'multi': [
            'nrf2','nuclear factor','nf-kb','nfkb','tnf','interleukin','il-6','il-1',
            'sirt1','sirt3','ampk','autophagy','beclin','mitochondri','ros',
            'oxidative stress','neuroinflammation','apoptosis','bcl-2','caspase-3',
            'pi3k','akt','mtor','bdnf','ngf','trkb','neurotrophin','neuroprotect',
        ],
    }
    target_text = ' '.join(chembl_targets).lower()
    for disease, keywords in TARGET_KW_DISEASE.items():
        for kw in keywords:
            if kw in target_text:
                if disease == 'multi':
                    for d in disease_scores: disease_scores[d] += 1
                else:
                    disease_scores[disease] += 2

    # ── Compound name / class keyword scoring ──
    name_lower = compound_name.lower()
    kw_map = {
        'alzheimers': ['alzheimer','donepezil','memantine','rivastigmine','galantamine',
                       'aducanumab','lecanemab','amyloid','bace','ache inhibit',
                       'quercetin','fisetin','curcumin','resveratrol','egcg',
                       'berberine','apigenin','luteolin','kaempferol','naringenin'],
        'parkinsons': ['parkinson','levodopa','carbidopa','selegiline','rasagiline',
                       'dopamine','maob','synuclein','lrrk2','coenzyme q','coq10',
                       'pqq','ubiquinol','quercetin','fisetin','curcumin','resveratrol'],
        'als':        ['als','riluzole','edaravone','sod1','tdp-43','fus','motor neuron',
                       'quercetin','curcumin','resveratrol','bacopa','lion mane'],
        'huntingtons':['huntington','tetrabenazine','deutetrabenazine','htt','striatum',
                       'spermidine','rapamycin','fisetin','cystamine','coenzyme q'],
    }
    for disease, keywords in kw_map.items():
        for kw in keywords:
            if kw in name_lower: disease_scores[disease] += 2

    # ── Compound CLASS broad rules ──
    compound_classes = {
        'flavonoid':     {'alzheimers':3,'parkinsons':3,'als':2,'huntingtons':2},
        'polyphenol':    {'alzheimers':3,'parkinsons':3,'als':2,'huntingtons':2},
        'stilbenoid':    {'alzheimers':3,'parkinsons':2,'als':2,'huntingtons':2},
        'catechin':      {'alzheimers':3,'parkinsons':3,'als':2,'huntingtons':1},
        'curcumin':      {'alzheimers':3,'parkinsons':3,'als':3,'huntingtons':2},
        'anthocyanin':   {'alzheimers':3,'parkinsons':2,'als':1,'huntingtons':1},
        'terpene':       {'alzheimers':2,'parkinsons':2,'als':1,'huntingtons':1},
        'alkaloid':      {'alzheimers':2,'parkinsons':3,'als':1,'huntingtons':1},
        'polyamine':     {'alzheimers':3,'parkinsons':2,'als':2,'huntingtons':2},
        'coenzyme':      {'alzheimers':2,'parkinsons':3,'als':2,'huntingtons':3},
        'omega':         {'alzheimers':2,'parkinsons':2,'als':2,'huntingtons':1},
        'urolithin':     {'alzheimers':2,'parkinsons':2,'als':2,'huntingtons':2},
        'spermidine':    {'alzheimers':3,'parkinsons':2,'als':2,'huntingtons':2},
        'bacopa':        {'alzheimers':3,'parkinsons':2,'als':2,'huntingtons':1},
        'ashwagandha':   {'alzheimers':3,'parkinsons':2,'als':2,'huntingtons':2},
        'lion':          {'alzheimers':3,'parkinsons':2,'als':2,'huntingtons':1},
        'saffron':       {'alzheimers':3,'parkinsons':2,'als':1,'huntingtons':1},
        'gotu':          {'alzheimers':3,'parkinsons':2,'als':1,'huntingtons':1},
    }
    for cls_kw, scores in compound_classes.items():
        if cls_kw in name_lower:
            for d, s in scores.items():
                disease_scores[d] += s

    # ── Convert scores to High/Med/Low ──
    result = {}
    for d, score in disease_scores.items():
        if score >= 4:   result[d] = 'High'
        elif score >= 2: result[d] = 'Med'
        else:            result[d] = 'Low'
    return result

def predict_pathways(compound_name, chembl_targets=[]):
    """Predict relevant NDD pathways from compound name + target keywords."""
    name_lower = compound_name.lower()
    target_text = ' '.join(chembl_targets).lower()
    combined = name_lower + ' ' + target_text

    matched = []
    for pathway, keywords in NDD_PATHWAYS.items():
        for kw in keywords:
            if kw in combined:
                matched.append(pathway)
                break

    # Broad-spectrum defaults for well-known compound classes
    broad_kw = {
        'Nrf2/GSH':     ['polyphenol','flavon','curcumin','egcg','resveratrol','quercetin'],
        'NF-kB':        ['anti-inflam','ibuprofen','aspirin','curcumin','flavon'],
        'PI3K/Akt/mTOR':['rapamycin','spermidine','berberine','metformin','resveratrol'],
        'BDNF/TrkB':    ['bacopa','lion mane','ashwagandha','ngf','neurotrophin'],
        'Autophagy':    ['spermidine','rapamycin','metformin','resveratrol','urolithin'],
    }
    for pathway, keywords in broad_kw.items():
        if pathway not in matched:
            for kw in keywords:
                if kw in name_lower and pathway not in matched:
                    matched.append(pathway)

    # Class-based pathway defaults for common compound types
    class_pathways = {
        'quercetin':   ['Nrf2/GSH','NF-kB','PI3K/Akt/mTOR','BDNF/TrkB','Tau/Aggregation'],
        'flavon':      ['Nrf2/GSH','NF-kB','PI3K/Akt/mTOR'],
        'curcumin':    ['Nrf2/GSH','NF-kB','PI3K/Akt/mTOR','Autophagy','Tau/Aggregation'],
        'resveratrol': ['Nrf2/GSH','Sirtuin/AMPK','Autophagy','PI3K/Akt/mTOR'],
        'egcg':        ['Nrf2/GSH','NF-kB','PI3K/Akt/mTOR','Tau/Aggregation','Alpha-syn'],
        'spermidine':  ['Autophagy','PI3K/Akt/mTOR','Sirtuin/AMPK'],
        'urolithin':   ['Mitophagy/PINK1','Autophagy','Sirtuin/AMPK','GSH/NAD+/ATP'],
        'bacopa':      ['BDNF/TrkB','Sirtuin/AMPK','Nrf2/GSH'],
        'lion':        ['BDNF/TrkB','Neurogenesis','Nrf2/GSH'],
        'coenzyme q':  ['GSH/NAD+/ATP','Mitophagy/PINK1'],
        'pqq':         ['GSH/NAD+/ATP','Mitophagy/PINK1','BDNF/TrkB'],
        'omega':       ['NF-kB','BDNF/TrkB','Nrf2/GSH'],
        'berberine':   ['Autophagy','PI3K/Akt/mTOR','NF-kB','AMPK'],
        'fisetin':     ['Nrf2/GSH','NF-kB','PI3K/Akt/mTOR','Tau/Aggregation'],
        'apigenin':    ['Nrf2/GSH','NF-kB','BDNF/TrkB','Caspase/Apoptosis'],
    }
    for cls, paths in class_pathways.items():
        if cls in name_lower:
            for p in paths:
                if p not in matched:
                    matched.append(p)

    return matched if matched else ['NF-kB', 'Nrf2/GSH']  # minimum default

# ── Train v3 Model ────────────────────────────────────────────────────────────
def train_v3_model(compounds_json_path, output_dir='models_v3'):
    """Train RF+Ridge ensemble with Morgan fingerprints."""
    os.makedirs(output_dir, exist_ok=True)

    with open(compounds_json_path) as f:
        data = json.load(f)
    items = list(data.items()) if isinstance(data, dict) else list(enumerate(data))

    NPS_DIMS = ['antioxidant','anti_inflammatory','mitochondrial_support',
                'aggregation_modulation','cognitive_enhancement','neurogenesis','synaptic_plasticity']

    X, y_nps = [], []
    print("Fetching SMILES for training compounds...")
    for i, (name, item) in enumerate(items):
        compound_name = item.get('name', str(name))
        smiles, props = get_smiles_pubchem(compound_name)
        if smiles:
            feats = smiles_to_features(smiles)
            if feats is not None:
                dim_vals = [float(item.get(d, 0) or 0) for d in NPS_DIMS]
                nps = np.mean(dim_vals) * 10
                X.append(feats)
                y_nps.append(nps)
                if (i+1) % 10 == 0:
                    print(f"  Processed {i+1}/{len(items)} — {compound_name} ✓")
        time.sleep(0.3)  # PubChem rate limit

    X = np.array(X)
    y = np.array(y_nps)
    print(f"\nTraining set: {len(X)} compounds with SMILES features")
    print(f"Feature dims: {X.shape[1]} (2048 Morgan + 8 physicochemical)")

    # Train RF regressor
    rf = RandomForestRegressor(n_estimators=200, max_depth=8,
                                min_samples_leaf=3, random_state=42, n_jobs=-1)
    rf_cv = cross_val_score(rf, X, y, cv=5, scoring='r2')
    print(f"RF v3 CV R²: {rf_cv.mean():.3f} ± {rf_cv.std():.3f}")
    rf.fit(X, y)

    # Train Ridge
    ridge = Pipeline([('scaler', StandardScaler()), ('ridge', Ridge(alpha=10))])
    ridge_cv = cross_val_score(ridge, X, y, cv=5, scoring='r2')
    print(f"Ridge v3 CV R²: {ridge_cv.mean():.3f} ± {ridge_cv.std():.3f}")
    ridge.fit(X, y)

    joblib.dump(rf, f'{output_dir}/rf_v3.pkl')
    joblib.dump(ridge, f'{output_dir}/ridge_v3.pkl')
    print(f"\n✅ Models saved to {output_dir}/")
    return rf, ridge, X.shape[1]

# ── Main Prediction Function ──────────────────────────────────────────────────
def predict_new_compound(compound_name, models_dir='models_v3'):
    """Full prediction pipeline for any new compound."""
    print(f"\n{'='*55}")
    print(f"  Predicting: {compound_name}")
    print(f"{'='*55}")

    # 1. Fetch SMILES + properties
    print("1. Fetching SMILES from PubChem...")
    smiles, props = get_smiles_pubchem(compound_name)
    if not smiles:
        return {"error": f"SMILES not found for {compound_name}"}
    print(f"   SMILES: {smiles[:60]}...")
    print(f"   MW={props.get('MolecularWeight','?')}, LogP={props.get('XLogP','?')}, TPSA={props.get('TPSA','?')}")

    # 2. Compute features
    feats = smiles_to_features(smiles)
    if feats is None:
        return {"error": "RDKit feature extraction failed"}

    # 3. Predict NPS
    nps_pred = None
    mol = Chem.MolFromSmiles(smiles)
    phys = np.array([
        float(props.get('MolecularWeight', 400)),
        float(props.get('XLogP', 2)),
        float(props.get('TPSA', 80)),
        rdMolDescriptors.CalcNumHBD(mol) if mol else 2,
        rdMolDescriptors.CalcNumHBA(mol) if mol else 4,
        rdMolDescriptors.CalcNumRotatableBonds(mol) if mol else 4,
        rdMolDescriptors.CalcNumAromaticRings(mol) if mol else 1,
        Descriptors.NumHeteroatoms(mol) if mol else 3,
    ]).reshape(1,-1)

    # Track A: NPS from physicochemical features
    if os.path.exists(f'{models_dir}/rf_nps.pkl'):
        rf_nps    = joblib.load(f'{models_dir}/rf_nps.pkl')
        ridge_nps = joblib.load(f'{models_dir}/ridge_nps.pkl')
        rf_pred   = rf_nps.predict(phys)[0]
        rdg_pred  = ridge_nps.predict(phys)[0]
        nps_pred  = round(max(0, min(100, (rf_pred + rdg_pred) / 2)), 1)
        print(f"2. NPS Prediction: {nps_pred}/100  (RF={rf_pred:.1f}, Ridge={rdg_pred:.1f})")

    # Track C: Similarity search
    if os.path.exists(f'{models_dir}/similarity_index.pkl') and feats is not None:
        idx = joblib.load(f'{models_dir}/similarity_index.pkl')
        query_fp = feats[:1024]
        sims = [np.dot(query_fp, f) / (np.linalg.norm(query_fp)*np.linalg.norm(f)+1e-9)
                for f in idx['fps']]
        top_i  = np.argmax(sims)
        print(f"   Nearest compound: {idx['names'][top_i]} (similarity={sims[top_i]:.2f})")
    else:
        # v3 not trained yet — use rule-based NPS estimate from physicochemical + disease
        mw_tmp   = float(props.get('MolecularWeight', 500))
        logp_tmp = float(props.get('XLogP', 3))
        tpsa_tmp = float(props.get('TPSA', 90))
        # Heuristic: good drug-like + multi-disease → higher NPS
        base = 50.0
        if mw_tmp < 400: base += 5
        if 1 <= logp_tmp <= 4: base += 5
        if tpsa_tmp < 90: base += 5
        nps_pred = min(round(base, 1), 75.0)
        print(f"2. NPS Estimate (rule-based, pre-training): {nps_pred}/100")
        print("   ⚠️  Train v3 model for accurate NPS: python3 ml_v3_engine.py train")

    # 4. BBB from physicochemical rules
    mw   = float(props.get('MolecularWeight', 500))
    logp = float(props.get('XLogP', 3))
    tpsa = float(props.get('TPSA', 80))
    if mw < 400 and 1 <= logp <= 4 and tpsa < 75:
        bbb = 'High'
    elif mw < 500 and logp <= 5 and tpsa < 90:
        bbb = 'Medium'
    elif mw < 600 and tpsa < 120:
        bbb = 'Low-Med'
    else:
        bbb = 'Low'
    print(f"3. BBB Class: {bbb}")

    # 5. ChEMBL targets
    print("4. Fetching ChEMBL targets...")
    targets = get_chembl_targets(compound_name)
    print(f"   Found {len(targets)} target annotations")

    # 6. Disease relevance
    diseases = predict_disease_relevance(compound_name, smiles, targets)
    print(f"5. Disease Relevance:")
    for d, rel in diseases.items():
        print(f"   {d:15s}: {rel}")

    # 7. Pathway annotation
    pathways = predict_pathways(compound_name, targets)
    print(f"6. Predicted Pathways: {', '.join(pathways)}")

    result = {
        'compound':  compound_name,
        'smiles':    smiles,
        'nps':       nps_pred,
        'bbb':       bbb,
        'mw':        mw,
        'logp':      logp,
        'tpsa':      tpsa,
        'diseases':  diseases,
        'pathways':  pathways,
        'targets_n': len(targets),
    }
    return result

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'train':
        train_v3_model('compounds.json')
    elif len(sys.argv) > 1:
        name = ' '.join(sys.argv[1:])
        result = predict_new_compound(name)
        print(json.dumps(result, indent=2))
    else:
        # Quick test without trained model
        result = predict_new_compound('Quercetin')
        print(json.dumps(result, indent=2))
