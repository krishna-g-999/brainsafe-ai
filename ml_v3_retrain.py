"""
BrainSafe AI v3 — Fixed Architecture
Track A: NPS regression    → 8 physicochemical features only (proven CV R²=0.20)
Track B: BBB classification → Morgan FP + PCA(50) (fingerprints work for classification)
Track C: Similarity search  → Morgan FP Tanimoto (best use of fingerprints at n=116)
"""
import json, numpy as np, joblib, os, time, requests, urllib.parse
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.linear_model import Ridge
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import cross_val_score
os.makedirs('models_v3', exist_ok=True)

NPS_DIMS = ['antioxidant','anti_inflammatory','mitochondrial_support',
            'aggregation_modulation','cognitive_enhancement','neurogenesis','synaptic_plasticity']

def get_smiles(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(name)}/property/IsomericSMILES,MolecularWeight,XLogP,TPSA/JSON"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            props = r.json()['PropertyTable']['Properties'][0]
            smiles = props.get('IsomericSMILES') or props.get('SMILES','')
            return smiles, props
    except: pass
    return None, {}

def physchem_features(props):
    """8 physicochemical features — proven for NPS regression."""
    return np.array([
        float(props.get('MolecularWeight', 400)),
        float(props.get('XLogP', 2)),
        float(props.get('TPSA', 80)),
        0, 0, 0, 0, 0  # HBD, HBA, RotBonds, ArRings, Heteroatoms from RDKit
    ])

def physchem_from_mol(mol, props):
    """Full 8 features including RDKit-computed ones."""
    try:
        return np.array([
            float(props.get('MolecularWeight', Descriptors.MolWt(mol))),
            float(props.get('XLogP', Descriptors.MolLogP(mol))),
            float(props.get('TPSA', Descriptors.TPSA(mol))),
            rdMolDescriptors.CalcNumHBD(mol),
            rdMolDescriptors.CalcNumHBA(mol),
            rdMolDescriptors.CalcNumRotatableBonds(mol),
            rdMolDescriptors.CalcNumAromaticRings(mol),
            Descriptors.NumHeteroatoms(mol),
        ])
    except:
        return physchem_features(props)

def morgan_fp(smiles, nbits=1024):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            gen = AllChem.GetMorganGenerator(radius=2, fpSize=nbits)
            return np.array(gen.GetFingerprint(mol)), mol
    except: pass
    return None, None

# ── Load training data ────────────────────────────────────────────────────────
print("Loading compounds.json...")
with open('compounds.json') as f:
    data = json.load(f)
items = list(data.items()) if isinstance(data, dict) else list(enumerate(data))

X_phys, X_fp, y_nps, y_bbb, fps_store, names_store = [], [], [], [], [], []
BBB_ENCODE = {'High':3,'Medium':2,'Low-Med':1,'Low':0,'Moderate':2}

print(f"Fetching SMILES for {len(items)} compounds...")
for i, (name, item) in enumerate(items):
    cname = item.get('name', str(name))
    smiles, props = get_smiles(cname)
    if not smiles: continue

    fp_arr, mol = morgan_fp(smiles)
    if mol is None: continue

    phys = physchem_from_mol(mol, props)
    dim_vals = [float(item.get(d, 0) or 0) for d in NPS_DIMS]
    nps = np.mean(dim_vals) * 10
    bbb = item.get('bbb','Low')
    if bbb == 'Moderate': bbb = 'Medium'

    X_phys.append(phys)
    X_fp.append(fp_arr)
    y_nps.append(nps)
    y_bbb.append(BBB_ENCODE.get(bbb, 0))
    fps_store.append(fp_arr)
    names_store.append(cname)
    if (i+1) % 20 == 0:
        print(f"  {i+1}/{len(items)} — {cname} ✓")
    time.sleep(0.3)

X_phys = np.array(X_phys)
X_fp   = np.array(X_fp)
y_nps  = np.array(y_nps)
y_bbb  = np.array(y_bbb)
print(f"\nTraining set: {len(X_phys)} compounds")

# ── Track A: NPS Regression (physicochemical only) ───────────────────────────
print("\n── Track A: NPS Regression (8 physicochemical features) ──")
rf_nps = RandomForestRegressor(n_estimators=300, max_depth=6,
                                min_samples_leaf=4, random_state=42)
ridge_nps = Pipeline([('scaler', StandardScaler()), ('ridge', Ridge(alpha=10))])

rf_cv  = cross_val_score(rf_nps,    X_phys, y_nps, cv=5, scoring='r2')
rdg_cv = cross_val_score(ridge_nps, X_phys, y_nps, cv=5, scoring='r2')
print(f"RF    CV R²: {rf_cv.mean():.3f} ± {rf_cv.std():.3f}")
print(f"Ridge CV R²: {rdg_cv.mean():.3f} ± {rdg_cv.std():.3f}")

rf_nps.fit(X_phys, y_nps)
ridge_nps.fit(X_phys, y_nps)
joblib.dump(rf_nps,    'models_v3/rf_nps.pkl')
joblib.dump(ridge_nps, 'models_v3/ridge_nps.pkl')
print("✅ NPS models saved")

# ── Track B: BBB Classification (Morgan FP + PCA 50) ────────────────────────
print("\n── Track B: BBB Classification (Morgan FP + PCA 50) ──")
bbb_pipe = Pipeline([
    ('pca',   PCA(n_components=50)),
    ('scaler',StandardScaler()),
    ('rf',    RandomForestClassifier(n_estimators=200, max_depth=8,
                                      min_samples_leaf=2, random_state=42))
])
bbb_cv = cross_val_score(bbb_pipe, X_fp, y_bbb, cv=5, scoring='f1_weighted')
print(f"BBB Classifier CV F1: {bbb_cv.mean():.3f} ± {bbb_cv.std():.3f}")
bbb_pipe.fit(X_fp, y_bbb)
joblib.dump(bbb_pipe, 'models_v3/bbb_classifier.pkl')
joblib.dump({'0':'Low','1':'Low-Med','2':'Medium','3':'High'}, 'models_v3/bbb_labels.pkl')
print("✅ BBB classifier saved")

# ── Track C: Similarity index ────────────────────────────────────────────────
print("\n── Track C: Building similarity index ──")
joblib.dump({'fps': fps_store, 'names': names_store}, 'models_v3/similarity_index.pkl')
print(f"✅ Similarity index saved ({len(fps_store)} compounds)")

print("\n" + "="*50)
print("  v3 TRAINING COMPLETE")
print("="*50)
print(f"  NPS RF CV R²  : {rf_cv.mean():.3f}")
print(f"  NPS Ridge CV R²: {rdg_cv.mean():.3f}")
print(f"  BBB CV F1     : {bbb_cv.mean():.3f}")
print(f"  Similarity idx: {len(fps_store)} compounds")
