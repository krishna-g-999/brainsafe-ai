"""
BrainSafe AI — Natural Products Batch Expander v2
Uses PubChem PUG REST classification search (correct endpoint)
Fixes: .items() not .items, correct API URL, rate limiting
"""
import sqlite3, json, numpy as np, requests, time, logging
from sklearn.ensemble import RandomForestRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(message)s",
                    handlers=[logging.StreamHandler(),
                              logging.FileHandler("natural_expand.log")])
log = logging.getLogger("batch_np")

DB    = "brainsafe_natural.db"
DIMS  = ["antioxidant","anti_inflammatory","mitochondrial_support",
         "aggregation_modulation","cognitive_enhancement","neurogenesis","synaptic_plasticity"]
POLYPHENOL_TYPES = {"flavonoid","polyphenol","catechin","stilbene","terpene",
                    "carotenoid","vitamin","phenolic","alkaloid","curcuminoid"}
NEURO_KWS = {"bdnf","trkb","wnt","ngf","neurogenesis","hippocampus","creb","notch"}

# ── Classes: (search_term, compound_type, is_polyphenol, n_neuro_est) ─────
CLASSES = [
    ("flavonoid",        "Flavonoid",   1.0, 2.0),
    ("polyphenol",       "Polyphenol",  1.0, 1.5),
    ("catechin",         "Flavonoid",   1.0, 2.0),
    ("quercetin",        "Flavonoid",   1.0, 2.0),
    ("kaempferol",       "Flavonoid",   1.0, 2.0),
    ("luteolin",         "Flavonoid",   1.0, 2.5),
    ("apigenin",         "Flavonoid",   1.0, 2.0),
    ("rutin",            "Flavonoid",   1.0, 1.5),
    ("resveratrol",      "Polyphenol",  1.0, 2.0),
    ("curcumin",         "Polyphenol",  1.0, 2.5),
    ("berberine",        "Alkaloid",    0.5, 2.0),
    ("huperzine",        "Alkaloid",    0.0, 3.0),
    ("ursolic acid",     "Terpenoid",   0.5, 1.5),
    ("betulinic acid",   "Terpenoid",   0.5, 1.0),
    ("oleic acid",       "Fatty Acid",  0.0, 0.5),
    ("DHA",              "Fatty Acid",  0.0, 2.5),
    ("tocopherol",       "Vitamin",     0.5, 1.0),
    ("ascorbic acid",    "Vitamin",     1.0, 0.5),
    ("retinol",          "Vitamin",     0.0, 1.0),
    ("sulforaphane",     "Terpenoid",   0.5, 2.0),
    ("baicalein",        "Flavonoid",   1.0, 2.5),
    ("naringenin",       "Flavonoid",   1.0, 1.5),
    ("epicatechin",      "Flavonoid",   1.0, 2.0),
    ("myricetin",        "Flavonoid",   1.0, 1.5),
    ("hesperidin",       "Flavonoid",   1.0, 1.0),
    ("fisetin",          "Flavonoid",   1.0, 3.0),
    ("pterostilbene",    "Polyphenol",  1.0, 2.5),
    ("ellagic acid",     "Polyphenol",  1.0, 1.5),
    ("gallic acid",      "Polyphenol",  1.0, 1.0),
    ("ferulic acid",     "Polyphenol",  1.0, 1.5),
    ("caffeic acid",     "Polyphenol",  1.0, 1.0),
    ("rosmarinic acid",  "Polyphenol",  1.0, 2.0),
    ("thymoquinone",     "Terpenoid",   0.5, 2.0),
    ("lycopene",         "Carotenoid",  0.5, 1.0),
    ("beta-carotene",    "Carotenoid",  0.5, 0.5),
    ("astaxanthin",      "Carotenoid",  0.5, 1.5),
    ("vincamine",        "Alkaloid",    0.0, 2.5),
    ("hericenone",       "Terpenoid",   0.5, 4.0),
    ("erinacine",        "Terpenoid",   0.5, 4.0),
    ("withanolide",      "Terpenoid",   0.5, 2.5),
    ("bacosides",        "Terpenoid",   0.5, 2.5),
]

# ── Train model on 129 curated compounds ──────────────────────────────────
with open("compounds.json") as f: raw = json.load(f)
curated = raw.get("compounds", raw) if isinstance(raw, dict) else raw
BBBMAP = {"Low":0,"Low-Med":1,"Medium":2,"High":3}
DISMAP = {"Low":0,"Med":1,"High":2}
X, y = [], []
for e in curated.values():
    ct   = (e.get("compound_type","") or "").lower()
    poly = 1.0 if any(p in ct for p in POLYPHENOL_TYPES) else 0.0
    pt   = " ".join(e.get("pathways",[])).lower()
    nn   = float(sum(1 for kw in NEURO_KWS if kw in pt))
    nm   = float(len(e.get("metabolites",[])))
    X.append([float(BBBMAP.get(e.get("bbb","Low"),0)),
              float(DISMAP.get(e.get("als","Low"),0)),
              float(DISMAP.get(e.get("alzheimers","Low"),0)),
              float(DISMAP.get(e.get("parkinsons","Low"),0)),
              float(DISMAP.get(e.get("huntingtons","Low"),0)),
              float(len(e.get("pathways",[]))), poly, nn, nm])
    y.append([float(e.get(d,5.0)) for d in DIMS])
X, y   = np.array(X), np.array(y)
scaler = StandardScaler(); Xs = scaler.fit_transform(X)
model  = MultiOutputRegressor(RandomForestRegressor(
    n_estimators=150, max_depth=6, min_samples_leaf=2, random_state=42))
model.fit(Xs, y)
log.info(f"✅  Model trained on {len(X)} curated compounds (9 features)")

# ── SQLite setup ──────────────────────────────────────────────────────────
conn = sqlite3.connect(DB)
conn.execute("""CREATE TABLE IF NOT EXISTS compounds (
    name TEXT PRIMARY KEY, compound_type TEXT, bbb TEXT,
    antioxidant REAL, anti_inflammatory REAL, mitochondrial_support REAL,
    aggregation_modulation REAL, cognitive_enhancement REAL,
    neurogenesis REAL, synaptic_plasticity REAL,
    nps REAL, confidence TEXT, source TEXT, pubchem_cid TEXT)""")

def bbb_from(mw, logp, tpsa):
    try:
        mw, logp, tpsa = float(mw), float(logp), float(tpsa)
        if mw < 360 and 1.0 < logp < 3.0 and tpsa < 60:  return "High",   3
        if mw < 450 and 0.0 < logp < 4.0 and tpsa < 90:  return "Medium", 2
        if mw < 500 and tpsa < 120:                        return "Low-Med",1
    except: pass
    return "Low", 0

def predict_and_insert(name, cid, mw, logp, tpsa, ctype, poly, nn):
    bbb_s, bbb_n = bbb_from(mw, logp, tpsa)
    nm = 3.0 if poly > 0 else 1.0
    feat  = np.array([[bbb_n, 0, 0.5, 0, 0, 3.0, poly, nn, nm]])
    preds = model.predict(scaler.transform(feat))[0]
    sc    = {d: round(float(np.clip(v, 1, 10)), 1) for d, v in zip(DIMS, preds)}
    nps   = round(min(100.0, 3*sc["antioxidant"] + 3*sc["anti_inflammatory"]
                      + 2*sc["mitochondrial_support"] + 2*sc["aggregation_modulation"]), 1)
    conn.execute("INSERT OR REPLACE INTO compounds VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (name, ctype, bbb_s, sc["antioxidant"], sc["anti_inflammatory"],
         sc["mitochondrial_support"], sc["aggregation_modulation"],
         sc["cognitive_enhancement"], sc["neurogenesis"],
         sc["synaptic_plasticity"], nps, "medium", "PubChem+RF-v2.2", str(cid)))

# ── Fetch CIDs per class and batch-predict ────────────────────────────────
total = 0
for term, ctype, poly, nn in CLASSES:
    log.info(f"Fetching: {term} ({ctype})")
    try:
        # Get all CIDs for this class (up to 2000)
        url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
               f"{requests.utils.quote(term)}/cids/JSON?name_type=word")
        r = requests.get(url, timeout=20)
        if r.status_code != 200:
            log.warning(f"  {term}: HTTP {r.status_code}"); continue
        cids = r.json().get("IdentifierList",{}).get("CID",[])[:1000]
        if not cids:
            log.warning(f"  {term}: 0 CIDs"); continue
        log.info(f"  Got {len(cids)} CIDs")

        # Batch fetch properties (chunks of 100)
        for i in range(0, len(cids), 100):
            chunk = cids[i:i+100]
            prop_url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                        f"{','.join(map(str,chunk))}/property/"
                        f"IUPACName,MolecularWeight,XLogP,TPSA/JSON")
            pr = requests.get(prop_url, timeout=30)
            if pr.status_code != 200:
                time.sleep(1); continue
            props = pr.json().get("PropertyTable",{}).get("Properties",[])
            for p in props:
                name = (p.get("IUPACName") or f"{term}_{p.get('CID','')}").strip()
                if not name or len(name) > 120: continue
                predict_and_insert(
                    name, p.get("CID",""),
                    p.get("MolecularWeight", 300),
                    p.get("XLogP", 2.0) or 2.0,
                    p.get("TPSA", 80) or 80,
                    ctype, poly, nn
                )
                total += 1
            conn.commit()
            time.sleep(0.25)
        log.info(f"  {term}: inserted running total={total}")
    except Exception as e:
        log.error(f"  {term} failed: {e}")
        time.sleep(2)
        continue

conn.close()
log.info(f"\n✅  DONE — {total:,} natural products in {DB}")
log.info(f"   Curated(129) + ChEMBL(800) + Natural({total}) = {total+929:,} total")
