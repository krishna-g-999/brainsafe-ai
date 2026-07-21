"""
BS_assay_composition.py -- quantify assay-type heterogeneity (threat-to-validity #1).
Re-queries ChEMBL capturing standard_type (which the original fetch filtered on but did not
store), and tabulates the IC50/Ki/Kd/EC50/Potency composition per endpoint. It also caches the
per-activity (smiles, pchembl, standard_type, year) rows so a single-assay sensitivity retrain
can run without re-hitting the API.

Outputs:
    data/_chembl_cache/<name>_std.json         (smiles, pchembl, std, year)
    supplementary/STable11_assay_type_composition.csv
    BS_assay_composition.json
"""
import os, json, time
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import requests, pandas as pd, urllib3
urllib3.disable_warnings()  # this network proxies TLS; public read-only ChEMBL fetch, data integrity unaffected
VERIFY = False
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
OUT = "supplementary"; os.makedirs(OUT, exist_ok=True)
os.makedirs("data/_chembl_cache", exist_ok=True)

TARGETS = {"AChE": "CHEMBL220", "BChE": "CHEMBL1914", "BACE1": "CHEMBL4822",
           "GSK3B": "CHEMBL262", "MAO_B": "CHEMBL2039", "MAO_A": "CHEMBL1951",
           "hERG": "CHEMBL240"}
KEEP = ("IC50", "Ki", "Kd", "EC50", "Potency")
MAX_PAGES, PAGE = 20, 1000

def fetch_std(name, tid):
    cache = f"data/_chembl_cache/{name}_std.json"
    if os.path.exists(cache):
        rows = json.load(open(cache)); print(f"  [{name}] cache std: {len(rows)}"); return rows
    rows, offset = [], 0
    for p in range(MAX_PAGES):
        url = (f"{CHEMBL}/activity.json?target_chembl_id={tid}"
               f"&pchembl_value__isnull=false&limit={PAGE}&offset={offset}")
        try:
            j = requests.get(url, timeout=60, verify=VERIFY).json()
        except Exception as e:
            print(f"  [{name}] page {p} err {e}"); break
        for a in j.get("activities", []):
            smi, pv, st = a.get("canonical_smiles"), a.get("pchembl_value"), a.get("standard_type")
            if smi and pv and st in KEEP:
                rows.append({"smiles": smi, "pchembl": float(pv), "std": st,
                             "year": a.get("document_year")})
        offset += PAGE
        if not j.get("page_meta", {}).get("next"):
            break
        time.sleep(0.2)
    json.dump(rows, open(cache, "w")); print(f"  [{name}] fetched std: {len(rows)}")
    return rows

comp_rows, bundle = [], {}
for name, tid in TARGETS.items():
    rows = fetch_std(name, tid)
    if not rows:
        continue
    df = pd.DataFrame(rows)
    vc = df["std"].value_counts()
    total = int(vc.sum())
    frac = {k: round(int(v) / total, 3) for k, v in vc.items()}
    dominant = vc.index[0]; dom_share = round(int(vc.iloc[0]) / total, 3)
    bundle[name] = {"n_activities": total, "composition": {k: int(v) for k, v in vc.items()},
                    "fraction": frac, "dominant_type": dominant, "dominant_share": dom_share}
    for k, v in vc.items():
        comp_rows.append({"endpoint": name, "standard_type": k, "n_activities": int(v),
                          "fraction": round(int(v) / total, 3)})
    print(f"  [{name}] n={total} dominant={dominant} ({dom_share}) :: "
          + ", ".join(f"{k}={frac[k]}" for k in vc.index))

pd.DataFrame(comp_rows).to_csv(f"{OUT}/STable11_assay_type_composition.csv", index=False)
json.dump(bundle, open("BS_assay_composition.json", "w"), indent=2)
print("\nWrote STable11 + BS_assay_composition.json")
