"""
BS_fetch_antioxidant.py — assemble a MEASURED antioxidant dataset to replace the weak
curated labels. Pulls DPPH radical-scavenging assays from ChEMBL, keeps IC50/EC50 with
concentration units, converts to pIC50 = -log10(M), aggregates per compound (median).
Out: data/endpoints_reg/antioxidant_dpph.csv (smiles, y, year)
"""
import os, json, time, math
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import requests, pandas as pd
B = "https://www.ebi.ac.uk/chembl/api/data"
os.makedirs("data/endpoints_reg", exist_ok=True)
CACHE = "data/_chembl_cache/antioxidant_dpph_raw.json"
MAX_ASSAYS = 1600
UNIT_TO_M = {"nM": 1e-9, "uM": 1e-6, "µM": 1e-6, "mM": 1e-3, "M": 1.0, "pM": 1e-12}


def assay_ids(term="DPPH"):
    ids, offset = [], 0
    while len(ids) < MAX_ASSAYS:
        u = f"{B}/assay.json?description__icontains={term}&limit=1000&offset={offset}"
        j = requests.get(u, timeout=60).json()
        for a in j.get("assays", []):
            ids.append(a["assay_chembl_id"])
        if not j.get("page_meta", {}).get("next"):
            break
        offset += 1000
    return ids[:MAX_ASSAYS]


def fetch_activities(ids):
    rows = []
    for i in range(0, len(ids), 25):
        batch = ",".join(ids[i:i + 25]); offset = 0
        while True:
            u = (f"{B}/activity.json?assay_chembl_id__in={batch}"
                 f"&standard_type__in=IC50,EC50&limit=1000&offset={offset}")
            try:
                j = requests.get(u, timeout=60).json()
            except Exception:
                break
            for a in j.get("activities", []):
                smi, val, unit = a.get("canonical_smiles"), a.get("standard_value"), a.get("standard_units")
                if smi and val and unit in UNIT_TO_M:
                    try:
                        molar = float(val) * UNIT_TO_M[unit]
                        if molar > 0:
                            rows.append({"smiles": smi, "pIC50": -math.log10(molar),
                                         "year": a.get("document_year")})
                    except Exception:
                        pass
            if not j.get("page_meta", {}).get("next"):
                break
            offset += 1000
        time.sleep(0.15)
    return rows


def main():
    if os.path.exists(CACHE):
        rows = json.load(open(CACHE)); print("cache:", len(rows))
    else:
        print("collecting DPPH assay ids..."); ids = assay_ids("DPPH")
        print(f"  {len(ids)} DPPH assays; fetching IC50/EC50 activities...")
        rows = fetch_activities(ids)
        json.dump(rows, open(CACHE, "w"))
        print("  raw activities:", len(rows))
    df = pd.DataFrame(rows)
    df = df[(df.pIC50 > 2) & (df.pIC50 < 10)]               # sane radical-scavenging range
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    g = df.groupby("smiles").agg(y=("pIC50", "median"), year=("year", "min")).reset_index()
    g.to_csv("data/endpoints_reg/antioxidant_dpph.csv", index=False)
    print(f"MEASURED antioxidant (DPPH) dataset: {len(g)} unique compounds, "
          f"pIC50 {g.y.min():.1f}-{g.y.max():.1f} (median {g.y.median():.1f}, std {g.y.std():.2f})")


if __name__ == "__main__":
    main()
