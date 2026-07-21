"""
BS_fetch_endpoints.py — fetch REAL measured bioactivity for brain-relevant endpoints.
Now captures document_year (for TEMPORAL validation) and an expanded CNS target panel.

Labelling: active = pChEMBL>=6 (<=1uM), inactive = pChEMBL<5 (>10uM); grey zone dropped;
per-compound median pChEMBL; per-compound EARLIEST document_year kept for time-splits.
Sources: ChEMBL REST API (targets) + B3DB (BBB).
Outputs: data/endpoints/<endpoint>.csv (smiles,label,pchembl,year)
"""
import os, json, time, io
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import requests, pandas as pd

OUT = "data/endpoints"; os.makedirs(OUT, exist_ok=True)
os.makedirs("data/_chembl_cache", exist_ok=True)
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

TARGETS = {
    # core CNS targets
    "AChE":  ("CHEMBL220",  "Alzheimer's / cognition"),
    "BChE":  ("CHEMBL1914", "Alzheimer's / cholinergic"),
    "BACE1": ("CHEMBL4822", "Alzheimer's / amyloid"),
    "GSK3B": ("CHEMBL262",  "tau / neuroprotection"),
    "MAO_B": ("CHEMBL2039", "Parkinson's / dopamine"),
    "MAO_A": ("CHEMBL1951", "mood / depression"),
    # expanded panel
    "D2":    ("CHEMBL217",  "Parkinson's / psychosis (dopamine D2)"),
    "A2A":   ("CHEMBL251",  "Parkinson's (adenosine A2A)"),
    "HT2A":  ("CHEMBL224",  "mood / psychosis (5-HT2A)"),
    "SERT":  ("CHEMBL228",  "depression (serotonin transporter)"),
    # safety anti-target
    "hERG":  ("CHEMBL240",  "SAFETY: cardiotoxicity liability (hERG block)"),
}
MAX_PAGES, PAGE = 16, 1000


def fetch_target_activities(name, tid):
    cache = f"data/_chembl_cache/{name}_y.json"
    if os.path.exists(cache):
        rows = json.load(open(cache)); print(f"  [{name}] cache: {len(rows)}"); return rows
    rows, offset = [], 0
    for p in range(MAX_PAGES):
        url = (f"{CHEMBL}/activity.json?target_chembl_id={tid}"
               f"&pchembl_value__isnull=false&limit={PAGE}&offset={offset}")
        try:
            j = requests.get(url, timeout=45).json()
        except Exception as e:
            print(f"  [{name}] page {p} err {e}"); break
        for a in j.get("activities", []):
            smi, pv, st = a.get("canonical_smiles"), a.get("pchembl_value"), a.get("standard_type")
            if smi and pv and st in ("IC50", "Ki", "Kd", "EC50", "Potency"):
                rows.append({"smiles": smi, "pchembl": float(pv),
                             "year": a.get("document_year")})
        offset += PAGE
        if not j.get("page_meta", {}).get("next"):
            break
        time.sleep(0.25)
    json.dump(rows, open(cache, "w")); print(f"  [{name}] fetched: {len(rows)}")
    return rows


def build(name, rows):
    if not rows:
        return None
    df = pd.DataFrame(rows)
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    g = df.groupby("smiles").agg(pchembl=("pchembl", "median"), year=("year", "min")).reset_index()
    g["label"] = g["pchembl"].apply(lambda v: 1 if v >= 6.0 else (0 if v < 5.0 else -1))
    g = g[g["label"] >= 0].reset_index(drop=True)
    return g[["smiles", "label", "pchembl", "year"]]


def fetch_b3db():
    url = "https://raw.githubusercontent.com/theochem/B3DB/main/B3DB/B3DB_classification.tsv"
    try:
        df = pd.read_csv(io.StringIO(requests.get(url, timeout=60).text), sep="\t")
        col = "BBB+/BBB-" if "BBB+/BBB-" in df.columns else [c for c in df.columns if "BBB" in c][0]
        smi = "SMILES" if "SMILES" in df.columns else [c for c in df.columns if c.lower() == "smiles"][0]
        o = pd.DataFrame({"smiles": df[smi],
                          "label": df[col].astype(str).str.upper().str.contains("BBB\\+").astype(int)})
        o = o.dropna().drop_duplicates("smiles"); o.to_csv(f"{OUT}/BBB.csv", index=False)
        print(f"  [BBB] {len(o)} (pos={int(o.label.sum())})")
    except Exception as e:
        print(f"  [BBB] error {e}")


def main():
    print("B3DB (BBB)..."); fetch_b3db()
    print("\nChEMBL targets (with document_year)...")
    summary = {}
    for name, (tid, meaning) in TARGETS.items():
        try:
            rows = fetch_target_activities(name, tid)
        except Exception as e:
            print(f"  [{name}] err {e}"); continue
        g = build(name, rows)
        if g is not None and len(g) >= 100 and g["label"].nunique() == 2:
            g.to_csv(f"{OUT}/{name}.csv", index=False)
            yr = g["year"].dropna()
            summary[name] = {"n": len(g), "pos": int(g.label.sum()), "meaning": meaning,
                             "year_min": int(yr.min()) if len(yr) else None,
                             "year_max": int(yr.max()) if len(yr) else None,
                             "year_known": int(yr.notna().sum())}
            print(f"  [{name}] saved n={len(g)} pos={int(g.label.sum())} "
                  f"years {summary[name]['year_min']}-{summary[name]['year_max']} ({summary[name]['year_known']} dated)")
        else:
            print(f"  [{name}] insufficient -> skip")
    json.dump(summary, open(f"{OUT}/_summary.json", "w"), indent=2)
    print("\nDone.")


if __name__ == "__main__":
    main()
