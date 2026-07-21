"""
BS_fetch_clinical.py — gather REAL clinical/translational data (ChEMBL) to address
'engagement != efficacy': nervous-system (ATC level-1 N) molecules that reached a clinical
phase, with SMILES, max clinical phase, name, and ATC-derived disease class.
This is a measured clinical-precedent reference (NOT an efficacy prediction).
Out: data/clinical_cns_reference.csv
"""
import os, time
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import requests, pandas as pd
B = "https://www.ebi.ac.uk/chembl/api/data"

# ATC subcode -> disease class relevant to our endpoints
ATC_DISEASE = {
    "N06D": "Alzheimer's / dementia", "N04": "Parkinson's",
    "N06A": "Depression", "N05A": "Psychosis", "N06B": "Cognition/ADHD",
    "N05B": "Anxiety", "N03": "Epilepsy", "N07": "Other CNS",
}


def disease_for(codes):
    for c in codes:
        for k, v in ATC_DISEASE.items():
            if c.startswith(k):
                return v
    return "CNS (other)"


def main():
    rows, offset = [], 0
    while True:
        u = (f"{B}/molecule.json?atc_classifications__level1=N&max_phase__gte=1"
             f"&limit=1000&offset={offset}")
        try:
            j = requests.get(u, timeout=60).json()
        except Exception as e:
            print("err", e); break
        for m in j.get("molecules", []):
            smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if not smi:
                continue
            atc = [a for a in (m.get("atc_classifications") or []) if isinstance(a, str)]
            rows.append({"name": m.get("pref_name"), "smiles": smi,
                         "max_phase": m.get("max_phase"),
                         "atc": ";".join(atc),
                         "disease": disease_for(atc)})
        if not j.get("page_meta", {}).get("next"):
            break
        offset += 1000; time.sleep(0.2)
    df = pd.DataFrame(rows).dropna(subset=["smiles"]).drop_duplicates("smiles")
    df.to_csv("data/clinical_cns_reference.csv", index=False)
    print(f"Clinical CNS reference: {len(df)} molecules with phase+structure")
    print(df["disease"].value_counts().to_string())
    print("phase 4 (approved):", int((df.max_phase == 4).sum()))


if __name__ == "__main__":
    main()
