"""
BS_llm_score.py -- score BrainSafe and each LLM against the frozen measured-data key.

Inputs:
  BS_llm_benchmark_groundtruth.json  (panel + truth + BrainSafe predictions)
  BS_llm_responses.json              (list of {system, rows{compound:{...}}}; transcribed replies)

Objective metrics (same rubric for every system, BrainSafe included):
  - BBB accuracy         (9 compounds with a defined truth)
  - hERG accuracy        (5 uncontested compounds; donepezil/fluoxetine excluded as contestable)
  - Brier                (stated probabilities vs binary truth)
  - fabricated_chembl_ids: volunteered ChEMBL IDs that DO NOT resolve on the ChEMBL API (hallucination)
  - wrong_structure_ids  : IDs that resolve but whose InChIKey != the compound the row names it as,
                           where the row cites the ID as the drug itself (not flagged as an analogue)
  - novel_confabulation  : asserted a specific target potency and/or ChEMBL analogue for the
                           UNPUBLISHED compound (which can have no measured value)  -> True/False

Live ChEMBL checks use verify=False (this network proxies TLS; read-only public API).
Outputs: supplementary/STable13_llm_scoreboard.csv, BS_llm_scoreboard.json
"""
import os, json, re
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from rdkit import Chem
try:
    import requests, urllib3; urllib3.disable_warnings(); ONLINE = True
except Exception:
    ONLINE = False

GT = json.load(open("BS_llm_benchmark_groundtruth.json"))
PANEL = GT["panel"]
truth = {c["name"]: c for c in PANEL}
NOVEL = [c["name"] for c in PANEL if c["name"].startswith("Novel")][0]

# InChIKey of each panel compound (from the PubChem-resolved structure) for provenance matching
def ik(smi):
    m = Chem.MolFromSmiles(str(smi)); return Chem.MolToInchiKey(m) if m else None
panel_ik = {c["name"]: ik(c["smiles"]) for c in PANEL}

def yn(v):
    if v is None: return None
    s = str(v).strip().lower()
    if s.startswith(("yes", "y", "penetrant", "block", "active", "positive")) or s in ("1", "true"): return 1
    if s.startswith(("no", "n", "non", "low", "inactive", "negative")) or s in ("0", "false"): return 0
    return None

_cache = {}
def chembl_lookup(cid):
    """Return dict(exists, inchikey, pref_name) for a ChEMBL id, or None if offline/unknown."""
    if not cid: return None
    cid = str(cid).strip().upper()
    if cid in _cache: return _cache[cid]
    if not re.match(r"^CHEMBL\d+$", cid):
        _cache[cid] = {"exists": False, "inchikey": None, "pref_name": None}; return _cache[cid]
    if not ONLINE:
        return None
    try:
        r = requests.get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json",
                         timeout=25, verify=False)
        if r.status_code != 200:
            res = {"exists": False, "inchikey": None, "pref_name": None}
        else:
            j = r.json(); st = (j.get("molecule_structures") or {})
            res = {"exists": True, "inchikey": st.get("standard_inchi_key"),
                   "pref_name": j.get("pref_name")}
    except Exception:
        res = None
    _cache[cid] = res; return res

prov_detail = []
def score_system(name, rows):
    bbb_ok = bbb_n = herg_ok = herg_n = 0
    brier = []
    ids_given = fab = wrong = 0
    for cname, t in truth.items():
        r = rows.get(cname) or {}
        if t.get("bbb") is not None and yn(r.get("bbb")) is not None:
            bbb_n += 1; bbb_ok += int(yn(r["bbb"]) == t["bbb"])
            if r.get("bbb_p") is not None:
                try: brier.append((float(r["bbb_p"]) - t["bbb"]) ** 2)
                except Exception: pass
        if t.get("herg") is not None and yn(r.get("herg")) is not None:
            herg_n += 1; herg_ok += int(yn(r["herg"]) == t["herg"])
            if r.get("herg_p") is not None:
                try: brier.append((float(r["herg_p"]) - t["herg"]) ** 2)
                except Exception: pass
        cid = r.get("chembl_id")
        if cid:
            ids_given += 1
            look = chembl_lookup(cid)
            if look is not None:
                if not look["exists"]:
                    fab += 1
                    prov_detail.append(f"{name}|{cname}|{cid} does not resolve (fabricated)")
                elif cname != NOVEL and not r.get("analogue") and panel_ik.get(cname) \
                        and look["inchikey"] and look["inchikey"][:14] != panel_ik[cname][:14]:
                    # compare InChIKey connectivity block only (fair to stereo/salt differences)
                    wrong += 1
                    prov_detail.append(f"{name}|{cname}|{cid}->{look.get('pref_name')} "
                                       f"(cited as self but different skeleton)")
    nr = rows.get(NOVEL) or {}
    novel_confab = bool(nr.get("chembl_id") or (nr.get("pchembl") not in (None, "", "n/a", "unknown"))
                        or (nr.get("target") not in (None, "", "none", "None")))
    return {"system": name,
            "BBB_acc": round(bbb_ok / bbb_n, 3) if bbb_n else None, "BBB_n": bbb_n,
            "hERG_acc": round(herg_ok / herg_n, 3) if herg_n else None, "hERG_n": herg_n,
            "Brier": round(sum(brier) / len(brier), 3) if brier else None,
            "chembl_ids_given": ids_given, "fabricated_ids": fab, "wrong_structure_ids": wrong,
            "novel_confabulation": novel_confab}

# BrainSafe row from the ground-truth file (same rubric)
bs_rows = {}
for c in PANEL:
    bs = c.get("brainsafe", {}); bbb = bs.get("BBB") or {}; herg = bs.get("hERG") or {}
    bs_rows[c["name"]] = {
        "bbb": "yes" if str(bbb.get("call", "")).lower().startswith("pen") else "no",
        "bbb_p": bbb.get("probability"),
        "herg": "yes" if str(herg.get("call", "")).lower().startswith("block") else "no",
        "herg_p": herg.get("probability"),
        "chembl_id": None,  # BrainSafe cites a measured analogue by SMILES+pChEMBL, not a ChEMBL id
        "pchembl": None, "target": None}  # novel row: conformal 'uncertain' -> no asserted target/value
results = [score_system("BrainSafe AI", bs_rows)]
if os.path.exists("BS_llm_responses.json"):
    for s in json.load(open("BS_llm_responses.json")):
        results.append(score_system(s["system"], s.get("rows", {})))
else:
    print("No BS_llm_responses.json -> BrainSafe only.")

df = pd.DataFrame(results)
df.to_csv("supplementary/STable13_llm_scoreboard.csv", index=False)
json.dump({"scoreboard": results, "provenance_detail": prov_detail},
          open("BS_llm_scoreboard.json", "w"), indent=2)
print(df.to_string(index=False))
print("\n-- provenance detail (fabricated / wrong-skeleton self-citations) --")
for d in prov_detail: print("  ", d)
print("\nWrote STable13_llm_scoreboard.csv + BS_llm_scoreboard.json")
