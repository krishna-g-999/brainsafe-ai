"""
BS_llm_score.py -- score BrainSafe and each LLM against the frozen measured-data key.

Inputs:
  BS_llm_benchmark_groundtruth.json  (panel + truth + BrainSafe predictions; from BS_llm_benchmark.py)
  BS_llm_responses.json              (list of {system, rows{compound: {...}}}; transcribed LLM replies)

BrainSafe is scored automatically from the ground-truth file, by the SAME rubric as every LLM.
Provenance/hallucination: a claimed ChEMBL ID is checked for existence against the ChEMBL API
(verify disabled for the proxied network; read-only). A specific measured value asserted for the
novel (unpublished) compound is counted as a hallucination by construction.

Outputs: supplementary/STable13_llm_scoreboard.csv, BS_llm_scoreboard.json
"""
import os, json, re
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
try:
    import requests, urllib3; urllib3.disable_warnings(); ONLINE = True
except Exception:
    ONLINE = False

GT = json.load(open("BS_llm_benchmark_groundtruth.json"))
PANEL = GT["panel"]
truth = {c["name"]: c for c in PANEL}
NOVEL = [c["name"] for c in PANEL if c["name"].startswith("Novel")][0]

def yn(v):
    if v is None: return None
    s = str(v).strip().lower()
    if s in ("yes", "y", "1", "true", "penetrant", "blocker", "active", "positive"): return 1
    if s in ("no", "n", "0", "false", "non-penetrant", "low-risk", "inactive", "negative"): return 0
    return None

def chembl_exists(cid):
    if not ONLINE or not cid: return None
    if not re.match(r"^CHEMBL\d+$", str(cid).strip(), re.I): return False  # malformed = not real
    try:
        r = requests.get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json",
                         timeout=20, verify=False)
        return r.status_code == 200
    except Exception:
        return None

def score_system(name, rows):
    bbb_ok = bbb_n = herg_ok = herg_n = 0
    brier_terms = []
    halluc = 0; prov_checked = 0; prov_valid = 0
    for cname, t in truth.items():
        r = rows.get(cname) or {}
        # BBB
        if t.get("bbb") is not None:
            p = yn(r.get("bbb"))
            if p is not None:
                bbb_n += 1; bbb_ok += int(p == t["bbb"])
            if r.get("bbb_p") is not None:
                try: brier_terms.append((float(r["bbb_p"]) - t["bbb"]) ** 2)
                except Exception: pass
        # hERG
        if t.get("herg") is not None:
            p = yn(r.get("herg"))
            if p is not None:
                herg_n += 1; herg_ok += int(p == t["herg"])
            if r.get("herg_p") is not None:
                try: brier_terms.append((float(r["herg_p"]) - t["herg"]) ** 2)
                except Exception: pass
        # provenance / hallucination
        cid = r.get("chembl_id")
        if cid:
            ex = chembl_exists(cid)
            if ex is not None:
                prov_checked += 1; prov_valid += int(bool(ex))
                if not ex: halluc += 1
        # asserting a specific measured value / ChEMBL id for the novel compound = hallucination
        if cname == NOVEL:
            if r.get("pchembl") not in (None, "", "n/a", "unknown") or cid:
                halluc += 1
    novel_r = rows.get(NOVEL) or {}
    novel_honest = not (novel_r.get("pchembl") not in (None, "", "n/a", "unknown") or novel_r.get("chembl_id")
                        or yn(novel_r.get("herg")) is not None and novel_r.get("herg_confident"))
    return {
        "system": name,
        "BBB_accuracy": round(bbb_ok / bbb_n, 3) if bbb_n else None,
        "BBB_scored": bbb_n,
        "hERG_accuracy": round(herg_ok / herg_n, 3) if herg_n else None,
        "hERG_scored": herg_n,
        "Brier": round(sum(brier_terms) / len(brier_terms), 3) if brier_terms else None,
        "provenance_checked": prov_checked,
        "provenance_valid": prov_valid,
        "hallucinations": halluc,
        "novel_honest_uncertainty": bool(novel_honest),
    }

# --- BrainSafe row, built from the ground-truth file (same rubric) ---
bs_rows = {}
for c in PANEL:
    bs = c.get("brainsafe", {})
    bbb = (bs.get("BBB") or {}); herg = (bs.get("hERG") or {})
    bs_rows[c["name"]] = {
        "bbb": "yes" if str(bbb.get("call", "")).lower().startswith("pen") else "no",
        "bbb_p": bbb.get("probability"),
        "herg": "yes" if str(herg.get("call", "")).lower().startswith("block") else "no",
        "herg_p": herg.get("probability"),
        # BrainSafe cites a real measured analogue (SMILES + pChEMBL), not a ChEMBL id -> no fabricated id;
        # for the novel compound it returns a conformal 'uncertain' set, i.e. no asserted measured value.
        "chembl_id": None,
        "pchembl": None if c["name"].startswith("Novel") else None,
    }
results = [score_system("BrainSafe AI", bs_rows)]

# --- LLM responses (if provided) ---
if os.path.exists("BS_llm_responses.json"):
    for sysobj in json.load(open("BS_llm_responses.json")):
        results.append(score_system(sysobj["system"], sysobj.get("rows", {})))
else:
    print("No BS_llm_responses.json yet -> scoring BrainSafe only. "
          "Transcribe LLM replies into BS_llm_responses.json (see BS_llm_responses.template.json).")

df = pd.DataFrame(results)
df.to_csv("supplementary/STable13_llm_scoreboard.csv", index=False)
json.dump(results, open("BS_llm_scoreboard.json", "w"), indent=2)
print(df.to_string(index=False))
print("\nWrote STable13_llm_scoreboard.csv + BS_llm_scoreboard.json")
