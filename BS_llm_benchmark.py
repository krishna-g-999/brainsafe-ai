"""
BS_llm_benchmark.py -- ground-truth-anchored head-to-head between BrainSafe and general LLMs.

We fix a panel of compounds whose relevant pharmacology is uncontested (drug labels / primary
literature), record BrainSafe's calibrated predictions + measured-analogue provenance, and emit:
  - BS_llm_benchmark_groundtruth.json  (panel, established truth, BrainSafe predictions)
The user runs the fixed prompt (BS_LLM_benchmark_protocol.md) on named LLMs and pastes the
responses; BS_llm_score.py then scores every system against the same key.

Ground truth is limited to items that are not in scientific dispute; ambiguous items are marked
"not_scored" so the benchmark stays objective.
"""
import os, json
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import requests, urllib3
urllib3.disable_warnings()
from rdkit import Chem
import BS_brain_predict as B

def resolve(name):
    """Canonical SMILES from PubChem (authoritative); returns None on failure."""
    try:
        u = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
             f"{requests.utils.quote(name)}/property/CanonicalSMILES/TXT")
        r = requests.get(u, timeout=30, verify=False)
        if r.ok and r.text.strip():
            smi = r.text.strip().splitlines()[0].strip()
            if Chem.MolFromSmiles(smi):
                return smi
    except Exception:
        pass
    return None

# canonical, uncontested pharmacology. bbb: 1 penetrant / 0 non-penetrant; herg: 1 clinically
# relevant block / 0 no relevant block / None not-scored; target: primary CNS target in our panel.
PANEL = [
 {"name":"Donepezil","smiles":"O=C1CC(CC2CCN(Cc3ccccc3)CC2)c2cc(OC)c(OC)cc21",
  "bbb":1,"herg":None,"target":"AChE","note":"AChE inhibitor (Alzheimer's); CNS-penetrant. hERG NOT scored: signal (IC50 ~5-35 uM) is weak relative to therapeutic Cmax, so clinical relevance is contestable"},
 {"name":"Galantamine","smiles":"CO[C@@H]1[C@@H]2CC3=C4C1=C(O)C(OC)=CC4=CC=C3C[N@@]2CC=C[C@@H]1O", # placeholder-corrected below
  "bbb":1,"herg":0,"target":"AChE","note":"AChE inhibitor; CNS-penetrant; no major hERG liability"},
 {"name":"Rivastigmine","smiles":"CCN(C)C(=O)Oc1cccc([C@@H](C)N(C)C)c1",
  "bbb":1,"herg":0,"target":"AChE","note":"AChE/BChE inhibitor; CNS-penetrant"},
 {"name":"Rasagiline","smiles":"C#CC[NH2+][C@H]1CCc2ccccc21", # will be neutralised on canon
  "bbb":1,"herg":None,"target":"MAO_B","note":"MAO-B inhibitor (Parkinson's); CNS-penetrant"},
 {"name":"Selegiline","smiles":"C#CCN(C)[C@@H](C)Cc1ccccc1",
  "bbb":1,"herg":None,"target":"MAO_B","note":"MAO-B inhibitor (Parkinson's); CNS-penetrant"},
 {"name":"Fluoxetine","smiles":"CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
  "bbb":1,"herg":None,"target":"SERT","note":"SERT inhibitor (depression); CNS-penetrant"},
 {"name":"Terfenadine","smiles":"CC(C)(C)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1",
  "bbb":0,"herg":1,"target":None,"note":"peripheral H1 antihistamine; classic hERG blocker (withdrawn); poor CNS entry"},
 {"name":"Astemizole","smiles":"COc1ccc(CCN2CCC(Nc3nc4ccccc4n3Cc3ccc(F)cc3)CC2)cc1",
  "bbb":0,"herg":1,"target":None,"note":"peripheral H1 antihistamine; potent hERG blocker (withdrawn)"},
 {"name":"Quercetin","smiles":"O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
  "bbb":0,"herg":0,"target":None,"note":"flavonoid antioxidant; poor BBB penetration; no relevant hERG block"},
 {"name":"Novel-arylpiperazine (unpublished)","smiles":"O=C1CCc2ccccc2N1CCN1CCN(c2ccccc2)CC1",
  "bbb":None,"herg":None,"target":None,"note":"hypothetical unpublished scaffold; no established ground truth -> tests confident hallucination vs honest uncertainty"},
]

def bs_row(smiles):
    try:
        prof = B.predict_brain_profile(smiles)
    except Exception as e:
        return {"error": str(e)}
    eps = {e["endpoint"]: e for e in prof.get("endpoints", [])}
    def cap(ep):
        e = eps.get(ep)
        if not e: return None
        ev = (e.get("evidence") or [{}])[0]
        return {"probability": e.get("probability"), "call": e.get("call"),
                "in_domain": e.get("in_domain"),
                "conformal_set": (e.get("conformal") or {}).get("set"),
                "nearest_analogue_tanimoto": ev.get("similarity"),
                "nearest_analogue_pchembl": ev.get("pchembl")}
    # primary CNS target = highest-probability target endpoint
    tgts = {k: v for k, v in eps.items() if v.get("kind") == "target"}
    prim = max(tgts.items(), key=lambda kv: kv[1].get("probability", 0)) if tgts else (None, {})
    return {"BBB": cap("BBB"), "hERG": cap("hERG"),
            "predicted_primary_target": prim[0],
            "predicted_primary_target_prob": prim[1].get("probability") if prim[1] else None}

out = []
for c in PANEL:
    row = dict(c)
    # prefer authoritative PubChem structure for named drugs; keep given SMILES for the novel one
    resolved = None if c["name"].startswith("Novel") else resolve(c["name"])
    row["smiles_source"] = "PubChem" if resolved else "curated"
    row["smiles"] = resolved or c["smiles"]
    row["brainsafe"] = bs_row(row["smiles"])
    out.append(row)
    bs = row["brainsafe"]
    bbb = (bs.get("BBB") or {}).get("call"); herg = (bs.get("hERG") or {}).get("call")
    print(f"{c['name']:32} [{row['smiles_source']:7}] BBB={bbb} hERG={herg} primary={bs.get('predicted_primary_target')}")

json.dump({"panel": out,
           "truth_provenance": "Established pharmacology (approved-drug labels and primary literature); "
                               "only uncontested items scored, ambiguous items marked null.",
           "scoring_dimensions": ["BBB_call_accuracy","hERG_call_accuracy","primary_target_accuracy",
                                  "provenance_validity(does cited ChEMBL id / pChEMBL exist)",
                                  "calibration(Brier if probabilities given)",
                                  "honest_uncertainty_on_novel_compound"]},
          open("BS_llm_benchmark_groundtruth.json", "w"), indent=2)
print("\nWrote BS_llm_benchmark_groundtruth.json")
