"""
BrainSafe AI — Live External Data Client
Queries PubChem REST API and ChEMBL REST API for any compound by name.
Computes BBB penetration prediction using simplified CNS-MPO-inspired rules.

References:
- PubChem REST API: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- ChEMBL REST API: https://www.ebi.ac.uk/chembl/api/data/docs
- Wager TT et al. (2010). CNS MPO approach. ACS Chem Neurosci 1(6):435-449.
- Lipinski CA et al. (2001). Rule of Five. Adv Drug Deliv Rev 46:3-26.
"""

import requests
from functools import lru_cache

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL_BASE  = "https://www.ebi.ac.uk/chembl/api/data"
TIMEOUT      = 10


@lru_cache(maxsize=256)
def fetch_pubchem(compound_name: str) -> dict:
    """
    Fetch compound properties from PubChem by compound name.
    Returns a dict with CID, physicochemical properties, synonyms, and BBB prediction.
    On failure, returns {"error": <message>}.
    """
    result = {}
    try:
        name_enc = requests.utils.quote(compound_name)

        r = requests.get(
            f"{PUBCHEM_BASE}/compound/name/{name_enc}/JSON",
            timeout=TIMEOUT
        )
        if r.status_code == 404:
            return {"error": f"'{compound_name}' was not found in the PubChem database."}
        if r.status_code != 200:
            return {"error": f"PubChem returned status {r.status_code}. Try again later."}

        compounds = r.json().get("PC_Compounds", [])
        if not compounds:
            return {"error": "PubChem returned no compound data."}

        cid = compounds[0].get("id", {}).get("id", {}).get("cid")
        if not cid:
            return {"error": "Could not extract CID from PubChem response."}

        result["cid"]         = cid
        result["pubchem_url"] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"

        props_url = (
            f"{PUBCHEM_BASE}/compound/cid/{cid}/property/"
            "MolecularFormula,MolecularWeight,XLogP,TPSA,"
            "HBondDonorCount,HBondAcceptorCount,RotatableBondCount,"
            "HeavyAtomCount,IsomericSMILES,InChIKey/JSON"
        )
        pr = requests.get(props_url, timeout=TIMEOUT)
        if pr.status_code == 200:
            props = pr.json().get("PropertyTable", {}).get("Properties", [{}])[0]
            result.update({
                "formula":          props.get("MolecularFormula", "N/A"),
                "mw":               props.get("MolecularWeight",  "N/A"),
                "xlogp":            props.get("XLogP",            "N/A"),
                "tpsa":             props.get("TPSA",             "N/A"),
                "hbd":              props.get("HBondDonorCount",  "N/A"),
                "hba":              props.get("HBondAcceptorCount","N/A"),
                "rotatable_bonds":  props.get("RotatableBondCount","N/A"),
                "heavy_atoms":      props.get("HeavyAtomCount",   "N/A"),
                "smiles":           props.get("IsomericSMILES",   "N/A"),
                "inchikey":         props.get("InChIKey",         "N/A"),
            })

        syn_r = requests.get(
            f"{PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON",
            timeout=TIMEOUT
        )
        if syn_r.status_code == 200:
            syns = syn_r.json().get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
            result["synonyms"] = syns[:8]

        result["bbb_prediction"] = _predict_bbb(result)

    except requests.exceptions.Timeout:
        return {"error": "PubChem API timed out. Please try again."}
    except Exception as exc:
        return {"error": f"PubChem query failed: {exc}"}

    return result


@lru_cache(maxsize=256)
def fetch_chembl(compound_name: str) -> dict:
    """
    Fetch compound data from ChEMBL by preferred name.
    Returns molecule properties, drug approval phase, Ro5 compliance,
    QED score, and mechanism-of-action entries.
    """
    result = {}
    headers = {"Accept": "application/json"}

    try:
        name_enc = requests.utils.quote(compound_name)
        r = requests.get(
            f"{CHEMBL_BASE}/molecule?pref_name__iexact={name_enc}&format=json&limit=1",
            timeout=TIMEOUT, headers=headers
        )
        molecules = []
        if r.status_code == 200:
            molecules = r.json().get("molecules", [])

        if not molecules:
            r2 = requests.get(
                f"{CHEMBL_BASE}/molecule?molecule_synonyms__synonym__iexact={name_enc}&format=json&limit=1",
                timeout=TIMEOUT, headers=headers
            )
            if r2.status_code == 200:
                molecules = r2.json().get("molecules", [])

        if not molecules:
            return {"error": f"'{compound_name}' was not found in ChEMBL."}

        mol       = molecules[0]
        chembl_id = mol.get("molecule_chembl_id", "")
        props     = mol.get("molecule_properties") or {}

        result.update({
            "chembl_id":    chembl_id,
            "chembl_url":   f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/",
            "mol_type":     mol.get("molecule_type", "N/A"),
            "max_phase":    mol.get("max_phase"),
            "oral":         mol.get("oral"),
            "alogp":        props.get("alogp",              "N/A"),
            "hba":          props.get("hba",                "N/A"),
            "hbd":          props.get("hbd",                "N/A"),
            "mw":           props.get("mw_freebase",        "N/A"),
            "tpsa":         props.get("psa",                "N/A"),
            "ro5_violations": props.get("num_ro5_violations","N/A"),
            "qed":          props.get("qed_weighted",        "N/A"),
            "ro5_pass":     (props.get("num_ro5_violations") or 99) <= 1,
        })

        if chembl_id:
            mr = requests.get(
                f"{CHEMBL_BASE}/mechanism?molecule_chembl_id={chembl_id}&format=json&limit=10",
                timeout=TIMEOUT, headers=headers
            )
            if mr.status_code == 200:
                result["mechanisms"] = [
                    {
                        "target": m.get("target_name", "Unknown target"),
                        "action": m.get("action_type",  "N/A"),
                    }
                    for m in mr.json().get("mechanisms", [])[:8]
                ]

    except requests.exceptions.Timeout:
        return {"error": "ChEMBL API timed out. Please try again."}
    except Exception as exc:
        return {"error": f"ChEMBL query failed: {exc}"}

    return result


def _predict_bbb(props: dict) -> dict:
    """
    Predict BBB penetration using a simplified CNS-MPO-inspired method.

    CNS-MPO reference: Wager TT et al. (2010).
    ACS Chem Neurosci 1(6):435-449. DOI:10.1021/cn100008c

    Lipinski Ro5 reference: Lipinski CA et al. (2001).
    Adv Drug Deliv Rev 46:3-26. DOI:10.1016/S0169-409X(00)00129-0
    """
    try:
        mw   = float(props.get("mw",   999))
        logp = float(props.get("xlogp", 5)) if props.get("xlogp") not in (None, "N/A") else 5.0
        tpsa = float(props.get("tpsa", 999)) if props.get("tpsa")  not in (None, "N/A") else 999.0
        hbd  = int(props.get("hbd",    99))  if props.get("hbd")   not in (None, "N/A") else 99
        hba  = int(props.get("hba",    99))  if props.get("hba")   not in (None, "N/A") else 99
    except (TypeError, ValueError):
        return {
            "level": "Unknown",
            "rationale": "Insufficient physicochemical data for BBB prediction.",
            "ref": ""
        }

    ro5_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

    if mw <= 360 and 1.0 <= logp <= 3.0 and tpsa <= 60 and hbd <= 1:
        level = "High"
        colour = "#14532D"
        rationale = (
            f"MW {mw:.0f} Da, XLogP {logp:.1f}, TPSA {tpsa:.0f} Å², HBD {hbd} — "
            f"all within the CNS-MPO optimal range. Consistent with high passive BBB penetration."
        )
    elif mw <= 450 and 0.0 <= logp <= 4.0 and tpsa <= 90 and hbd <= 3:
        level = "Medium"
        colour = "#92710C"
        rationale = (
            f"MW {mw:.0f} Da, XLogP {logp:.1f}, TPSA {tpsa:.0f} Å², HBD {hbd} — "
            f"moderate CNS penetration expected. Partial passive diffusion; active transport may contribute."
        )
    elif mw <= 500 and tpsa <= 120:
        level = "Low-Med"
        colour = "#7C3009"
        rationale = (
            f"MW {mw:.0f} Da, TPSA {tpsa:.0f} Å² — borderline CNS penetration. "
            f"Likely limited passive diffusion; carrier-mediated transport may be required."
        )
    else:
        level = "Low"
        colour = "#7F1D1D"
        rationale = (
            f"MW {mw:.0f} Da, XLogP {logp:.1f}, TPSA {tpsa:.0f} Å², HBD {hbd} — "
            f"poor predicted CNS penetration due to high molecular weight, polarity, or hydrogen-bonding capacity."
        )

    return {
        "level":        level,
        "colour":       colour,
        "rationale":    rationale,
        "ro5_violations": ro5_violations,
        "ref":          (
            "Wager TT et al. (2010). CNS MPO approach. "
            "ACS Chem Neurosci 1(6):435–449. | "
            "Lipinski CA et al. (2001). Ro5. Adv Drug Deliv Rev 46:3–26."
        )
    }


@lru_cache(maxsize=256)
def fetch_kegg_pathways(compound_name: str) -> dict:
    """
    Fetch KEGG compound ID and associated pathway links for a compound.

    Pipeline:
    1. Search KEGG compound database by name → extract KEGG compound ID (e.g., C10443)
    2. Retrieve pathway associations via KEGG LINK API
    3. Batch-fetch pathway names via KEGG LIST API

    Returns a dict with:
      kegg_id   — KEGG compound identifier
      kegg_url  — link to KEGG compound page
      pathways  — list of {id, name, url, category} dicts
      error     — present if lookup failed

    Reference: Kanehisa M & Goto S (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes.
    Nucleic Acids Res 28(1):27–30. DOI:10.1093/nar/28.1.27
    """
    KEGG_BASE = "https://rest.kegg.jp"
    result    = {}
    try:
        name_enc = requests.utils.quote(compound_name)

        r = requests.get(f"{KEGG_BASE}/find/compound/{name_enc}", timeout=TIMEOUT)
        if r.status_code != 200 or not r.text.strip():
            return {"error": f"'{compound_name}' not found in KEGG compound database."}

        first_line = r.text.strip().split("\n")[0]
        kegg_id    = first_line.split("\t")[0].replace("cpd:", "")
        kegg_name  = first_line.split("\t")[1] if "\t" in first_line else compound_name

        result["kegg_id"]   = kegg_id
        result["kegg_name"] = kegg_name.split(";")[0].strip()
        result["kegg_url"]  = f"https://www.kegg.jp/compound/{kegg_id}"

        r2 = requests.get(f"{KEGG_BASE}/link/pathway/cpd:{kegg_id}", timeout=TIMEOUT)
        if r2.status_code != 200 or not r2.text.strip():
            result["pathways"] = []
            return result

        pathway_ids = []
        for line in r2.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                pathway_ids.append(parts[1].replace("path:", ""))

        if not pathway_ids:
            result["pathways"] = []
            return result

        human_ids  = [p for p in pathway_ids if p.startswith("hsa")][:12]
        map_ids    = [p for p in pathway_ids if p.startswith("map")][:6]
        fetch_ids  = (human_ids + map_ids) if human_ids else map_ids[:15]

        pathways = []
        if fetch_ids:
            batch = "+".join(fetch_ids)
            r3 = requests.get(f"{KEGG_BASE}/list/{batch}", timeout=TIMEOUT)
            if r3.status_code == 200:
                for line in r3.text.strip().split("\n"):
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        pid   = parts[0].replace("path:", "")
                        pname = parts[1].split(" - Homo sapiens")[0].strip()
                        cat   = "Human (HSA)" if pid.startswith("hsa") else "Reference Map"
                        pathways.append({
                            "id":       pid,
                            "name":     pname,
                            "url":      f"https://www.kegg.jp/pathway/{pid}",
                            "category": cat,
                        })

        result["pathways"] = pathways

    except requests.exceptions.Timeout:
        return {"error": "KEGG API timed out. Please try again."}
    except Exception as exc:
        return {"error": f"KEGG query failed: {exc}"}

    return result
