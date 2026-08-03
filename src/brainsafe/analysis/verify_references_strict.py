"""Resolve citations by EXACT TITLE against CrossRef and Europe PMC, with a strict match test.

A first attempt matched on a keyword appearing anywhere in the title and produced confidently wrong
results: a search for the random forest paper returned a 2026 epidemiology article, and a search for
scikit-learn returned an unrelated statistics-package validation. Relevance ranking is not identity.

This version therefore states the exact title of each intended work, queries CrossRef (which covers
the machine-learning venues that Europe PMC does not) and falls back to Europe PMC, and accepts a
record only when the returned title matches the requested one with a normalised similarity of at
least ACCEPT. Anything below that is reported unresolved and is not cited. Software without a
citable paper is recorded as a software entry rather than being attached to an unrelated article.

Output: manuscript/references_verified.json, manuscript/references.md
"""
from __future__ import annotations

import difflib
import json
import re
import time
from pathlib import Path

import requests

requests.packages.urllib3.disable_warnings()
ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "manuscript"
CROSSREF = "https://api.crossref.org/works"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
ACCEPT = 0.82
HEADERS = {"User-Agent": "BrainSafeAI/1.0 (research; mailto:krishnasalinig@sssihl.edu.in)"}

# key -> exact title of the intended work
TITLES = {
    "chembl": "The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods",
    "bindingdb": "BindingDB in 2015: A public database for medicinal chemistry, computational chemistry and systems pharmacology",
    "b3db": "B3DB: a curated blood-brain barrier dataset with unified binary and numerical experimental values",
    "ecfp": "Extended-Connectivity Fingerprints",
    "random_forest": "Random Forests",
    "sklearn": "Scikit-learn: Machine Learning in Python",
    "tdc": "Therapeutics Data Commons: Machine Learning Datasets and Tasks for Drug Discovery and Development",
    "moleculenet": "MoleculeNet: a benchmark for molecular machine learning",
    "kegg": "KEGG: kyoto encyclopedia of genes and genomes",
    "reactome": "The Reactome Pathway Knowledgebase 2024",
    "iuphar": "The IUPHAR/BPS Guide to PHARMACOLOGY in 2024",
    "bemis_murcko": "The properties of known drugs. 1. Molecular frameworks",
    "dude": "Directory of useful decoys, enhanced (DUD-E): better ligands and decoys for better benchmarking",
    "cns_mpo": "Moving beyond rules: the development of a central nervous system multiparameter optimization (CNS MPO) approach to enable alignment of druglike properties",
    "kpuu": "The use of drug transporter data in drug discovery",
    "xgboost": "XGBoost: A Scalable Tree Boosting System",
    "gin_gnn": "How Powerful are Graph Neural Networks?",
    "calibration": "Predicting good probabilities with supervised learning",
    "platt": "Probabilistic Outputs for Support Vector Machines and Comparisons to Regularized Likelihood Methods",
    "conformal": "Conformal Prediction in Drug Discovery",
    "ad_qsar": "QSAR applicability domain estimation by projection of the training set in descriptor space: a review",
    "herg_pred": "Deep learning based prediction of hERG blockers",
    "bbb_ml": "Prediction of blood-brain barrier penetration using machine learning",
    "delong": "Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach",
    "wilson_ci": "Probable Inference, the Law of Succession, and Statistical Inference",
    "cns_attrition": "Central nervous system drug discovery: challenges and opportunities",
    "ache_ad": "Acetylcholinesterase inhibitors for Alzheimer's disease",
    "bace1_fail": "BACE1 inhibitors: attractive therapeutics for Alzheimer's disease",
    "mao_b_pd": "Monoamine oxidase B inhibitors in Parkinson's disease",
    "lrrk2_pd": "LRRK2 inhibitors for Parkinson's disease",
    "orexin_insomnia": "Orexin receptor antagonists for the treatment of insomnia",
    "nlrp3_neuro": "NLRP3 inflammasome in neurodegenerative diseases",
    "nrf2_neuro": "The Keap1-Nrf2 pathway in neurodegenerative diseases",
    "hdac_hd": "Histone deacetylase inhibitors in Huntington's disease",
    "riluzole_als": "Riluzole for amyotrophic lateral sclerosis",
}
SOFTWARE = {
    "rdkit": {"kind": "software", "text": "RDKit: Open-source cheminformatics", "url": "https://www.rdkit.org"},
    "streamlit": {"kind": "software", "text": "Streamlit: an open-source app framework", "url": "https://streamlit.io"},
}


def _safe(x):
    return str(x).encode('ascii','replace').decode('ascii')


def norm(s):
    s = re.sub(r"<[^>]+>", "", s or "").lower()
    s = re.sub(r"[^a-z0-9 ]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def sim(a, b):
    return difflib.SequenceMatcher(None, norm(a), norm(b)).ratio()


def try_crossref(title):
    try:
        r = requests.get(CROSSREF, params={"query.bibliographic": title, "rows": 6},
                         timeout=15, headers=HEADERS, verify=False)
        if r.status_code != 200:
            return None
        best = None
        for it in r.json()["message"]["items"]:
            t = (it.get("title") or [""])[0]
            s = sim(t, title)
            if best is None or s > best[0]:
                auth = it.get("author") or []
                names = ", ".join(f"{a.get('family','')} {a.get('given','')[:1]}".strip()
                                  for a in auth[:3] if a.get("family"))
                if len(auth) > 3:
                    names += " et al."
                yr = ""
                for f in ("published-print", "published-online", "issued", "created"):
                    if it.get(f, {}).get("date-parts", [[None]])[0][0]:
                        yr = str(it[f]["date-parts"][0][0]); break
                best = (s, {"title": t.rstrip("."), "authors": names,
                            "journal": (it.get("container-title") or [""])[0],
                            "year": yr, "doi": it.get("DOI", ""), "source": "CrossRef"})
        return best
    except Exception:
        return None


def try_epmc(title):
    try:
        r = requests.get(EPMC, params={"query": f'TITLE:"{title}"', "format": "json",
                                       "pageSize": 6, "resultType": "core"},
                         timeout=15, verify=False)
        if r.status_code != 200:
            return None
        best = None
        for res in r.json().get("resultList", {}).get("result", []):
            t = res.get("title") or ""
            s = sim(t, title)
            if best is None or s > best[0]:
                best = (s, {"title": t.rstrip("."), "authors": res.get("authorString", ""),
                            "journal": (res.get("journalInfo", {}) or {}).get("journal", {}).get("title", ""),
                            "year": res.get("pubYear", ""), "doi": res.get("doi", ""),
                            "pmid": res.get("pmid", ""), "source": "EuropePMC"})
        return best
    except Exception:
        return None


def main():
    verified, unresolved = {}, []
    for key, title in TITLES.items():
        cands = [c for c in (try_crossref(title), try_epmc(title)) if c]
        best = max(cands, key=lambda c: c[0]) if cands else None
        if best and best[0] >= ACCEPT:
            rec = best[1]; rec["match_score"] = round(best[0], 3); rec["requested_title"] = title
            verified[key] = rec
            print(_safe(f"OK   {key:16} {best[0]:.2f}  {rec['year']}  {rec['title'][:58]}"), flush=True)
        else:
            unresolved.append(key)
            got = best[1]["title"][:52] if best else "no result"
            print(_safe(f"MISS {key:16} {(best[0] if best else 0):.2f}  best was: {got}"), flush=True)
        time.sleep(0.05)

    payload = {"papers": verified, "software": SOFTWARE, "unresolved": unresolved,
               "accept_threshold": ACCEPT}
    (OUT / "references_verified.json").write_text(json.dumps(payload, indent=2))

    lines = ["# References", "",
             f"Each entry was resolved by exact-title query against CrossRef or Europe PMC and "
             f"accepted only above a normalised title-similarity of {ACCEPT}. The requested title, "
             f"the matched title and the similarity score are recorded in "
             f"`references_verified.json`, so every entry can be re-checked mechanically. "
             f"None is written from memory.", ""]
    for i, (k, v) in enumerate(sorted(verified.items(), key=lambda kv: (kv[1]["year"] or "0")), 1):
        bits = [x for x in [v["authors"].rstrip(". ") if v["authors"] else "",
                            v["title"], v["journal"], v["year"]] if x]
        ref = ". ".join(bits)
        if v.get("doi"):
            ref += f". doi:{v['doi']}"
        elif v.get("pmid"):
            ref += f". PMID:{v['pmid']}"
        lines.append(f"{i}. {ref}")
    lines += ["", "## Software", ""]
    for k, v in SOFTWARE.items():
        lines.append(f"- {v['text']}. {v['url']}")
    if unresolved:
        lines += ["", f"Requested but not resolved above the similarity threshold, and therefore "
                      f"not cited: {', '.join(unresolved)}."]
    (OUT / "references.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"\nverified {len(verified)}/{len(TITLES)}; unresolved: {unresolved}")


if __name__ == "__main__":
    main()
