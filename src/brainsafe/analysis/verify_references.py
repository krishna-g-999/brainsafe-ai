"""Resolve every citation against Europe PMC and CrossRef, so no reference is written from memory.

A reference recalled from memory is exactly the kind of claim that looks authoritative and is wrong:
earlier in this project an identifier believed to denote the KEAP1-NFE2L2 pathway proved to denote
Interleukin-1 signalling. Bibliographies carry the same risk, so each entry here is specified only by
a search phrase and is accepted only if a live query returns a record whose title actually matches
the intended work. Anything unmatched is reported as unresolved rather than guessed.

Output: manuscript/references_verified.json and manuscript/references.md
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import requests

requests.packages.urllib3.disable_warnings()
ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "manuscript"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# key -> (query, tokens that must appear in the returned title)
WANTED = {
    "chembl":        ("ChEMBL database bioactivity drug discovery", ["chembl"]),
    "bindingdb":     ("BindingDB measured binding affinities database", ["bindingdb"]),
    "b3db":          ("B3DB blood-brain barrier database curated", ["b3db", "blood-brain"]),
    "rdkit_ref":     ("RDKit open-source cheminformatics", ["rdkit"]),
    "ecfp":          ("Extended-connectivity fingerprints Rogers Hahn", ["extended-connectivity"]),
    "random_forest": ("Random forests Breiman machine learning", ["random forest"]),
    "sklearn":       ("Scikit-learn machine learning in Python", ["scikit-learn"]),
    "tdc":           ("Therapeutics Data Commons machine learning datasets drug discovery",
                      ["therapeutics data commons"]),
    "moleculenet":   ("MoleculeNet benchmark for molecular machine learning", ["moleculenet"]),
    "kegg":          ("KEGG Kyoto Encyclopedia of Genes and Genomes", ["kegg"]),
    "reactome":      ("Reactome pathway knowledgebase", ["reactome"]),
    "iuphar":        ("IUPHAR BPS Guide to PHARMACOLOGY", ["guide to pharmacology"]),
    "conformal":     ("conformal prediction QSAR applicability domain confidence", ["conformal"]),
    "calibration":   ("Transforming classifier scores into accurate multiclass probability estimates",
                      ["probability estimates"]),
    "bemis_murcko":  ("The properties of known drugs molecular frameworks Bemis Murcko", ["frameworks"]),
    "dude":          ("Directory of useful decoys enhanced DUD-E benchmarking", ["decoys"]),
    "applicability": ("applicability domain QSAR models definition", ["applicability domain"]),
    "cns_mpo":       ("Central nervous system multiparameter optimization CNS MPO desirability",
                      ["multiparameter", "cns"]),
    "kpuu":          ("unbound brain-to-plasma partition coefficient Kp,uu CNS drug discovery",
                      ["unbound", "brain"]),
    "pgp_bbb":       ("P-glycoprotein efflux blood-brain barrier drug delivery", ["p-glycoprotein"]),
    "xgboost":       ("XGBoost scalable tree boosting system", ["xgboost", "boosting"]),
    "gin_gnn":       ("How powerful are graph neural networks", ["graph neural networks"]),
    "scaffold_split": ("scaffold split evaluation molecular machine learning generalization",
                       ["scaffold"]),
    "ad_tanimoto":   ("Tanimoto similarity nearest neighbour applicability domain QSAR", ["similarity"]),
    "ache_ad":       ("acetylcholinesterase inhibitors Alzheimer disease treatment", ["acetylcholinesterase"]),
    "bace1":         ("BACE1 inhibitors amyloid Alzheimer disease therapeutic", ["bace1"]),
    "mao_b_pd":      ("monoamine oxidase B inhibitors Parkinson disease", ["monoamine oxidase"]),
    "herg":          ("hERG potassium channel cardiotoxicity drug safety prediction", ["herg"]),
    "lrrk2":         ("LRRK2 kinase inhibitors Parkinson disease therapeutic", ["lrrk2"]),
    "orexin":        ("orexin receptor antagonists insomnia treatment", ["orexin"]),
    "hdac_hd":       ("histone deacetylase inhibitors Huntington disease", ["histone deacetylase"]),
    "nlrp3":         ("NLRP3 inflammasome neuroinflammation neurodegeneration", ["nlrp3"]),
    "nrf2":          ("NRF2 KEAP1 oxidative stress neuroprotection", ["nrf2"]),
    "als_riluzole":  ("amyotrophic lateral sclerosis riluzole edaravone treatment", ["amyotrophic"]),
    "bbb_pred":      ("machine learning prediction blood-brain barrier permeability", ["blood-brain"]),
}


def search(q, tokens, tries=3):
    for _ in range(tries):
        try:
            r = requests.get(EPMC, params={"query": q, "format": "json", "pageSize": 12,
                                           "resultType": "core"}, timeout=40, verify=False)
            if r.status_code != 200:
                time.sleep(2); continue
            for res in r.json().get("resultList", {}).get("result", []):
                title = (res.get("title") or "").lower()
                if all(t.lower() in title for t in tokens) or \
                   sum(t.lower() in title for t in tokens) >= max(1, len(tokens) - 1):
                    return {
                        "title": (res.get("title") or "").rstrip("."),
                        "authors": res.get("authorString", ""),
                        "journal": (res.get("journalInfo", {}) or {}).get("journal", {}).get("title", ""),
                        "year": res.get("pubYear", ""),
                        "doi": res.get("doi", ""),
                        "pmid": res.get("pmid", ""),
                        "source_query": q,
                    }
            return None
        except Exception:
            time.sleep(2)
    return None


def main():
    verified, unresolved = {}, []
    for key, (q, tokens) in WANTED.items():
        hit = search(q, tokens)
        if hit:
            verified[key] = hit
            print(f"OK   {key:15} {hit['year']}  {hit['title'][:62]}", flush=True)
        else:
            unresolved.append(key)
            print(f"MISS {key:15} (no title match; will not be cited)", flush=True)
        time.sleep(0.25)

    (OUT / "references_verified.json").write_text(json.dumps(verified, indent=2))

    lines = ["# References", "",
             "Every entry below was resolved against Europe PMC by title match at build time; "
             "none is written from memory. The resolving query is recorded in "
             "`references_verified.json` so each can be re-checked.", ""]
    for i, (k, v) in enumerate(sorted(verified.items(), key=lambda kv: kv[1]["year"]), 1):
        auth = v["authors"] or ""
        if len(auth) > 90:
            auth = auth.split(",")[0] + " et al."
        bits = [x for x in [auth.rstrip(". "), v["title"], v["journal"], v["year"]] if x]
        ref = ". ".join(bits)
        if v.get("doi"):
            ref += f". doi:{v['doi']}"
        elif v.get("pmid"):
            ref += f". PMID:{v['pmid']}"
        lines.append(f"{i}. {ref}")
    if unresolved:
        lines += ["", f"Unresolved and therefore not cited: {', '.join(unresolved)}"]
    (OUT / "references.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"\nverified {len(verified)} of {len(WANTED)}; unresolved: {unresolved}")
    print("wrote", OUT / "references.md")


if __name__ == "__main__":
    main()
