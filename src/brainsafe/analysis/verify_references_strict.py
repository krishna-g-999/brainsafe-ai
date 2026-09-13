"""Resolve citations by EXACT TITLE against CrossRef and Europe PMC, with a strict match test.

A first attempt matched on a keyword appearing anywhere in the title and produced confidently wrong
results: a search for the random forest paper returned a 2026 epidemiology article, and a search for
scikit-learn returned an unrelated statistics-package validation. Relevance ranking is not identity.

This version therefore states the exact title of each intended work, queries CrossRef (which covers
the machine-learning venues that Europe PMC does not) and falls back to Europe PMC, and accepts a
record only when the returned title matches the requested one with a normalised similarity of at
least ACCEPT. Anything below that is reported unresolved and is not cited. Software without a
citable paper is recorded as a software entry rather than being attached to an unrelated article.

A live search is not deterministic run to run: a near-duplicate title (a book chapter that paraphrases
a journal article's title, most often) can occasionally rank above the intended work on one run and
below it on the next, both scoring above ACCEPT. An early version of this script overwrote
references_verified.json unconditionally on every run, and once did exactly that: it replaced a
correct 2017 journal match, with a PMID, for "Monoamine Oxidase B Inhibitors in Parkinson's Disease"
with an unrelated 2003 book chapter of similar title, no PMID, and a mojibake-corrupted title, scoring
0.981 against the correct match's 1.0. Nothing caught it until a manuscript rebuild put the garbled
title in a reference list. Existing entries are now kept unless a fresh search returns a strictly
higher match score, so a later run can improve a record but cannot silently replace a better match
with a merely-acceptable one.

This script no longer writes manuscript/references.md. cite.py owns that file, generating it from
this one in the manuscript's own citation order (order of first appearance); this script writing it
too, in year order, is the reason references.md once disagreed with every number the manuscript
actually printed, and this run would revert that fix on every use if it still did.

Output: manuscript/references_verified.json only. Run cite.py afterwards to carry any change into
manuscript/references.md and the built manuscript.
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
    # The B3DB release paper is titled by what it contains, not by the acronym; searching the acronym
    # returns nothing. The abstract of this record names B3DB explicitly.
    "b3db": "A curated diverse molecular database of blood-brain barrier permeability with chemical descriptors",
    "ecfp": "Extended-Connectivity Fingerprints",
    "random_forest": "Random Forests",
    "tdc": "Artificial intelligence foundation for therapeutic science",
    "moleculenet": "MoleculeNet: a benchmark for molecular machine learning",
    "kegg": "KEGG: kyoto encyclopedia of genes and genomes",
    "reactome": "The Reactome Pathway Knowledgebase 2024",
    "iuphar": "The IUPHAR/BPS Guide to PHARMACOLOGY in 2024",
    "bemis_murcko": "The properties of known drugs. 1. Molecular frameworks",
    "dude": "Directory of useful decoys, enhanced (DUD-E): better ligands and decoys for better benchmarking",
    "cns_mpo": "Moving beyond rules: the development of a central nervous system multiparameter optimization (CNS MPO) approach to enable alignment of druglike properties",
    "kpuu": "Unbound Brain-to-Plasma Partition Coefficient, Kp,uu,brain-a Game Changing Parameter for CNS Drug Discovery and Development",
    "xgboost": "XGBoost: A Scalable Tree Boosting System",
    "calibration": "Predicting good probabilities with supervised learning",
    # Platt's 1999 chapter has no registered DOI. The algorithm scikit-learn actually implements is
    # the corrected form from this paper, so this is the citation that describes the running code.
    "platt": "A note on Platt's probabilistic outputs for support vector machines",
    "conformal": "Introducing Conformal Prediction in Predictive Modeling. A Transparent and Flexible Alternative to Applicability Domain Determination",
    "ad_qsar": "QSAR applicability domain estimation by projection of the training set in descriptor space: a review",
    "herg_pred": "Cardiac safety, drug-induced QT prolongation and torsade de pointes (TdP)",
    # No DeLong test is run anywhere in this pipeline. An earlier manuscript draft claimed one; the
    # claim was removed rather than the test added, so the reference goes with it. A bibliography
    # entry for a method that is not used is padding.
    "wilson_ci": "Probable Inference, the Law of Succession, and Statistical Inference",
    "shap_trees": "From local explanations to global understanding with explainable AI for trees",
    "cns_attrition": "Drug metabolism and pharmacokinetics, the blood-brain barrier, and central nervous system drug discovery",
    "ache_ad": "Acetylcholinesterase inhibitors for Alzheimer's disease",
    "bace1_fail": "BACE1 inhibitors: attractive therapeutics for Alzheimer's disease",
    "mao_b_pd": "Monoamine oxidase B inhibitors in Parkinson's disease",
    "lrrk2_pd": "Achieving neuroprotection with LRRK2 kinase inhibitors in Parkinson disease",
    "orexin_insomnia": "Orexin receptor antagonists for the treatment of insomnia",
    "nlrp3_neuro": "Role of the NLRP3 inflammasome in neurodegenerative diseases and therapeutic implications",
    "nrf2_neuro": "The Keap1-Nrf2 pathway in neurodegenerative diseases",
    # No review of HDAC *inhibitors* in Huntington's disease resolved above threshold. This is a
    # primary study, and the manuscript claim is written to say only what it demonstrates: removing
    # histone deacetylases modifies Huntington's pathology in mice.
    "hdac_hd": "Histone deacetylase knockouts modify transcription, CAG instability and nuclear pathology in Huntington disease mice",
    "riluzole_als": "Riluzole for amyotrophic lateral sclerosis (ALS)/motor neuron disease (MND)",
}

# Works whose identity is known by DOI. A title query ranks by relevance, and for conference
# proceedings that ranking is unreliable: the query for the XGBoost paper returned an RFID article
# at 0.53 while the record itself was registered and correct. Resolving the DOI and then checking
# that the returned title matches is a stronger test than searching, not a weaker one, because it
# fixes the record first and asks the registry to confirm it.
#
# key -> (doi, intended title, first author's family name). ACM registered the XGBoost paper under
# the short title "XGBoost", so a similarity test against the full title scores 0.30 on a record
# that is unambiguously right. The test applied to these entries is therefore: the registered title
# must be a leading fragment of the intended one, and the first author must match. Both conditions
# come from the registry, neither from this file.
#
# Two entries also sit here because a title match is not a work match, and both had resolved to the
# wrong work at similarity 1.00. "Extended-Connectivity Fingerprints" returned the IUPAC Gold Book
# terminology entry rather than the method paper that defines the descriptor this server computes,
# and "Random Forests" returned a 2020 textbook chapter rather than the algorithm's source. Naming
# the DOI fixes the work; the registry then confirms the title and the author.
DOIS = {
    "xgboost": ("10.1145/2939672.2939785", "XGBoost: A Scalable Tree Boosting System", "Chen"),
    "ecfp": ("10.1021/ci100050t", "Extended-Connectivity Fingerprints", "Rogers"),
    # TreeExplainer, which is the exact-for-trees estimator actually used, rather than the 2017
    # NeurIPS paper that introduced SHAP in general and has no registered DOI.
    "shap_trees": ("10.1038/s42256-019-0138-9",
                   "From local explanations to global understanding with explainable AI for trees",
                   "Lundberg"),
    "random_forest": ("10.1023/a:1010933404324", "Random Forests", "Breiman"),
}
# Works with no registered DOI, so no title query can resolve them. Each is recorded with the page
# that establishes the citation, rather than being attached to a different paper that happens to
# mention the same words.
SOFTWARE = {
    "rdkit": {"kind": "software", "text": "RDKit: Open-source cheminformatics", "url": "https://www.rdkit.org"},
    "streamlit": {"kind": "software", "text": "Streamlit: an open-source app framework", "url": "https://streamlit.io"},
    "sklearn": {"kind": "software",
                "text": "Pedregosa F, Varoquaux G, Gramfort A et al. Scikit-learn: Machine Learning "
                        "in Python. Journal of Machine Learning Research. 2011;12:2825-2830",
                "url": "https://www.jmlr.org/papers/v12/pedregosa11a.html"},
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


def try_doi(doi, expect_title, expect_first_author):
    """Resolve a DOI and confirm the registry agrees on the title and the first author."""
    try:
        r = requests.get(f"{CROSSREF}/{doi}", timeout=15, headers=HEADERS, verify=False)
        if r.status_code != 200:
            return None
        it = r.json()["message"]
        t = (it.get("title") or [""])[0]
        auth = it.get("author") or []
        names = ", ".join(f"{a.get('family','')} {a.get('given','')[:1]}".strip()
                          for a in auth[:3] if a.get("family"))
        if len(auth) > 3:
            names += " et al."
        yr = ""
        for f in ("published-print", "published-online", "issued", "created"):
            if it.get(f, {}).get("date-parts", [[None]])[0][0]:
                yr = str(it[f]["date-parts"][0][0]); break
        got_first = (auth[0].get("family", "") if auth else "").lower()
        title_ok = norm(t) and norm(expect_title).startswith(norm(t))
        author_ok = got_first == expect_first_author.lower()
        if not (title_ok and author_ok):
            return None
        # The registry confirmed both conditions, so the record is identified; report it as such
        # rather than as a fuzzy title score that would understate a DOI-exact match.
        return (1.0,
                {"title": expect_title, "registered_title": t.rstrip("."), "authors": names,
                 "journal": (it.get("container-title") or [""])[0],
                 "year": yr, "doi": it.get("DOI", ""), "source": "CrossRef (by DOI)"})
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
    out_path = OUT / "references_verified.json"
    existing = json.loads(out_path.read_text(encoding="utf-8")) if out_path.exists() else {}
    existing_papers = existing.get("papers", {})

    verified, unresolved, kept, improved = {}, [], [], []
    for key, title in TITLES.items():
        cands = [c for c in (try_crossref(title), try_epmc(title)) if c]
        if key in DOIS:
            cands = [c for c in (try_doi(*DOIS[key]),) if c] or cands
        best = max(cands, key=lambda c: c[0]) if cands else None
        prior = existing_papers.get(key)
        prior_score = prior.get("match_score", 0) if prior else -1

        if best and best[0] >= ACCEPT:
            if prior and prior_score >= best[0]:
                # A fresh search is not guaranteed to reproduce or improve on an already-accepted
                # match; keep the better (or equal) one already on record rather than replace it
                # with a merely-acceptable result that happened to rank highest this run.
                verified[key] = prior
                kept.append(key)
                print(_safe(f"OK   {key:16} kept existing {prior_score:.3f} "
                            f"(fresh search: {best[0]:.2f})"), flush=True)
            else:
                rec = best[1]; rec["match_score"] = round(best[0], 3); rec["requested_title"] = title
                verified[key] = rec
                if prior:
                    improved.append(key)
                print(_safe(f"OK   {key:16} {best[0]:.2f}  {rec['year']}  {rec['title'][:58]}"),
                      flush=True)
        elif prior:
            # Nothing above threshold this run; do not discard a record that was previously
            # verified, since a transient API issue should not un-cite a work.
            verified[key] = prior
            kept.append(key)
            print(_safe(f"OK   {key:16} kept existing {prior_score:.3f} "
                        f"(no acceptable match this run)"), flush=True)
        else:
            unresolved.append(key)
            got = best[1]["title"][:52] if best else "no result"
            print(_safe(f"MISS {key:16} {(best[0] if best else 0):.2f}  best was: {got}"), flush=True)
        time.sleep(0.05)

    payload = {"papers": verified, "software": SOFTWARE, "unresolved": unresolved,
               "accept_threshold": ACCEPT}
    out_path.write_text(json.dumps(payload, indent=2))

    print(f"\nverified {len(verified)}/{len(TITLES)}; unresolved: {unresolved}")
    if kept:
        print(f"kept the existing match rather than a fresh one for: {', '.join(kept)}")
    if improved:
        print(f"replaced with a strictly better match for: {', '.join(improved)}")
    print("\nrun src/brainsafe/analysis/cite.py next to carry any change into "
          "manuscript/references.md and the built manuscript")


if __name__ == "__main__":
    main()
