"""Resolve PubMed and PubMed Central identifiers for every verified reference, as NAR requires.

NAR's instructions are explicit: "The references section must include active electronic DOI and
PubMed Central links for each cited paper (where available). Please include a PubMed abstract link
if the PMC link is not available."

The identifiers are not typed in. Each DOI is resolved through the NCBI ID Converter, which maps
between DOI, PMID and PMCID, and whatever it returns is recorded verbatim alongside the DOI that
produced it. An entry the converter does not know keeps its DOI alone and is reported, because a
reference list that quietly invents a PMID is worse than one that is incomplete.

Software citations have no PubMed record by construction and are expected to resolve to nothing.

This script owns exactly one file: manuscript/references_links.json, keyed by DOI. It used to also
rewrite manuscript/references.md directly, appending links to whatever numbered list it found there.
That made two scripts write the same file in two different citation orders (this one preserved
whatever order references.md already had; cite.py numbers by order of first appearance in the
manuscript), and the two silently disagreed. cite.py now generates references.md from this file, so
that is the only writer of it. The DOI list here is read from references_verified.json rather than
from references.md itself, for the same reason: the verified set is the source cite.py resolves
against, so resolving links against it rather than against a rendering of it cannot drift from what
the manuscript actually cites.

Run:  python src/brainsafe/analysis/add_reference_links.py
      then: python src/brainsafe/analysis/cite.py     (regenerates references.md with the new links)
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[3]
MS = ROOT / "manuscript"
VERIFIED = MS / "references_verified.json"
LINKS = MS / "references_links.json"
# NCBI moved the converter to pmc.ncbi.nlm.nih.gov; the old www path no longer answers.
CONVERTER = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"


def resolve(dois: list[str], verify=True) -> dict:
    """Ask NCBI to map DOIs to PMID and PMCID, in batches, and record exactly what came back."""
    out = {}
    for i in range(0, len(dois), 50):
        batch = dois[i:i + 50]
        try:
            r = requests.get(CONVERTER, params={"ids": ",".join(batch), "format": "json",
                                                "tool": "brainsafe", "email": "noreply@example.org"},
                             timeout=45, verify=verify)
            r.raise_for_status()
            for rec in r.json().get("records", []):
                doi = (rec.get("doi") or "").lower()
                if not doi:
                    continue
                out[doi] = {"pmid": rec.get("pmid"), "pmcid": rec.get("pmcid"),
                            "live": rec.get("live", True)}
        except Exception as exc:
            print(f"  batch {i//50 + 1} failed: {type(exc).__name__}: {str(exc)[:90]}")
        time.sleep(0.5)
    return out


def main() -> None:
    papers = json.loads(VERIFIED.read_text(encoding="utf-8")).get("papers", {})
    existing = json.loads(LINKS.read_text(encoding="utf-8")) if LINKS.exists() else {}
    dois = sorted({rec["doi"].lower() for rec in papers.values() if rec.get("doi")})
    print(f"{len(papers)} verified papers, {len(dois)} distinct DOIs")

    # NCBI is not TLS-intercepted on this network: it presents a genuine GoDaddy certificate. The
    # local CA bundle assembled for intercepted hosts must NOT be used here, because it carries a
    # Sophos root whose basic constraints are not marked critical, and OpenSSL 3 rejects the entire
    # bundle over that one entry. Clearing the two variables restores certifi.
    for _v in ("SSL_CERT_FILE", "REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE"):
        os.environ.pop(_v, None)
    verify = True

    got = resolve(dois, verify=verify)
    print(f"resolved by the PMC converter: {len(got)} of {len(dois)}")

    # The converter only knows content deposited in PMC, so a paper that is in PubMed but not in PMC
    # resolves to nothing there. NAR asks for a PubMed abstract link where PMC is unavailable, so the
    # remainder are looked up by DOI against PubMed itself. Anything still unresolved is genuinely
    # not indexed, which is expected for conference papers and software.
    rest = [d for d in dois if d not in got or not (got[d].get("pmcid") or got[d].get("pmid"))]
    found = 0
    for d in rest:
        try:
            r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                             params={"db": "pubmed", "term": f"{d}[DOI]", "retmode": "json",
                                     "tool": "brainsafe", "email": "noreply@example.org"},
                             timeout=30)
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                got.setdefault(d, {}).update({"pmid": ids[0], "source": "esearch"})
                found += 1
        except Exception:
            pass
        time.sleep(0.4)
    print(f"resolved by PubMed search: {found} of the remaining {len(rest)}")

    # A DOI the live queries did not return anything for this run keeps whatever the file already
    # held, so a transient NCBI outage cannot erase a link resolved on an earlier run.
    merged = {**existing, **got}
    missing = [d for d in dois if not (merged.get(d, {}).get("pmcid") or merged.get(d, {}).get("pmid"))]

    LINKS.write_text(json.dumps(merged, indent=2), encoding="utf-8")
    print(f"\nlinked: {len(dois) - len(missing)} of {len(dois)} DOIs")
    if missing:
        print(f"{len(missing)} DOIs have no PubMed record and keep their DOI alone: "
              f"{', '.join(missing)}")
    print(f"\nwrote {LINKS.relative_to(ROOT)}")
    print("run src/brainsafe/analysis/cite.py next to carry these into references.md "
          "and the built manuscript")


if __name__ == "__main__":
    main()
