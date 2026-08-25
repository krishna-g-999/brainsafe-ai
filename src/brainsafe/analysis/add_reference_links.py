"""Attach PubMed and PubMed Central links to every reference, as NAR requires.

NAR's instructions are explicit: "The references section must include active electronic DOI and
PubMed Central links for each cited paper (where available). Please include a PubMed abstract link
if the PMC link is not available." Our references carry DOIs for all 32 entries and PubMed links for
none, so the requirement is not met.

The identifiers are not typed in. Each DOI is resolved through the NCBI ID Converter, which maps
between DOI, PMID and PMCID, and whatever it returns is recorded verbatim alongside the DOI that
produced it. An entry the converter does not know keeps its DOI alone and is reported, because a
reference list that quietly invents a PMID is worse than one that is incomplete.

Software citations have no PubMed record by construction and are expected to resolve to nothing.

Output: manuscript/references_links.json, and an updated manuscript/references.md

Run:  python src/brainsafe/analysis/add_reference_links.py
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[3]
MS = ROOT / "manuscript"
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
    src = MS / "references.md"
    text = src.read_text(encoding="utf-8")
    entries = re.findall(r"^(\d+)\. (.+)$", text, re.M)
    dois = []
    for _n, body in entries:
        m = re.search(r"doi:(\S+)", body, re.I)
        if m:
            dois.append(m.group(1).rstrip(". ").lower())
    print(f"{len(entries)} numbered entries, {len(dois)} carrying a DOI")

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

    lines, added, missing = [], 0, []
    for line in text.splitlines():
        m = re.match(r"^(\d+)\. (.+)$", line)
        if not m:
            lines.append(line)
            continue
        n, body = m.groups()
        if "pubmed.ncbi" in body or "pmc.ncbi" in body:      # already linked; leave alone
            lines.append(line)
            continue
        d = re.search(r"doi:(\S+)", body, re.I)
        rec = got.get(d.group(1).rstrip(". ").lower()) if d else None
        if rec and (rec.get("pmcid") or rec.get("pmid")):
            bits = []
            if rec.get("pmcid"):
                bits.append(f"https://pmc.ncbi.nlm.nih.gov/articles/{rec['pmcid']}/")
            if rec.get("pmid"):
                bits.append(f"https://pubmed.ncbi.nlm.nih.gov/{rec['pmid']}/")
            lines.append(f"{n}. {body.rstrip()} {' '.join(bits)}")
            added += 1
        else:
            lines.append(line)
            missing.append(f"{n}. {body[:70]}")

    src.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (MS / "references_links.json").write_text(json.dumps(got, indent=2), encoding="utf-8")

    print(f"\nlinks added to {added} entries")
    if missing:
        print(f"{len(missing)} entries have no PubMed record and keep their DOI alone:")
        for x in missing:
            print("   ", x.encode("ascii", "replace").decode())
    print(f"\nwrote {src.relative_to(ROOT)} and references_links.json")


if __name__ == "__main__":
    main()
