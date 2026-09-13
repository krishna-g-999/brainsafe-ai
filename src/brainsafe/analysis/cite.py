"""Resolve [@key] citation tokens into NAR-style numbers, and emit the reference list.

NAR numbers references in order of first citation. Doing that by hand is how bibliographies drift:
a paragraph moves, every number after it is wrong, and nothing complains. Here the manuscript is
written with keys, and numbering is derived at build time from the order the keys first appear.

Only keys present in manuscript/references_verified.json resolve. A key that is not there is a
citation to a work no live query returned, so the build fails rather than emitting a number pointing
at nothing. Keys that are verified but never cited are reported too, since an uncited entry in a
reference list is padding.

NAR's instructions for a Web Server paper are explicit: "The references section must include active
electronic DOI and PubMed Central links for each cited paper (where available). Please include a
PubMed abstract link if the PMC link is not available." references_links.json, built by
add_reference_links.py against the NCBI ID Converter and PubMed esearch, carries that data by DOI.
This was previously wired only into a standalone copy of the reference list, manuscript/references.md,
which add_reference_links.py rewrote directly and which was numbered by year rather than by citation
order: the exact bug this module's own docstring warns a hand-numbered list invites. That copy could
silently drift from the one actually built into the manuscript, and did. Loading the links here, in
the same function that produces the manuscript's own reference list, is what makes the two the same
list rather than two lists that happen to start out matching.

Used by build_manuscript.py; runnable alone to check the citation state of the draft, and to
regenerate manuscript/references.md as a byproduct so it cannot re-diverge from the built manuscript.

Run:  python src/brainsafe/analysis/cite.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
MS = ROOT / "manuscript"
VERIFIED = MS / "references_verified.json"
LINKS = MS / "references_links.json"
TOKEN = re.compile(r"\[@([a-z0-9_]+)\]")


def load() -> tuple[dict, dict]:
    d = json.loads(VERIFIED.read_text(encoding="utf-8"))
    return d.get("papers", {}), d.get("software", {})


def load_links() -> dict:
    if not LINKS.exists():
        return {}
    return {k.lower(): v for k, v in json.loads(LINKS.read_text(encoding="utf-8")).items()}


def format_paper(rec: dict, links: dict | None = None) -> str:
    # "et al." already ends in a full stop, so the separator must not add a second one
    authors = rec.get("authors", "").strip().rstrip(".")
    bits = [authors, rec.get("title", "").strip().rstrip(".")]
    tail = ". ".join(b for b in bits if b)
    journal, year, doi = rec.get("journal", ""), rec.get("year", ""), rec.get("doi", "")
    out = tail
    if journal:
        out += f". {journal}"
    if year:
        out += f". {year}"
    if doi:
        out += f". doi:{doi}"
    out += "."
    # PMC link if available; otherwise a PubMed abstract link if available; neither if the work is
    # not indexed there (expected for conference proceedings and some society journals).
    link = (links or {}).get((doi or "").rstrip(". ").lower())
    if link:
        bits = []
        if link.get("pmcid"):
            bits.append(f"https://pmc.ncbi.nlm.nih.gov/articles/{link['pmcid']}/")
        if link.get("pmid"):
            bits.append(f"https://pubmed.ncbi.nlm.nih.gov/{link['pmid']}/")
        if bits:
            out += " " + " ".join(bits)
    return out


def format_software(rec: dict) -> str:
    text = rec.get("text", "").strip().rstrip(".")
    url = rec.get("url", "")
    return f"{text}. {url}" if url else f"{text}."


def resolve(text: str) -> tuple[str, list[str], list[str], list[str]]:
    """Replace tokens with numbers assigned by first appearance.

    Returns the rewritten text, the ordered keys, the keys that could not be resolved, and the
    verified keys that were never cited.
    """
    papers, software = load()
    known = {**papers, **software}
    order: list[str] = []
    unknown: list[str] = []

    for m in TOKEN.finditer(text):
        key = m.group(1)
        if key not in known:
            if key not in unknown:
                unknown.append(key)
            continue
        if key not in order:
            order.append(key)

    number = {k: i + 1 for i, k in enumerate(order)}

    def sub(m):
        key = m.group(1)
        return f"({number[key]})" if key in number else f"[@{key}: UNRESOLVED]"

    uncited = [k for k in known if k not in order]
    return TOKEN.sub(sub, text), order, unknown, uncited


def reference_section(order: list[str], header: str = "## References") -> str:
    papers, software = load()
    links = load_links()
    lines = [header, "",
             "Every entry was resolved by a live query against CrossRef or Europe PMC and accepted "
             "only on a title match, or, where the identity is known and the registered title is a "
             "short form, by resolving the DOI and confirming the title and first author. The "
             "requested title, the matched title and the score are recorded in "
             "`manuscript/references_verified.json`, so the list can be re-checked mechanically. "
             "None is written from memory. A PubMed Central or PubMed abstract link is given where "
             "NCBI indexes the work (`manuscript/references_links.json`); neither exists for a "
             "work outside PubMed's coverage, which is expected for some conference proceedings and "
             "for the software citations.", ""]
    for i, key in enumerate(order, 1):
        if key in papers:
            lines.append(f"{i}. {format_paper(papers[key], links)}")
        else:
            lines.append(f"{i}. {format_software(software[key])}")
    return "\n".join(lines)


def main() -> None:
    src = MS / "NAR_WebServer_BrainSafe_draft.md"
    text = src.read_text(encoding="utf-8")
    _out, order, unknown, uncited = resolve(text)
    print(f"{len(order)} citations resolved, in order of first appearance")
    for i, k in enumerate(order, 1):
        print(f"  {i:2d}  {k}")
    if unknown:
        print(f"\nUNRESOLVED ({len(unknown)}): {', '.join(unknown)}")
        print("  these are cited but not in references_verified.json; the build will fail")
    if uncited:
        print(f"\nverified but never cited ({len(uncited)}): {', '.join(sorted(uncited))}")

    # manuscript/references.md used to be a second, hand-maintained copy of this list, numbered by
    # year rather than by citation order, and add_reference_links.py wrote PubMed/PMC links into
    # that copy alone. The two could disagree, and did: every in-text number pointed at the wrong
    # entry in that file. Writing it here, from the same order and the same formatter the built
    # manuscript uses, makes that impossible rather than merely unlikely.
    if not unknown:
        out_md = MS / "references.md"
        out_md.write_text(reference_section(order, header="# References") + "\n", encoding="utf-8")
        print(f"\nwrote {out_md.relative_to(ROOT)} ({len(order)} entries, matching citation order)")


if __name__ == "__main__":
    main()
