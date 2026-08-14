"""Resolve [@key] citation tokens into NAR-style numbers, and emit the reference list.

NAR numbers references in order of first citation. Doing that by hand is how bibliographies drift:
a paragraph moves, every number after it is wrong, and nothing complains. Here the manuscript is
written with keys, and numbering is derived at build time from the order the keys first appear.

Only keys present in manuscript/references_verified.json resolve. A key that is not there is a
citation to a work no live query returned, so the build fails rather than emitting a number pointing
at nothing. Keys that are verified but never cited are reported too, since an uncited entry in a
reference list is padding.

Used by build_manuscript.py; runnable alone to check the citation state of the draft.

Run:  python src/brainsafe/analysis/cite.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
MS = ROOT / "manuscript"
VERIFIED = MS / "references_verified.json"
TOKEN = re.compile(r"\[@([a-z0-9_]+)\]")


def load() -> tuple[dict, dict]:
    d = json.loads(VERIFIED.read_text(encoding="utf-8"))
    return d.get("papers", {}), d.get("software", {})


def format_paper(rec: dict) -> str:
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
    return out + "."


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


def reference_section(order: list[str]) -> str:
    papers, software = load()
    lines = ["## References", "",
             "Every entry was resolved by a live query against CrossRef or Europe PMC and accepted "
             "only on a title match, or, where the identity is known and the registered title is a "
             "short form, by resolving the DOI and confirming the title and first author. The "
             "requested title, the matched title and the score are recorded in "
             "`manuscript/references_verified.json`, so the list can be re-checked mechanically. "
             "None is written from memory.", ""]
    for i, key in enumerate(order, 1):
        if key in papers:
            lines.append(f"{i}. {format_paper(papers[key])}")
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


if __name__ == "__main__":
    main()
