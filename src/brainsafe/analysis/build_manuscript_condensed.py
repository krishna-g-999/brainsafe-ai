"""Assemble the condensed, NAR-length manuscript: resolve citations, then convert to .docx.

NAR's Web Server Issue instructions ask for four to five printed pages, and the full manuscript
(build_manuscript.py) runs well past that once every validation is described in full. This is the
same paper at that length, covering the same claims with the same numbers, for the length-constrained
submission itself. It shares build_manuscript.py's citation machinery (cite.py) rather than a second,
hand-rolled numbering scheme, for the same reason that machinery exists at all: two documents citing
the same works by two different hand-assigned number sequences is exactly the kind of drift this
project's own audit discipline has repeatedly found and fixed elsewhere (most recently in this
manuscript's own reference list, which briefly disagreed with itself the same way).

The condensed draft resolves its own citation order independently of the full draft: a work cited
third here and ninth there gets two different numbers, one per document, which is correct NAR style
and not an inconsistency, since each paper's reference list is numbered by its own order of first
appearance in that paper.

Run:  python src/brainsafe/analysis/build_manuscript_condensed.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pypandoc

import cite

ROOT = Path(__file__).resolve().parents[3]
MS = ROOT / "manuscript"
SRC = MS / "NAR_condensed_draft.md"
BUILT = MS / "NAR_condensed_built.md"
DOCX = MS / "NAR_WebServer_BrainSafe_condensed.docx"
REFERENCE_DOCX = ROOT / "docs" / "nar_manuscript_reference.docx"


def main():
    text = SRC.read_text(encoding="utf-8")

    text, order, unknown, uncited = cite.resolve(text)
    if unknown:
        raise SystemExit("cited but not verified, so no number can be assigned: "
                         + ", ".join(unknown))
    if "<!-- REFERENCES -->" not in text:
        raise SystemExit("marker <!-- REFERENCES --> not found in the condensed draft")
    text = text.replace("<!-- REFERENCES -->", cite.reference_section(order))
    BUILT.write_text(text, encoding="utf-8")
    print(f"citations resolved: {len(order)}")
    if uncited:
        print(f"verified but uncited ({len(uncited)}): {', '.join(sorted(uncited))}")

    pypandoc.convert_file(str(BUILT), "docx", outputfile=str(DOCX),
                          extra_args=[f"--resource-path={MS}",
                                      f"--reference-doc={REFERENCE_DOCX}"])
    print("wrote", BUILT.name, f"({len(text):,} chars)")
    print("wrote", DOCX.name, f"({DOCX.stat().st_size/1024:.0f} KB)")

    body_end = text.find("## Data availability")
    body_start = text.find("## Introduction")
    if body_start != -1 and body_end != -1:
        n_words = len(text[body_start:body_end].split())
        print(f"main text (Introduction through Discussion): {n_words:,} words")


if __name__ == "__main__":
    main()
