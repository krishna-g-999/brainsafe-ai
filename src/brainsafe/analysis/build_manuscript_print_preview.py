"""Build a two-column print preview of the condensed manuscript.

NAR's "4 to 5 printed pages" limit for the Web Server Issue means their own two-column typeset
page, not a page of a single-column Word document. A single-spaced, single-column Word page holds
roughly a third of what a two-column journal page holds, so judging length from
NAR_WebServer_BrainSafe_condensed.docx alone overstates the manuscript by about that factor.

This script re-flows the same built markdown (produced by build_manuscript_condensed.py, which
must run first) into a two-column layout at a print-realistic type size, so the page count it
reports is comparable to what NAR's production system would actually produce. It is a preview for
judging length, not a submission file: NAR still wants a plain single-column .doc/.docx/.pdf at
submission, and typesets the two-column version itself after acceptance.

Run:  python src/brainsafe/analysis/build_manuscript_print_preview.py
      (after build_manuscript_condensed.py)
"""
from __future__ import annotations

import copy
from pathlib import Path

import docx
import pypandoc
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Cm

ROOT = Path(__file__).resolve().parents[3]
MS = ROOT / "manuscript"
BUILT = MS / "NAR_condensed_built.md"
DOCX = MS / "NAR_WebServer_BrainSafe_condensed_print_preview.docx"
REFERENCE_DOCX = ROOT / "docs" / "nar_twocolumn_preview_reference.docx"

COLUMN_WIDTH = Cm(8.85)   # A4 minus 1.6cm margins and a 0.75cm gutter, split two ways


def _make_full_width_until_introduction(path: Path) -> None:
    """Title, key points and abstract span the full page; the body below is two columns.

    Matches how NAR (and every two-column journal) actually typesets the title block, and is
    needed for the page count this script reports to mean anything: a title block squeezed into
    one column reads nothing like the real thing.
    """
    d = docx.Document(path)
    paras = d.paragraphs
    intro_idx = next((i for i, p in enumerate(paras)
                      if p.style.name == "Heading 2" and p.text.strip() == "Introduction"), None)
    if intro_idx is None:
        raise SystemExit("Introduction heading not found; has the section title changed?")
    prev_p = paras[intro_idx - 1]

    body_sectPr = d.element.body.find(qn('w:sectPr'))
    one_col = copy.deepcopy(body_sectPr)
    cols = one_col.find(qn('w:cols'))
    if cols is not None:
        one_col.remove(cols)
    cols = OxmlElement('w:cols')
    cols.set(qn('w:num'), '1')
    one_col.append(cols)
    sect_type = one_col.find(qn('w:type')) or OxmlElement('w:type')
    if sect_type.getparent() is None:
        one_col.insert(0, sect_type)
    sect_type.set(qn('w:val'), 'continuous')
    prev_p._p.get_or_add_pPr().append(one_col)

    body_type = body_sectPr.find(qn('w:type')) or OxmlElement('w:type')
    if body_type.getparent() is None:
        body_sectPr.insert(0, body_type)
    body_type.set(qn('w:val'), 'continuous')
    d.save(path)


def _shrink_images_to_column_width(path: Path) -> None:
    """Figures sized for a single full-width Word page overflow a print column and get clipped.

    Scaling every inline image to the column width, preserving aspect ratio, is a preview
    compromise (real NAR figures are drawn at column size from the start) but keeps the page
    count this script reports honest rather than corrupted by a rendering artefact.
    """
    d = docx.Document(path)
    for shape in d.inline_shapes:
        if shape.width > COLUMN_WIDTH:
            ratio = COLUMN_WIDTH / shape.width
            shape.height = int(shape.height * ratio)
            shape.width = COLUMN_WIDTH
    d.save(path)


def main() -> None:
    if not BUILT.exists():
        raise SystemExit(f"{BUILT.name} not found; run build_manuscript_condensed.py first")
    pypandoc.convert_file(str(BUILT), "docx", outputfile=str(DOCX),
                          extra_args=[f"--resource-path={MS}",
                                      f"--reference-doc={REFERENCE_DOCX}"])
    _make_full_width_until_introduction(DOCX)
    _shrink_images_to_column_width(DOCX)
    print("wrote", DOCX.name, f"({DOCX.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
