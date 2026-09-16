"""Build a reading copy of the manuscript with the graphical abstract embedded in the body.

NAR's Web Server Issue asks for the graphical abstract as a separate submission file, not as a
figure inside the manuscript text, and NAR_WebServer_BrainSafe_built.md and its .docx stay exactly
that way so the submission itself is not second-guessed. This script produces a companion copy, for
circulation and review rather than submission, with the same graphical abstract shown inline right
after the title block, since a reader following a link or opening a shared file benefits from seeing
it without opening a second attachment.

It runs after build_manuscript.py and does not duplicate its citation or table logic: it reads the
built markdown that script already produced and inserts one image line.

Output: manuscript/NAR_WebServer_BrainSafe_with_figure.docx

Run:  python src/brainsafe/analysis/build_manuscript_with_figure.py
"""
from __future__ import annotations

from pathlib import Path

import pypandoc

ROOT = Path(__file__).resolve().parents[3]
MS = ROOT / "manuscript"
BUILT = MS / "NAR_WebServer_BrainSafe_built.md"
DOCX = MS / "NAR_WebServer_BrainSafe_with_figure.docx"
REFERENCE_DOCX = ROOT / "docs" / "nar_manuscript_reference.docx"
MARKER = "**Graphical abstract:** submitted as a separate file, `manuscript/figures/GraphicalAbstract.png`\n(and `.pdf`)."


def main() -> None:
    if not BUILT.exists():
        raise SystemExit(f"{BUILT.name} not found; run build_manuscript.py first")
    text = BUILT.read_text(encoding="utf-8")
    if MARKER not in text:
        raise SystemExit("graphical-abstract marker not found in the built manuscript; "
                         "has its wording changed?")
    text = text.replace(
        MARKER,
        MARKER + "\n\n![Graphical abstract](figures/GraphicalAbstract.png)\n\n"
        "*Reading copy only: the submission itself carries the graphical abstract as the separate "
        "file NAR requires, not inline. This copy embeds it here for circulation.*")

    tmp = MS / "_with_figure_tmp.md"
    tmp.write_text(text, encoding="utf-8")
    try:
        pypandoc.convert_file(str(tmp), "docx", outputfile=str(DOCX),
                              extra_args=[f"--resource-path={MS}",
                                          f"--reference-doc={REFERENCE_DOCX}"])
    finally:
        tmp.unlink(missing_ok=True)
    print("wrote", DOCX.name, f"({DOCX.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
