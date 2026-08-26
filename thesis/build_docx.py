"""Convert every thesis chapter to Word, using the project's existing reference styling.

Chapters are written in Markdown so that they stay diffable and greppable, and so that a number in a
chapter can be checked against an artefact by the same tools that check the manuscript. Word is the
delivery format, not the source: edit the .md, re-run this, and the .docx follows.

Run:  brainsafe_env/Scripts/python.exe thesis/build_docx.py
      brainsafe_env/Scripts/python.exe thesis/build_docx.py --chapter 02
Out:  thesis/docx/chapterNN_*.docx  and  thesis/docx/BrainSafe_thesis.docx (all chapters, one file)
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pypandoc

ROOT = Path(__file__).resolve().parents[1]
THESIS = ROOT / "thesis"
OUT = THESIS / "docx"
REFERENCE = ROOT / "docs" / "report_reference.docx"


def chapters() -> list[Path]:
    return sorted(p for p in THESIS.glob("chapter*.md"))


def strip_session_notes(text: str) -> str:
    """Remove the provenance block-quote at the top and the outstanding-items section at the foot.

    Both exist for the supervisor and neither belongs in a submitted chapter. They are kept in the
    Markdown, which is the working copy, and dropped from the Word file, which is the deliverable.
    """
    lines = text.splitlines()
    out, i = [], 0
    # drop a leading block quote that follows the H1
    while i < len(lines):
        if lines[i].startswith("> ") or (out and out[-1].startswith("> ") and not lines[i].strip()):
            i += 1
            continue
        out.append(lines[i])
        i += 1
    text = "\n".join(out)
    text = re.split(r"\n---\n+## Outstanding items", text)[0]
    return text.rstrip() + "\n"


def convert(src: Path, dest: Path, keep_notes: bool) -> None:
    text = src.read_text(encoding="utf-8")
    if not keep_notes:
        text = strip_session_notes(text)
    args = ["--from", "markdown+tex_math_dollars+pipe_tables", "--to", "docx",
            f"--resource-path={ROOT}"]
    if REFERENCE.exists():
        args += [f"--reference-doc={REFERENCE}"]
    pypandoc.convert_text(text, "docx", format="markdown+tex_math_dollars+pipe_tables",
                          outputfile=str(dest),
                          extra_args=[a for a in args if not a.startswith("--from")
                                      and not a.startswith("--to") and a not in ("markdown+tex_math_dollars+pipe_tables", "docx")])
    print(f"  {dest.relative_to(ROOT).as_posix():<52} {dest.stat().st_size / 1024:6.0f} KB")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Build Word versions of the thesis chapters.")
    ap.add_argument("--chapter", help="build one chapter, by its number, e.g. 02")
    ap.add_argument("--keep-notes", action="store_true",
                    help="keep the provenance header and outstanding-items footer")
    args = ap.parse_args(argv)

    chs = chapters()
    if args.chapter:
        chs = [c for c in chs if c.name.startswith(f"chapter{args.chapter}")]
        if not chs:
            raise SystemExit(f"no chapter matching {args.chapter!r} in {THESIS}")
    if not chs:
        raise SystemExit(f"no chapter*.md found in {THESIS}")

    OUT.mkdir(parents=True, exist_ok=True)
    print(f"reference styling: "
          f"{REFERENCE.relative_to(ROOT).as_posix() if REFERENCE.exists() else 'pandoc default'}")
    print(f"{len(chs)} chapter(s):")
    for c in chs:
        convert(c, OUT / f"{c.stem}.docx", args.keep_notes)

    # the combined document, only when building everything
    if not args.chapter:
        merged = "\n\n\\newpage\n\n".join(
            strip_session_notes(c.read_text(encoding="utf-8")) for c in chapters())
        dest = OUT / "BrainSafe_thesis.docx"
        extra = [f"--resource-path={ROOT}", "--toc", "--toc-depth=2"]
        if REFERENCE.exists():
            extra.append(f"--reference-doc={REFERENCE}")
        pypandoc.convert_text(merged, "docx", format="markdown+tex_math_dollars+pipe_tables",
                              outputfile=str(dest), extra_args=extra)
        print(f"  {dest.relative_to(ROOT).as_posix():<52} {dest.stat().st_size / 1024:6.0f} KB"
              f"   ({len(chapters())} chapters, with contents)")


if __name__ == "__main__":
    main()
