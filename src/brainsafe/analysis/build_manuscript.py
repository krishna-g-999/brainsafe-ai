"""Assemble the NAR manuscript: inject the generated tables, then convert to .docx.

The manuscript source carries a <!-- TABLES --> marker; the tables are injected from
manuscript/tables_generated.md, which is itself computed from the saved cross-validation
predictions. This guarantees that no number in the manuscript can drift from the results.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pypandoc

import cite

ROOT = Path(__file__).resolve().parents[3]
MS = ROOT / "manuscript"
SRC = MS / "NAR_WebServer_BrainSafe_draft.md"
TABLES = MS / "tables_generated.md"
BUILT = MS / "NAR_WebServer_BrainSafe_built.md"
DOCX = MS / "NAR_WebServer_BrainSafe.docx"


def main():
    text = SRC.read_text(encoding="utf-8")
    tbl = TABLES.read_text(encoding="utf-8")
    # drop the generated file's own H1 header, keep the tables
    tbl = "\n".join(l for l in tbl.splitlines() if not l.startswith("# "))
    if "<!-- TABLES -->" not in text:
        raise SystemExit("marker <!-- TABLES --> not found in manuscript source")
    text = text.replace("<!-- TABLES -->", tbl.strip())

    # Citations are keys in the source and numbers only in the built file, so moving a paragraph
    # renumbers the bibliography instead of silently invalidating it.
    text, order, unknown, uncited = cite.resolve(text)
    if unknown:
        raise SystemExit("cited but not verified, so no number can be assigned: "
                         + ", ".join(unknown))
    if "<!-- REFERENCES -->" not in text:
        raise SystemExit("marker <!-- REFERENCES --> not found in manuscript source")
    text = text.replace("<!-- REFERENCES -->", cite.reference_section(order))
    BUILT.write_text(text, encoding="utf-8")
    print(f"citations resolved: {len(order)}")
    if uncited:
        print(f"verified but uncited ({len(uncited)}): {', '.join(sorted(uncited))}")

    pypandoc.convert_file(str(BUILT), "docx", outputfile=str(DOCX),
                          extra_args=[f"--resource-path={MS}", "--toc", "--toc-depth=2"])
    print("wrote", BUILT.name, f"({len(text):,} chars)")
    print("wrote", DOCX.name, f"({DOCX.stat().st_size/1024:.0f} KB)")
    n_tables = text.count("\n| Endpoint |") + text.count("\n| Endpoint|")
    print("tables embedded:", text.count("## Table"))
    print("figures referenced:", text.count("](figures/"))


if __name__ == "__main__":
    main()
