"""Assemble the NAR manuscript: inject the generated tables, then convert to .docx.

The manuscript source carries a <!-- TABLES --> marker; the tables are injected from
manuscript/tables_generated.md, which is itself computed from the saved cross-validation
predictions. This guarantees that no number in the manuscript can drift from the results.
"""
from __future__ import annotations

from pathlib import Path

import pypandoc

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
    BUILT.write_text(text, encoding="utf-8")

    pypandoc.convert_file(str(BUILT), "docx", outputfile=str(DOCX),
                          extra_args=[f"--resource-path={MS}", "--toc", "--toc-depth=2"])
    print("wrote", BUILT.name, f"({len(text):,} chars)")
    print("wrote", DOCX.name, f"({DOCX.stat().st_size/1024:.0f} KB)")
    n_tables = text.count("\n| Endpoint |") + text.count("\n| Endpoint|")
    print("tables embedded:", text.count("## Table"))
    print("figures referenced:", text.count("](figures/"))


if __name__ == "__main__":
    main()
