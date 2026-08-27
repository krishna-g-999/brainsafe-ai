"""Shared chrome for the thesis defence decks.

The palette, the type scale, the numbered-disc motif and the title auto-sizing live here so that the
chapter decks are one visual system rather than several that resemble each other. A change to the
palette or to how a long title is fitted reaches every deck that imports this.

Colour is used as a code and not as decoration:
    DEEP     structure, the default for headings and data
    TEAL     the secondary voice, and anything about uncertainty
    AMBER    caveats, refutations, and results that do not flatter the system
    CRIMSON  withdrawal and outright failure
Nothing else should be introduced without a reason of the same kind.
"""
from __future__ import annotations

import csv
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results" / "tables"
INV = ROOT / "inversion" / "results"
PKG = ROOT / "submission_package"
FIG = ROOT / "manuscript" / "figures"

INK = RGBColor(0x0F, 0x1A, 0x2E)
DEEP = RGBColor(0x1F, 0x3A, 0x5F)
TEAL = RGBColor(0x1C, 0x72, 0x93)
AMBER = RGBColor(0xD9, 0x83, 0x24)
CRIMSON = RGBColor(0x9E, 0x2A, 0x2B)
MUTED = RGBColor(0x5F, 0x6B, 0x7A)
TINT = RGBColor(0xEE, 0xF2, 0xF7)
PAPER = RGBColor(0xFF, 0xFF, 0xFF)
CHALK = RGBColor(0xD7, 0xDF, 0xEA)
GRID = RGBColor(0xE2, 0xE8, 0xF0)

HEAD = "Cambria"
BODY = "Calibri"

W, H = 13.333, 7.5
M = 0.62


def rows(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def fmt(v, k=3):
    return f"{v:.{k}f}"


class Deck:
    """A 16:9 deck with the shared chrome. Every method returns the shape it made."""

    grid = GRID

    def __init__(self):
        self.prs = Presentation()
        self.prs.slide_width, self.prs.slide_height = Inches(W), Inches(H)
        self._blank = self.prs.slide_layouts[6]

    # -- slides ------------------------------------------------------------
    def _slide(self, dark: bool):
        s = self.prs.slides.add_slide(self._blank)
        f = s.background.fill
        f.solid()
        f.fore_color.rgb = INK if dark else PAPER
        return s

    def dark(self):
        return self._slide(True)

    def light(self):
        return self._slide(False)

    # -- primitives --------------------------------------------------------
    def text(self, s, x, y, w, h, runs, size=15, color=None, font=BODY, bold=False,
             italic=False, align=PP_ALIGN.LEFT, space=6, line=None, anchor=MSO_ANCHOR.TOP):
        tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
        tf = tb.text_frame
        tf.word_wrap = True
        tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
        tf.vertical_anchor = anchor
        items = runs if isinstance(runs, list) else [runs]
        for i, item in enumerate(items):
            p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
            p.alignment = align
            p.space_after = Pt(space)
            if line:
                p.line_spacing = line
            body, over = (item if isinstance(item, tuple) else (item, {}))
            r = p.add_run()
            r.text = body
            fo = r.font
            fo.name = over.get("font", font)
            fo.size = Pt(over.get("size", size))
            fo.bold = over.get("bold", bold)
            fo.italic = over.get("italic", italic)
            fo.color.rgb = over.get("color", color or INK)
        return tb

    def card(self, s, x, y, w, h, fill=TINT):
        sh = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y),
                                Inches(w), Inches(h))
        sh.adjustments[0] = 0.05
        sh.fill.solid()
        sh.fill.fore_color.rgb = fill
        sh.line.fill.background()
        sh.shadow.inherit = False
        sh.text_frame.text = ""
        return sh

    def dot(self, s, x, y, label, fill=TEAL, fg=PAPER, dia=0.42):
        sh = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(x), Inches(y), Inches(dia), Inches(dia))
        sh.fill.solid()
        sh.fill.fore_color.rgb = fill
        sh.line.fill.background()
        sh.shadow.inherit = False
        tf = sh.text_frame
        tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
        tf.vertical_anchor = MSO_ANCHOR.MIDDLE
        p = tf.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        r = p.add_run()
        r.text = label
        r.font.name = BODY
        r.font.size = Pt(15 if dia >= 0.40 else 13)
        r.font.bold = True
        r.font.color.rgb = fg
        return sh

    def head(self, s, num_label, title, sub=None):
        """Auto-sized so a long title stays on one line instead of colliding with the subtitle."""
        self.dot(s, M, 0.52, num_label)
        ln = len(title)
        size = 30 if ln <= 50 else (27 if ln <= 58 else 24)
        self.text(s, M + 0.62, 0.40, W - 2 * M - 0.62, 0.74, title, size=size, bold=True,
                  font=HEAD, color=INK, anchor=MSO_ANCHOR.MIDDLE)
        if sub:
            self.text(s, M + 0.62, 1.16, W - 2 * M - 0.62, 0.42, sub, size=14, color=MUTED,
                      italic=True)

    def source(self, s, note):
        self.text(s, M, H - 0.52, W - 2 * M, 0.32, note, size=10, color=MUTED)

    def notes(self, s, txt):
        s.notes_slide.notes_text_frame.text = txt

    def stat(self, s, x, y, w, value, label, color=TEAL, vsize=40):
        if len(value) > 5:
            vsize = vsize * (5.0 / len(value)) ** 0.55
        self.text(s, x, y, w, 0.66, value, size=vsize, bold=True, font=HEAD, color=color)
        self.text(s, x, y + 0.70, w, 0.72, label, size=12, color=MUTED, line=1.15)

    def equation(self, s, x, y, w, parts, size=27, color=DEEP):
        """parts: list of (text, is_subscript), rendered as true subscripts."""
        tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(0.8))
        tf = tb.text_frame
        tf.word_wrap = True
        tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
        p = tf.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        for body, is_sub in parts:
            r = p.add_run()
            r.text = body
            r.font.name = HEAD
            r.font.bold = True
            r.font.size = Pt(size * 0.72 if is_sub else size)
            r.font.color.rgb = color
            if is_sub:
                r.font._rPr.set("baseline", "-25000")
        return tb

    def save(self, out: Path):
        self.prs.save(out)
        n = len(self.prs.slides._sldIdLst)
        print(f"wrote {out.relative_to(ROOT).as_posix()}  ({n} slides)")
        return n
