"""Check a defence deck for the layout defects that a build script cannot see.

Every one of these checks exists because the fault it looks for was actually shipped in one of these
decks and found only by rendering the slide to an image and looking at it. Rendering is slow and
needs PowerPoint; these checks are fast and need neither, so they run first and catch the same faults
earlier.

  off-slide          a shape whose box leaves the 13.333 by 7.5 inch canvas. Caused by a column
                     positioned from the left when it should have been sized to the remaining width
  invisible text     a run whose colour matches the fill immediately behind it. Caused by using a
                     dark palette colour for a card title on a dark card
  overset text       a text box whose content cannot fit at the requested size. Caused by writing a
                     longer sentence into a box sized for the previous one
  collision          two runs of text that overlap. A card behind text is fine and is not flagged,
                     and neither is a box whose declared height exceeds its content, since these
                     decks size boxes generously on purpose. Only the estimated ink is compared
  no notes           a slide with no speaker notes. Every slide in these decks is meant to be
                     defensible aloud

Overset detection is an estimate, not a layout engine. It assumes an average glyph width of 0.50 em
for the body face and 0.52 for the heading face, which is close enough to catch a box that is 30 per
cent too small and deliberately will not flag one that is 5 per cent too small.

Run:  brainsafe_env/Scripts/python.exe thesis/presentations/validate_deck.py
      brainsafe_env/Scripts/python.exe thesis/presentations/validate_deck.py chapter09_defence.pptx
Exit: 0 if every deck passes, 1 otherwise.
"""
from __future__ import annotations

import sys
from pathlib import Path

from pptx import Presentation
from pptx.util import Emu

HERE = Path(__file__).resolve().parent
EMU_IN = 914400.0
W_IN, H_IN = 13.333, 7.5
TOL = 0.02                 # inches of slack at the slide edge

GLYPH = {"Cambria": 0.52, "Calibri": 0.50}
LINE_FACTOR = 1.22         # line height as a multiple of point size, before explicit line_spacing


def _rgb(color):
    try:
        return tuple(color.rgb)
    except Exception:                                      # theme colour or inherited
        return None


def _fill_rgb(shape):
    try:
        f = shape.fill
        if f.type is not None and f.type == 1:             # MSO_FILL.SOLID
            return _rgb(f.fore_color)
    except Exception:
        pass
    return None


def _bg_rgb(slide):
    try:
        f = slide.background.fill
        if f.type == 1:
            return _rgb(f.fore_color)
    except Exception:
        pass
    return None


def _box(sh):
    return (sh.left / EMU_IN, sh.top / EMU_IN,
            (sh.left + sh.width) / EMU_IN, (sh.top + sh.height) / EMU_IN)


def _overlap(a, b):
    return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])


def _near(c1, c2, tol=48):
    """Are two colours close enough that text in one is unreadable on the other?"""
    if c1 is None or c2 is None:
        return False
    return sum(abs(x - y) for x, y in zip(c1, c2)) < tol


def _estimate_height(sh) -> float:
    """Inches of text the frame needs, summed over paragraphs."""
    tf = sh.text_frame
    width_in = sh.width / EMU_IN
    total = 0.0
    for p in tf.paragraphs:
        runs = [r for r in p.runs if r.text]
        if not runs:
            total += 0.10
            continue
        size = max((r.font.size.pt if r.font.size else 12) for r in runs)
        face = runs[0].font.name or "Calibri"
        chars = sum(len(r.text) for r in runs)
        per_line = max(1, int(width_in * 72 / (size * GLYPH.get(face, 0.50))))
        # An explicit newline forces a break, but the text either side of it is short, so the
        # wrapped estimate and the break count are two views of the same lines rather than two
        # sets to be added. Take the larger.
        breaks = sum(r.text.count("\n") for r in runs)
        lines = max(1, -(-chars // per_line), breaks + 1)
        spacing = float(p.line_spacing) if isinstance(p.line_spacing, float) else 1.0
        after = (p.space_after.pt if p.space_after else 0) / 72.0
        total += lines * size * LINE_FACTOR * spacing / 72.0 + after
    return total


def validate(path: Path) -> list[str]:
    prs = Presentation(str(path))
    problems: list[str] = []
    for i, slide in enumerate(prs.slides, 1):
        bg = _bg_rgb(slide)
        cards = []                                          # (box, fill) of non-text shapes
        texts = []                                          # (box, shape) of shapes carrying text

        for sh in slide.shapes:
            if sh.left is None or sh.top is None:
                continue
            box = _box(sh)
            if (box[0] < -TOL or box[1] < -TOL
                    or box[2] > W_IN + TOL or box[3] > H_IN + TOL):
                problems.append(
                    f"slide {i}: off-slide shape at "
                    f"({box[0]:.2f}, {box[1]:.2f}) to ({box[2]:.2f}, {box[3]:.2f})")
            has_text = sh.has_text_frame and sh.text_frame.text.strip()
            fill = _fill_rgb(sh)
            if has_text:
                texts.append((box, sh))
            if fill is not None:
                cards.append((box, fill))

        for box, sh in texts:
            behind = bg
            for cbox, cfill in cards:
                if cbox != box and _overlap(box, cbox):
                    behind = cfill
            own = _fill_rgb(sh)
            if own is not None:
                behind = own
            for p in sh.text_frame.paragraphs:
                for r in p.runs:
                    if not r.text.strip():
                        continue
                    if _near(_rgb(r.font.color), behind):
                        problems.append(
                            f"slide {i}: text \"{r.text[:38]}\" is the colour of what is behind it")
                        break

            need = _estimate_height(sh)
            have = box[3] - box[1]
            if need > have * 1.30 and need - have > 0.12:
                problems.append(
                    f"slide {i}: text box may overset, needs about {need:.2f} in and has "
                    f"{have:.2f} in: \"{sh.text_frame.text[:44].strip()}\"")

        # Compare the ink, not the container. These decks give a text box more height than its
        # content needs so that a later edit has room, and the box below it starts where the text
        # actually ends. Comparing declared boxes flags every one of those as a collision.
        ink = []
        for box, sh in texts:
            need = min(_estimate_height(sh), box[3] - box[1])
            top = box[1]
            if sh.text_frame.vertical_anchor is not None and str(
                    sh.text_frame.vertical_anchor).startswith("MIDDLE"):
                top = box[1] + max(0.0, (box[3] - box[1] - need) / 2)
            ink.append(((box[0], top, box[2], top + need), sh))

        for a in range(len(ink)):
            for b in range(a + 1, len(ink)):
                box_a, sh_a = ink[a]
                box_b, sh_b = ink[b]
                if not _overlap(box_a, box_b):
                    continue
                ov = (min(box_a[2], box_b[2]) - max(box_a[0], box_b[0])) * \
                     (min(box_a[3], box_b[3]) - max(box_a[1], box_b[1]))
                if ov > 0.06:
                    problems.append(
                        f"slide {i}: text overlaps by {ov:.2f} sq in: "
                        f"\"{sh_a.text_frame.text[:26].strip()}\" and "
                        f"\"{sh_b.text_frame.text[:26].strip()}\"")

        if not slide.has_notes_slide or not slide.notes_slide.notes_text_frame.text.strip():
            problems.append(f"slide {i}: no speaker notes")

    return problems


def main() -> int:
    names = sys.argv[1:] or sorted(p.name for p in HERE.glob("chapter*_defence.pptx"))
    bad = 0
    for name in names:
        path = HERE / name
        problems = validate(path)
        print(f"\n{name}: {len(Presentation(str(path)).slides._sldIdLst)} slides, "
              f"{len(problems)} problem(s)")
        for p in problems:
            print(f"  {p}")
        bad += len(problems)
    print(f"\n{'PASS' if not bad else 'FAIL'}: {bad} problem(s) across {len(names)} deck(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
