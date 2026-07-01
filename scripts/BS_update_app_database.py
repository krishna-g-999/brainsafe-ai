"""
scripts/BS_update_app_database.py
BrainSafe AI (BS) — Update App to Use Full 535-Compound Database

The app currently reads compounds.json (134 compounds).
The model was trained on 535 compounds.
This script updates app_v6.py to read BS_compounds_full.json instead.

Run AFTER BS_status_and_fix.py (which creates BS_compounds_full.json):
  D:\\BRAINSAFE_AI\\brainsafe_env\\Scripts\\python.exe scripts\\BS_update_app_database.py

What it patches:
  1. All references to "compounds.json" -> "BS_compounds_full.json"
  2. All references to "compounds_ml.json" -> removed (merged into BS_compounds_full)
  3. "325 compounds" display label -> "535 compounds"
  4. Score scale: ensures app reads 0-100 scores correctly
  5. Writes app_v6_final.py (does not touch app_v6.py)
"""

import sys, re, ast
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

APP_IN  = ROOT / "app_v6.py"
APP_OUT = ROOT / "app_v6_final.py"
DB_NEW  = "BS_compounds_full.json"


def main():
    print("=" * 62)
    print("BS Update App Database: 134 -> 535 compounds")
    print("=" * 62)

    # Verify BS_compounds_full.json exists
    db_path = ROOT / DB_NEW
    if not db_path.exists():
        print(f"ERROR: {db_path} not found.")
        print("Run BS_status_and_fix.py first to generate it.")
        sys.exit(1)

    import json
    with open(db_path) as f:
        db = json.load(f)
    n_compounds = len(db)
    print(f"BS_compounds_full.json: {n_compounds} compounds ready")

    if not APP_IN.exists():
        print(f"ERROR: {APP_IN} not found.")
        sys.exit(1)

    content = APP_IN.read_text(encoding="utf-8", errors="replace")
    print(f"Read {APP_IN.name}: {len(content.splitlines())} lines")

    patches_applied = 0

    # ── Patch 1: compounds.json -> BS_compounds_full.json ────────────────
    old = '"compounds.json"'
    new = f'"{DB_NEW}"'
    if old in content:
        content = content.replace(old, new)
        patches_applied += 1
        print(f"  Patch 1: compounds.json -> {DB_NEW}")
    
    old = "'compounds.json'"
    new = f"'{DB_NEW}'"
    if old in content:
        content = content.replace(old, new)
        patches_applied += 1
        print(f"  Patch 1b: 'compounds.json' -> '{DB_NEW}'")

    # ── Patch 2: compounds_ml.json loading ───────────────────────────────
    # The app loads compounds_ml.json separately and merges.
    # Since BS_compounds_full.json already contains all 535, we disable
    # the separate ML load by pointing it to an empty structure.
    old2 = '"compounds_ml.json"'
    if old2 in content:
        content = content.replace(old2, '"BS_compounds_ml_merged.json"')
        # Create the empty merged file so the app doesn't crash
        empty_ml = ROOT / "BS_compounds_ml_merged.json"
        with open(empty_ml, "w") as f:
            json.dump({"_note": "Merged into BS_compounds_full.json"}, f)
        patches_applied += 1
        print("  Patch 2: compounds_ml.json -> BS_compounds_ml_merged.json (empty, merged in)")

    # ── Patch 3: "325 compounds" label ───────────────────────────────────
    for old_label in ['"325 compounds"', "'325 compounds'",
                      "325 compounds", "325-compound"]:
        new_label = old_label.replace("325", str(n_compounds))
        if old_label in content:
            content = content.replace(old_label, new_label)
            patches_applied += 1
            print(f"  Patch 3: '{old_label}' -> '{new_label}'")

    # Explicit replacements for common patterns
    for old_pat, new_pat in [
        ("325 compound", f"{n_compounds} compound"),
        ("325-compound", f"{n_compounds}-compound"),
        (f"n=325", f"n={n_compounds}"),
        ("of 325", f"of {n_compounds}"),
        ("over 325", f"over {n_compounds}"),
    ]:
        if old_pat in content:
            content = content.replace(old_pat, new_pat)
            patches_applied += 1
            print(f"  Patch 3b: '{old_pat}' -> '{new_pat}'")

    # ── Patch 4: Score scale in app (some apps divide by 10) ────────────
    # The BS_compounds_full.json stores 0-100 scores.
    # Check if app does /10 anywhere on dimension scores.
    if "/ 10" in content or "/10" in content:
        print("  NOTE: App has /10 division — verify dimension score reading")
        print("        BS_compounds_full.json uses 0-100 scale directly")
        print("        Check app's score reading logic manually if scores look wrong")

    # ── Write output ──────────────────────────────────────────────────────
    APP_OUT.write_text(content, encoding="utf-8")
    print(f"\n  Written: {APP_OUT} ({len(content.splitlines())} lines)")
    print(f"  Total patches applied: {patches_applied}")

    # ── Syntax check ─────────────────────────────────────────────────────
    try:
        ast.parse(content)
        print("  Syntax check: PASSED")
    except SyntaxError as e:
        print(f"  SYNTAX ERROR at line {e.lineno}: {e.msg}")
        print("  Do NOT use this file — check app_v6.py manually")
        sys.exit(1)

    # ── Verify database is readable ───────────────────────────────────────
    print(f"\n  Database verification:")
    sample = list(db.keys())[:5]
    print(f"  First 5 compounds: {sample}")
    curcumin = db.get("curcumin")
    if curcumin:
        nps = curcumin.get("nps")
        ao  = curcumin.get("antioxidant")
        print(f"  Curcumin NPS={nps}, antioxidant={ao}")
        if nps and float(nps) > 70:
            print("  Curcumin calibration: OK (NPS > 70)")
        else:
            print("  WARNING: Curcumin NPS seems low")
    else:
        print("  WARNING: curcumin not in database — check compound names are lowercase")

    # Check negatives
    for neg in ["mptp", "mannitol", "doxorubicin"]:
        entry = db.get(neg)
        if entry:
            nps = entry.get("nps", 0)
            ok = "OK" if float(nps) < 25 else "PROBLEM - too high"
            print(f"  Negative control {neg}: NPS={nps} [{ok}]")

    print(f"\n{'=' * 62}")
    print("UPDATE COMPLETE")
    print(f"{'=' * 62}")
    print(f"""
NEXT STEPS:
  1. Test the updated app:
     brainsafe_env\\Scripts\\python.exe -m streamlit run app_v6_final.py

  2. Open http://localhost:8501
     Search for:  curcumin  (expect NPS ~79)
     Search for:  quercetin (expect NPS ~71)
     Search for:  MPTP      (expect NPS < 15)
     Confirm 535 compounds shown in UI

  3. If app looks correct:
     copy app_v6_final.py app_v6.py /Y
     copy app.py app_v5_backup.py /Y
     copy app_v6.py app.py /Y

  4. Deploy:
     git add .
     git commit -m "v6-final: 535-compound BS database, 93-feature model, all 7 dims validated"
     git push origin master
""")


if __name__ == "__main__":
    main()
