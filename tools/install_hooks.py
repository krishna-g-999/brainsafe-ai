"""Install the git hooks that make the freshness check unskippable.

The freshness check only prevents the incidents it was written for if it runs every time, and a
check a person has to remember is a check that eventually does not run. Four stale-artefact
incidents in this project happened while a careful person was paying attention.

The pre-commit hook runs tools/check_freshness.py and refuses the commit if any declared artefact is
older than its inputs. It is deliberately overridable, because a legitimate case exists: committing
the very fix that will bring an artefact back into date. Overriding is explicit and leaves a trace,
which is the point.

    git commit --no-verify        skip the check, deliberately and visibly

Run:  python tools/install_hooks.py
      python tools/install_hooks.py --uninstall
"""
from __future__ import annotations

import argparse
import os
import stat
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

HOOK = r"""#!/bin/sh
# Installed by tools/install_hooks.py. Refuses a commit whose declared artefacts are stale.
# Override deliberately with: git commit --no-verify

PY="./brainsafe_env/Scripts/python.exe"
[ -x "$PY" ] || PY="python"

echo "[pre-commit] checking artefact freshness ..."
"$PY" tools/check_freshness.py --quiet
STATUS=$?

if [ $STATUS -ne 0 ]; then
    echo ""
    echo "[pre-commit] BLOCKED: an artefact is older than something it was derived from."
    echo "[pre-commit] Rebuild it, or accept the current state if you have verified it:"
    echo "[pre-commit]     python tools/check_freshness.py --accept"
    echo "[pre-commit] To commit anyway, deliberately:"
    echo "[pre-commit]     git commit --no-verify"
    exit 1
fi

echo "[pre-commit] artefacts are consistent with their inputs"
exit 0
"""


def hooks_dir() -> Path:
    out = subprocess.run(["git", "rev-parse", "--git-path", "hooks"], cwd=ROOT,
                         capture_output=True, text=True).stdout.strip()
    return (ROOT / out) if out else (ROOT / ".git" / "hooks")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Install or remove the pre-commit freshness hook.")
    ap.add_argument("--uninstall", action="store_true")
    args = ap.parse_args(argv)

    hd = hooks_dir()
    hd.mkdir(parents=True, exist_ok=True)
    path = hd / "pre-commit"

    if args.uninstall:
        if path.exists():
            path.unlink()
            print(f"removed {path}")
        else:
            print("no pre-commit hook installed")
        return

    if path.exists() and "check_freshness" not in path.read_text(encoding="utf-8", errors="ignore"):
        backup = path.with_suffix(".before-freshness")
        backup.write_text(path.read_text(encoding="utf-8", errors="ignore"), encoding="utf-8")
        print(f"an unrelated pre-commit hook was already here; kept a copy at {backup.name}")

    path.write_text(HOOK, encoding="utf-8", newline="\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    print(f"installed {path}")
    print("every commit now runs tools/check_freshness.py first")
    print("override deliberately with: git commit --no-verify")


if __name__ == "__main__":
    main(sys.argv[1:])
