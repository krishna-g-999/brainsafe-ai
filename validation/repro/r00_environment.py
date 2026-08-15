"""Capture everything needed to say what produced a number, before producing any.

A reproduced metric is worth nothing without the state that produced it. This records the commit,
whether the tree was dirty at the time, the interpreter, the installed packages, the hardware, and
the seeds declared in the source. Every row of the reproduction ledger cites the commit recorded
here, so a reader can tell whether two rows came from the same state.

The code under test must be committed. A reproduction run against uncommitted edits to the pipeline
cannot be repeated by anyone else, so this fails if any tracked file is modified. Untracked files are
recorded but are not fatal: the reproduction scripts themselves are new files, and adding a script
that reads the pipeline does not change the pipeline it reads.

Output: validation/repro/environment.json, validation/repro/pip_freeze.txt

Run:  python validation/repro/r00_environment.py
"""
from __future__ import annotations

import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "validation" / "repro"


def sh(*args) -> str:
    try:
        return subprocess.run(args, cwd=ROOT, capture_output=True, text=True,
                              timeout=120).stdout.strip()
    except Exception as exc:                                   # pragma: no cover
        return f"<failed: {exc}>"


def declared_seeds() -> dict[str, list[str]]:
    """Every literal seed in the source, so a claim of determinism can be checked."""
    found: dict[str, list[str]] = {}
    for p in sorted((ROOT / "src").rglob("*.py")):
        hits = []
        for i, line in enumerate(p.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
            if re.search(r"\b(SEED|random_state|random\.seed|default_rng)\b", line):
                hits.append(f"{i}: {line.strip()[:110]}")
        if hits:
            found[p.relative_to(ROOT).as_posix()] = hits
    return found


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    porcelain = [ln for ln in sh("git", "status", "--porcelain").splitlines() if ln.strip()]
    modified = [ln for ln in porcelain if not ln.startswith("??")]
    untracked = [ln for ln in porcelain if ln.startswith("??")]

    freeze = sh(sys.executable, "-m", "pip", "freeze")
    (OUT / "pip_freeze.txt").write_text(freeze + "\n", encoding="utf-8")

    key = {}
    for line in freeze.splitlines():
        name = re.split(r"[=<>@ ]", line, maxsplit=1)[0].lower()
        if name in {"scikit-learn", "numpy", "pandas", "rdkit", "scipy", "xgboost", "joblib",
                    "matplotlib", "streamlit"}:
            key[name] = line.strip()

    env = {
        "captured_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "commit": sh("git", "rev-parse", "HEAD"),
        "commit_subject": sh("git", "log", "-1", "--pretty=%s"),
        "branch": sh("git", "rev-parse", "--abbrev-ref", "HEAD"),
        "tracked_files_clean": not modified,
        "modified_tracked_files": modified,
        "untracked_files": untracked,
        "python": sys.version.replace("\n", " "),
        "executable": sys.executable,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "machine": platform.machine(),
        "cpu_count": __import__("os").cpu_count(),
        "key_packages": key,
        "declared_seeds": declared_seeds(),
    }
    (OUT / "environment.json").write_text(json.dumps(env, indent=2), encoding="utf-8")

    print(f"commit      {env['commit']}  ({env['commit_subject'][:60]})")
    print(f"tracked     {'clean' if not modified else 'MODIFIED'}"
          f"{'' if not untracked else f', {len(untracked)} untracked path(s)'}")
    if modified:
        print("  tracked files are modified, so this run cannot be repeated from the commit:")
        for d in modified[:20]:
            print(f"    {d}")
    print(f"python      {sys.version.split()[0]}")
    for k, v in sorted(key.items()):
        print(f"  {v}")
    print(f"cpu         {env['cpu_count']} logical, {env['platform']}")
    print(f"seed sites  {sum(len(v) for v in env['declared_seeds'].values())} across "
          f"{len(env['declared_seeds'])} files")
    print(f"\nwrote {(OUT / 'environment.json').relative_to(ROOT).as_posix()} and pip_freeze.txt")
    if modified:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
