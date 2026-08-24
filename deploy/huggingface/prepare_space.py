"""Assemble a Hugging Face Space directory containing exactly what the server needs.

Copying the repository wholesale is the usual way this goes wrong. It carries the 155 MB of
scaffold-hold-out twins, the reviewer package, the manuscript figures and the validation artefacts,
none of which serve a prediction, and the upload is then several times larger than it has to be for
no benefit.

This copies the runtime set, writes the Space card and the git-LFS rules, and reports the size. It
does not push: pushing is a credentialled act and belongs to whoever owns the Space.

Run:  python deploy/huggingface/prepare_space.py --out ../brainsafe-space
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

# What the server reads at runtime. Anything not named here is not copied, which is the point.
# requirements.txt is deliberately absent: the root one is the training environment and the Space
# gets the runtime subset from this directory instead. See deploy/huggingface/requirements.txt.
FILES = ["app.py", "api.py", "serve.py", "model_fetch.py", "models_manifest.json",
         "CITATION.cff", "LICENSE"]
DIRS = ["src", "assets", "results", "docs", ".streamlit"]

# data/ is 582 MB and the server reads 18 MB of it. The rest is the raw pulls, the API caches and
# the external sets that the training and validation scripts consume, none of which is opened to
# answer a query: app.py reads data/endpoints only, and does not import the pools module that
# reaches for the external drug list. Copying data/ wholesale made the upload three times larger
# than the thing being uploaded.
DATA_SUBDIRS = ["endpoints", "endpoints_reg", "adme", "readacross"]

# models_rf/ is copied separately so holdout/ can be left behind: it holds the scaffold-split twins
# used for validation and is never consulted to answer a query.
MODELS_SKIP = {"holdout"}

LFS = """*.joblib filter=lfs diff=lfs merge=lfs -text
*.pkl filter=lfs diff=lfs merge=lfs -text
*.tar.gz filter=lfs diff=lfs merge=lfs -text
*.xlsx filter=lfs diff=lfs merge=lfs -text
"""


def copy_tree(src: Path, dst: Path, skip: set[str] = frozenset()) -> int:
    total = 0
    for p in src.rglob("*"):
        if not p.is_file():
            continue
        rel = p.relative_to(src)
        if rel.parts and rel.parts[0] in skip:
            continue
        if "__pycache__" in rel.parts or p.suffix == ".pyc":
            continue
        out = dst / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, out)
        total += p.stat().st_size
    return total


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Assemble a Space directory.")
    ap.add_argument("--out", required=True, help="target directory (the cloned Space)")
    ap.add_argument("--with-holdout", action="store_true",
                    help="include models_rf/holdout, which the server never reads")
    args = ap.parse_args(argv)

    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    print(f"assembling into {out}\n")

    total = 0
    for name in FILES:
        src = ROOT / name
        if not src.exists():
            print(f"  missing, skipped   {name}")
            continue
        shutil.copy2(src, out / name)
        total += src.stat().st_size
        print(f"  {src.stat().st_size/1e6:8.2f} MB  {name}")

    for name in DIRS:
        src = ROOT / name
        if not src.exists():
            print(f"  missing, skipped   {name}/")
            continue
        n = copy_tree(src, out / name)
        total += n
        print(f"  {n/1e6:8.2f} MB  {name}/")

    for sub in DATA_SUBDIRS:
        src = ROOT / "data" / sub
        if not src.exists():
            print(f"  missing, skipped   data/{sub}/")
            continue
        n = copy_tree(src, out / "data" / sub)
        total += n
        print(f"  {n/1e6:8.2f} MB  data/{sub}/")

    skip = set() if args.with_holdout else MODELS_SKIP
    n = copy_tree(ROOT / "models_rf", out / "models_rf", skip=skip)
    total += n
    print(f"  {n/1e6:8.2f} MB  models_rf/" + ("" if args.with_holdout else "  (holdout excluded)"))

    shutil.copy2(HERE / "README.md", out / "README.md")       # the Space card, with its YAML block
    shutil.copy2(HERE / "requirements.txt", out / "requirements.txt")   # runtime subset, not root
    shutil.copy2(HERE / "Dockerfile", out / "Dockerfile")               # the SDK no longer builds it
    print(f"  {'':8s}    requirements.txt  (runtime subset)")
    print(f"  {'':8s}    Dockerfile")
    (out / ".gitattributes").write_text(LFS, encoding="utf-8")
    print(f"  {'':8s}    README.md  (Space card)")
    print(f"  {'':8s}    .gitattributes  (git-LFS rules)")

    print(f"\n  total {total/1e9:.2f} GB")
    print("\nNext, from inside that directory:")
    print("    git lfs install")
    print("    git add .gitattributes && git commit -m 'Track model files with LFS'")
    print("    git add -A && git commit -m 'BrainSafe AI'")
    print("    git push")
    print("\nTrack LFS before adding the models, or the push is rejected for file size.")


if __name__ == "__main__":
    main(sys.argv[1:])
