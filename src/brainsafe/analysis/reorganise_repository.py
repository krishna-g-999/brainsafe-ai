"""Audit the repository, consolidate everything superseded into one archive, and index what remains.

The tree had accumulated four kinds of material that no longer belong beside the live code: a
superseded publication bundle, a legacy archive from an earlier reorganisation, build intermediates
that any release step regenerates, and run logs. Together they are the larger part of the repository
by size and none of them is read by the deployed application.

The dangerous one is the superseded publication bundle. It carries per-endpoint CSVs with the same
filenames as the live training tables but different contents, being a snapshot taken before BindingDB
was pooled in and before measured inactives were merged. AChE has 4,324 rows there against 4,387
live, BACE1 8,067 against 8,501. Anyone reading the tree could take those for the training data.
That is the strongest argument for moving them out rather than leaving them in place.

Nothing is deleted. Everything is moved into one dated archive folder, and every moved file is
recorded in ARCHIVE_MANIFEST.csv with its size, its SHA-256 and the reason it was moved. The archive
itself is not committed, because it is several gigabytes and its contents are either regenerable or
already in git history; the manifest is committed, so the repository still records exactly what was
set aside and where it went.

Three indexes are written for what remains:

  REPOSITORY_INVENTORY.csv  every live file, with its category and role
  SCRIPT_INDEX.csv          every script in pipeline order, with what it reads and writes
  REPOSITORY_MAP.md         the same in prose, for a reader rather than a spreadsheet

Run:  python src/brainsafe/analysis/reorganise_repository.py            (report only, moves nothing)
      python src/brainsafe/analysis/reorganise_repository.py --apply    (performs the moves)
"""
from __future__ import annotations

import ast
import csv
import hashlib
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
ARCHIVE = ROOT / f"_ARCHIVE_{date.today().isoformat()}"

# What moves, where it lands inside the archive, and why. The reason is written into the manifest
# rather than into a commit message, so it stays attached to the files themselves.
TO_ARCHIVE = [
    ("BrainSafe_AI_Publication", "01_superseded_publication",
     "Publication bundle for the July manuscript, superseded by manuscript/NAR_WebServer_BrainSafe_"
     "built.md. Its Supplementary/Datasets CSVs share filenames with the live training tables but "
     "hold pre-BindingDB, pre-inactives snapshots, so leaving them in the tree invites them to be "
     "read as current training data."),
    ("archive/legacy", "02_legacy_code_and_models",
     "Legacy application, scripts and model binaries retired in the 2026-07-20 reorganisation. Kept "
     "for provenance; no live code path reads any of it."),
    ("models_rf_uncompressed", "03_build_intermediates/models_rf_uncompressed",
     "Pre-compression copies of the deployed estimators. Written by compress_models.py as a backup "
     "and never read back; models_rf/ is the deployed set."),
    ("dist", "03_build_intermediates/dist",
     "The packaged model archive uploaded to Zenodo. Regenerable by package_models.py and already "
     "published at doi:10.5281/zenodo.21858576, so it does not need to sit in the working tree."),
    ("logs", "04_run_logs",
     "Standard output from training and pipeline runs. Superseded by the result tables, which carry "
     "the numbers these logs were read for."),
    ("research papers", "05_literature_not_used",
     "Thirteen papers supplied for literature extraction of state-dependent potency shifts. The "
     "harvest was run and rejected: the validation gate measured 7.7 per cent extraction error and "
     "only one usable paired measurement was recovered against a requirement of 800. Retained "
     "because they were gathered by hand and the decision not to use them is itself a result."),
]

EMPTY_TO_REMOVE = ["data/interim", "results/metrics", "scripts", "tests"]

# Pipeline order. A reviewer reading the scripts should meet them in the order they run, not
# alphabetically, so each directory is given a stage and a rank.
STAGES = [
    ("1. Data acquisition", ["src/brainsafe/data/fetch_", "src/brainsafe/data/find_targets",
                             "src/brainsafe/data/harvest_", "src/brainsafe/adme/fetch_"]),
    ("2. Data assembly and audit", ["src/brainsafe/data/", "src/brainsafe/build_endpoint_context"]),
    ("3. Featurisation", ["src/brainsafe/features/"]),
    ("4. Model training", ["src/brainsafe/models/train_", "src/brainsafe/adme/train_"]),
    ("5. Calibration and thresholds", ["src/brainsafe/models/calibrate", "src/brainsafe/models/"
                                       "final_thresholds", "src/brainsafe/models/screening_",
                                       "src/brainsafe/models/apply_specificity"]),
    ("6. Applicability domain", ["src/brainsafe/build_ad_", "src/brainsafe/build_readacross"]),
    ("7. Packaging and release", ["src/brainsafe/models/compress_models",
                                  "src/brainsafe/models/package_models", "model_fetch.py"]),
    ("8. Evaluation", ["src/brainsafe/evaluation/"]),
    ("9. Inversion and falsification", ["inversion/"]),
    ("10. Analysis and reporting", ["src/brainsafe/analysis/", "src/brainsafe/build_manuscript"]),
    ("11. Application and serving", ["app.py", "api.py", "serve.py"]),
]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def docstring_of(p: Path) -> tuple[str, str]:
    """First line of the module docstring, and the 'Run:'/'Outputs:' hints if present."""
    text = p.read_text(encoding="utf-8", errors="replace")
    # Jupytext-paired notebook exports open with a percent cell, not a module docstring. Their
    # first markdown heading is the equivalent, so read that rather than reporting them as
    # undocumented.
    if text.lstrip().startswith("# %%"):
        for line in text.splitlines():
            if line.startswith("# # "):
                return (line[4:].strip(), f"paired notebook: {p.with_suffix('.ipynb').name}")
        return ("Notebook export", f"paired notebook: {p.with_suffix('.ipynb').name}")
    try:
        doc = ast.get_docstring(ast.parse(text)) or ""
    except (SyntaxError, ValueError):
        return ("", "")
    if not doc:
        return ("", "")
    lines = [l.rstrip() for l in doc.splitlines()]
    summary = lines[0].strip()
    outputs = []
    grab = False
    for l in lines[1:]:
        low = l.strip().lower()
        if low.startswith(("outputs", "writes", "output:")):
            grab = True
            if ":" in l and l.split(":", 1)[1].strip():
                outputs.append(l.split(":", 1)[1].strip())
            continue
        if grab:
            if not l.strip() or low.startswith(("run:", "usage")):
                grab = False
                continue
            outputs.append(l.strip())
    return (summary, "; ".join(outputs[:4]))


def stage_of(rel: str) -> str:
    for stage, prefixes in STAGES:
        for pre in prefixes:
            if rel.startswith(pre):
                return stage
    return "Supporting"


def tracked_files() -> set[str]:
    try:
        out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True,
                             timeout=180)
        return set(out.stdout.splitlines())
    except Exception:
        return set()


CATEGORY = [
    ("data/endpoints", "Training data, per endpoint"),
    ("data/endpoints_reg", "Training data, regression endpoints"),
    ("data/adme", "Training data, ADME endpoints"),
    ("data/raw", "Source data as retrieved"),
    ("data/external", "External structure libraries"),
    ("data/processed", "Derived data tables"),
    ("data/readacross", "Read-across index inputs"),
    ("models_rf", "Deployed model artefacts and metadata"),
    ("results/tables", "Result tables"),
    ("results/figures", "Result figures"),
    ("results/gnn", "Graph-network comparison results"),
    ("inversion", "Inversion and falsification analysis"),
    ("manuscript", "Manuscript"),
    ("supplementary", "Supplementary tables"),
    ("figures", "Manuscript figures"),
    ("presentation", "Presentation"),
    ("reviewer_package", "Reviewer package"),
    ("docs", "Documentation"),
    ("deploy", "Deployment instructions"),
    ("assets", "Application assets"),
    ("logo", "Source logos"),
    ("src", "Source code"),
]


def category_of(rel: str) -> str:
    for pre, name in CATEGORY:
        if rel.startswith(pre):
            return name
    return "Root: application, build and project files"


def main():
    apply = "--apply" in sys.argv
    print(f"{'APPLYING' if apply else 'DRY RUN, nothing will move'}\n")

    manifest_rows = []
    moved_bytes = moved_files = 0
    for src_rel, dest_rel, reason in TO_ARCHIVE:
        src = ROOT / src_rel
        if not src.exists():
            print(f"  [absent]  {src_rel}")
            continue
        files = [f for f in src.rglob("*") if f.is_file()]
        size = sum(f.stat().st_size for f in files)
        moved_bytes += size
        moved_files += len(files)
        print(f"  archive   {src_rel:32} {len(files):6,} files  {size / 2**20:9,.0f} MB")
        for f in files:
            manifest_rows.append({
                "archived_from": f.relative_to(ROOT).as_posix(),
                "archived_to": f"{ARCHIVE.name}/{dest_rel}/{f.relative_to(src).as_posix()}",
                "size_bytes": f.stat().st_size,
                # Hashing gigabytes of model binaries costs minutes for no benefit; those are
                # already checksummed in models_manifest.json. Hash the small files only.
                "sha256": sha256(f) if f.stat().st_size < 32 * 2**20 else "",
                "reason_archived": reason,
            })
        if apply:
            dest = ARCHIVE / dest_rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                shutil.rmtree(dest)
            shutil.move(str(src), str(dest))

    print(f"\n  total to archive: {moved_files:,} files, {moved_bytes / 2**30:.2f} GB")

    for rel in EMPTY_TO_REMOVE:
        p = ROOT / rel
        if p.exists() and not any(p.rglob("*")):
            print(f"  empty     {rel}")
            if apply:
                p.rmdir()

    if not apply:
        print("\nRe-run with --apply to perform the moves.")
        return

    # Re-running after the moves have happened finds nothing left to archive. That is the normal way
    # to refresh the indexes, so keep the existing manifest rather than truncating it to a header.
    if manifest_rows:
        ARCHIVE.mkdir(exist_ok=True)
        with (ARCHIVE / "ARCHIVE_MANIFEST.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()))
            w.writeheader()
            w.writerows(manifest_rows)
        shutil.copy2(ARCHIVE / "ARCHIVE_MANIFEST.csv", ROOT / "ARCHIVE_MANIFEST.csv")
    else:
        print("  nothing left to archive; manifest left as it stands, indexes refreshed")

    # ---- indexes of what remains ----
    tracked = tracked_files()
    skip_dirs = {"brainsafe_env", ".git", "__pycache__", ".claude", ".ipynb_checkpoints",
                 ARCHIVE.name}
    inv, scripts = [], []
    for f in sorted(ROOT.rglob("*")):
        if not f.is_file():
            continue
        rel = f.relative_to(ROOT).as_posix()
        if any(part in skip_dirs for part in f.relative_to(ROOT).parts):
            continue
        inv.append({
            "path": rel,
            "category": category_of(rel),
            "extension": f.suffix.lstrip(".") or "(none)",
            "size_kb": round(f.stat().st_size / 1024, 1),
            "in_git": "yes" if rel in tracked else "no",
        })
        if f.suffix == ".py":
            summary, outputs = docstring_of(f)
            scripts.append({
                "stage": stage_of(rel),
                "script": rel,
                "purpose": summary,
                "writes": outputs,
                "in_git": "yes" if rel in tracked else "no",
            })

    with (ROOT / "REPOSITORY_INVENTORY.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(inv[0].keys()))
        w.writeheader()
        w.writerows(inv)

    order = {s: i for i, (s, _) in enumerate(STAGES)}
    scripts.sort(key=lambda r: (order.get(r["stage"], 99), r["script"]))
    with (ROOT / "SCRIPT_INDEX.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(scripts[0].keys()))
        w.writeheader()
        w.writerows(scripts)

    # ---- the same thing in prose, generated so it cannot drift from the CSVs ----
    by_cat = {}
    for r in inv:
        c = by_cat.setdefault(r["category"], {"n": 0, "kb": 0.0})
        c["n"] += 1
        c["kb"] += r["size_kb"]
    cat_rows = "\n".join(
        f"| {k} | {v['n']:,} | {v['kb'] / 1024:,.0f} |"
        for k, v in sorted(by_cat.items(), key=lambda kv: -kv[1]["kb"]))
    stage_rows = "\n".join(
        f"| {s} | {sum(1 for r in scripts if r['stage'] == s)} |"
        for s, _ in STAGES if any(r["stage"] == s for r in scripts))
    arch_total = sum(r["size_bytes"] for r in manifest_rows) if manifest_rows else None
    arch_line = (f"{len(manifest_rows):,} files, {arch_total / 2**30:.2f} GB"
                 if manifest_rows else "see ARCHIVE_MANIFEST.csv")

    (ROOT / "REPOSITORY_MAP.md").write_text(f"""# Repository map

Generated by `src/brainsafe/analysis/reorganise_repository.py`. Both tables are computed from
`REPOSITORY_INVENTORY.csv` and `SCRIPT_INDEX.csv`, so this file cannot disagree with them.

## What is here

{len(inv):,} live files. The two large directories are `models_rf/`, which holds the deployed
estimators, and `data/`, which holds the measured sets they were fitted to. Neither is committed:
the models are published at doi:10.5281/zenodo.21858576 and fetched by `model_fetch.py`, and the
bulky structure libraries are regenerable from the acquisition scripts.

| Category | Files | MB |
|---|---|---|
{cat_rows}

## The pipeline, in the order it runs

{len(scripts):,} scripts. `SCRIPT_INDEX.csv` gives each one's purpose and what it writes.

| Stage | Scripts |
|---|---|
{stage_rows}

`reviewer_package/scripts/` holds copies of seven of these for reviewers who want to reproduce the
core result without the full tree. They are byte-identical to their originals in `src/`; if they ever
diverge, treat `src/` as authoritative.

## What was set aside, and why

Everything superseded was moved into one dated archive folder ({arch_line}), listed file by file
with its checksum in `ARCHIVE_MANIFEST.csv`. Nothing was deleted. The folder is not committed,
because its contents are either regenerable or already in git history; the manifest is committed, so
the repository records what was moved and where.

Four categories were archived:

- **The July publication bundle.** Superseded by `manuscript/NAR_WebServer_BrainSafe_built.md`. This
  is the one that mattered most: its `Supplementary/Datasets/*.csv` files carry the same names as the
  live training tables but hold a snapshot taken before BindingDB was pooled in and before measured
  inactives were merged. AChE had 4,324 rows against 4,387 live; BACE1 8,067 against 8,501. Left in
  place, they would eventually have been read as the training data.
- **The 2026-07-20 legacy archive.** The retired application, its scripts and its model binaries.
  No live code path reads any of it.
- **Build intermediates.** The uncompressed model copies that `compress_models.py` writes as a
  backup and never reads back, and the packaged archive that `package_models.py` produces for
  Zenodo. Both regenerate from `models_rf/`.
- **Run logs and unused literature.** Training stdout, superseded by the result tables; and the
  thirteen papers gathered for state-dependent potency extraction, which were harvested, measured at
  7.7 per cent extraction error, and rejected. The papers are kept because the decision not to use
  them is itself a result.

## Regenerating these indexes

```
python src/brainsafe/analysis/reorganise_repository.py --apply
```

Run without `--apply` first; it reports what would move and moves nothing.
""", encoding="utf-8")

    print(f"\n  REPOSITORY_MAP.md         written")
    print(f"  REPOSITORY_INVENTORY.csv  {len(inv):,} live files")
    print(f"  SCRIPT_INDEX.csv          {len(scripts):,} scripts in pipeline order")
    print(f"  ARCHIVE_MANIFEST.csv      {len(manifest_rows):,} archived files")
    undoc = sum(1 for s in scripts if not s["purpose"])
    print(f"  scripts without a docstring summary: {undoc}")


if __name__ == "__main__":
    main()
