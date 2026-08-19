"""Is any derived artefact older than something it was derived from?

Every stale-artefact incident in this project has had the same shape. A file is produced from an
input, the input is later regenerated, and nothing notices that the output is now describing
something that no longer exists. The failure is silent by construction: the stale file is internally
consistent, so every content check passes, and every command exits 0.

Content checks cannot catch this. `check_manuscript_numbers.py` compares numbers against the tables,
and it is satisfied when the manuscript agrees with a table that is itself six days out of date.
Freshness is a different property and needs a different check.

This declares the dependency graph explicitly and compares modification times. If an output is older
than any of its inputs it is STALE, the command that rebuilds it is printed, and the exit code is
non-zero so a pipeline or a hook can refuse to proceed.

What it deliberately does NOT do: rebuild anything. Deciding what to regenerate is a scientific
judgement, because rebuilding an artefact changes what published numbers refer to. This reports and
stops.

Staleness is judged on content, not on timestamps. Timestamps alone produce false alarms that train
a reader to ignore the tool: re-running the threshold sequence rewrote binder_modes.json with
byte-identical thresholds, and a timestamp check called nine downstream artefacts stale when nothing
had changed. A check that cries wolf is worse than no check.

Each input therefore gets a fingerprint. Small text artefacts are fingerprinted by SHA-256 of their
content. Large binary trees, where hashing 0.78 GB on every run would make the check too slow to use,
are fingerprinted by file count and total size, which is stated here because it is weaker: a retrain
that produced exactly the same number of files at exactly the same total size would not be caught.

`--accept` records the current fingerprints in tools/freshness_state.json, and means "I have checked
that everything is consistent as of now". Later runs compare against that record, so an artefact is
stale only when an input's fingerprint has actually moved since it was last accepted.

models_manifest.json is verified by recomputing every SHA-256 it records, which is stronger still.

Run:  python tools/check_freshness.py
      python tools/check_freshness.py --accept    record the current state as consistent
      python tools/check_freshness.py --json      machine-readable
      python tools/check_freshness.py --quiet     only the stale ones
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Four steps that all mutate models_rf/binder_modes.json, in the only order that is correct. See the
# ORDERING HAZARD note in GRAPH below.
THRESHOLD_SEQUENCE = (
    "python src/brainsafe/models/final_thresholds.py && "
    "python src/brainsafe/models/screening_thresholds.py && "
    "python src/brainsafe/models/apply_specificity_decisions.py && "
    "python src/brainsafe/models/calibrate_background_specificity.py"
)

# output -> (inputs, command that rebuilds it)
# Inputs may be files, directories or globs. A directory is represented by its newest file, which is
# what "this directory changed" means for these purposes.
GRAPH: list[tuple[str, list[str], str]] = [
    # ---- models are derived from the endpoint tables -----------------------------------------
    ("models_rf/BBB.joblib", ["data/endpoints/*.csv"],
     "python src/brainsafe/models/train_rf.py"),
    ("results/tables/rf_cv_summary.csv", ["models_rf/BBB.joblib"],
     "python src/brainsafe/models/train_rf.py"),
    ("results/tables/rf_cv_folds.csv", ["models_rf/BBB.joblib"],
     "python src/brainsafe/models/train_rf.py"),

    # ---- calibration sits on the core models --------------------------------------------------
    ("results/tables/calibration.csv", ["models_rf/BBB.joblib"],
     "python src/brainsafe/models/calibrate.py"),

    # ---- the binder panel ----------------------------------------------------------------------
    ("models_rf/binder_modes.json", ["data/endpoints/*.csv"],
     "python src/brainsafe/models/train_binders_hybrid.py"),
    # ORDERING HAZARD. These four steps all write models_rf/binder_modes.json and each depends on
    # the one before it. Running any of them alone silently reverts the later ones: re-running
    # final_thresholds.py by itself undid the background-specificity tightening for HT1A, HT2A and
    # Nav1_5, and running it without apply_specificity_decisions.py afterwards would re-deploy two
    # endpoints that were withdrawn for firing on glucose and atenolol. The whole sequence is
    # therefore quoted as the rebuild command for every member of it.
    # These three depend on the fitted binder models, NOT on binder_modes.json, even though every
    # step of the sequence reads it. binder_modes.json is a co-output of this same sequence: steps
    # three and four rewrite it after step one has written final_thresholds.csv, so an edge from it
    # to these tables can never be satisfied. A correct, complete run left them reported stale, and
    # re-running to clear the report re-ran step one alone, which is precisely the ordering hazard
    # above. A check that cannot be satisfied by doing the right thing teaches the reader to ignore
    # it. The binder .joblib files carry the same signal without the cycle: they change when the
    # panel is retrained, which is when these thresholds genuinely must be re-derived, and the
    # threshold sequence never writes them.
    ("results/tables/final_thresholds.csv", ["models_rf/*_binder.joblib"], THRESHOLD_SEQUENCE),
    ("results/tables/screening_thresholds.csv", ["models_rf/*_binder.joblib"],
     THRESHOLD_SEQUENCE),
    ("results/tables/background_specificity.csv", ["models_rf/*_binder.joblib"],
     THRESHOLD_SEQUENCE),

    # ---- the evaluation layer, all of it downstream of the models ------------------------------
    ("results/tables/external_bbb_validation.csv", ["models_rf/BBB.joblib"],
     "python src/brainsafe/evaluation/external_validation.py"),
    ("results/tables/rf_conformal.csv", ["models_rf/BBB.joblib"],
     "python src/brainsafe/evaluation/rf_conformal_temporal.py"),
    ("results/tables/inversion_validation.csv", ["models_rf/BBB.joblib"],
     "python src/brainsafe/evaluation/validate_inversion.py"),
    ("results/tables/integrity_audit.csv", ["models_rf/BBB.joblib"],
     "python src/brainsafe/evaluation/integrity_audit.py"),
    ("results/tables/noncns_specificity_summary.csv",
     ["models_rf/BBB.joblib", "models_rf/binder_modes.json"],
     "python src/brainsafe/evaluation/noncns_specificity.py"),
    ("results/tables/binder_cv_summary.csv", ["models_rf/binder_modes.json"],
     "python src/brainsafe/evaluation/binder_cv_per_fold.py"),

    # The two-step hold-out. This pair is the incident that motivated the whole file: the panel was
    # re-run and the report was not, so the inputs were current and the output was a day stale.
    ("models_rf/holdout", ["models_rf/binder_modes.json"],
     "python src/brainsafe/evaluation/scaffold_holdout_panel.py"),
    ("results/tables/scaffold_holdout_results.csv", ["models_rf/holdout"],
     "python src/brainsafe/evaluation/scaffold_holdout_report.py"),

    # ---- the reviewer package -------------------------------------------------------------------
    ("reviewer_package/model_outputs/TRAINING_MATRIX.csv", ["data/endpoints/*.csv"],
     "python src/brainsafe/analysis/build_reviewer_matrix.py"),
    ("reviewer_package/model_outputs/BrainSafe_models_inputs_and_folds.xlsx",
     ["reviewer_package/model_outputs/TRAINING_MATRIX.csv", "results/tables/rf_cv_folds.csv"],
     "python src/brainsafe/analysis/build_reviewer_workbook_full.py"),

    # ---- figures read the tables ----------------------------------------------------------------
    ("manuscript/figures/Figure3_cv_design.png", ["results/tables/rf_cv_summary.csv"],
     "python src/brainsafe/figures/fig03_cv_design.py"),
    ("manuscript/figures/Figure4_pools_and_thresholds.png", ["models_rf/binder_modes.json"],
     "python src/brainsafe/figures/fig04_pools_and_thresholds.py"),
    ("manuscript/figures/Figure5_negative_class.png", ["results/tables/rf_cv_summary.csv"],
     "python src/brainsafe/figures/fig05_negative_class.py"),
    ("manuscript/figures/Figure6_validation.png",
     ["results/tables/calibration.csv", "results/tables/scaffold_holdout_results.csv",
      "results/tables/noncns_specificity_summary.csv",
      "results/tables/external_bbb_validation.csv", "results/tables/inversion_validation.csv"],
     "python src/brainsafe/figures/fig06_validation.py"),
    ("manuscript/figures/Figure7_binder_panel.png", ["models_rf/binder_modes.json"],
     "python src/brainsafe/figures/fig07_binder_panel.py"),

    # ---- the manuscript is downstream of everything it quotes -----------------------------------
    ("manuscript/tables_generated.md",
     ["results/tables/rf_cv_summary.csv", "models_rf/binder_modes.json"],
     "python src/brainsafe/analysis/manuscript_tables.py"),
    ("manuscript/NAR_WebServer_BrainSafe_built.md",
     ["manuscript/NAR_WebServer_BrainSafe_draft.md", "manuscript/tables_generated.md",
      "manuscript/references_verified.json", "manuscript/figures/Figure6_validation.png"],
     "python src/brainsafe/analysis/build_manuscript.py"),
]

# Checked by content rather than by timestamp, because it carries checksums of what it describes.
CHECKSUMMED = "models_manifest.json"
STATE = ROOT / "tools" / "freshness_state.json"
HASH_LIMIT = 32 * 1024 * 1024        # above this, fingerprint by size rather than content


def fingerprint(pattern: str) -> str | None:
    """A content fingerprint for a file, glob or directory.

    Small files are hashed. Large binary trees are described by file count and total size, which is
    weaker and is documented as such: it would miss a retrain that happened to produce the same
    number of files at the same total size.
    """
    p = ROOT / pattern
    if "*" in pattern:
        hits = sorted(ROOT.glob(pattern))
    elif p.is_dir():
        hits = sorted(q for q in p.rglob("*") if q.is_file())
    elif p.exists():
        hits = [p]
    else:
        return None
    total = sum(q.stat().st_size for q in hits)
    if len(hits) == 1 and total <= HASH_LIMIT:
        h = hashlib.sha256()
        with open(hits[0], "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return "sha256:" + h.hexdigest()
    if total <= HASH_LIMIT:
        h = hashlib.sha256()
        for q in hits:
            h.update(q.relative_to(ROOT).as_posix().encode())
            h.update(q.read_bytes())
        return "sha256:" + h.hexdigest()
    return f"count-size:{len(hits)}:{total}"


def newest(pattern: str):
    """The newest file matching a path, glob or directory, and which file it was."""
    p = ROOT / pattern
    if "*" in pattern:
        hits = sorted(ROOT.glob(pattern))
    elif p.is_dir():
        hits = [q for q in p.rglob("*") if q.is_file()]
    elif p.exists():
        hits = [p]
    else:
        hits = []
    if not hits:
        return None, None
    best = max(hits, key=lambda q: q.stat().st_mtime)
    return best.stat().st_mtime, best


def stamp(ts) -> str:
    return dt.datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M") if ts else "missing"


def check_manifest() -> dict:
    """Verify every checksum the manifest records, which is stronger than any timestamp."""
    path = ROOT / CHECKSUMMED
    if not path.exists():
        return {"artefact": CHECKSUMMED, "state": "MISSING", "detail": "no manifest",
                "rebuild": "python src/brainsafe/models/package_models.py 1.1"}
    man = json.loads(path.read_text(encoding="utf-8"))
    files = man.get("files", {})
    ok = bad = missing = 0
    for rel, rec in files.items():
        q = ROOT / rel
        if not q.exists():
            missing += 1
            continue
        want = rec["sha256"] if isinstance(rec, dict) else rec
        h = hashlib.sha256()
        with open(q, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        ok += 1 if h.hexdigest() == want else 0
        bad += 0 if h.hexdigest() == want else 1
    state = "OK" if (bad == 0 and missing == 0) else "STALE"
    return {"artefact": CHECKSUMMED, "state": state,
            "detail": f"{ok} verified, {bad} checksum mismatch, {missing} missing, "
                      f"of {len(files)} entries",
            "rebuild": "python src/brainsafe/models/package_models.py 1.1"}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Report artefacts older than their inputs.")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--quiet", action="store_true", help="print only stale artefacts")
    ap.add_argument("--accept", action="store_true",
                    help="record the current fingerprints as a consistent state")
    args = ap.parse_args(argv)

    recorded = json.loads(STATE.read_text(encoding="utf-8")) if STATE.exists() else {}
    current: dict[str, str] = {}
    for _out, ins, _cmd in GRAPH:
        for src in ins:
            if src not in current:
                fp = fingerprint(src)
                if fp:
                    current[src] = fp

    results = []
    for out, ins, cmd in GRAPH:
        out_ts, _ = newest(out)
        if out_ts is None:
            results.append({"artefact": out, "state": "MISSING", "detail": "not built",
                            "rebuild": cmd})
            continue
        stale_by = []
        for src in ins:
            if src not in current:
                stale_by.append(f"{src} (missing)")
                continue
            # Staleness needs both conditions. An input touched after the output was built is
            # suspicious, but only actually stale if its content moved: re-running a step that
            # rewrites a file identically must not condemn everything downstream. Conversely a
            # content change that predates the rebuild is already incorporated, which is the case
            # every time an input is edited and its output regenerated straight afterwards.
            in_ts, in_file = newest(src)
            touched_after_build = bool(in_ts and in_ts > out_ts)
            if not touched_after_build:
                continue
            was = recorded.get(src)
            if was is None:
                stale_by.append(f"{in_file.relative_to(ROOT).as_posix()} newer "
                                f"@ {stamp(in_ts)} (no accepted baseline)")
            elif was != current[src]:
                stale_by.append(f"{src} changed after this was built")
        results.append({
            "artefact": out,
            "state": "STALE" if stale_by else "OK",
            "detail": (f"built {stamp(out_ts)}; " + "; ".join(stale_by)) if stale_by
                      else f"built {stamp(out_ts)}",
            "rebuild": cmd,
        })
    results.append(check_manifest())

    bad = [r for r in results if r["state"] != "OK"]
    if args.accept:
        STATE.write_text(json.dumps(current, indent=2, sort_keys=True), encoding="utf-8")
        print(f"recorded {len(current)} input fingerprints in "
              f"{STATE.relative_to(ROOT).as_posix()}")
        print("later runs will call an artefact stale only when an input's content has moved")
        raise SystemExit(0)
    if args.json:
        print(json.dumps({"stale": len(bad), "results": results}, indent=2))
    else:
        width = max(len(r["artefact"]) for r in results)
        for r in results:
            if args.quiet and r["state"] == "OK":
                continue
            mark = {"OK": "  ok  ", "STALE": " STALE", "MISSING": "MISSING"}[r["state"]]
            print(f"[{mark}] {r['artefact']:<{width}}  {r['detail']}")
        print()
        if bad:
            print(f"{len(bad)} artefact(s) are stale or missing. Rebuild, in this order:")
            seen = set()
            for r in bad:
                if r["rebuild"] not in seen:
                    print(f"    {r['rebuild']}")
                    seen.add(r["rebuild"])
        else:
            print("every declared artefact is newer than its inputs")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main(sys.argv[1:])
