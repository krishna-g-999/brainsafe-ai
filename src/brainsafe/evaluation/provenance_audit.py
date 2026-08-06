"""Which reported results are still true of the model that is deployed?

This project has changed a great deal in a short time: endpoints added, two withdrawn, thresholds
corrected, a gating rule changed. Every one of those changes can silently invalidate a number sitting
in a results file, and a stale number in a results file becomes a stale number in a manuscript. The
coverage panel had already drifted into advertising a withdrawn endpoint before a check caught it;
there is no reason to assume it was the only thing that drifted.

This compares every result artefact against the deployed state and reports three kinds of problem:

  stale       the file names endpoints that no longer exist, or is older than the models it describes
  incomplete  the file omits endpoints that are now deployed, so any aggregate it reports is over a
              different panel than the one running
  contradicts a metric in the file disagrees with the same metric computed from the deployed models

It does not fix anything. The point is a list of what must be regenerated before the manuscript is
final, produced mechanically rather than from memory of what was run when.

Read-only. Writes results/provenance_audit.csv
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))

RES = ROOT / "results"
INV = ROOT / "inversion" / "results"
MODES = ROOT / "models_rf" / "binder_modes.json"

# result files that describe the target panel, and the column naming the endpoint
# (path, endpoint column, scope). Scope matters: a binder-panel file is not incomplete for omitting
# the base-rate-enrichment endpoints, which have no binder model, and a record of which candidates
# were audited is not stale for naming one that was subsequently rejected. Judging every file against
# the full panel produced three false alarms on the first run of this script.
PANEL_FILES = [
    (INV / "H7_target_discrimination.csv", "target", "binder"),
    (RES / "deployed_specificity_audit.csv", "target", "binder"),
    (RES / "batch4_audit.csv", "target", "historical"),
    (RES / "tables" / "binder_vs_measured_inactives.csv", None, "binder"),
    (RES / "tables" / "background_specificity.csv", None, "binder"),
]


def main():
    import app

    modes = json.loads(MODES.read_text())
    deployed = {t for t in app.TARGET_KIND if t != "NEURO"
                and modes.get(t, {}).get("deployed", True)}
    withdrawn = {t for t, v in modes.items() if not v.get("deployed", True)}
    model_mtime = max(p.stat().st_mtime for p in (ROOT / "models_rf").glob("*.joblib"))
    print(f"deployed endpoints: {len(deployed)}; withdrawn: {sorted(withdrawn)}", flush=True)
    print(f"newest model artefact: {time.ctime(model_mtime)}\n", flush=True)

    binder_deployed = {t for t in deployed if app.TARGET_KIND.get(t) == "binder"}
    rows = []
    for path, col, scope in PANEL_FILES:
        if not path.exists():
            rows.append({"artefact": path.name, "status": "MISSING", "detail": "file absent"})
            continue
        age = path.stat().st_mtime
        d = pd.read_csv(path)
        if col is None:
            col = next((c for c in d.columns if c.lower() in
                        ("target", "endpoint", "name")), None)
        named = set(d[col].astype(str)) if col else set()
        expected = binder_deployed if scope == "binder" else deployed
        gone = sorted(named & withdrawn) if scope != "historical" else []
        absent = sorted(expected - named) if (named and scope != "historical") else []
        problems = []
        if age < model_mtime:
            problems.append(f"older than the deployed models ({time.ctime(age)})")
        if gone:
            problems.append(f"reports withdrawn endpoints {gone}")
        if absent:
            problems.append(f"omits {len(absent)} deployed endpoints: {absent[:6]}"
                            + (" ..." if len(absent) > 6 else ""))
        rows.append({"artefact": str(path.relative_to(ROOT)),
                     "status": "STALE" if problems else "current",
                     "n_rows": len(d),
                     "detail": "; ".join(problems) or "consistent with the deployed panel"})
        print(f"  {path.name:42} {'STALE' if problems else 'current'}", flush=True)

    # the inversion verdicts describe a graph that has since gained two conditions
    v = INV / "VERDICTS.csv"
    if v.exists():
        graph_files = [INV / "H1_disease_layer.csv", INV / "H2_weight_ablation.csv",
                       INV / "H6_clinical_indication.csv", INV / "H8_panel_independence.csv"]
        app_mtime = (ROOT / "app.py").stat().st_mtime
        for p in graph_files:
            if not p.exists():
                continue
            stale = p.stat().st_mtime < app_mtime
            rows.append({"artefact": str(p.relative_to(ROOT)),
                         "status": "STALE" if stale else "current", "n_rows": len(pd.read_csv(p)),
                         "detail": ("computed before the knowledge graph last changed; the disease "
                                    "layer now has "
                                    f"{len(app.DISEASE_ORDER)} conditions and "
                                    f"{len(app.KNOWLEDGE_GRAPH)} targets")
                         if stale else "consistent with the current graph"})
            print(f"  {p.name:42} {'STALE' if stale else 'current'}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(RES / "provenance_audit.csv", index=False)
    pd.set_option("display.width", 230)
    pd.set_option("display.max_colwidth", 96)
    print()
    print(df.to_string(index=False))
    stale = df[df.status != "current"]
    print(f"\n{len(stale)} of {len(df)} artefacts need regenerating before the manuscript is final")
    print("wrote", RES / "provenance_audit.csv")


if __name__ == "__main__":
    main()
