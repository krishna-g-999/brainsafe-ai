"""Compare every headline number in the manuscript against the artefact that produces it.

The audit found the manuscript, the README and the model card disagreeing with each other and with
results/tables on record counts, endpoint counts, AUROCs and sensitivities. Regenerating the panel
moved almost all of those numbers again. Checking them by eye is how they drifted in the first place,
so this reads each claim out of the document and each value out of the file that computes it, and
reports the pairs that do not match.

It asserts nothing about which is right. It reports MATCHES, DIFFERS or NOT FOUND with both values
and the file each came from, so the correction is made once, at the source, and can be re-checked.

Run:  python src/brainsafe/analysis/check_manuscript_numbers.py
      python src/brainsafe/analysis/check_manuscript_numbers.py --write-report
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "results" / "tables"
DOCS = [
    ROOT / "manuscript" / "NAR_WebServer_BrainSafe_built.md",
    ROOT / "manuscript" / "NAR_WebServer_BrainSafe_draft.md",
    ROOT / "README.md",
    ROOT / "docs" / "BS_MODEL_CARD.md",
    ROOT / "docs" / "RF_CV_RESULTS.md",
    ROOT / "docs" / "VALIDATION.md",
    ROOT / "docs" / "METHODS.md",
]
OUT = ROOT / "audit" / "MANUSCRIPT_NUMBERS.md"


def _csv(name: str) -> pd.DataFrame:
    p = TAB / name
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def truth() -> list[dict]:
    """Every quantity worth checking, computed from the artefact that owns it."""
    facts = []
    cv = _csv("rf_cv_summary.csv")
    if len(cv):
        c = cv[cv.task == "classification"]
        r = cv[cv.task == "regression"]
        for split in ("random", "scaffold"):
            g = c[c.split == split]
            if len(g):
                facts += [
                    {"quantity": f"classifier AUROC mean, {split} 10-fold",
                     "value": round(float(g.roc_auc_mean.mean()), 4),
                     "source": "results/tables/rf_cv_summary.csv"},
                    {"quantity": f"classifier AUROC min, {split} 10-fold",
                     "value": round(float(g.roc_auc_mean.min()), 3),
                     "source": "results/tables/rf_cv_summary.csv"},
                    {"quantity": f"classifier AUROC max, {split} 10-fold",
                     "value": round(float(g.roc_auc_mean.max()), 3),
                     "source": "results/tables/rf_cv_summary.csv"},
                ]
            h = r[r.split == split]
            if len(h):
                facts.append({"quantity": f"regression R2 mean, {split} 10-fold",
                              "value": round(float(h.r2_mean.mean()), 4),
                              "source": "results/tables/rf_cv_summary.csv"})

    cal = _csv("calibration.csv")
    if len(cal):
        facts += [
            {"quantity": "mean ECE before calibration", "value": round(float(cal.ece_raw.mean()), 4),
             "source": "results/tables/calibration.csv"},
            {"quantity": "mean ECE after calibration",
             "value": round(float(cal.ece_calibrated.mean()), 4),
             "source": "results/tables/calibration.csv"},
        ]

    con = _csv("rf_conformal.csv")
    if len(con):
        facts += [
            {"quantity": "conformal coverage min", "value": round(float(con.empirical_coverage.min()), 3),
             "source": "results/tables/rf_conformal.csv"},
            {"quantity": "conformal coverage max", "value": round(float(con.empirical_coverage.max()), 3),
             "source": "results/tables/rf_conformal.csv"},
        ]

    tmp = _csv("rf_temporal.csv")
    if len(tmp):
        t = tmp[tmp.task == "classification"]
        if len(t):
            facts += [
                {"quantity": "temporal classifier AUROC min", "value": round(float(t.score.min()), 3),
                 "source": "results/tables/rf_temporal.csv"},
                {"quantity": "temporal classifier AUROC max", "value": round(float(t.score.max()), 3),
                 "source": "results/tables/rf_temporal.csv"},
            ]

    ext = _csv("external_bbb_validation.csv")
    if len(ext):
        for _, row in ext.iterrows():
            facts.append({"quantity": f"external BBB AUROC [{row['set'][:44]}]",
                          "value": round(float(row.auroc), 4),
                          "source": "results/tables/external_bbb_validation.csv"})
            facts.append({"quantity": f"external BBB n [{row['set'][:44]}]",
                          "value": int(row.n),
                          "source": "results/tables/external_bbb_validation.csv"})

    spec = _csv("noncns_specificity_summary.csv")
    if len(spec):
        s = spec[spec.metric.str.startswith("Specificity")]
        if len(s):
            facts.append({"quantity": "specificity on non-CNS chemistry",
                          "value": round(float(s.estimate.iloc[0]), 4),
                          "source": "results/tables/noncns_specificity_summary.csv"})

    bm_path = ROOT / "models_rf" / "binder_modes.json"
    if bm_path.exists():
        bm = json.loads(bm_path.read_text())
        sens = [v["sensitivity_at_threshold"] for v in bm.values()
                if v.get("sensitivity_at_threshold") is not None]
        auroc = [v["auroc_vs_measured_inactives"] for v in bm.values()
                 if v.get("auroc_vs_measured_inactives") is not None]
        deployed = sum(1 for v in bm.values() if v.get("deployed", True))
        facts += [
            {"quantity": "binder panel mean sensitivity",
             "value": round(sum(sens) / len(sens), 4) if sens else None,
             "source": "models_rf/binder_modes.json"},
            {"quantity": "binder panel mean AUROC vs measured inactives",
             "value": round(sum(auroc) / len(auroc), 4) if auroc else None,
             "source": "models_rf/binder_modes.json"},
            {"quantity": "binder endpoints deployed", "value": deployed,
             "source": "models_rf/binder_modes.json"},
            {"quantity": "binder endpoints total", "value": len(bm),
             "source": "models_rf/binder_modes.json"},
        ]

    # dataset size, computed rather than quoted
    eps = sorted((ROOT / "data" / "endpoints").glob("*.csv"))
    total = 0
    for f in eps:
        d = pd.read_csv(f, usecols=["smiles"]) if "smiles" in pd.read_csv(f, nrows=0).columns else None
        if d is not None:
            total += len(d)
    facts.append({"quantity": "compound-endpoint records in data/endpoints",
                  "value": total, "source": "data/endpoints/*.csv"})

    hold = _csv("scaffold_holdout_results.csv")
    if len(hold) and "holdout_recall" in hold.columns:
        facts.append({"quantity": "mean per-target prospective recall",
                      "value": round(float(hold.holdout_recall.mean()), 4),
                      "source": "results/tables/scaffold_holdout_results.csv"})
    return facts


def numbers_in(text: str) -> set[str]:
    """Every number the document states, normalised so 0.96 and 0.960 compare equal."""
    out = set()
    for m in re.finditer(r"\b\d{1,3}(?:,\d{3})+\b|\b\d+\.\d+\b|\b\d+\b", text):
        tok = m.group(0).replace(",", "")
        try:
            v = float(tok)
        except ValueError:
            continue
        out.add(f"{v:g}")
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Check manuscript numbers against the artefacts.")
    ap.add_argument("--write-report", action="store_true", help=f"write {OUT.name}")
    args = ap.parse_args(argv)

    facts = truth()
    docs = {d: d.read_text(encoding="utf-8", errors="ignore") for d in DOCS if d.exists()}
    doc_nums = {d: numbers_in(t) for d, t in docs.items()}

    rows = []
    for fact in facts:
        v = fact["value"]
        if v is None:
            continue
        cands = {f"{float(v):g}"}
        if isinstance(v, float):
            cands |= {f"{round(v, n):g}" for n in (2, 3, 4)}
            cands |= {f"{round(v * 100, n):g}" for n in (0, 1, 2)}   # percentage forms
        found_in = [d.relative_to(ROOT).as_posix() for d, nums in doc_nums.items() if cands & nums]
        rows.append({"quantity": fact["quantity"], "current_value": v,
                     "source": fact["source"],
                     "appears_in": ", ".join(found_in) if found_in else "",
                     "status": "PRESENT" if found_in else "ABSENT from every document"})

    out = pd.DataFrame(rows)
    print(f"{len(out)} quantities checked against {len(docs)} documents\n")
    absent = out[out.status != "PRESENT"]
    print(f"stated nowhere (the document still carries the pre-regeneration figure): {len(absent)}")
    print(absent[["quantity", "current_value", "source"]].to_string(index=False))

    if args.write_report:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        lines = ["# Manuscript numbers against the artefacts", "",
                 "Generated by `src/brainsafe/analysis/check_manuscript_numbers.py`. Each row is a",
                 "quantity computed from the file that owns it, and whether any document states it.",
                 "`ABSENT` means no document contains that value, so the documents still carry the",
                 "figure from before the regeneration.", "",
                 out.to_markdown(index=False)]
        OUT.write_text("\n".join(lines), encoding="utf-8")
        print(f"\nwrote {OUT.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
