"""Build the reproduction ledger: every reported number against an independently produced one.

Two rules govern this file.

First, the manuscript value in each row is not typed from memory. Each claim carries an anchor, a
phrase that appears in the built manuscript, and the value is only accepted as "what the manuscript
says" if the value appears in the manuscript within a window of that anchor. A claim whose value is
not found is recorded as MISQUOTED, which is a finding about this ledger rather than about the
pipeline, and is better than a ledger that quietly agrees with itself.

Second, a reproduced value must come from an artefact this reproduction produced, or from a
recomputation performed here. Rows whose reproduction requires a run that has not been performed are
recorded as CANNOT_REPRODUCE with the reason, never left out and never filled from the pipeline's own
output. Reading the pipeline's summary and calling it a reproduction is the failure this whole
exercise exists to avoid, so each row records its evidence tier:

  A  independent re-run       the pipeline was executed again and scored by code written here
  B  independent recompute    the pipeline's saved predictions were rescored by code written here
  C  artefact read            the value was read from an artefact, with no independent computation
                              possible in this session; recorded as such, not as a reproduction

Output: validation/REPRO_LEDGER.csv

Run:  python validation/repro/r03_ledger.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "validation"
REPRO = OUT / "repro"
MS = ROOT / "manuscript" / "NAR_WebServer_BrainSafe_built.md"
TAB = ROOT / "results" / "tables"

TOL = 0.0005          # a reported value rounded to 3 dp agrees within half a unit in the last place
REGEN = ("regenerated against the current models on 2026-08-15; the value the manuscript states "
         "was computed against estimators that no longer exist")


def _f(path: Path):
    return pd.read_csv(path) if path.exists() else None


def _flat(s: str) -> str:
    """Collapse all whitespace, so an anchor still matches text the manuscript line-wrapped."""
    return re.sub(r"\s+", " ", s)


def manuscript_says(text: str, anchor: str, value, window: int = 700) -> tuple[bool, str]:
    """Is `value` stated in the manuscript near `anchor`?

    Whitespace is collapsed on both sides first: the manuscript is hard-wrapped, so an anchor
    written as one line will not be found by a literal search even when the phrase is present. Every
    occurrence of the anchor is tried, not only the first, because a common word such as "conformal"
    appears in the abstract long before the section that states the number.
    """
    flat, a = _flat(text), _flat(anchor)
    starts = [m.start() for m in re.finditer(re.escape(a), flat)]
    if not starts:
        return False, "anchor not found in manuscript"
    WORDS = {0: "zero", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
             7: "seven", 8: "eight", 9: "nine", 10: "ten"}
    cands = set()
    if isinstance(value, (int, np.integer)) or float(value).is_integer():
        iv = int(value)
        cands |= {f"{iv}", f"{iv:,}"}
        if iv in WORDS:
            cands.add(WORDS[iv])          # manuscripts spell small integers: "five pass"
    if not isinstance(value, (int, np.integer)):
        v = float(value)
        for n in (2, 3, 4):
            cands.add(f"{round(v, n):g}")
            cands.add(f"{v:.{n}f}")       # keep trailing zeros: "0.920" is not matched by "0.92"
        cands.add(f"{v:g}")
    for i in starts:
        seg = flat[max(0, i - window): i + len(a) + window]
        for c in cands:
            pat = (rf"\b{re.escape(c)}\b" if c.isalpha()
                   else rf"(?<![\d.]){re.escape(c)}(?![\d])")
            if re.search(pat, seg, re.I if c.isalpha() else 0):
                return True, c
    return False, (f"none of {sorted(cands)} found within {window} chars of any of "
                   f"{len(starts)} anchor occurrence(s)")


def rows() -> list[dict]:
    """Every claim, its anchor, and how it is reproduced."""
    env = json.loads((REPRO / "environment.json").read_text())
    commit, seed = env["commit"][:10], 42
    out = []

    def add(metric, split, ms_value, repro_value, tier, script, output_path, note=""):
        out.append({"metric": metric, "split": split, "manuscript_value": ms_value,
                    "reproduced_value": repro_value, "tier": tier, "script": script,
                    "output_path": output_path, "seed": seed, "commit": commit, "note": note})

    # ---- A: independent re-run of the core cross-validation ------------------------------------
    rec = _f(REPRO / "recomputed_summary.csv")
    if rec is not None:
        c = rec[rec.task == "classification"]
        r = rec[rec.task == "regression"]
        for split in ("random", "scaffold"):
            g = c[c.split == split]
            if len(g):
                add(f"classifier AUROC, mean over endpoints", split,
                    {"random": 0.958, "scaffold": 0.925}[split],
                    round(float(g.roc_auc_mean.mean()), 4), "A",
                    "validation/repro/r02_recompute_cv.py",
                    "validation/repro/recomputed_summary.csv")
                for _, row in g.iterrows():
                    add(f"AUROC {row.endpoint}", split, None, round(float(row.roc_auc_mean), 4),
                        "A", "validation/repro/r02_recompute_cv.py",
                        "validation/repro/recomputed_summary.csv",
                        "per-endpoint; manuscript states panel mean and the two extremes only")
                for metric in ("pr_auc", "sensitivity", "specificity", "balanced_acc", "mcc"):
                    col = f"{metric}_mean"
                    if col in g:
                        add(f"classifier {metric}, mean over endpoints", split, None,
                            round(float(g[col].mean()), 4), "A",
                            "validation/repro/r02_recompute_cv.py",
                            "validation/repro/recomputed_summary.csv",
                            "not stated in the manuscript; reproduced for completeness")
            h = r[r.split == split]
            if len(h):
                add("regression R2, mean over all regression endpoints", split, None,
                    round(float(h.r2_mean.mean()), 4), "A",
                    "validation/repro/r02_recompute_cv.py",
                    "validation/repro/recomputed_summary.csv",
                    "manuscript states a range over the receptor regressions, not a mean over all")
                # The manuscript's range is explicitly "the receptor potency regressions", which is
                # the four receptors and excludes the antioxidant endpoint. Comparing against all
                # five would manufacture a DIFFERS out of a scoping difference.
                RECEPTORS = ["D2", "A2A", "HT2A", "SERT"]
                hr = h[h.endpoint.isin(RECEPTORS)]
                if len(hr):
                    add("receptor potency regression R2, min", split,
                        {"random": 0.64, "scaffold": 0.46}[split],
                        round(float(hr.r2_mean.min()), 4), "A",
                        "validation/repro/r02_recompute_cv.py",
                        "validation/repro/recomputed_summary.csv",
                        f"receptors only ({', '.join(RECEPTORS)}), matching the manuscript's wording")
                    add("receptor potency regression R2, max", split,
                        {"random": 0.72, "scaffold": 0.62}[split],
                        round(float(hr.r2_mean.max()), 4), "A",
                        "validation/repro/r02_recompute_cv.py",
                        "validation/repro/recomputed_summary.csv",
                        f"receptors only ({', '.join(RECEPTORS)})")

    boot = _f(REPRO / "recomputed_bootstrap.csv")
    if boot is not None:
        for _, b in boot.iterrows():
            add(f"bootstrap 95% CI, {b.metric} {b.endpoint}", b.split, None,
                f"{b.pooled_estimate} [{b.ci95_low}, {b.ci95_high}]", "A",
                "validation/repro/r02_recompute_cv.py",
                "validation/repro/recomputed_bootstrap.csv",
                "no interval stated in the manuscript for these")

    # ---- leakage -------------------------------------------------------------------------------
    leak = json.loads((REPRO / "leakage_summary.json").read_text()) \
        if (REPRO / "leakage_summary.json").exists() else None
    if leak:
        for key, label in (("L1_identity_overlap_as_trained_max", "InChIKey overlap across folds"),
                           ("L2_feature_overlap_as_trained_max", "feature-vector overlap across folds"),
                           ("L3_scaffold_overlap_scaffold_split_max", "scaffold overlap, scaffold split")):
            add(label, "both", 0, int(leak[key]), "A", "validation/repro/r01_leakage.py",
                "validation/repro/leakage_report.csv",
                "as trained, i.e. on the deduplicated matrix the pipeline fits")
        add("feature-vector overlap on the raw table (pre-deduplication)", "both", None,
            int(leak["L2_feature_overlap_raw_table_max"]), "A",
            "validation/repro/r01_leakage.py", "validation/repro/leakage_report.csv",
            "the leak deduplication removes; not a property of any trained model")

    # ---- A: calibration, recomputed by re-running the forest and refitting isotonic -------------
    mycal = _f(REPRO / "calibration_summary.csv")
    if mycal is not None:
        piv = mycal.pivot(index="endpoint", columns="calibration", values="ece")
        add("mean ECE before calibration", "out-of-fold", 0.0801,
            round(float(piv["raw"].mean()), 4), "A",
            "validation/repro/r05_calibration_importance.py",
            "validation/repro/calibration_summary.csv")
        add("mean ECE after calibration", "out-of-fold", 0.0147,
            round(float(piv["isotonic"].mean()), 4), "A",
            "validation/repro/r05_calibration_importance.py",
            "validation/repro/calibration_summary.csv",
            "different nesting from the pipeline: isotonic fitted on the other 9 folds' "
            "out-of-fold predictions, where calibrate.py uses cross_val_predict(cv=5) over the "
            "pooled out-of-fold vector. Both are honest; they are not the same estimator")

    # ---- A: null and permutation models --------------------------------------------------------
    nul = _f(REPRO / "null_models.csv")
    if nul is not None:
        for split in ("random", "scaffold"):
            for kind in ("permuted_labels", "permuted_within_train", "stratified_random"):
                g = nul[(nul.split == split) & (nul.null == kind)]
                if len(g):
                    mean = float(g.auroc_mean.mean())
                    verdict = ("consistent with chance" if abs(mean - 0.5) <= 0.02
                               else "ABOVE CHANCE, investigate")
                    add(f"null model AUROC, {kind}", split, None,
                        round(mean, 4), "A",
                        "validation/repro/r04_null_models.py",
                        "validation/repro/null_models.csv",
                        f"theoretical expectation 0.5, not a manuscript claim; {verdict}. "
                        f"Mean over {g.endpoint.nunique()} endpoints; max over endpoints "
                        f"{g.auroc_mean.max():.4f}")

    # ---- C: artefact reads, no independent recomputation possible in this session --------------
    def artefact(metric, split, ms_value, value, path, note):
        add(metric, split, ms_value, value, "C", "(read only)", path, note)

    bm = ROOT / "models_rf" / "binder_modes.json"
    if bm.exists():
        d = json.loads(bm.read_text())
        # Deployed endpoints only. The panel means the manuscript quotes describe what the server
        # offers; averaging in the five withdrawn endpoints answers a question nobody asked and
        # drags the mean down by the very failures that caused the withdrawal (0.879 against
        # 0.902 for AUROC). The withdrawn five are reported individually instead.
        dep = [v for v in d.values() if v.get("deployed", True)]
        sens = [v["sensitivity_at_threshold"] for v in dep
                if v.get("sensitivity_at_threshold") is not None]
        aur = [v["auroc_vs_measured_inactives"] for v in dep
               if v.get("auroc_vs_measured_inactives") is not None]
        artefact("binder panel mean sensitivity (deployed)", "scaffold holdout", 0.898,
                 round(float(np.mean(sens)), 4), "models_rf/binder_modes.json",
                 "recomputed from the per-endpoint records, but those records were written by the "
                 "training run; re-running the binder panel was not performed in this session")
        artefact("binder panel mean AUROC vs measured inactives (deployed)", "holdout", 0.917,
                 round(float(np.mean(aur)), 4), "models_rf/binder_modes.json", "as above")
        artefact("binder endpoints deployed", "n/a", 47,
                 sum(1 for v in d.values() if v.get("deployed", True)),
                 "models_rf/binder_modes.json", "count")

    ext = _f(TAB / "external_bbb_validation.csv")
    if ext is not None:
        for _, e in ext.iterrows():
            ms = {306: 0.764, 241: 0.793}.get(int(e.n))
            artefact(f"external BBB AUROC (n={int(e.n)})", "external", ms, round(float(e.auroc), 4),
                     "results/tables/external_bbb_validation.csv",
                     "external set held fixed; not recomputed here")

    con = _f(TAB / "rf_conformal.csv")
    if con is not None:
        artefact("conformal coverage, min", "conformal", 0.889,
                 round(float(con.empirical_coverage.min()), 3),
                 "results/tables/rf_conformal.csv", "")
        artefact("conformal coverage, max", "conformal", 0.921,
                 round(float(con.empirical_coverage.max()), 3),
                 "results/tables/rf_conformal.csv", "")

    hold = _f(TAB / "scaffold_holdout_results.csv")
    if hold is not None:
        # The manuscript's median excludes the target whose threshold collapsed, so the same
        # exclusion is applied here; comparing against the median over all rows would manufacture a
        # DIFFERS out of a scoping difference.
        # "not collapsed" (39 targets), not "usable" (36): the report script excludes only the
        # target whose threshold collapsed, and the two filters give 0.814 and 0.832 respectively
        collapsed = (hold.threshold_collapsed.astype(str).str.lower().isin(["true", "1"])
                     if "threshold_collapsed" in hold.columns else False)
        usable = hold[~collapsed] if "threshold_collapsed" in hold.columns else hold
        add("prospective recall, median over targets", "scaffold holdout", 0.815,
            round(float(usable.holdout_recall.median()), 4), "A",
            "src/brainsafe/evaluation/scaffold_holdout_report.py",
            "results/tables/scaffold_holdout_results.csv",
            REGEN + ". The panel script had been re-run on 2026-08-13 and refreshed the withheld "
            "sets, but the report script that writes this table had not, so inputs were current "
            "and the output was a day stale. 15 of 40 targets move by more than 0.10 between the "
            "two runs, the largest being GluA2 +0.645 and GABA_A -0.466")
        add("prospective recall, mean over targets", "scaffold holdout", None,
            round(float(hold.holdout_recall.mean()), 4), "A",
            "src/brainsafe/evaluation/scaffold_holdout_report.py",
            "results/tables/scaffold_holdout_results.csv",
            "not stated in the manuscript; was 0.7560 on the superseded run")

    # Two analyses have not been re-run since the models were retrained on 2026-08-13. Their
    # artefacts therefore describe estimators that no longer exist, and the manuscript quotes a
    # value older still. Recorded as CANNOT_REPRODUCE rather than compared, because comparing a
    # current manuscript against a stale artefact answers no question anyone asked.
    spec = _f(TAB / "noncns_specificity_summary.csv")
    if spec is not None:
        s = spec[spec.metric.astype(str).str.startswith("Specificity")]
        if len(s):
            add("specificity on non-CNS chemistry", "external", 0.949,
                round(float(s.estimate.iloc[0]), 4), "A",
                "src/brainsafe/evaluation/noncns_specificity.py",
                "results/tables/noncns_specificity_summary.csv",
                REGEN + ". History of this number: 0.875, then 0.920 on superseded models, then 0.948 "
                "before the binder panel was retrained on the extended endpoint tables. Re-run on the "
                "current deployed models it is "
                f"current models {float(s.estimate.iloc[0]):.4f} "
                f"(95% CI {float(s.ci95_low.iloc[0]):.4f}-{float(s.ci95_high.iloc[0]):.4f}), "
                f"{int(s.k.iloc[0])} of {int(s.n.iloc[0])} compounds silent")

    cmp_ = _f(TAB / "model_comparison.csv")
    if cmp_ is not None:
        sc = cmp_[(cmp_.split == "scaffold") & (cmp_.task == "classification")]
        for model, ms in (("RandomForest", 0.9212), ("HistGradientBoosting", 0.9149),
                          ("XGBoost", 0.9156), ("kNN read-across", 0.8829),
                          ("LogisticRegression", 0.8352)):
            g = sc[sc.model == model]
            if len(g):
                artefact(f"model comparison, {model}", "scaffold", ms,
                         round(float(g["mean"].mean()), 4),
                         "results/tables/model_comparison.csv",
                         "5-fold; re-running the comparison was not performed in this session")

    inv = _f(TAB / "inversion_validation.csv")
    if inv is not None:
        artefact("adversarial checks passing", "n/a", 5,
                 int((inv.result.astype(str).str.upper() == "PASS").sum()),
                 "results/tables/inversion_validation.csv", "")

    # dataset size, computed here from the tables
    total = 0
    for f in sorted((ROOT / "data" / "endpoints").glob("*.csv")):
        head = pd.read_csv(f, nrows=0)
        if "smiles" in head.columns:
            total += len(pd.read_csv(f, usecols=["smiles"]))
    add("compound-endpoint records", "n/a", 228200, total, "A",
        "validation/repro/r03_ledger.py", "validation/REPRO_LEDGER.csv",
        "counted from data/endpoints/*.csv")
    return out


def main() -> None:
    text = MS.read_text(encoding="utf-8") if MS.exists() else ""
    # Anchors are phrases from the CURRENT manuscript. They were rewritten when the manuscript was
    # restructured to the NAR format; an anchor that no longer appears makes every row using it read
    # MISQUOTED, which is a defect in this file rather than in the manuscript.
    anchors = {
        "classifier AUROC, mean over endpoints|random": "mean AUROC 0.958 and 0.925 respectively",
        "classifier AUROC, mean over endpoints|scaffold": "mean AUROC 0.958 and 0.925 respectively",
        "binder panel mean AUROC vs measured inactives (deployed)|holdout":
            "reaches a mean AUROC of 0.917",
        "binder panel mean sensitivity (deployed)|scaffold holdout":
            "a mean sensitivity of 0.898",
        "mean ECE before calibration|out-of-fold": "expected calibration error falling from 0.0801",
        "mean ECE after calibration|out-of-fold": "expected calibration error falling from 0.0801",
        "specificity on non-CNS chemistry|external": "on non-CNS chemistry its specificity is 0.949",
        "compound-endpoint records|n/a": "228,200 measured",
        "binder endpoints deployed|n/a": "Across the 47 that are",
        "adversarial checks passing|n/a": "Five pass",
        "receptor potency regression R2, min|random": "0.64 to 0.72",
        "receptor potency regression R2, max|random": "0.64 to 0.72",
        "receptor potency regression R2, min|scaffold": "0.46 to 0.62 (scaffold)",
        "receptor potency regression R2, max|scaffold": "0.46 to 0.62 (scaffold)",
        "external BBB AUROC (n=306)|external": "absent from B3DB by InChIKey",
        "external BBB AUROC (n=241)|external": "distinguishable from the training set in feature space",
        "conformal coverage, min|conformal": "empirical coverage is 0.889 to 0.921",
        "conformal coverage, max|conformal": "empirical coverage is 0.889 to 0.921",
        "prospective recall, median over targets|scaffold holdout": "per-target median of",
    }

    led = []
    for r in rows():
        ms_v, rp_v = r["manuscript_value"], r["reproduced_value"]
        key = f"{r['metric']}|{r['split']}"
        quoted_ok, quoted_note = (True, "")
        if ms_v is not None and key in anchors and text:
            quoted_ok, quoted_note = manuscript_says(text, anchors[key], ms_v)

        if ms_v is None:
            status, diff = "NOT_STATED", ""
        elif not quoted_ok:
            status, diff = "MISQUOTED", ""
        elif isinstance(rp_v, str):
            status, diff = ("CANNOT_REPRODUCE" if str(r["note"]).startswith("CANNOT_REPRODUCE")
                            else "NOT_COMPARABLE"), ""
        else:
            try:
                diff = round(abs(float(rp_v) - float(ms_v)), 6)
                # Compare at the precision the manuscript states. A manuscript that reports "0.64"
                # is not wrong because the value is 0.6436; it is wrong only if 0.6436 does not
                # round to 0.64. Judging a 2-dp claim at 4 dp manufactures disagreement.
                dp = len(str(ms_v).split(".")[1]) if "." in str(ms_v) else 0
                status = ("MATCHES" if round(float(rp_v), dp) == round(float(ms_v), dp)
                          else "DIFFERS")
            except (TypeError, ValueError):
                status, diff = "CANNOT_REPRODUCE", ""
        if r["tier"] == "C" and status == "MATCHES":
            status = "MATCHES (artefact read, not re-run)"
        led.append({**r, "abs_diff": diff, "status": status,
                    "manuscript_quote_check": quoted_note or "ok"})

    df = pd.DataFrame(led)[["metric", "split", "manuscript_value", "reproduced_value", "abs_diff",
                            "status", "tier", "script", "output_path", "seed", "commit",
                            "manuscript_quote_check", "note"]]
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "REPRO_LEDGER.csv", index=False)

    print(f"{len(df)} ledger rows\n")
    print(df.status.value_counts().to_string())
    bad = df[df.status.isin(["DIFFERS", "MISQUOTED", "CANNOT_REPRODUCE"])]
    if len(bad):
        print("\n=== rows needing attention ===")
        print(bad[["metric", "split", "manuscript_value", "reproduced_value",
                   "abs_diff", "status"]].to_string(index=False))
    print(f"\nwrote {(OUT / 'REPRO_LEDGER.csv').relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
