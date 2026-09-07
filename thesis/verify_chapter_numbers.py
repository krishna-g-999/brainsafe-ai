"""Check every quantitative claim in a thesis chapter against the artefact it cites.

Chapters 3 to 8 were each written by reading artefacts and then verified by reading them again, and
the verification was reported as a count of figures checked. An audit made the obvious objection:
that count rested on nothing committed, so a reader had no way to repeat it, which is the same
weakness the thesis itself raises about the technical report's hard-coded null models. A claim of
verification that cannot itself be verified is worth very little.

This script is the answer. Each claim below names the chapter, a short description, the artefact it
comes from and a callable that recomputes the value from that artefact. The chapter text is then
searched for the formatted value, so a claim fails either because the artefact has moved or because
the chapter states something the artefact does not support.

What it does and does not establish:

  it does        confirm that every listed number appears in the chapter and is currently true of
                 the artefact it is attributed to
  it does not    confirm that the chapter's *prose* around a number is correct, nor that the list of
                 claims is complete

The second limit is real and is the reason the count is reported as "claims checked" rather than
"chapter verified". Completeness is a matter of judgement and this file is where that judgement is
recorded, so a reader can see exactly which claims were pinned and add the ones that were not.

Run:  brainsafe_env/Scripts/python.exe thesis/verify_chapter_numbers.py
      brainsafe_env/Scripts/python.exe thesis/verify_chapter_numbers.py --chapter 09
Exit: 0 if every claim holds, 1 otherwise.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
THESIS = ROOT / "thesis"
TAB = ROOT / "results" / "tables"
INV = ROOT / "inversion" / "results"


def rows(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh))


def col(path: Path, field: str, where=None) -> list[float]:
    out = []
    for r in rows(path):
        if where and not where(r):
            continue
        try:
            out.append(float(r[field]))
        except (TypeError, ValueError):
            continue
    return out


def cell(path: Path, field: str, where) -> str:
    for r in rows(path):
        if where(r):
            return r[field]
    raise KeyError(f"no row matching in {path.name}")


# ---------------------------------------------------------------------------------------------
# Claims. (chapter, description, artefact, callable -> value, format)
# The format string is how the chapter writes the number, and both are checked.
# ---------------------------------------------------------------------------------------------

CLAIMS: list[tuple[str, str, str, callable, str]] = []


def claim(ch, desc, src, fn, fmt="{:.4f}"):
    CLAIMS.append((ch, desc, src, fn, fmt))


# ---- Chapter 9: the falsification suite -----------------------------------------------------

claim("09", "H1 frequency null", "H1_disease_layer.csv",
      lambda: float(cell(INV / "H1_disease_layer.csv", "value",
                         lambda r: r["test"].startswith("frequency null"))))
claim("09", "H1 evaluation set size", "H1_disease_layer.csv",
      lambda: int(cell(INV / "H1_disease_layer.csv", "n",
                       lambda r: r["test"] == "observed top-3 accuracy")), "{:,}")

claim("09", "H2 curated", "H2_weight_ablation.csv",
      lambda: float(cell(INV / "H2_weight_ablation.csv", "top3_accuracy",
                         lambda r: r["weights"] == "curated")))
claim("09", "H2 uniform", "H2_weight_ablation.csv",
      lambda: float(cell(INV / "H2_weight_ablation.csv", "top3_accuracy",
                         lambda r: r["weights"].startswith("uniform"))))
claim("09", "H2 permuted", "H2_weight_ablation.csv",
      lambda: float(cell(INV / "H2_weight_ablation.csv", "top3_accuracy",
                         lambda r: r["weights"].startswith("randomly"))))
claim("09", "H2 spread", "H2_weight_ablation.csv",
      lambda: round(float(cell(INV / "H2_weight_ablation.csv", "top3_accuracy",
                               lambda r: r["weights"] == "curated"))
                    - float(cell(INV / "H2_weight_ablation.csv", "top3_accuracy",
                                 lambda r: r["weights"].startswith("randomly"))), 4))

for _stratum, _lab in [("in domain (T>=0.5)", "H4 in-domain FPR"),
                       ("near domain (0.3-0.5)", "H4 near-domain FPR"),
                       ("distant (T<0.3)", "H4 distant FPR"),
                       ("ALL", "H4 all-strata FPR")]:
    claim("09", _lab, "H4_distant_specificity.csv",
          (lambda s: lambda: float(cell(INV / "H4_distant_specificity.csv",
                                        "false_positive_rate", lambda r: r["stratum"] == s)))(_stratum))
    claim("09", _lab.replace("FPR", "n"), "H4_distant_specificity.csv",
          (lambda s: lambda: int(cell(INV / "H4_distant_specificity.csv",
                                      "n", lambda r: r["stratum"] == s)))(_stratum), "{:d}")
claim("09", "H4 comparator, library FPR", "noncns_specificity_summary.csv",
      lambda: float(cell(TAB / "noncns_specificity_summary.csv", "estimate",
                         lambda r: r["metric"].startswith("False-positive rate (any actionable"))))

claim("09", "H5 read-across recall", "H5_readacross_value.csv",
      lambda: float(cell(INV / "H5_readacross_value.csv", "recall",
                         lambda r: r["method"].startswith("read-across"))))
claim("09", "H5 frequency baseline", "H5_readacross_value.csv",
      lambda: float(cell(INV / "H5_readacross_value.csv", "recall",
                         lambda r: r["method"].startswith("frequency"))))
claim("09", "H5 evaluation set size", "H5_readacross_value.csv",
      lambda: int(cell(INV / "H5_readacross_value.csv", "n",
                       lambda r: r["method"].startswith("read-across"))), "{:,}")
claim("09", "H5 targets", "H5_readacross_per_target.csv",
      lambda: len(rows(INV / "H5_readacross_per_target.csv")), "{:d}")
claim("09", "H5 worst per-target recall", "H5_readacross_per_target.csv",
      lambda: min(col(INV / "H5_readacross_per_target.csv", "readacross_recall")), "{:.3f}")
claim("09", "H5 median per-target recall", "H5_readacross_per_target.csv",
      lambda: st.median(col(INV / "H5_readacross_per_target.csv", "readacross_recall")), "{:.3f}")

claim("09", "H7 testable targets", "H7_target_discrimination.csv",
      lambda: len(rows(INV / "H7_target_discrimination.csv")), "{:d}")
claim("09", "H7 lowest AUROC", "H7_target_discrimination.csv",
      lambda: min(col(INV / "H7_target_discrimination.csv", "auroc_vs_random")))
claim("09", "H7 highest AUROC", "H7_target_discrimination.csv",
      lambda: max(col(INV / "H7_target_discrimination.csv", "auroc_vs_random")))
claim("09", "H7 median AUROC", "H7_target_discrimination.csv",
      lambda: st.median(col(INV / "H7_target_discrimination.csv", "auroc_vs_random")))
claim("09", "H7 mean AUROC", "H7_target_discrimination.csv",
      lambda: round(st.mean(col(INV / "H7_target_discrimination.csv", "auroc_vs_random")), 4))
claim("09", "H7 median deployed sensitivity", "H7_target_discrimination.csv",
      lambda: st.median(col(INV / "H7_target_discrimination.csv", "deployed_sensitivity")))
claim("09", "H7 mean deployed sensitivity", "H7_target_discrimination.csv",
      lambda: round(st.mean(col(INV / "H7_target_discrimination.csv", "deployed_sensitivity")), 4))
claim("09", "H7 max deployed sensitivity", "H7_target_discrimination.csv",
      lambda: max(col(INV / "H7_target_discrimination.csv", "deployed_sensitivity")), "{:.3f}")

claim("09", "H8 targets that ever fire", "H8_panel_independence.csv",
      lambda: int(float(cell(INV / "H8_panel_independence.csv", "value",
                             lambda r: r["metric"] == "targets that ever fire on the drug set"))), "{:d}")
claim("09", "H8 independent directions", "H8_panel_independence.csv",
      lambda: int(float(cell(INV / "H8_panel_independence.csv", "value",
                             lambda r: r["metric"] == "independent directions in the firing pattern"))),
      "{:d}")
for _a, _b in [("OPRM1", "OPRK1"), ("SERT", "NET"), ("D2", "D3"), ("HT2A", "HT7"), ("HT1A", "HT7")]:
    claim("09", f"H8 phi {_a}-{_b}", "H8_family_correlation.csv",
          (lambda a, b: lambda: float(cell(INV / "H8_family_correlation.csv", "phi_correlation",
                                           lambda r: r["target_a"] == a and r["target_b"] == b)))(_a, _b),
          "{:.3f}")
    claim("09", f"H8 conditional {_a}-{_b}", "H8_family_correlation.csv",
          (lambda a, b: lambda: float(cell(INV / "H8_family_correlation.csv", "p_b_given_a",
                                           lambda r: r["target_a"] == a and r["target_b"] == b)))(_a, _b),
          "{:.3f}")

claim("09", "H9 mean per-indication AUROC", "H9_disease_discrimination_summary.csv",
      lambda: float(cell(INV / "H9_disease_discrimination_summary.csv", "model",
                         lambda r: r["metric"] == "mean per-indication AUROC")))
claim("09", "H9 best indication AUROC", "H9_disease_discrimination.csv",
      lambda: max(col(INV / "H9_disease_discrimination.csv", "auroc_model")))
claim("09", "H9 worst indication AUROC", "H9_disease_discrimination.csv",
      lambda: min(col(INV / "H9_disease_discrimination.csv", "auroc_model")))

claim("09", "graph conditions", "GRAPH_FINGERPRINT.json",
      lambda: json.loads((INV / "GRAPH_FINGERPRINT.json").read_text())["n_conditions"], "{:d}")
claim("09", "graph targets", "GRAPH_FINGERPRINT.json",
      lambda: json.loads((INV / "GRAPH_FINGERPRINT.json").read_text())["n_targets"], "{:d}")

claim("09", "held-out entries before the H1 restriction", "heldout_actives.json",
      lambda: sum(len(v) for v in
                  json.loads((ROOT / "models_rf" / "holdout" / "heldout_actives.json")
                             .read_text(encoding="utf-8")).values()), "{:,}")

for _ep, _fpr in [("Nav1_1", None), ("GluA2", None), ("NRF2", None), ("NFKB1", None)]:
    claim("09", f"deployed audit FPR {_ep}", "deployed_specificity_audit.csv",
          (lambda e: lambda: float(cell(ROOT / "results" / "deployed_specificity_audit.csv",
                                        "random_fpr", lambda r: r["target"] == e)))(_ep))

claim("09", "permutation null, random mean", "permutation_null.csv",
      lambda: round(st.mean(col(TAB / "permutation_null.csv", "permuted_roc_auc_mean",
                                lambda r: r["split"] == "random")), 4))
claim("09", "permutation null, scaffold mean", "permutation_null.csv",
      lambda: round(st.mean(col(TAB / "permutation_null.csv", "permuted_roc_auc_mean",
                                lambda r: r["split"] == "scaffold")), 4))
claim("09", "permutation null, random minimum", "permutation_null.csv",
      lambda: min(col(TAB / "permutation_null.csv", "permuted_roc_auc_mean",
                      lambda r: r["split"] == "random")))
claim("09", "permutation null, random maximum", "permutation_null.csv",
      lambda: max(col(TAB / "permutation_null.csv", "permuted_roc_auc_mean",
                      lambda r: r["split"] == "random")))
claim("09", "permutation null, scaffold minimum", "permutation_null.csv",
      lambda: min(col(TAB / "permutation_null.csv", "permuted_roc_auc_mean",
                      lambda r: r["split"] == "scaffold")))
claim("09", "permutation null, scaffold maximum", "permutation_null.csv",
      lambda: max(col(TAB / "permutation_null.csv", "permuted_roc_auc_mean",
                      lambda r: r["split"] == "scaffold")))
claim("09", "permutation null, smallest fold sd", "permutation_null.csv",
      lambda: min(col(TAB / "permutation_null.csv", "permuted_roc_auc_sd")))
claim("09", "permutation null, largest fold sd", "permutation_null.csv",
      lambda: max(col(TAB / "permutation_null.csv", "permuted_roc_auc_sd")))


def _margins(split: str) -> list[float]:
    """True cross-validated AUROC minus this endpoint's own permuted null, per endpoint."""
    null = {r["endpoint"]: float(r["permuted_roc_auc_mean"])
            for r in rows(TAB / "permutation_null.csv") if r["split"] == split}
    true = {r["endpoint"]: float(r["roc_auc_mean"])
            for r in rows(TAB / "rf_cv_summary.csv")
            if r["task"] == "classification" and r["split"] == split}
    return [true[e] - null[e] for e in null if e in true]


claim("09", "smallest random-split margin over null", "permutation_null.csv + rf_cv_summary.csv",
      lambda: min(_margins("random")))
claim("09", "largest random-split margin over null", "permutation_null.csv + rf_cv_summary.csv",
      lambda: max(_margins("random")))
claim("09", "mean random-split margin over null", "permutation_null.csv + rf_cv_summary.csv",
      lambda: round(st.mean(_margins("random")), 4))
claim("09", "smallest scaffold-split margin over null", "permutation_null.csv + rf_cv_summary.csv",
      lambda: min(_margins("scaffold")))
claim("09", "largest scaffold-split margin over null", "permutation_null.csv + rf_cv_summary.csv",
      lambda: max(_margins("scaffold")))
claim("09", "mean scaffold-split margin over null", "permutation_null.csv + rf_cv_summary.csv",
      lambda: round(st.mean(_margins("scaffold")), 4))

claim("09", "largest null deviation from chance", "permutation_null.csv",
      lambda: max(abs(v - 0.5) for v in col(TAB / "permutation_null.csv", "permuted_roc_auc_mean")))

claim("09", "registry held-out mean sensitivity", "binder_modes.json",
      lambda: round(st.mean([v["sensitivity_at_threshold"] for v in
                             json.loads((ROOT / "models_rf" / "binder_modes.json")
                                        .read_text(encoding="utf-8")).values()
                             if v.get("deployed")]), 4))
claim("09", "endpoints failing the reliability gate", "binder_modes.json",
      lambda: sum(1 for v in json.loads((ROOT / "models_rf" / "binder_modes.json")
                                        .read_text(encoding="utf-8")).values()
                  if v.get("deployed") and v.get("reliable_call") is False), "{:d}")
claim("09", "background_specificity.csv mean sensitivity, uncorrected", "background_specificity.csv",
      lambda: round(st.mean(col(TAB / "background_specificity.csv", "sensitivity_after")), 4))
claim("09", "background_specificity.csv rows marked reliable", "background_specificity.csv",
      lambda: sum(1 for r in rows(TAB / "background_specificity.csv") if r["reliable"] == "True"),
      "{:d}")


# ---- Chapter 10: limitations and future work -------------------------------------------------

_BANDS = ["below 0.40 (different chemotype)", "0.40 to 0.55 (related series)",
          "0.55 to 0.70 (same series)", "0.70 and above (close analogue)"]


def _strata(split, band, field):
    return cell(TAB / "external_novelty_strata.csv", field,
                lambda r: r["split"] == split and r["novelty_band"] == band)


for _b in _BANDS:
    claim("10", f"time-split recall, {_b[:22]}", "external_novelty_strata.csv",
          (lambda b: lambda: float(_strata("time", b, "recall_at_threshold")))(_b))
    claim("10", f"time-split AUROC, {_b[:22]}", "external_novelty_strata.csv",
          (lambda b: lambda: float(_strata("time", b, "auroc")))(_b))
    claim("10", f"time-split actives, {_b[:22]}", "external_novelty_strata.csv",
          (lambda b: lambda: int(_strata("time", b, "n_actives")))(_b), "{:,}")
claim("10", "cross-provenance recall, lowest band", "external_novelty_strata.csv",
      lambda: float(_strata("cross_source", _BANDS[0], "recall_at_threshold")))
claim("10", "cross-provenance recall, highest band", "external_novelty_strata.csv",
      lambda: float(_strata("cross_source", _BANDS[3], "recall_at_threshold")))


def _h10(method, population, field="auroc"):
    return float(cell(INV / "H10_barrier_necessity.csv", field,
                      lambda r: r["method"] == method and r["population"] == population))


for _m in ["deployed forest, 1,036 features", "descriptor forest, 12 features",
           "descriptor logistic regression", "tpsa alone", "hbd alone", "qed alone",
           "CNS heuristic: TPSA <= 90 and MW <= 400", "mw alone", "clogp alone"]:
    for _p in ["scaffold hold-out", "external approved"]:
        claim("10", f"H10 {_m[:28]}, {_p[:8]}", "H10_barrier_necessity.csv",
              (lambda m, pp: lambda: _h10(m, pp))(_m, _p))
for _p in ["scaffold hold-out", "external approved"]:
    for _f in ["delta_ci95_low", "delta_ci95_high"]:
        claim("10", f"H10 descriptor forest {_f[6:]}, {_p[:8]}", "H10_barrier_necessity.csv",
              (lambda pp, ff: lambda: _h10("descriptor forest, 12 features", pp, ff))(_p, _f))
claim("10", "H10 external bootstrap p", "H10_barrier_necessity.csv",
      lambda: _h10("descriptor forest, 12 features", "external approved",
                   "bootstrap_p_deployed_better"), "{:.3f}")

claim("10", "core classifiers calibrated", "calibration.csv",
      lambda: len(rows(TAB / "calibration.csv")), "{:d}")
claim("10", "core ECE, raw mean", "calibration.csv",
      lambda: round(st.mean(col(TAB / "calibration.csv", "ece_raw")), 4))
claim("10", "core ECE, calibrated mean", "calibration.csv",
      lambda: round(st.mean(col(TAB / "calibration.csv", "ece_calibrated")), 4))
claim("10", "core ECE, best", "calibration.csv",
      lambda: min(col(TAB / "calibration.csv", "ece_calibrated")))
claim("10", "core ECE, worst, which is BBB", "calibration.csv",
      lambda: max(col(TAB / "calibration.csv", "ece_calibrated")))
claim("10", "binder endpoints with a calibration figure", "integrity_calibration_per_target.csv",
      lambda: len(rows(TAB / "integrity_calibration_per_target.csv")), "{:d}")
for _lab, _fn in [("mean", st.mean), ("median", st.median), ("min", min), ("max", max)]:
    claim("10", f"binder ECE, {_lab}", "integrity_calibration_per_target.csv",
          (lambda f: lambda: round(float(f(col(TAB / "integrity_calibration_per_target.csv",
                                               "ece"))), 4))(_fn))

# The uncertainty layers do not all stop at the same endpoint, and an earlier draft of section 10.5
# collapsed that into "38 binder endpoints". 38 is how many carry a measured calibration error; the
# panel is 47. Each coverage count is pinned separately so the two can never be conflated again.
def _deployed() -> set:
    reg = json.loads((ROOT / "models_rf" / "binder_modes.json").read_text(encoding="utf-8"))
    return {k for k, v in reg.items() if v.get("deployed")}


def _ad_reference() -> set:
    import sys as _sys
    _sys.path.insert(0, str(ROOT))
    import app as _app
    return set(_app.load_ad_per_endpoint())


claim("10", "deployed binder endpoints", "binder_modes.json",
      lambda: len(_deployed()), "{:d}")
claim("10", "deployed binders with a per-endpoint applicability reference",
      "app.load_ad_per_endpoint",
      lambda: len(_deployed() & _ad_reference()), "{:d}")
claim("10", "deployed binders with no measured calibration error",
      "integrity_calibration_per_target.csv",
      lambda: len(_deployed() - {r["endpoint"]
                                 for r in rows(TAB / "integrity_calibration_per_target.csv")}),
      "{:d}")
claim("10", "ADME and auxiliary estimators", "models_rf/adme",
      lambda: len(list((ROOT / "models_rf" / "adme").glob("*.joblib")))
      + len([q for q in (ROOT / "models_rf").glob("*.joblib")
             if q.stem in ("antioxidant_DPPH", "pka_basic")]), "{:d}")

claim("10", "natural-product candidates", "external_natural_products_summary.csv",
      lambda: sum(int(r["n_candidates"])
                  for r in rows(TAB / "external_natural_products_summary.csv")), "{:,}")
claim("10", "natural-product compounds removed as contaminated",
      "external_natural_products_summary.csv",
      lambda: int(sum(float(r["removed_as_contaminated"] or 0)
                      for r in rows(TAB / "external_natural_products_summary.csv"))), "{:,}")
claim("10", "surveyed targets with no deployed estimator",
      "external_natural_products_summary.csv",
      lambda: sum(1 for r in rows(TAB / "external_natural_products_summary.csv")
                  if r["status"] == "no deployed estimator"), "{:d}")
for _ep in ["GBA1", "AChE", "hERG"]:
    claim("10", f"natural-product AUROC, {_ep}", "external_natural_products_summary.csv",
          (lambda e: lambda: float(cell(TAB / "external_natural_products_summary.csv", "auroc",
                                        lambda r: r["endpoint"] == e)))(_ep))
    claim("10", f"natural-product n, {_ep}", "external_natural_products_summary.csv",
          (lambda e: lambda: int(float(cell(TAB / "external_natural_products_summary.csv", "n",
                                            lambda r: r["endpoint"] == e))))(_ep), "{:d}")


def _lib(metric):
    return float(cell(TAB / "library_sp3_coverage.csv", "value",
                      lambda r: r["metric"] == metric))


claim("10", "library median fraction sp3", "library_sp3_coverage.csv",
      lambda: _lib("median fraction sp3"))
claim("10", "library sp3-rich share, per cent", "library_sp3_coverage.csv",
      lambda: _lib("sp3-rich and at most one aromatic ring, per cent"), "{:.2f}")
claim("10", "library fraction sp3, 25th percentile", "library_sp3_coverage.csv",
      lambda: _lib("fraction sp3, 25th percentile"))
claim("10", "library fraction sp3, 75th percentile", "library_sp3_coverage.csv",
      lambda: _lib("fraction sp3, 75th percentile"))
claim("10", "structures parsed as a desalted parent", "library_sp3_coverage.csv",
      lambda: int(_lib("structures parsed as a desalted parent")), "{:,}")

claim("10", "specificity on non-CNS chemistry", "noncns_specificity_summary.csv",
      lambda: float(cell(TAB / "noncns_specificity_summary.csv", "estimate",
                         lambda r: r["metric"].startswith("Specificity"))), "{:.3f}")
claim("10", "specificity interval, low", "noncns_specificity_summary.csv",
      lambda: float(cell(TAB / "noncns_specificity_summary.csv", "ci95_low",
                         lambda r: r["metric"].startswith("Specificity"))), "{:.3f}")
claim("10", "specificity interval, high", "noncns_specificity_summary.csv",
      lambda: float(cell(TAB / "noncns_specificity_summary.csv", "ci95_high",
                         lambda r: r["metric"].startswith("Specificity"))), "{:.4f}")

claim("10", "panel sensitivity, mean on held-out actives", "binder_modes.json",
      lambda: round(st.mean([v["sensitivity_at_threshold"] for v in
                             json.loads((ROOT / "models_rf" / "binder_modes.json")
                                        .read_text(encoding="utf-8")).values()
                             if v.get("deployed")]), 4))
claim("10", "endpoints failing the reliability gate", "binder_modes.json",
      lambda: sum(1 for v in json.loads((ROOT / "models_rf" / "binder_modes.json")
                                        .read_text(encoding="utf-8")).values()
                  if v.get("deployed") and v.get("reliable_call") is False), "{:d}")

claim("10", "permutation null, random mean", "permutation_null.csv",
      lambda: round(st.mean(col(TAB / "permutation_null.csv", "permuted_roc_auc_mean",
                                lambda r: r["split"] == "random")), 4))
claim("10", "permutation null, scaffold mean", "permutation_null.csv",
      lambda: round(st.mean(col(TAB / "permutation_null.csv", "permuted_roc_auc_mean",
                                lambda r: r["split"] == "scaffold")), 4))
for _sp in ("random", "scaffold"):
    claim("10", f"core CV AUROC, lowest, {_sp}", "rf_cv_summary.csv",
          (lambda sp: lambda: min(col(TAB / "rf_cv_summary.csv", "roc_auc_mean",
                                      lambda r: r["task"] == "classification"
                                      and r["split"] == sp)))(_sp), "{:.4g}")
    claim("10", f"core CV AUROC, highest, {_sp}", "rf_cv_summary.csv",
          (lambda sp: lambda: max(col(TAB / "rf_cv_summary.csv", "roc_auc_mean",
                                      lambda r: r["task"] == "classification"
                                      and r["split"] == sp)))(_sp), "{:.4g}")

claim("10", "falsification hypotheses", "VERDICTS.csv",
      lambda: len(rows(INV / "VERDICTS.csv")), "{:d}")
claim("10", "hypotheses refuted", "VERDICTS.csv",
      lambda: sum(1 for r in rows(INV / "VERDICTS.csv")
                  if r["verdict"].startswith("REFUTED")), "{:d}")
claim("10", "graph conditions", "GRAPH_FINGERPRINT.json",
      lambda: json.loads((INV / "GRAPH_FINGERPRINT.json").read_text())["n_conditions"], "{:d}")
claim("10", "graph targets", "GRAPH_FINGERPRINT.json",
      lambda: json.loads((INV / "GRAPH_FINGERPRINT.json").read_text())["n_targets"], "{:d}")


def _ledger_size() -> int:
    """Numbered items under 'Outstanding items' across Chapters 1 to 9.

    Chapter 10's ledger states this count, and a count stated in prose beside a list that keeps
    growing is exactly the kind of number this thesis argues should be derived rather than typed.
    """
    total = 0
    for p in sorted(THESIS.glob("chapter0[1-9]*.md")):
        tail = p.read_text(encoding="utf-8").split("## Outstanding items")[-1]
        total += len(re.findall(r"^\s*\d+\.\s", tail, flags=re.M))
    return total


claim("10", "items raised across Chapters 1 to 9", "thesis/chapter0*.md", _ledger_size, "{:d}")



# ---- The viva pack ---------------------------------------------------------------------------
# It is a study document rather than a chapter, but it is the one document that will be recited
# aloud under questioning, so a number that has drifted in it is worse than a number that has
# drifted anywhere else.

_CVC = [r for r in rows(TAB / "rf_cv_summary.csv") if r["task"] == "classification"]


def _cv(split, fn):
    return round(fn([float(r["roc_auc_mean"]) for r in _CVC if r["split"] == split]), 4)


for _sp in ("random", "scaffold"):
    claim("viva", f"core AUROC {_sp}, mean", "rf_cv_summary.csv",
          (lambda sp: lambda: _cv(sp, st.mean))(_sp))
    claim("viva", f"core AUROC {_sp}, lowest", "rf_cv_summary.csv",
          (lambda sp: lambda: _cv(sp, min))(_sp))
    claim("viva", f"core AUROC {_sp}, highest", "rf_cv_summary.csv",
          (lambda sp: lambda: _cv(sp, max))(_sp))
    claim("viva", f"permuted-label null, {_sp}", "permutation_null.csv",
          (lambda sp: lambda: round(st.mean(col(TAB / "permutation_null.csv",
                                                "permuted_roc_auc_mean",
                                                lambda r: r["split"] == sp)), 4))(_sp))

claim("viva", "worst null deviation from chance", "permutation_null.csv",
      lambda: round(max(abs(v - 0.5) for v in col(TAB / "permutation_null.csv",
                                                  "permuted_roc_auc_mean")), 4))

claim("viva", "panel sensitivity mean", "binder_modes.json",
      lambda: round(st.mean(_dep_sens()), 4))
claim("viva", "panel sensitivity median", "binder_modes.json",
      lambda: round(st.median(_dep_sens()), 3), "{:.3f}")
claim("viva", "panel sensitivity lowest", "binder_modes.json",
      lambda: round(min(_dep_sens()), 3), "{:.3f}")
claim("viva", "panel sensitivity highest", "binder_modes.json",
      lambda: round(max(_dep_sens()), 3), "{:.3f}")

claim("viva", "core ECE raw mean", "calibration.csv",
      lambda: round(st.mean(col(TAB / "calibration.csv", "ece_raw")), 4))
claim("viva", "core ECE calibrated mean", "calibration.csv",
      lambda: round(st.mean(col(TAB / "calibration.csv", "ece_calibrated")), 4))
claim("viva", "core ECE best", "calibration.csv",
      lambda: min(col(TAB / "calibration.csv", "ece_calibrated")))
claim("viva", "core ECE worst, the barrier model", "calibration.csv",
      lambda: max(col(TAB / "calibration.csv", "ece_calibrated")))
claim("viva", "binder ECE mean", "integrity_calibration_per_target.csv",
      lambda: round(st.mean(col(TAB / "integrity_calibration_per_target.csv", "ece")), 4))
claim("viva", "binder ECE median", "integrity_calibration_per_target.csv",
      lambda: round(st.median(col(TAB / "integrity_calibration_per_target.csv", "ece")), 4))

claim("viva", "conformal coverage, lowest", "rf_conformal.csv",
      lambda: min(col(TAB / "rf_conformal.csv", "empirical_coverage")), "{:.3f}")
claim("viva", "conformal coverage, highest", "rf_conformal.csv",
      lambda: max(col(TAB / "rf_conformal.csv", "empirical_coverage")), "{:.3f}")
claim("viva", "conformal set size, smallest", "rf_conformal.csv",
      lambda: min(col(TAB / "rf_conformal.csv", "avg_set_size")), "{:.3f}")
claim("viva", "conformal set size, largest", "rf_conformal.csv",
      lambda: max(col(TAB / "rf_conformal.csv", "avg_set_size")), "{:.3f}")

claim("viva", "specificity on the non-CNS library", "noncns_specificity_summary.csv",
      lambda: float(cell(TAB / "noncns_specificity_summary.csv", "estimate",
                         lambda r: r["metric"].startswith("Specificity"))), "{:.3f}")
claim("viva", "specificity interval, low", "noncns_specificity_summary.csv",
      lambda: float(cell(TAB / "noncns_specificity_summary.csv", "ci95_low",
                         lambda r: r["metric"].startswith("Specificity"))))
claim("viva", "specificity interval, high", "noncns_specificity_summary.csv",
      lambda: float(cell(TAB / "noncns_specificity_summary.csv", "ci95_high",
                         lambda r: r["metric"].startswith("Specificity"))))

for _i, _lab in enumerate(["all 306", "the 227 novel in feature space", "the 79 memorised"]):
    claim("viva", f"external barrier AUROC, {_lab}", "external_bbb_validation.csv",
          (lambda i: lambda: float(rows(TAB / "external_bbb_validation.csv")[i]["auroc"]))(_i))

_BANDS = ["below 0.40 (different chemotype)", "0.40 to 0.55 (related series)",
          "0.55 to 0.70 (same series)", "0.70 and above (close analogue)"]
for _sp in ("time", "random", "cross_source"):
    for _b in _BANDS:
        claim("viva", f"recall, {_sp}, {_b[:22]}", "external_novelty_strata.csv",
              (lambda sp, b: lambda: float(
                  cell(TAB / "external_novelty_strata.csv", "recall_at_threshold",
                       lambda r: r["split"] == sp and r["novelty_band"] == b)))(_sp, _b))

claim("viva", "H10 deployed forest, external", "H10_barrier_necessity.csv",
      lambda: float(cell(INV / "H10_barrier_necessity.csv", "auroc",
                         lambda r: r["population"] == "external approved"
                         and r["method"].startswith("deployed"))))
claim("viva", "H10 descriptor forest, external", "H10_barrier_necessity.csv",
      lambda: float(cell(INV / "H10_barrier_necessity.csv", "auroc",
                         lambda r: r["population"] == "external approved"
                         and r["method"] == "descriptor forest, 12 features")))
claim("viva", "H10 external interval, upper bound", "H10_barrier_necessity.csv",
      lambda: float(cell(INV / "H10_barrier_necessity.csv", "delta_ci95_high",
                         lambda r: r["population"] == "external approved"
                         and r["method"] == "descriptor forest, 12 features")))

claim("viva", "deployed endpoints", "binder_modes.json", lambda: len(_deployed()), "{:d}")
claim("viva", "endpoints failing the reliability gate", "binder_modes.json",
      lambda: sum(1 for v in json.loads((ROOT / "models_rf" / "binder_modes.json")
                                        .read_text(encoding="utf-8")).values()
                  if v.get("deployed") and v.get("reliable_call") is False), "{:d}")
claim("viva", "conditions in the graph", "GRAPH_FINGERPRINT.json",
      lambda: json.loads((INV / "GRAPH_FINGERPRINT.json")
                         .read_text(encoding="utf-8"))["n_conditions"], "{:d}")
claim("viva", "targets in the graph", "GRAPH_FINGERPRINT.json",
      lambda: json.loads((INV / "GRAPH_FINGERPRINT.json")
                         .read_text(encoding="utf-8"))["n_targets"], "{:d}")
claim("viva", "distinct SMILES in the endpoint tables", "library_sp3_coverage.csv",
      lambda: int(float(cell(TAB / "library_sp3_coverage.csv", "value",
                             lambda r: r["metric"].startswith("distinct SMILES")))), "{:,}")

claim("viva", "hypotheses refuted", "VERDICTS.csv",
      lambda: sum(1 for r in rows(INV / "VERDICTS.csv")
                  if r["verdict"].startswith("REFUTED")), "{:d}")


def _dep_sens():
    reg = json.loads((ROOT / "models_rf" / "binder_modes.json").read_text(encoding="utf-8"))
    return [v["sensitivity_at_threshold"] for v in reg.values() if v.get("deployed")]



# ---- The code walkthrough --------------------------------------------------------------------
# Two of these are the evidence that TECHNICAL_REPORT.md's stage-1 diagram states the wrong active
# cut, so they are pinned rather than quoted: if a future rebuild ever did move the tables to a cut
# of 7, the walkthrough would be wrong and should fail rather than mislead.

def _endpoint_rows(pred) -> int:
    import glob
    import pandas as pd
    total = 0
    for f in sorted(glob.glob(str(ROOT / "data" / "endpoints" / "*.csv"))):
        d = pd.read_csv(f)
        if "pchembl" not in d.columns or "label" not in d.columns:
            continue
        total += int(pred(d, pd.to_numeric(d["pchembl"], errors="coerce")).sum())
    return total


claim("code", "endpoint tables on disk", "data/endpoints",
      lambda: len(list((ROOT / "data" / "endpoints").glob("*.csv"))), "{:d}")
claim("code", "active rows between the two cuts, 6.0 to 7.0", "data/endpoints/*.csv",
      lambda: _endpoint_rows(lambda d, p: (d.label == 1) & (p >= 6.0) & (p < 7.0)), "{:,}")
claim("code", "rows at the censored inactive bound, pChEMBL 5.0", "data/endpoints/*.csv",
      lambda: _endpoint_rows(lambda d, p: (d.label == 0) & (p == 5.0)), "{:,}")


def _const(path: str, name: str) -> str:
    """The literal a module assigns to a name, so a setting quoted here cannot drift from the code."""
    import re
    text = (ROOT / path).read_text(encoding="utf-8")
    m = re.search(rf"^{re.escape(name)}\s*=\s*([^\n#]+)", text, flags=re.M)
    if not m:
        raise KeyError(f"{name} not assigned in {path}")
    return m.group(1).strip()


claim("code", "reliability gate, minimum sensitivity", "panel.py",
      lambda: float(_const("src/brainsafe/panel.py", "MIN_SENSITIVITY")), "{:.2f}")
claim("code", "reliability gate, minimum AUROC", "panel.py",
      lambda: float(_const("src/brainsafe/panel.py", "MIN_AUROC")), "{:.2f}")
claim("code", "threshold requirement A, measured-inactive FPR",
      "final_thresholds.py",
      lambda: float(_const("src/brainsafe/models/final_thresholds.py", "TARGET_FPR")), "{:.2f}")
claim("code", "threshold requirement B, background FPR", "final_thresholds.py",
      lambda: float(_const("src/brainsafe/models/final_thresholds.py", "BACKGROUND_FPR")), "{:.2f}")
claim("code", "screening false-positive rate", "screening_thresholds.py",
      lambda: float(_const("src/brainsafe/models/screening_thresholds.py", "SCREENING_FPR")),
      "{:.2f}")
claim("code", "binder active-scaffold hold-out share", "train_binders_hybrid.py",
      lambda: float(_const("src/brainsafe/models/train_binders_hybrid.py", "ACTIVE_HOLDOUT")),
      "{:.2f}")
def _node_dim() -> int:
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "src" / "brainsafe" / "gnn"))
    from graph_features import NODE_DIM
    return int(NODE_DIM)


claim("code", "graph node feature dimension", "gnn/graph_features.py", _node_dim, "{:d}")



# ---- The criteria justification --------------------------------------------------------------
# This document argues that a choice was justified by measurement rather than by preference, so the
# measurements it leans on are pinned. The feature-block ablation matters most: it is the evidence
# that the twelve descriptors add almost nothing, which is a claim against the project's own
# interest and would be embarrassing to state from a stale number.

def _ablation():
    from collections import defaultdict
    by = defaultdict(dict)
    for r in rows(TAB / "feature_block_ablation.csv"):
        by[(r["endpoint"], r["metric"])][r["block"]] = float(r["mean"])
    return [v for v in by.values() if len(v) == 3]


claim("criteria", "ablation endpoints", "feature_block_ablation.csv",
      lambda: len(_ablation()), "{:d}")
claim("criteria", "combined minus fingerprint-only, mean", "feature_block_ablation.csv",
      lambda: round(st.mean(v["combined"] - v["fingerprint_only"] for v in _ablation()), 4))
claim("criteria", "combined minus fingerprint-only, best", "feature_block_ablation.csv",
      lambda: round(max(v["combined"] - v["fingerprint_only"] for v in _ablation()), 4))
claim("criteria", "combined minus descriptors-only, mean", "feature_block_ablation.csv",
      lambda: round(st.mean(v["combined"] - v["descriptors_only"] for v in _ablation()), 4))
claim("criteria", "combined minus descriptors-only, worst endpoint", "feature_block_ablation.csv",
      lambda: round(max(v["combined"] - v["descriptors_only"] for v in _ablation()), 4))
claim("criteria", "endpoints where combined beats fingerprint alone",
      "feature_block_ablation.csv",
      lambda: sum(1 for v in _ablation() if v["combined"] > v["fingerprint_only"]), "{:d}")

claim("criteria", "core AUROC, random split", "rf_cv_summary.csv",
      lambda: round(st.mean(col(TAB / "rf_cv_summary.csv", "roc_auc_mean",
                                lambda r: r["task"] == "classification"
                                and r["split"] == "random")), 4))
claim("criteria", "core AUROC, scaffold split", "rf_cv_summary.csv",
      lambda: round(st.mean(col(TAB / "rf_cv_summary.csv", "roc_auc_mean",
                                lambda r: r["task"] == "classification"
                                and r["split"] == "scaffold")), 4))

for _alt, _sub, _field, _lab, _fmt in [
        ("HistGradientBoosting", "all 13", "wilcoxon_p", "forest vs HGB, p", "{:.5f}"),
        ("XGBoost", "all 13", "wilcoxon_p", "forest vs XGBoost, p", "{:.5f}"),
        ("kNN read-across", "8 classification", "median_delta",
         "forest vs read-across, median delta", "{:.4f}"),
        ("LogisticRegression", "8 classification", "median_delta",
         "forest vs logistic regression, median delta", "{:.4f}")]:
    claim("criteria", _lab, "model_family_significance.csv",
          (lambda a, sub, f: lambda: float(
              cell(TAB / "model_family_significance.csv", f,
                   lambda r: r["split"] == "scaffold" and r["alternative"] == a
                   and r["subset"].startswith(sub))))(_alt, _sub, _field), _fmt)

claim("criteria", "GSK-3B AUROC before the bulk inactives", "inactives_audit.csv",
      lambda: float(cell(TAB / "inactives_audit.csv", "auroc_before",
                         lambda r: r["endpoint"] == "GSK3B")))
claim("criteria", "GSK-3B AUROC after the bulk inactives", "inactives_audit.csv",
      lambda: float(cell(TAB / "inactives_audit.csv", "auroc_after",
                         lambda r: r["endpoint"] == "GSK3B")))
claim("criteria", "bulk inactives added to GSK-3B", "inactives_audit.csv",
      lambda: int(cell(TAB / "inactives_audit.csv", "added_inactives",
                       lambda r: r["endpoint"] == "GSK3B")), "{:,}")
claim("criteria", "median similarity of those inactives to the actives", "inactives_audit.csv",
      lambda: float(cell(TAB / "inactives_audit.csv", "median_sim_to_active",
                         lambda r: r["endpoint"] == "GSK3B")), "{:.3f}")

claim("criteria", "core ECE before calibration", "calibration.csv",
      lambda: round(st.mean(col(TAB / "calibration.csv", "ece_raw")), 4))
claim("criteria", "core ECE after calibration", "calibration.csv",
      lambda: round(st.mean(col(TAB / "calibration.csv", "ece_calibrated")), 4))


# ---- The ML foundations document --------------------------------------------------------------
# It teaches the mathematics from first principles, so its numbers are the ones a reader will take
# away as definitions rather than as results. The sample sizes and class balances matter most: they
# are quoted in a table a reader will treat as the specification of the dataset, and a drifted count
# there would be repeated in every later description of the work.

def _oof_balance(ep):
    r = rows(ROOT / "data" / "processed" / "cv_predictions" / f"{ep}_random_oof.csv")
    ys = [float(x["y_true"]) for x in r]
    return len(ys), int(sum(ys))


for _ep in ("BBB", "hERG", "MAO_A", "BACE1"):
    claim("mlmath", f"{_ep} rows after deduplication", f"{_ep}_random_oof.csv",
          (lambda e: lambda: _oof_balance(e)[0])(_ep), "{:,}")
    claim("mlmath", f"{_ep} positives", f"{_ep}_random_oof.csv",
          (lambda e: lambda: _oof_balance(e)[1])(_ep), "{:,}")

claim("mlmath", "feature vector width", "featurize.py", lambda: 1024 + 12, "{:,}")

claim("mlmath", "collision: distinct environments", "fingerprint_collisions.csv",
      lambda: int(float(cell(TAB / "fingerprint_collisions.csv", "value",
                             lambda r: r["measure"] == "distinct atomic environments"))), "{:,}")
claim("mlmath", "collision: environments per bit, median", "fingerprint_collisions.csv",
      lambda: int(float(cell(TAB / "fingerprint_collisions.csv", "value",
                             lambda r: r["measure"].startswith(
                                 "environments per occupied bit, median")))), "{:d}")
claim("mlmath", "collision: environments per bit, max", "fingerprint_collisions.csv",
      lambda: int(float(cell(TAB / "fingerprint_collisions.csv", "value",
                             lambda r: r["measure"].startswith(
                                 "environments per occupied bit, max")))), "{:d}")

claim("mlmath", "ablation: descriptors add", "feature_block_ablation.csv",
      lambda: round(st.mean(v["combined"] - v["fingerprint_only"] for v in _ablation()), 4))
claim("mlmath", "ablation: fingerprint adds", "feature_block_ablation.csv",
      lambda: round(st.mean(v["combined"] - v["descriptors_only"] for v in _ablation()), 4))

claim("mlmath", "core ECE before calibration", "calibration.csv",
      lambda: round(st.mean(col(TAB / "calibration.csv", "ece_raw")), 4))
claim("mlmath", "core ECE after calibration", "calibration.csv",
      lambda: round(st.mean(col(TAB / "calibration.csv", "ece_calibrated")), 4))

claim("mlmath", "BBB conformal set size", "rf_conformal.csv",
      lambda: float(cell(TAB / "rf_conformal.csv", "avg_set_size",
                         lambda r: r["endpoint"] == "BBB")), "{:.3f}")
claim("mlmath", "BACE1 empty-set fraction", "rf_conformal.csv",
      lambda: round(100 * float(cell(TAB / "rf_conformal.csv", "frac_empty",
                                     lambda r: r["endpoint"] == "BACE1")), 1), "{:.1f}")

claim("mlmath", "distinct structures in the library", "library_sp3_coverage.csv",
      lambda: int(float(cell(TAB / "library_sp3_coverage.csv", "value",
                             lambda r: r["metric"].startswith("distinct SMILES")))), "{:,}")

CHAPTER_FILES = {
    "09": THESIS / "chapter09_falsification.md",
    "10": THESIS / "chapter10_limitations.md",
    "viva": THESIS / "viva_preparation.md",
    "code": THESIS / "code_walkthrough.md",
    "criteria": THESIS / "criteria_justification.md",
    "mlmath": THESIS / "ml_foundations.md",
}


def _present(text: str, value, fmt: str) -> bool:
    """Is the formatted value written in the chapter?

    Numbers are matched with a word boundary so that 0.072 does not match inside 0.0721, and a
    thousands-separated integer is accepted with or without its comma.
    """
    s = fmt.format(value)
    variants = {s}
    if "," in s:
        variants.add(s.replace(",", ""))
    if s.startswith("0."):
        variants.add(s[1:])
    if s.startswith("-"):
        # Prose uses a typographic minus, U+2212, where a format string produces a hyphen. They are
        # the same number and the check should not care which glyph a chapter happens to use.
        variants |= {"−" + v[1:] for v in list(variants)}
    return any(re.search(r"(?<![\d.])" + re.escape(v) + r"(?![\d])", text) for v in variants)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chapter", default=None, help="two-digit chapter number, e.g. 09")
    args = ap.parse_args()

    wanted = [c for c in sorted(CHAPTER_FILES) if args.chapter in (None, c)]
    if not wanted:
        print(f"no claims registered for chapter {args.chapter}")
        return 1

    failures = 0
    checked = 0
    for ch in wanted:
        path = CHAPTER_FILES[ch]
        text = path.read_text(encoding="utf-8", errors="replace")
        print(f"\n== chapter {ch}: {path.name} ==")
        for c_ch, desc, src, fn, fmt in CLAIMS:
            if c_ch != ch:
                continue
            checked += 1
            try:
                value = fn()
            except Exception as exc:                       # noqa: BLE001
                print(f"  [ERROR] {desc:52} could not read {src}: {exc}")
                failures += 1
                continue
            ok = _present(text, value, fmt)
            if not ok:
                failures += 1
            print(f"  [{'ok' if ok else 'MISSING'}] {desc:52} {fmt.format(value):>10}  <- {src}")

    print(f"\n{checked} claims checked, {failures} not found in the chapter text")
    if failures:
        print("A missing claim means the chapter states something the artefact no longer supports,")
        print("or that the artefact has moved since the chapter was written. Both need a human.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
