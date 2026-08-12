"""One workbook holding every input and every per-fold result, for a reader who wants the mechanics.

Written after a reviewer objected, fairly, that a compound-by-endpoint matrix is mostly empty. It is
mostly empty because the evidence is: no compound has been measured against most targets. A sparse
matrix is the wrong shape for that, so the measurements are given here in long form, one row per
measurement, with no empty cells at all.

The other objection was that the summary did not show the mechanics. This workbook answers, for every
model: what the feature vector is, column by column; how the model was fitted; and what every
individual cross-validation fold returned, rather than only the mean.

Sheets:
  01_READ_ME                 what each sheet is, and the counts that matter
  02_FEATURE_VECTOR          all 1,036 input columns defined, one row each
  03_MODEL_INDEX             every model: algorithm, hyper-parameters, training set, CV scheme
  04_FOLDS_core              per-fold results, 13 core endpoints x 2 splits x 10 folds
  05_FOLDS_binders           per-fold results, binder panel
  06_FOLDS_adme              per-fold results, 9 ADME endpoints x 2 splits x 10 folds
  07_TRAINING_COMPOSITION    per model: actives, inactives, decoys, what was withheld
  08_MEASUREMENTS_LONG       one row per measurement: compound, endpoint, label, value. No blanks.
  09_DATA_SOURCES            per endpoint: file, rows, label rule, class balance
  10_PROVENANCE              when each source table was written, so staleness is visible

Run:  python src/brainsafe/analysis/build_reviewer_workbook_full.py
"""
from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import MORGAN_BITS, MORGAN_RADIUS, _DESCRIPTORS, feature_names  # noqa: E402

M = ROOT / "models_rf"
TAB = ROOT / "results" / "tables"
OUT = ROOT / "reviewer_package" / "model_outputs"
XLSX = OUT / "BrainSafe_models_inputs_and_folds.xlsx"

DESCRIPTOR_NOTES = {
    "mw": "molecular weight, Descriptors.MolWt",
    "clogp": "Crippen logP, Crippen.MolLogP",
    "tpsa": "topological polar surface area, rdMolDescriptors.CalcTPSA",
    "hbd": "hydrogen-bond donors",
    "hba": "hydrogen-bond acceptors",
    "rotatable_bonds": "rotatable bond count",
    "aromatic_rings": "aromatic ring count",
    "fraction_csp3": "fraction of sp3 carbons",
    "ring_count": "ring count",
    "heavy_atoms": "heavy-atom count",
    "formal_charge": "net formal charge",
    "qed": "quantitative estimate of drug-likeness",
}


def _json(p: Path):
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _mtime(p: Path) -> str:
    return _dt.datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M") if p.exists() else ""


def sheet_feature_vector() -> pd.DataFrame:
    """Every input column the estimator receives, defined individually."""
    names = feature_names()
    rows = []
    for i, n in enumerate(names):
        if n.startswith("ecfp4_"):
            rows.append({
                "index": i, "column": n, "block": "fingerprint", "dtype": "float32 (0.0 or 1.0)",
                "definition": f"bit {n.split('_')[1]} of a folded Morgan/ECFP-4 fingerprint, "
                              f"radius {MORGAN_RADIUS}, {MORGAN_BITS} bits, chirality NOT included",
                "computed_on": "largest organic fragment after salt stripping, sanitised",
                "note": "folded, so several substructure environments share a bit; bit identity is "
                        "not a unique substructure",
            })
        else:
            rows.append({
                "index": i, "column": n, "block": "descriptor", "dtype": "float32",
                "definition": DESCRIPTOR_NOTES.get(n, n),
                "computed_on": "largest organic fragment after salt stripping, sanitised",
                "note": "unscaled; random forests are invariant to monotone rescaling, so no "
                        "scaler is fitted and none can leak across a split",
            })
    return pd.DataFrame(rows)


def sheet_model_index() -> pd.DataFrame:
    cat = OUT / "OUTPUT_CATALOGUE.csv"
    if cat.exists():
        return pd.read_csv(cat)
    raise SystemExit("run build_reviewer_matrix.py first; OUTPUT_CATALOGUE.csv is missing")


def sheet_training_composition() -> pd.DataFrame:
    rows = []
    for f in sorted(M.glob("*_meta.json")):
        meta = _json(f)
        if "task" not in meta:
            continue
        ded = meta.get("deduplication") or {}
        rows.append({
            "model": meta["endpoint"], "family": "core (train_rf.py)",
            "positives": meta.get("positives"), "total_after_dedup": meta.get("n_compounds"),
            "rows_before_dedup": ded.get("rows_in"),
            "duplicate_rows_removed": ded.get("duplicate_rows_removed"),
            "contradictory_groups_dropped": ded.get("conflicting_groups_dropped"),
            "n_scaffolds": meta.get("n_scaffolds"),
            "decoys": None, "actives_withheld": None, "inactives_withheld": None,
            "deployed_fit_on": "all remaining rows after deduplication",
        })
    for ep, v in sorted(_json(M / "binder_modes.json").items()):
        rows.append({
            "model": f"{ep}_binder", "family": f"binder ({v.get('mode')})",
            "positives": v.get("n_positive"), "total_after_dedup": None,
            "rows_before_dedup": None, "duplicate_rows_removed": None,
            "contradictory_groups_dropped": None, "n_scaffolds": None,
            "decoys": v.get("n_decoy"),
            "actives_withheld": v.get("n_active_holdout"),
            "inactives_withheld": v.get("n_measured_inactive_holdout"),
            "deployed_fit_on": "80% train split, sigmoid calibration on the remaining 20%",
        })
    for f in sorted((M / "adme").glob("*_meta.json")):
        meta = _json(f)
        rows.append({
            "model": f"adme_{meta.get('endpoint', f.stem)}", "family": "ADME (train_adme.py)",
            "positives": None, "total_after_dedup": meta.get("n_compounds"),
            "rows_before_dedup": None, "duplicate_rows_removed": None,
            "contradictory_groups_dropped": None, "n_scaffolds": None,
            "decoys": None, "actives_withheld": None, "inactives_withheld": None,
            "deployed_fit_on": "all rows",
        })
    return pd.DataFrame(rows)


def sheet_measurements_long() -> pd.DataFrame:
    """One row per measurement. No empty cells, because absent measurements are absent rows."""
    src = OUT / "TRAINING_MATRIX.csv"
    if not src.exists():
        raise SystemExit("run build_reviewer_matrix.py first; TRAINING_MATRIX.csv is missing")
    wide = pd.read_csv(src)
    id_cols = ["inchikey", "parent_smiles"]
    value_cols = [c for c in wide.columns if "__" in c]
    long = wide.melt(id_vars=id_cols, value_vars=value_cols,
                     var_name="_col", value_name="value").dropna(subset=["value"])
    long["endpoint"] = long["_col"].str.split("__").str[0]
    long["measurement"] = long["_col"].str.split("__").str[1]
    return long[["inchikey", "parent_smiles", "endpoint", "measurement", "value"]]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    folds = {
        "04_FOLDS_core": TAB / "rf_cv_folds.csv",
        "05_FOLDS_binders": TAB / "binder_cv_folds.csv",
        "06_FOLDS_adme": TAB / "adme_cv_folds.csv",
    }
    print("assembling sheets ...", flush=True)
    feat = sheet_feature_vector()
    idx = sheet_model_index()
    comp = sheet_training_composition()
    long = sheet_measurements_long()
    srcs = pd.read_csv(OUT / "TRAINING_SOURCES.csv") if (OUT / "TRAINING_SOURCES.csv").exists() \
        else pd.DataFrame()

    fold_frames, prov = {}, []
    for sheet, path in folds.items():
        fold_frames[sheet] = pd.read_csv(path) if path.exists() else pd.DataFrame()
        prov.append({"sheet": sheet, "source": path.relative_to(ROOT).as_posix(),
                     "written": _mtime(path), "rows": len(fold_frames[sheet])})
    for extra in ("rf_cv_summary.csv", "binder_cv_summary.csv", "adme_cv_summary.csv"):
        p = TAB / extra
        prov.append({"sheet": "(summary, not included)", "source": f"results/tables/{extra}",
                     "written": _mtime(p), "rows": len(pd.read_csv(p)) if p.exists() else 0})
    prov.append({"sheet": "03_MODEL_INDEX", "source": "models_rf/*_meta.json, binder_modes.json",
                 "written": _mtime(M / "binder_modes.json"), "rows": len(idx)})
    prov = pd.DataFrame(prov)

    total_folds = sum(len(v) for v in fold_frames.values())
    readme = pd.DataFrame([
        ("Purpose", "Every model input and every per-fold result, in one workbook."),
        ("Generated", _dt.datetime.now().strftime("%Y-%m-%d %H:%M")),
        ("Generated by", "src/brainsafe/analysis/build_reviewer_workbook_full.py"),
        ("", ""),
        ("Feature vector", f"{MORGAN_BITS} ECFP-4 bits (Morgan, radius {MORGAN_RADIUS}, chirality "
                          f"OFF) + {len(_DESCRIPTORS)} descriptors = {MORGAN_BITS + len(_DESCRIPTORS)} "
                          f"columns. Defined column by column in 02_FEATURE_VECTOR."),
        ("Computed on", "the largest organic fragment after salt stripping, sanitised."),
        ("Scaling", "none. No scaler is fitted anywhere, so none can leak across a split."),
        ("", ""),
        ("Models on disk", len(idx)),
        ("Models scored per compound", int((idx.get("deployed") != False).sum())
         if "deployed" in idx.columns else len(idx)),
        ("Withdrawn, never scored", "Nav1_1, Cav3_2"),
        ("Dual-modelled endpoints", "D2, A2A, HT2A, SERT each carry a potency regression AND a "
                                    "binder classifier, so they appear twice."),
        ("", ""),
        ("Cross-validation", "10-fold, run twice per endpoint: StratifiedKFold/KFold with "
                             "shuffle and seed 42 (random), and GroupKFold on the Bemis-Murcko "
                             "scaffold (scaffold). Seed 42 throughout."),
        ("Per-fold rows in this workbook", total_folds),
        ("What a fold row is", "one fitted model. n_train and n_test are its split sizes; the "
                               "metric columns are that model scored on its own held-out fold."),
        ("Deployed model", "a separate refit on all data after CV; it is not any single fold."),
        ("", ""),
        ("08_MEASUREMENTS_LONG", f"{len(long):,} rows, one per measurement. There are no empty "
                                 f"cells: a compound not measured against an endpoint has no row, "
                                 f"rather than a blank."),
        ("Why the earlier matrix was sparse", "it was compound x endpoint, and no compound is "
                                              "measured against most endpoints. Long form is the "
                                              "correct shape for that and removes every blank."),
    ], columns=["item", "value"])

    # xlsxwriter, not openpyxl: pandas 3.0.3 with openpyxl 3.1.5 raises "At least one sheet must be
    # visible" on any write, including a two-cell one, so the engine is pinned rather than left to
    # the default and left to fail on someone else's machine.
    with pd.ExcelWriter(XLSX, engine="xlsxwriter") as xl:
        readme.to_excel(xl, sheet_name="01_READ_ME", index=False)
        feat.to_excel(xl, sheet_name="02_FEATURE_VECTOR", index=False)
        idx.to_excel(xl, sheet_name="03_MODEL_INDEX", index=False)
        for sheet, frame in fold_frames.items():
            frame.to_excel(xl, sheet_name=sheet, index=False)
        comp.to_excel(xl, sheet_name="07_TRAINING_COMPOSITION", index=False)
        long.to_excel(xl, sheet_name="08_MEASUREMENTS_LONG", index=False)
        if len(srcs):
            srcs.to_excel(xl, sheet_name="09_DATA_SOURCES", index=False)
        prov.to_excel(xl, sheet_name="10_PROVENANCE", index=False)

    # CSV twins, for anyone who would rather not open a workbook
    feat.to_csv(OUT / "FEATURE_VECTOR.csv", index=False)
    comp.to_csv(OUT / "TRAINING_COMPOSITION.csv", index=False)
    long.to_csv(OUT / "MEASUREMENTS_LONG.csv", index=False)
    prov.to_csv(OUT / "PROVENANCE.csv", index=False)

    print(f"wrote {XLSX.relative_to(ROOT).as_posix()}")
    print(f"  feature columns {len(feat)}, models {len(idx)}, fold rows {total_folds:,}, "
          f"measurements {len(long):,}")


if __name__ == "__main__":
    main()
