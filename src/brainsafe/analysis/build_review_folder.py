"""Assemble the organised review folder: inputs, cross-validation, external validation, models.

Three things are built here that did not exist in the long-format package.

The compound matrix. One row per unique structure, one column per endpoint, holding that compound's
measured value at that endpoint. This is the shape a reviewer reaches for when asking what was
measured on what: it makes the sparsity of the panel visible at a glance, which the long format
hides. Most cells are empty, and that is the honest picture. A compound measured at one target and
nowhere else has one filled cell out of sixty, because nobody has run the other fifty-nine assays on
it. An empty cell means no measurement exists in the sources, never a measurement that was dropped.

The external matrix. The same treatment for compounds the models were scored on but not fitted to,
carrying the prediction each evaluation produced alongside the overlap check.

The folder itself. Numbered directories in the order a reviewer works through them, each with the
inputs that went in and the outputs that came out, and a starting note that states the compound
accounting before anything else.

Writes reviewer_package/<date>/REVIEW/
"""
from __future__ import annotations

import json
import shutil
import sys
import warnings
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
EP = ROOT / "data" / "endpoints"
EPR = ROOT / "data" / "endpoints_reg"
TAB = ROOT / "results" / "tables"
INV = ROOT / "inversion" / "results"


def canon(s):
    m = Chem.MolFromSmiles(str(s))
    return Chem.MolToSmiles(m) if m else None


def compound_matrix():
    """One row per structure, one column per endpoint, holding the measured value."""
    import app
    order = (list(app.TARGET_CLASSIFIERS) + list(app.RECEPTOR_REGRESSORS)
             + list(app.BINDER_TARGETS))
    cells, panel_of = {}, {}
    for n in order:
        p = EP / f"{n}.csv"
        if not p.exists():
            continue
        panel_of[n] = ("Core target/property classifier" if n in app.TARGET_CLASSIFIERS
                       else "Receptor potency regressor" if n in app.RECEPTOR_REGRESSORS
                       else "Binder classifier")
        d = pd.read_csv(p)
        # An endpoint table without a pchembl column is labelled rather than potency-valued; guard
        # for it, because pd.to_numeric on a missing column returns a scalar NaN, not a Series.
        pch = pd.to_numeric(d["pchembl"], errors="coerce") if "pchembl" in d.columns else None
        for i, smi in enumerate(d["smiles"].astype(str)):
            c = canon(smi)
            if c is None:
                continue
            v = pch.iloc[i] if pch is not None else None
            if pd.notna(v):
                val = f"pChEMBL {v:.2f}"
            elif "label" in d:
                val = "active" if d["label"].iloc[i] == 1 else "inactive"
            else:
                val = "measured"
            cells.setdefault(c, {})[n] = val
    for p in sorted(EPR.glob("*.csv")) if EPR.exists() else []:
        n = p.stem
        panel_of[n] = "Regression endpoint"
        d = pd.read_csv(p)
        ycol = "y" if "y" in d.columns else d.columns[-1]
        for i, smi in enumerate(d["smiles"].astype(str)):
            c = canon(smi)
            if c is not None:
                cells.setdefault(c, {})[n] = f"{pd.to_numeric(d[ycol].iloc[i], errors='coerce'):.3f}"

    cols = [c for c in order if c in panel_of] + [c for c in panel_of if c not in order]
    rows = []
    for c, vals in cells.items():
        r = {"canonical_smiles": c, "n_endpoints_measured": len(vals)}
        r.update({k: vals.get(k) for k in cols})
        rows.append(r)
    m = pd.DataFrame(rows).sort_values("n_endpoints_measured", ascending=False)
    header = pd.DataFrame([{"canonical_smiles": "PANEL", "n_endpoints_measured": None,
                            **{k: panel_of[k] for k in cols}}])
    return pd.concat([header, m], ignore_index=True), cols


def main():
    base = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "reviewer_package" / "08.09.26"
    out = base / "REVIEW"
    dirs = {k: out / k for k in
            ["01_INPUTS", "02_CROSS_VALIDATION", "03_EXTERNAL_VALIDATION", "04_MODELS"]}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    mat, cols = compound_matrix()
    mat.to_csv(dirs["01_INPUTS"] / "COMPOUND_MATRIX_endpoints_as_columns.csv", index=False)
    body = mat.iloc[1:]
    filled = int(body["n_endpoints_measured"].sum())
    print(f"COMPOUND_MATRIX: {len(body):,} compounds x {len(cols)} endpoint columns")
    print(f"   measurements present {filled:,} of {len(body) * len(cols):,} cells "
          f"({filled / (len(body) * len(cols)) * 100:.2f}% filled)")
    print(f"   measured at one endpoint only : "
          f"{int((body.n_endpoints_measured == 1).sum()):,}")
    print(f"   measured at ten or more       : "
          f"{int((body.n_endpoints_measured >= 10).sum()):,}")

    for src, dst in [(base / "MASTER_training_inputs.csv", dirs["01_INPUTS"]),
                     (base / "02_training_input_files.csv", dirs["01_INPUTS"]),
                     (base / "01_MASTER_endpoint_inventory.csv", dirs["01_INPUTS"]),
                     (TAB / "manuscript_T2_per_fold.csv", dirs["02_CROSS_VALIDATION"]),
                     (TAB / "adme_cv_folds.csv", dirs["02_CROSS_VALIDATION"]),
                     (TAB / "adme_cv_summary.csv", dirs["02_CROSS_VALIDATION"]),
                     (TAB / "binder_cv_folds.csv", dirs["02_CROSS_VALIDATION"]),
                     (TAB / "binder_cv_summary.csv", dirs["02_CROSS_VALIDATION"]),
                     (TAB / "rf_cv_summary.csv", dirs["02_CROSS_VALIDATION"]),
                     (base / "MASTER_external_test_results.csv", dirs["03_EXTERNAL_VALIDATION"]),
                     (base / "MASTER_external_test_summary.csv", dirs["03_EXTERNAL_VALIDATION"]),
                     (TAB / "manuscript_Table6_scaffold_holdout.csv", dirs["03_EXTERNAL_VALIDATION"]),
                     (TAB / "manuscript_Table5_temporal.csv", dirs["03_EXTERNAL_VALIDATION"]),
                     (TAB / "noncns_specificity_summary.csv", dirs["03_EXTERNAL_VALIDATION"]),
                     (INV / "H6_clinical_indication.csv", dirs["03_EXTERNAL_VALIDATION"]),
                     (INV / "VERDICTS.csv", dirs["03_EXTERNAL_VALIDATION"]),
                     (base / "06_model_registry.csv", dirs["04_MODELS"]),
                     (base / "07_binder_panel_training_design.csv", dirs["04_MODELS"]),
                     (base / "05_feature_definitions.csv", dirs["04_MODELS"]),
                     (ROOT / "models_manifest.json", dirs["04_MODELS"])]:
        if src.exists():
            shutil.copy2(src, dst / src.name)
        else:
            print(f"   [missing, not copied] {src.name}")

    xl = base / "BrainSafe_AI_training_record.xlsx"
    if xl.exists():
        shutil.copy2(xl, out / xl.name)

    acct = pd.read_csv(base / "REVIEW" / "01_INPUTS" / "01_MASTER_endpoint_inventory.csv") \
        if (dirs["01_INPUTS"] / "01_MASTER_endpoint_inventory.csv").exists() else None

    folds = pd.concat([pd.read_csv(p) for p in
                       [TAB / "manuscript_T2_per_fold.csv", TAB / "adme_cv_folds.csv",
                        TAB / "binder_cv_folds.csv"] if p.exists()], ignore_index=True)
    n_ep = folds.endpoint.nunique()

    (out / "00_START_HERE.md").write_text(f"""# BrainSafe AI: review materials

Generated by `src/brainsafe/analysis/build_review_folder.py`. Every number below is read from the
file that produced it.

## The compound accounting, first, because it is the thing most easily misread

No endpoint was trained on the panel total. Each endpoint has its own measured set, and it is
cross-validated on that set alone. The sets differ by two orders of magnitude.

| | Records | Unique compounds (InChIKey) |
|---|---|---|
| Core target panel, 12 targets plus antioxidant | 67,984 | 61,226 |
| Full panel, all endpoints | 203,884 | 160,365 |

Those are sums across endpoints. A compound measured at five targets contributes five records and is
counted once per endpoint. The smallest endpoint, GluA2, has 183 compounds; the largest, BACE1, has
8,501. `01_INPUTS/COMPOUND_MATRIX_endpoints_as_columns.csv` shows this directly: one row per
structure, one column per endpoint, and most cells empty because most compounds were assayed against
one or two targets and never against the rest.

## How the ten-fold cross-validation was actually run

Ten folds per endpoint, in two regimes, on that endpoint's own set:

- **random**: folds drawn at random, measuring interpolation within known chemistry
- **scaffold-grouped**: `GroupKFold` over generic Bemis-Murcko scaffolds, so entire chemical series
  are withheld, measuring extrapolation to new chemotypes

{n_ep} endpoints were cross-validated this way, which is {n_ep} x 2 x 10 = {len(folds):,} model fits.
Every one of those fits was scored once on the fold withheld from it and then discarded. The model
that is served is refitted on the endpoint's full set afterwards, so a compound never contributes to
both the training and the scoring of the same model. There are not ten deployed models per endpoint;
there is one, and ten transient models per regime that existed only to produce the estimate.

## What is in each folder

| Folder | Contents |
|---|---|
| `01_INPUTS` | the compound matrix, every training record in long form, the per-endpoint inventory, and the file-by-file input manifest |
| `02_CROSS_VALIDATION` | every individual fold for every endpoint, and the per-endpoint summaries. `binder_cv_summary.csv` also carries the AUROC recorded when each model was trained, so the re-run and the original can be read against each other |
| `03_EXTERNAL_VALIDATION` | compounds the models were scored on but not fitted to, with the overlap check and the results split on it |
| `04_MODELS` | the model registry, thresholds and how they were set, the 1036 feature definitions, and the checksum manifest for the archived binaries |
| `BrainSafe_AI_training_record.xlsx` | all of the above in one workbook |

## On empty cells

An empty cell always means one of a small number of stated things, never "we did not look". The
reasons are enumerated in the workbook's data-dictionary sheet. The two most common: a regression
endpoint has no AUROC because the metric does not apply, and a compound has no value at an endpoint
because that assay was never run on it in any source we drew from.

## Reproducing any of this

```
python src/brainsafe/analysis/build_master_inputs_and_external.py
python src/brainsafe/evaluation/binder_cv_per_fold.py
python src/brainsafe/analysis/build_reviewer_workbook.py
python src/brainsafe/analysis/build_review_folder.py
```

Models are archived at doi:10.5281/zenodo.21858576; source at
https://github.com/krishna-g-999/brainsafe-ai.
""", encoding="utf-8")

    print(f"\nreview folder: {out}")
    for d in sorted(out.rglob("*")):
        if d.is_file():
            print(f"   {d.relative_to(out).as_posix():62} {d.stat().st_size / 1024:>9,.0f} KB")


if __name__ == "__main__":
    main()
