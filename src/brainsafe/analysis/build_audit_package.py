"""The audit package: what the model is, what it was fitted to, and what runs a prediction.

Written for a reader who has to satisfy themselves that the tool does what the manuscript says,
without reading the source. It answers, in order, the questions such a reader actually asks:

  how does the model work            one estimator family, one fixed 1036-dimensional input, and a
                                     count of exactly how many model files exist per endpoint and
                                     what each is for
  what data was it fitted to         every endpoint's measured set, its source, the rule that turned
                                     a measurement into a label, and the threshold used
  what does it receive               the feature matrix itself, structures as rows and the 1036
                                     model inputs as columns, with the values
  how was it tested                  every fold of every endpoint, both split regimes
  was anyone trying to break it      the inversion analysis, kept separate because a falsification
                                     record is not a results table
  what runs at prediction time       the files the served application opens for a single query,
                                     which are not the same files that trained it

Every CSV opens with three header lines stating what it holds, what it is for, and what wrote it,
so a file read in isolation still explains itself. Programmatic readers should pass skiprows=3.

Not comment='#'. SMILES use '#' for a triple bond, so treating it as a comment character truncates
every nitrile and alkyne in the file: C#N silently becomes C. Skipping three lines by position is
the only safe way to read these.

Writes AUDIT_PACKAGE/
"""
from __future__ import annotations

import json
import shutil
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize, feature_names  # noqa: E402

M = ROOT / "models_rf"
EP = ROOT / "data" / "endpoints"
EPR = ROOT / "data" / "endpoints_reg"
TAB = ROOT / "results" / "tables"
INV = ROOT / "inversion" / "results"
OUT = ROOT / "AUDIT_PACKAGE"

INDEX: list[dict] = []


def write_csv(df: pd.DataFrame, path: Path, contains: str, used_for: str, source: str):
    """Write a CSV that explains itself in its first three lines."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# WHAT THIS FILE CONTAINS: {contains}\n")
        fh.write(f"# WHAT IT IS USED FOR: {used_for}\n")
        fh.write(f"# GENERATED FROM: {source}\n")
        df.to_csv(fh, index=False, lineterminator="\n")
    INDEX.append({
        "file": path.relative_to(OUT).as_posix(),
        "rows": len(df), "columns": len(df.columns),
        "size_kb": round(path.stat().st_size / 1024, 1),
        "what_it_contains": contains, "what_it_is_used_for": used_for,
        "generated_from": source,
    })
    print(f"  {path.relative_to(OUT).as_posix():62} {len(df):8,} rows")


def model_inventory():
    """How many model files exist per endpoint, and what each one is for."""
    import app
    modes = json.loads((M / "binder_modes.json").read_text())
    rows = []
    role = {
        "": "base estimator, fitted on the endpoint's full set after cross-validation fixed the "
            "estimate",
        "_calibrated": "probability calibration wrapped around the frozen base estimator, so the "
                       "reported number is a probability and not a forest vote share",
        "_binder": "binder classifier trained against property-matched decoys plus measured "
                   "inactives, used where the measured set is almost entirely active",
    }
    for ep in sorted(set(app.TARGET_CLASSIFIERS) | set(app.RECEPTOR_REGRESSORS)
                     | set(app.BINDER_TARGETS) | set(app.ADME) | {"antioxidant_DPPH", "pka_basic"}):
        files = []
        for suffix, what in role.items():
            for cand in (M / f"{ep}{suffix}.joblib", M / "adme" / f"{ep}{suffix}.joblib"):
                if cand.exists():
                    files.append((cand.relative_to(ROOT).as_posix(), what,
                                  round(cand.stat().st_size / 2**20, 2)))
                    break
        if not files:
            continue
        rec = modes.get(ep, {})
        rows.append({
            "endpoint": ep,
            "deployed_model_files": len(files),
            "model_files": " | ".join(f[0] for f in files),
            "what_each_file_is": " | ".join(f[1] for f in files),
            "total_size_mb": round(sum(f[2] for f in files), 2),
            "evaluation_fits_during_10_fold": 20,
            "evaluation_fits_note": "ten folds in each of two split regimes; every one was scored "
                                    "on the fold withheld from it and then discarded",
            "estimator_family": ("RandomForestRegressor" if ep in app.RECEPTOR_REGRESSORS
                                 or ep in ("antioxidant_DPPH", "pka_basic")
                                 else "RandomForestClassifier"),
            "decision_threshold": rec.get("threshold"),
            "threshold_basis": rec.get("threshold_basis"),
        })
    return pd.DataFrame(rows)


def training_and_labelling():
    """Per endpoint: the measured set, where it came from, and the rule that made a label."""
    import app
    modes = json.loads((M / "binder_modes.json").read_text())
    ctx = json.loads((M / "endpoint_context.json").read_text()) \
        if (M / "endpoint_context.json").exists() else {"classifiers": {}}
    rows = []
    for ep in sorted(set(app.TARGET_CLASSIFIERS) | set(app.RECEPTOR_REGRESSORS)
                     | set(app.BINDER_TARGETS)):
        p = EP / f"{ep}.csv"
        if not p.exists():
            continue
        d = pd.read_csv(p)
        rec = modes.get(ep, {})
        is_binder = ep in app.BINDER_TARGETS
        if is_binder and rec.get("mode") == "hybrid_decoys_plus_measured_inactives":
            rule = ("positive if pChEMBL >= 7.0 (100 nM or better). Negative class is two things "
                    "kept distinct: property-matched decoys drawn from background chemistry at 3:1 "
                    "to actives with Tanimoto < 0.35 to every active, and half the measured "
                    "inactives. The other half of the measured inactives never enters training and "
                    "is used to set the decision threshold.")
            n_train = ((rec.get("n_positive") or 0) + (rec.get("n_decoy") or 0)
                       + (rec.get("n_measured_inactive_train") or 0))
        elif is_binder:
            rule = ("measured labels only; the source reports both actives and inactives, so no "
                    "decoys were generated")
            n_train = (rec.get("n_positive") or 0) + (rec.get("n_measured_inactive_train") or 0)
        elif ep in app.RECEPTOR_REGRESSORS:
            rule = ("no label. The measured pChEMBL value is regressed directly, because these sets "
                    "are almost entirely active and a binary split would be ill-posed.")
            n_train = len(d)
        else:
            rule = ("measured dataset label as supplied by the source; see the data source column "
                    "for which curation it follows")
            n_train = len(d)
        src = (d["source"].value_counts().to_dict() if "source" in d.columns else {})
        rows.append({
            "endpoint": ep,
            "panel": ("Core target/property classifier" if ep in app.TARGET_CLASSIFIERS
                      else "Receptor potency regressor" if ep in app.RECEPTOR_REGRESSORS
                      else "Binder classifier"),
            "task": "regression" if ep in app.RECEPTOR_REGRESSORS else "classification",
            "measured_records": len(d),
            "measured_actives": int((d["label"] == 1).sum()) if "label" in d else None,
            "measured_inactives": int((d["label"] == 0).sum()) if "label" in d else None,
            "decoys_generated": rec.get("n_decoy", 0) if is_binder else 0,
            "compounds_the_model_was_fitted_to": n_train,
            "how_a_label_was_assigned": rule,
            "measurement_sources": "; ".join(f"{k} {v:,}" for k, v in src.items()) or "see "
                                                                                      "data/raw/measured_endpoints_SOURCE.md",
            "base_rate_of_activity": ctx.get("classifiers", {}).get(ep, {}).get("base_rate"),
            "cross_validation": "10-fold, in two regimes: random, and GroupKFold over generic "
                                "Bemis-Murcko scaffolds",
            "training_table": f"data/endpoints/{ep}.csv",
        })
    return pd.DataFrame(rows).sort_values("compounds_the_model_was_fitted_to", ascending=False)


def feature_matrices():
    """The model input itself: structures as rows, the 1036 inputs as columns."""
    names = feature_names()
    desc = names[-12:]

    # Every compound, the twelve interpretable descriptors. The 1024 fingerprint bits are omitted
    # here only because 160,000 x 1036 is not a file anyone opens; the full-width matrix below
    # carries them for the tested set.
    smis = set()
    for p in list(EP.glob("*.csv")) + (list(EPR.glob("*.csv")) if EPR.exists() else []):
        smis |= set(pd.read_csv(p, usecols=lambda c: c == "smiles")["smiles"].astype(str))
    smis = sorted(smis)
    X, mask = featurize(smis)
    kept = [s for s, k in zip(smis, mask) if k]
    d12 = pd.DataFrame(X[:, -12:], columns=desc)
    d12.insert(0, "smiles", kept)
    d12.insert(1, "ecfp4_bits_set", (X[:, :1024] > 0).sum(axis=1).astype(int))

    # The tested compounds at full width: every one of the 1036 values the model receives.
    ext = ROOT / "reviewer_package" / "08.09.26" / "MASTER_external_test_results.csv"
    tested = (pd.read_csv(ext)["smiles"].dropna().astype(str).drop_duplicates().tolist()
              if ext.exists() else [])
    Xt, mt = featurize(tested)
    keptt = [s for s, k in zip(tested, mt) if k]
    full = pd.DataFrame(Xt, columns=names)
    full.insert(0, "smiles", keptt)
    return d12, full, len(smis) - len(kept)


def runtime_dependencies():
    """Every file the served application opens to answer one query."""
    rows = [
        ("models_rf/<endpoint>.joblib", "model", "The fitted estimator for each endpoint. 71 files.",
         "produces the raw score for that endpoint"),
        ("models_rf/<endpoint>_calibrated.joblib", "model",
         "Calibration wrapped around the frozen estimator.",
         "converts the forest's vote share into a probability that means what it says"),
        ("models_rf/<endpoint>_binder.joblib", "model",
         "Binder classifier for targets whose measured set is almost entirely active.",
         "engagement call for the binder panel"),
        ("models_rf/adme/<endpoint>.joblib", "model", "The nine ADME and exposure models.",
         "absorption, distribution and unbound brain exposure"),
        ("models_rf/binder_modes.json", "runtime table",
         "Per endpoint: how it was trained, its decision threshold and the basis for that "
         "threshold.",
         "decides whether a probability is reported as engagement. Models trained against decoys "
         "and against measured inactives sit on different scales and must not share a cut"),
        ("models_rf/endpoint_context.json", "runtime table",
         "Per endpoint base rate of activity, and the distribution of each regressor's training "
         "values.",
         "turns a raw probability into enrichment over the base rate, which is what the interface "
         "reports"),
        ("models_rf/ad_reference.pkl", "runtime data",
         "158,890 background structures with fingerprints.",
         "applicability domain. A query far from all of them is flagged as outside the domain "
         "rather than scored confidently"),
        ("models_rf/ad_per_endpoint.pkl", "runtime data",
         "Per-endpoint fingerprints of that endpoint's own measured chemistry.",
         "the same test, endpoint by endpoint, since the panel's domains differ"),
        ("models_rf/readacross_index.pkl", "runtime data",
         "Measured actives with fingerprints and the targets they hit.",
         "shows the nearest measured compounds and what they are active at, so a prediction can be "
         "checked against real data"),
        ("models_rf/<endpoint>_meta.json", "metadata",
         "Training set size, positives, features, hyperparameters, CV summary.",
         "displayed beside each prediction as provenance"),
        ("models_manifest.json", "integrity",
         "SHA-256 of the archive and of all 195 extracted files.",
         "model_fetch.py refuses to start on any mismatch, so a corrupted download cannot reach a "
         "user"),
        ("results/tables/rf_cv_summary.csv", "display only",
         "Cross-validated scores.",
         "shown in the interface. Read for display; it takes no part in computing a prediction"),
    ]
    return pd.DataFrame(rows, columns=["file", "kind", "what_it_holds",
                                       "what_it_does_at_prediction_time"])


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()
    print(f"building {OUT.name}/\n")

    d1 = OUT / "01_HOW_THE_MODEL_WORKS"
    d2 = OUT / "02_TRAINING_DATA_AND_LABELS"
    d3 = OUT / "03_FEATURE_MATRIX_what_the_model_receives"
    d4 = OUT / "04_TRAINING_TESTING_VALIDATION"
    d5 = OUT / "05_INVERSION_TESTS"
    d6 = OUT / "06_WHAT_RUNS_A_PREDICTION"

    mi = model_inventory()
    write_csv(mi, d1 / "MODEL_INVENTORY_per_endpoint.csv",
              "Every endpoint with the model files deployed for it, what each file does, its "
              "estimator family and its decision threshold.",
              "Answering how many models exist per endpoint, and why an endpoint may have more "
              "than one file.",
              "models_rf/ and models_rf/binder_modes.json")

    names = feature_names()
    fd = pd.DataFrame({
        "position": range(1, len(names) + 1), "feature": names,
        "block": ["ECFP-4 fingerprint bit"] * 1024 + ["Physicochemical descriptor"] * 12,
        "type": ["binary, 0 or 1"] * 1024 + ["continuous"] * 12,
        "meaning": ["presence of a circular substructure of radius 2, hashed to this bit"] * 1024 + [
            "molecular weight", "calculated logP", "topological polar surface area",
            "hydrogen-bond donors", "hydrogen-bond acceptors", "rotatable bonds", "aromatic rings",
            "fraction of sp3 carbons", "ring count", "heavy atoms", "formal charge",
            "quantitative estimate of drug-likeness"],
    })
    write_csv(fd, d1 / "FEATURE_DEFINITIONS_all_1036.csv",
              "All 1036 model inputs in the exact order the estimator receives them: 1024 ECFP-4 "
              "fingerprint bits followed by 12 named physicochemical descriptors.",
              "Reading the feature matrix, and confirming that every endpoint shares one fixed "
              "representation.",
              "src/brainsafe/features/featurize.py")

    tl = training_and_labelling()
    write_csv(tl, d2 / "ENDPOINT_TRAINING_AND_LABELLING.csv",
              "Per endpoint: how many measurements, how many actives and inactives, how many "
              "decoys were generated, how many compounds the model was actually fitted to, and the "
              "exact rule that turned a measurement into a label.",
              "Auditing the labelling. The label rule differs by panel and this is where that is "
              "stated per endpoint rather than in prose.",
              "data/endpoints/, models_rf/binder_modes.json, models_rf/endpoint_context.json")

    d12, full, unparsed = feature_matrices()
    write_csv(d12, d3 / "FEATURE_MATRIX_descriptors_all_compounds.csv",
              f"Every one of {len(d12):,} distinct training structures as a row, with the 12 named "
              f"physicochemical descriptors as columns and a count of how many of the 1024 "
              f"fingerprint bits it sets. {unparsed} structures could not be parsed by RDKit and "
              f"are absent.",
              "Seeing the chemical space the models were fitted over, at a width that opens in a "
              "spreadsheet.",
              "data/endpoints/ featurised by src/brainsafe/features/featurize.py")
    write_csv(full, d3 / "FEATURE_MATRIX_all_1036_columns_tested_compounds.csv",
              f"The complete model input for the {len(full):,} externally tested structures: one "
              f"row per compound, all 1036 columns, holding the exact values passed to the "
              f"estimator.",
              "Confirming precisely what the model receives. This is the array that goes in; there "
              "is no further transformation between this file and the forest.",
              "the compounds in MASTER_external_test_results.csv, featurised identically")

    folds = pd.concat([pd.read_csv(p) for p in
                       [TAB / "manuscript_T2_per_fold.csv", TAB / "adme_cv_folds.csv",
                        TAB / "binder_cv_folds.csv"] if p.exists()], ignore_index=True)
    write_csv(folds, d4 / "CROSS_VALIDATION_every_fold.csv",
              f"Every individual fold of every endpoint: {len(folds):,} rows covering "
              f"{folds.endpoint.nunique()} endpoints in two split regimes, with training and test "
              f"sizes and the metrics for that fold alone.",
              "Checking the ten-fold claim directly. Each row is one model fitted on nine folds and "
              "scored on the tenth.",
              "results/tables/manuscript_T2_per_fold.csv, adme_cv_folds.csv, binder_cv_folds.csv")

    for src, doc in [
        (INV / "VERDICTS.csv",
         "The verdict on each falsification hypothesis H1 to H8: whether the attempt to break the "
         "tool succeeded, and the headline finding."),
        (INV / "H6_clinical_indication.csv",
         "H6: whether a drug's licensed indication is recovered, stratified by whether the drug's "
         "chemistry was seen in training."),
        (INV / "H7_target_discrimination.csv",
         "H7: per target, whether the model separates held-out actives from random chemistry."),
        (INV / "H8_panel_independence.csv",
         "H8: whether the panel's endpoints fire independently or are driven by a shared signal."),
        (INV / "H8_family_correlation.csv",
         "H8 detail: co-firing rates between homologous target pairs."),
        (INV / "H4_distant_specificity.csv",
         "H4: false-positive rate on chemistry increasingly unlike anything trained on."),
        (INV / "H5_readacross_value.csv",
         "H5: whether similarity read-across alone would do as well as the trained models."),
        (INV / "H1_disease_layer.csv",
         "H1: whether the disease-scoring layer adds anything over the endpoint scores."),
        (INV / "H2_weight_ablation.csv",
         "H2: sensitivity of disease ranking to the mechanism weights."),
        (INV / "H6_mesh_mapping_audit.csv",
         "The MeSH mapping behind H6, every panel condition to every matched heading, so a wrong "
         "match can be found rather than assumed absent."),
    ]:
        if src.exists():
            df = pd.read_csv(src)
            write_csv(df, d5 / src.name, doc,
                      "The falsification record. These analyses were designed to make the tool "
                      "look wrong; the ones that succeeded are reported as such.",
                      f"inversion/results/{src.name}")
    for extra in ["REPORT.md", "PLAN.md"]:
        s = ROOT / "inversion" / extra
        if s.exists():
            shutil.copy2(s, d5 / f"INVERSION_{extra}")

    rt = runtime_dependencies()
    write_csv(rt, d6 / "FILES_LOADED_AT_PREDICTION_TIME.csv",
              "Every file the served application opens to answer a single query, and what each one "
              "does at that moment.",
              "Separating what trained the model from what runs it. The training tables are not "
              "read at prediction time; the applicability-domain and calibration data are.",
              "the load functions in app.py")

    idx = pd.DataFrame(INDEX)
    with (OUT / "00_FILE_INDEX.csv").open("w", encoding="utf-8", newline="") as fh:
        fh.write("# WHAT THIS FILE CONTAINS: every file in AUDIT_PACKAGE, with what it holds and "
                 "what it is for.\n")
        fh.write("# WHAT IT IS USED FOR: finding the file that answers a given question.\n")
        fh.write("# GENERATED FROM: src/brainsafe/analysis/build_audit_package.py\n")
        idx.to_csv(fh, index=False, lineterminator="\n")

    n_models = int(mi.deployed_model_files.sum())
    (OUT / "00_READ_ME_FIRST.md").write_text(f"""# BrainSafe AI: audit package

Written for a reader who has to satisfy themselves the tool does what the manuscript says, without
reading the source. Every file here is generated by `src/brainsafe/analysis/build_audit_package.py`
from the deployed models and the tables that produced them. Every CSV opens with three comment lines
saying what it holds and what it is for. **Read every one of them with `skiprows=3`, not
`comment='#'`**: SMILES use `#` for a triple bond, and a comment character truncates every nitrile
and alkyne in the file, turning `C#N` into `C` without warning.

## The five things most likely to be misread, stated first

**One estimator family, one fixed input.** Every endpoint is a random forest over the same
1036-dimensional vector: 1024 ECFP-4 fingerprint bits and 12 named physicochemical descriptors.
There is no per-endpoint feature selection. `01_HOW_THE_MODEL_WORKS/FEATURE_DEFINITIONS_all_1036.csv`
lists all of them in the order the estimator receives them.

**{n_models} model files are deployed across {len(mi)} endpoints**: {(mi.deployed_model_files == 1).sum()}
endpoints carry one file and {(mi.deployed_model_files == 2).sum()} carry two. The second file is
never a competing model. It is either a calibration wrapper frozen around the base estimator, so the
reported number is a probability rather than a forest vote share, or a binder classifier used where
the measured set is almost entirely active and a plain classifier has no negative class to learn
from. `MODEL_INVENTORY_per_endpoint.csv` names the files and says which is which.

**Ten-fold cross-validation fitted models to measure, not to deploy.** Each endpoint was fitted 20
times during evaluation, ten in each split regime, each scored once on the fold withheld from it and
then discarded. The served model is refitted on the full set afterwards. Across
{folds.endpoint.nunique()} endpoints that is {len(folds):,} evaluation fits, every one of them a row
in `04_TRAINING_TESTING_VALIDATION/CROSS_VALIDATION_every_fold.csv`.

**No endpoint was trained on the panel total.** The totals in the manuscript are sums across
endpoints. Each model saw its own measured set alone, and those range over two orders of magnitude.
`02_TRAINING_DATA_AND_LABELS/ENDPOINT_TRAINING_AND_LABELLING.csv` gives the per-endpoint figure.

**The label rule is not uniform, and it should not be.** Core classifiers use the source's own
labels; receptor endpoints are regressed on pChEMBL directly because their sets are almost entirely
active and a binary split would be ill-posed; binder endpoints call actives at pChEMBL >= 7 against
a negative class of property-matched decoys plus half the measured inactives, with the other half
withheld to set the threshold. Each endpoint's rule is stated on its own row rather than inferred.

## What is in each folder

| Folder | Question it answers |
|---|---|
| `01_HOW_THE_MODEL_WORKS` | what the model is, how many files per endpoint, what all 1036 inputs are |
| `02_TRAINING_DATA_AND_LABELS` | what it was fitted to, from which source, and how a measurement became a label |
| `03_FEATURE_MATRIX_what_the_model_receives` | the input array itself, structures as rows and features as columns |
| `04_TRAINING_TESTING_VALIDATION` | every fold of every endpoint, both split regimes |
| `05_INVERSION_TESTS` | the record of trying to make the tool wrong, kept separate from the results |
| `06_WHAT_RUNS_A_PREDICTION` | the files opened to answer one query, which are not the files that trained it |

`00_FILE_INDEX.csv` lists every file with its purpose.

## Two notes on the feature matrix

`FEATURE_MATRIX_all_1036_columns_tested_compounds.csv` is the complete input for the externally
tested compounds: one row per structure, all 1036 columns, the exact values passed to the estimator.
Nothing happens between that file and the forest.

`FEATURE_MATRIX_descriptors_all_compounds.csv` covers every training structure but carries only the
12 named descriptors plus a count of fingerprint bits set. The full-width version of that file would
be 160,000 rows by 1036 columns, which is not a file anyone opens. The 12 descriptors are the
interpretable part; the fingerprint block is defined exhaustively in the feature definitions.

## On the inversion tests

These are kept in their own folder because a falsification record is a different kind of evidence
from a results table. Each hypothesis was an attempt to make the tool look wrong. Some succeeded:
the ones that did are reported as findings, not buried. `VERDICTS.csv` gives the verdict on each,
and `INVERSION_REPORT.md` the reasoning. The full analysis, including the scripts, is in
`inversion/` at the repository root.

## Provenance

Models: doi:10.5281/zenodo.21858576. Source: https://github.com/krishna-g-999/brainsafe-ai.
Regenerate this package with `python src/brainsafe/analysis/build_audit_package.py`.
""", encoding="utf-8")

    # Read every file back the way the package tells a reader to read it, and check the structures
    # survive. This exists because the first version of this package told readers to use
    # comment='#', which truncates every SMILES at its first triple bond: C#N reads back as C. The
    # data was right and the instruction was wrong, which is the harder failure to notice.
    problems = []
    for f in sorted(OUT.rglob("*.csv")):
        try:
            d = pd.read_csv(f, skiprows=3)
        except Exception as exc:
            problems.append(f"{f.name}: will not parse with skiprows=3 ({exc})")
            continue
        if len(d.columns) < 2:
            problems.append(f"{f.name}: parsed to one column, header line count is wrong")
        if "smiles" in d.columns:
            if d["smiles"].duplicated().any():
                problems.append(f"{f.name}: duplicate SMILES after round-trip")
            n_bad = sum(1 for s in d["smiles"].astype(str) if Chem.MolFromSmiles(s) is None)
            if n_bad:
                problems.append(f"{f.name}: {n_bad} SMILES no longer parse as molecules")
    print("\n  verification: " + ("every file re-reads cleanly and every structure survives"
                                  if not problems else f"{len(problems)} PROBLEMS"))
    for p in problems:
        print(f"    {p}")

    print(f"\n  {len(INDEX)} CSVs written, all with self-describing headers")
    print(f"  {n_models} deployed model files across {len(mi)} endpoints")
    print(f"  feature matrix: {len(full):,} tested compounds x 1036 columns")


if __name__ == "__main__":
    main()
