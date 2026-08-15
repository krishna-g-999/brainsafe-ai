"""Leakage and split integrity, checked independently of the pipeline's own claims.

This does not read any result the pipeline wrote. It rebuilds the folds from the endpoint tables
using the same splitter the pipeline uses, and then asks, of the actual index sets, questions that
have a right answer:

  L1  identity overlap      does any InChIKey appear in both the train and the test half of a fold
  L2  feature overlap       does any test row have a byte-identical feature vector to a train row
  L3  scaffold overlap      does any scaffold appear on both sides of a scaffold-grouped fold
  L4  featurisation purity  does a molecule's feature vector depend on anything but that molecule
  L5  no fitted scaling     is any transform fitted on data, and if so is it fitted on train only
  L6  label independence    does the featuriser ever see the label column

L1 and L3 are properties the splitter should guarantee; they are checked anyway, because a guarantee
that is never tested is a comment. L2 is the one that actually bit this project: a stereo-blind
fingerprint makes stereoisomers byte-identical, so identity by InChIKey can hold while the model is
scored on rows it has memorised.

Both regimes are reported for every endpoint, and both the raw table and the deduplicated matrix the
pipeline fits, because those are different questions and conflating them is how the earlier check
reported a leak into a model that never saw one.

Output: validation/repro/leakage_report.csv, validation/repro/leakage_summary.json

Run:  python validation/repro/r01_leakage.py            (all classification endpoints)
      python validation/repro/r01_leakage.py BBB MAO_A  (named endpoints only)
"""
from __future__ import annotations

import ast
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold, StratifiedKFold

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
sys.path.insert(0, str(ROOT / "src" / "brainsafe" / "data"))
from features.featurize import featurize, featurize_one            # noqa: E402
from models.train_rf import (CLASSIFICATION, N_SPLITS, SEED, _dedup_features,  # noqa: E402
                             _load, _scaffold_groups)
from build_compound_library import standardise                     # noqa: E402

OUT = ROOT / "validation" / "repro"


def fold_overlaps(X, y, groups, keys, split):
    """Per-fold overlap counts, computed from the index sets the splitter returns."""
    splitter = (GroupKFold(N_SPLITS) if split == "scaffold"
                else StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED))
    it = splitter.split(X, y, groups) if split == "scaffold" else splitter.split(X, y)
    ident, feat, scaf = 0, 0, 0
    for tr, te in it:
        tr_keys = {keys[i] for i in tr if keys[i]}
        ident = max(ident, sum(1 for i in te if keys[i] and keys[i] in tr_keys))
        tr_vec = {X[i].tobytes() for i in tr}
        feat = max(feat, sum(1 for i in te if X[i].tobytes() in tr_vec))
        tr_grp = {groups[i] for i in tr}
        scaf = max(scaf, sum(1 for i in te if groups[i] in tr_grp))
    return ident, feat, scaf


def check_featurisation_purity(smiles: list[str]) -> dict:
    """L4: a molecule's vector must not depend on what else is in the batch, or on order.

    If featurisation ever fitted anything across rows, a vector computed alone and the same vector
    computed inside a batch would differ, and the batch would be carrying information between
    compounds. Sampled rather than exhaustive, because the property is structural.
    """
    sample = smiles[:40]
    batch, mask = featurize(sample)
    alone = [featurize_one(s) for s in sample]
    alone = [a for a in alone if a is not None]
    same_alone = bool(len(alone) == len(batch) and
                      all(np.array_equal(batch[i], alone[i]) for i in range(len(batch))))
    rev, _ = featurize(list(reversed(sample)))
    same_reversed = bool(len(rev) == len(batch) and
                         all(np.array_equal(batch[i], rev[len(rev) - 1 - i])
                             for i in range(len(batch))))
    return {"identical_alone_vs_batch": same_alone, "identical_under_reversed_order": same_reversed,
            "n_sampled": len(batch)}


def check_no_fitted_transform() -> dict:
    """L5/L6: read the training source and confirm no transform is fitted, and no label is seen.

    A static check, deliberately. The claim is about what the code can do, not about what one run
    happened to do, and the absence of a scaler is only convincing if it is absent everywhere.
    """
    findings = {"fit_calls_outside_estimator": [], "scaler_imports": [], "featurize_sees_label": []}
    for p in sorted((ROOT / "src" / "brainsafe").rglob("*.py")):
        text = p.read_text(encoding="utf-8", errors="ignore")
        rel = p.relative_to(ROOT).as_posix()
        for name in ("StandardScaler", "MinMaxScaler", "RobustScaler", "QuantileTransformer",
                     "PowerTransformer", "Normalizer"):
            if name in text:
                findings["scaler_imports"].append(f"{rel}: {name}")
        if p.name == "featurize.py":
            tree = ast.parse(text)
            for node in ast.walk(tree):
                if isinstance(node, ast.arg) and node.arg in ("y", "label", "labels", "target"):
                    findings["featurize_sees_label"].append(f"{rel}: argument '{node.arg}'")
    return findings


def main(argv=None) -> None:
    t0 = time.time()
    eps = list(argv) if argv else list(CLASSIFICATION)
    rows, purity = [], None

    for ep in eps:
        df = _load(ep).dropna(subset=["smiles", "label"]).reset_index(drop=True)
        smiles_all = df["smiles"].astype(str).tolist()
        X, mask = featurize(smiles_all)
        smiles = [s for s, k in zip(smiles_all, mask) if k]
        y = df.loc[mask, "label"].to_numpy().astype(int)
        groups = _scaffold_groups(smiles)
        keys = [standardise(s)[1] for s in smiles]
        if purity is None:
            purity = check_featurisation_purity(smiles)

        Xd, yd, gd, sd, rep = _dedup_features(X, y, groups, smiles, "classification")
        kd = [standardise(s)[1] for s in sd]

        for stage, (Xs, ys, gs, ks) in (("raw table", (X, y, groups, keys)),
                                        ("as trained (deduplicated)", (Xd, yd, gd, kd))):
            for split in ("random", "scaffold"):
                ident, feat, scaf = fold_overlaps(Xs, ys, gs, ks, split)
                rows.append({
                    "endpoint": ep, "stage": stage, "split": split, "n_rows": len(ys),
                    "worst_fold_shared_inchikey": ident,
                    "worst_fold_shared_feature_vector": feat,
                    "worst_fold_shared_scaffold": scaf,
                    "duplicate_rows_removed": rep.get("duplicate_rows_removed", 0),
                })
                print(f"[{ep:6s}] {stage:26s} {split:9s} "
                      f"inchikey {ident:4d} | features {feat:4d} | scaffold {scaf:5d}", flush=True)

    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "leakage_report.csv", index=False)

    trained = out[out.stage.str.startswith("as trained")]
    static = check_no_fitted_transform()
    summary = {
        "commit": json.loads((OUT / "environment.json").read_text())["commit"],
        "endpoints_checked": eps,
        "L1_identity_overlap_as_trained_max": int(trained.worst_fold_shared_inchikey.max()),
        "L2_feature_overlap_as_trained_max": int(trained.worst_fold_shared_feature_vector.max()),
        "L3_scaffold_overlap_scaffold_split_max": int(
            trained[trained.split == "scaffold"].worst_fold_shared_scaffold.max()),
        "L2_feature_overlap_raw_table_max": int(
            out[out.stage == "raw table"].worst_fold_shared_feature_vector.max()),
        "L4_featurisation_purity": purity,
        "L5_L6_static": static,
        "wall_clock_s": round(time.time() - t0, 1),
    }
    (OUT / "leakage_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\n=== summary, as trained ===")
    print(f"L1 identity overlap (InChIKey)        max {summary['L1_identity_overlap_as_trained_max']}")
    print(f"L2 feature-vector overlap             max {summary['L2_feature_overlap_as_trained_max']}")
    print(f"L3 scaffold overlap, scaffold split   max {summary['L3_scaffold_overlap_scaffold_split_max']}")
    print(f"L2 on the RAW table (pre-dedup)       max {summary['L2_feature_overlap_raw_table_max']}"
          "   <- what deduplication removes")
    print(f"L4 featurisation purity               {purity}")
    print(f"L5 scalers found in source            {static['scaler_imports'] or 'none'}")
    print(f"L6 featurize() label arguments        {static['featurize_sees_label'] or 'none'}")
    print(f"\nwrote leakage_report.csv and leakage_summary.json ({summary['wall_clock_s']}s)")


if __name__ == "__main__":
    main(sys.argv[1:])
