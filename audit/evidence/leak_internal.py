"""Within-endpoint duplicate check.

For each endpoint training table, count how many rows are duplicates of another row in the
SAME table under three equivalence relations, and whether duplicates carry conflicting labels.
Rows that are duplicates under the relation the model actually perceives (identical ECFP-4 +
identical descriptors) cannot be separated by any split that does not group them, so under
random K-fold they appear in both train and test.
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
ROOT = Path(r"D:\BRAINSAFE_AI")
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
sys.path.insert(0, str(ROOT / "src" / "brainsafe" / "data"))
from features.featurize import featurize  # noqa: E402
from build_compound_library import standardise  # noqa: E402

EP_DIR = ROOT / "data" / "endpoints"
rows = []
for path in sorted(EP_DIR.glob("*.csv")):
    ep = path.stem
    df = pd.read_csv(path)
    if "smiles" not in df.columns or "label" not in df.columns:
        continue
    n = len(df)
    smis = df["smiles"].astype(str).tolist()
    labels = df["label"].tolist()

    # exact SMILES-string duplicates (what a naive dedup would catch)
    n_smiles_dup = n - df["smiles"].nunique()

    # InChIKey skeleton duplicates (same constitution, ignoring stereo/salt/protonation)
    skel = defaultdict(list)
    for s, y in zip(smis, labels):
        _, ik = standardise(s)
        if ik:
            skel[ik.split("-")[0]].append(y)
    n_skel_unique = len(skel)
    n_skel_dup = sum(len(v) - 1 for v in skel.values() if len(v) > 1)
    skel_conflict = sum(1 for v in skel.values() if len(set(v)) > 1)

    # feature-identical duplicates: exactly what the model receives
    X, mask = featurize(smis)
    kept_labels = [y for y, m in zip(labels, mask) if m]
    feat = defaultdict(list)
    for vec, y in zip(X, kept_labels):
        feat[vec.tobytes()].append(y)
    n_feat_unique = len(feat)
    n_feat_dup = sum(len(v) - 1 for v in feat.values() if len(v) > 1)
    feat_conflict = sum(1 for v in feat.values() if len(set(v)) > 1)
    feat_conflict_rows = sum(len(v) for v in feat.values() if len(set(v)) > 1)

    rows.append({
        "endpoint": ep, "rows": n,
        "dup_smiles_string": n_smiles_dup,
        "unique_skeleton": n_skel_unique, "dup_skeleton": n_skel_dup,
        "skeleton_groups_with_conflicting_labels": skel_conflict,
        "featurised": int(mask.sum()),
        "unique_feature_vectors": n_feat_unique,
        "dup_feature_vectors": n_feat_dup,
        "dup_feature_pct": round(100.0 * n_feat_dup / max(int(mask.sum()), 1), 2),
        "feature_groups_with_conflicting_labels": feat_conflict,
        "rows_in_conflicting_feature_groups": feat_conflict_rows,
    })
    r = rows[-1]
    print(f"{ep:16s} n={n:6d} featdup={n_feat_dup:5d} ({r['dup_feature_pct']:5.2f}%) "
          f"skeldup={n_skel_dup:5d} conflict_groups={feat_conflict:4d}", flush=True)

out = pd.DataFrame(rows)
out.to_csv(Path(__file__).with_name("duplicate_audit.csv"), index=False)
print()
print("TOTAL rows:", out["rows"].sum())
print("TOTAL feature-identical duplicate rows:", out["dup_feature_vectors"].sum())
print("TOTAL feature groups with conflicting labels:",
      out["feature_groups_with_conflicting_labels"].sum())
print("wrote duplicate_audit.csv")
