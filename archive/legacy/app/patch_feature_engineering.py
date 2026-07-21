"""
scripts/patch_feature_engineering.py
Updates the feature pipeline in ml_v5_training.py and ml_v5_engine.py
to include 6 new structural features, making the total 93.

The 6 new features (appended after the existing 87):
  phenolic_oh_count, catechol_flag, aromatic_ring_count,
  hbd_count, tpsa, rotatable_bonds

These are pre-computed in brainsafe_SCIENTIFIC_FIXED.csv — they are
NOT computed on-the-fly during training. The feature matrix builder
just reads them from the DataFrame.

Run BEFORE training:
  D:\\BRAINSAFE_AI\\brainsafe_env\\Scripts\\python.exe scripts\\patch_feature_engineering.py
"""

import sys, re
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

NEW_FEATURE_COLS = [
    "phenolic_oh_count", "catechol_flag", "aromatic_ring_count",
    "hbd_count", "tpsa", "rotatable_bonds",
]
N_NEW_FEATURES    = len(NEW_FEATURE_COLS)   # 6
N_OLD_FEATURES    = 87
N_TOTAL_FEATURES  = N_OLD_FEATURES + N_NEW_FEATURES   # 93


def patch_training_script():
    """Patch ml_v5_training.py to include structural features."""
    path = ROOT / "ml_v5_training.py"
    if not path.exists():
        print(f"  NOT FOUND: {path}")
        return False

    content = path.read_text(encoding="utf-8", errors="replace")

    # 1. Update N_FEATURE_TOTAL constant
    content = re.sub(
        r'N_FEATURE_TOTAL\s*=\s*87',
        f'N_FEATURE_TOTAL = {N_TOTAL_FEATURES}',
        content
    )

    # 2. Update the assertion message
    content = content.replace(
        'f"Feature count mismatch: expected {N_FEATURE_TOTAL}, got {X.shape[1]}"',
        f'f"Feature count mismatch: expected {N_TOTAL_FEATURES}, got {{X.shape[1]}}"'
    )

    # 3. Find the build_feature_matrix function and add structural features
    # We inject the 6 new features just before the final np.hstack
    old_hstack = (
        "    X = np.hstack([ecfp_reduced, chemberta_reduced, disease_feats, bbb_vec])\n"
        "    assert X.shape[1] == N_FEATURE_TOTAL"
    )
    new_hstack = (
        "    # 6 new structural features (pre-computed in CSV)\n"
        "    struct_feats = np.zeros((len(df), 6), dtype=np.float32)\n"
        "    _struct_cols = [\n"
        "        'phenolic_oh_count', 'catechol_flag', 'aromatic_ring_count',\n"
        "        'hbd_count', 'tpsa', 'rotatable_bonds',\n"
        "    ]\n"
        "    for _si, _sc in enumerate(_struct_cols):\n"
        "        if _sc in df.columns:\n"
        "            struct_feats[:, _si] = pd.to_numeric(\n"
        "                df[_sc], errors='coerce').fillna(0).values\n"
        "    # Normalise: phenolic_oh (0-10), catechol (0-1), arings (0-10),\n"
        "    #            hbd (0-15), tpsa (0-200), rotbonds (0-20)\n"
        "    _norm = np.array([10.0, 1.0, 10.0, 15.0, 200.0, 20.0])\n"
        "    struct_feats = struct_feats / _norm\n"
        "\n"
        "    X = np.hstack([ecfp_reduced, chemberta_reduced,\n"
        "                   disease_feats, bbb_vec, struct_feats])\n"
        "    assert X.shape[1] == N_FEATURE_TOTAL"
    )

    if old_hstack in content:
        content = content.replace(old_hstack, new_hstack)
        print("  ✓ Patched build_feature_matrix (added 6 structural features)")
    else:
        # Try alternative form
        old2 = "    X = np.hstack([ecfp_reduced, chemberta_reduced, disease_feats, bbb_vec])"
        new2 = (
            "    _struct_cols = ['phenolic_oh_count','catechol_flag',\n"
            "        'aromatic_ring_count','hbd_count','tpsa','rotatable_bonds']\n"
            "    struct_feats = np.zeros((len(df), 6), dtype=np.float32)\n"
            "    for _si,_sc in enumerate(_struct_cols):\n"
            "        if _sc in df.columns:\n"
            "            struct_feats[:,_si]=pd.to_numeric(df[_sc],errors='coerce').fillna(0).values\n"
            "    struct_feats = struct_feats / np.array([10.,1.,10.,15.,200.,20.])\n"
            "    X = np.hstack([ecfp_reduced, chemberta_reduced, disease_feats, bbb_vec, struct_feats])"
        )
        if old2 in content:
            content = content.replace(old2, new2)
            print("  ✓ Patched build_feature_matrix (alternative form)")
        else:
            print("  ⚠️  Could not find hstack line — check ml_v5_training.py manually")
            print("     Manually add struct_feats before the np.hstack call")
            print("     See patch comment in scripts/patch_feature_engineering.py")

    path.write_text(content, encoding="utf-8")
    print(f"  N_FEATURE_TOTAL updated: 87 → {N_TOTAL_FEATURES}")
    return True


def patch_model_config():
    """Update N_FEATURE_TOTAL in model_config.py."""
    path = ROOT / "model_config.py"
    if not path.exists():
        print(f"  NOT FOUND: {path}")
        return
    content = path.read_text(encoding="utf-8", errors="replace")
    content = re.sub(
        r'N_FEATURE_TOTAL\s*:\s*int\s*=\s*87',
        f'N_FEATURE_TOTAL: int = {N_TOTAL_FEATURES}',
        content
    )
    content = re.sub(r'N_FEATURE_TOTAL = 87', f'N_FEATURE_TOTAL = {N_TOTAL_FEATURES}', content)
    path.write_text(content, encoding="utf-8")
    print(f"  ✓ model_config.py: N_FEATURE_TOTAL → {N_TOTAL_FEATURES}")


def verify_patch():
    """Verify the patch applied correctly."""
    path = ROOT / "ml_v5_training.py"
    content = path.read_text(encoding="utf-8", errors="replace")
    if f"N_FEATURE_TOTAL = {N_TOTAL_FEATURES}" in content:
        print(f"  ✓ N_FEATURE_TOTAL = {N_TOTAL_FEATURES} confirmed in ml_v5_training.py")
    else:
        print(f"  ⚠️  N_FEATURE_TOTAL not updated to {N_TOTAL_FEATURES}")
    if "phenolic_oh_count" in content:
        print("  ✓ phenolic_oh_count injection confirmed")
    else:
        print("  ⚠️  phenolic_oh_count not found in ml_v5_training.py")
        print("     Add manually: include struct_feats in the np.hstack call")


if __name__ == "__main__":
    print("=" * 62)
    print("Patching feature engineering pipeline (87 → 93 features)")
    print("=" * 62)
    patch_training_script()
    patch_model_config()
    print()
    verify_patch()
    print()
    print("=" * 62)
    print("NEXT — run scientific_fix_v2.py, then retrain:")
    print()
    print("  1. python scripts\\scientific_fix_v2.py")
    print("     > logs\\scientific_fix.log 2>&1")
    print()
    print("  2. Start-Process (training on brainsafe_SCIENTIFIC_FIXED.csv)")
    print("=" * 62)
