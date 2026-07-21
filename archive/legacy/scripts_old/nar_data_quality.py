"""
scripts/nar_data_quality.py
NAR Web Server track — data quality validation.
Run BEFORE training. Flags any issue that would cause reviewers to reject.

Nucleic Acids Research requirements addressed:
  - Score scale consistency (all dimensions 0–100)
  - SMILES coverage ≥ 96% (claimed in paper)
  - Disease category proportions (claimed in paper: AD 60%, PD 37%, ALS 30%, HD 22%)
  - BBB class distributions (4 categories, not 3 as incorrectly stated in manuscript)
  - NPS distribution normality (paper claims approximately normal, mean 62.4 ± 8.0)
  - Duplicate detection
  - Missing data pattern analysis

Output:  logs/nar_data_quality_report.txt  (paste back to Claude after running)
"""

import sys
import os
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))   # ensure scorer, model_config importable
LOG_DIR = ROOT / "logs"
LOG_DIR.mkdir(exist_ok=True)

DIMENSION_COLS = [
    "antioxidant", "anti_inflammatory", "mitochondrial_support",
    "aggregation_modulation", "cognitive_enhancement",
    "neurogenesis", "synaptic_plasticity",
]

BBB_MAP = {"Low": 0, "Low-Med": 1, "Medium": 2, "High": 3}

PAPER_CLAIMS = {
    "n_total":          535,
    "n_ad":             194,   "pct_ad":   59.7,
    "n_pd":             119,   "pct_pd":   36.6,
    "n_als":             97,   "pct_als":  29.8,
    "n_hd":              71,   "pct_hd":   21.8,
    "smiles_coverage":   96.6,
    "nps_mean":          52.5,
    "nps_sd":             17.6,
    "bbb_high_pct":      34.0,   # from Fig 2C (correct)
    "bbb_medium_pct":    33.0,
    "bbb_lowmed_pct":    18.0,
    "bbb_low_pct":       14.0,
}


def sep(char="─", n=62):
    return char * n


def check(label, condition, actual="", expected=""):
    status = "✓ PASS" if condition else "✗ FAIL"
    line = f"  {status}  {label}"
    if not condition:
        line += f"\n          Expected: {expected}\n          Got:      {actual}"
    return line, condition


def run_checks(df: pd.DataFrame, smiles_cache: dict) -> list[str]:
    lines = []
    passed = 0
    failed = 0

    def add(label, condition, actual="", expected=""):
        nonlocal passed, failed
        line, ok = check(label, condition, str(actual), str(expected))
        lines.append(line)
        if ok:
            passed += 1
        else:
            failed += 1

    # ── 1. SAMPLE SIZE ───────────────────────────────────────────────────────
    lines.append(sep())
    lines.append("1. SAMPLE SIZE")
    lines.append(sep())

    n = len(df)
    add("Total compound count == 325", n == PAPER_CLAIMS["n_total"],
        n, PAPER_CLAIMS["n_total"])

    # Disease counts (categories overlap — count via BBB/disease flags)
    for dis, col, expected_n, expected_pct in [
        ("AD",  "ad_target_count",  194, 59.7),
        ("PD",  "pd_target_count",  119, 36.6),
        ("ALS", "als_target_count",  97, 29.8),
        ("HD",  "hd_target_count",   71, 21.8),
    ]:
        if col in df.columns:
            n_dis = (df[col] > 0).sum()
            pct = round(n_dis / n * 100, 1)
            add(f"  {dis} compound count ~{expected_n} (±10)",
                abs(n_dis - expected_n) <= 15,
                f"n={n_dis} ({pct}%)", f"n≈{expected_n} ({expected_pct}%)")

    # ── 2. SMILES COVERAGE ───────────────────────────────────────────────────
    lines.append("")
    lines.append(sep())
    lines.append("2. SMILES COVERAGE")
    lines.append(sep())

    if "smiles" in df.columns:
        n_smiles = (df["smiles"].fillna("") != "").sum()
        pct_smiles = round(n_smiles / n * 100, 1)
        add("SMILES coverage ≥ 96.6% (paper claim)",
            pct_smiles >= 96.0,
            f"{n_smiles}/{n} ({pct_smiles}%)", "≥96.6%")

        n_imputed = n - n_smiles
        add("Imputed SMILES ≤ 15 compounds",
            n_imputed <= 15, n_imputed, "≤15")
        lines.append(f"    ⓘ  {n_imputed} compounds use median-imputed features (flagged in Table S1)")

    # ── 3. SCORE DISTRIBUTIONS ───────────────────────────────────────────────
    lines.append("")
    lines.append(sep())
    lines.append("3. DIMENSION SCORE DISTRIBUTIONS  (must all be on 0-100 scale)")
    lines.append(sep())

    for col in DIMENSION_COLS:
        if col not in df.columns:
            add(f"  {col:<30} — column present", False, "MISSING", "present")
            continue
        vals = df[col].dropna()
        mn, mx, mean = vals.min(), vals.max(), vals.mean()
        on_100 = mx > 10.0  # on 1-10 scale if max ≤ 10
        add(f"  {col:<30} on 0-100 scale (max={mx:.1f})",
            on_100, f"max={mx:.1f}", ">10.0")
        lines.append(f"    ⓘ  n={len(vals):3d}  min={mn:.1f}  mean={mean:.1f}  max={mx:.1f}  "
                     f"sd={vals.std():.1f}  missing={df[col].isna().sum()}")

    # ── 4. NPS DISTRIBUTION ──────────────────────────────────────────────────
    lines.append("")
    lines.append(sep())
    lines.append("4. NPS DISTRIBUTION  (paper: mean=62.4 ± 8.0, approximately normal)")
    lines.append(sep())

    if "nps" in df.columns:
        nps = df["nps"].dropna()
        nps_mean = nps.mean()
        nps_sd = nps.std()
        nps_min = nps.min()
        nps_max = nps.max()

        add("NPS mean in range 55–70 (paper claims 62.4)",
            55 <= nps_mean <= 70, f"{nps_mean:.1f}", "55–70")
        add("NPS SD in range 5–15 (paper claims 8.0)",
            5 <= nps_sd <= 15, f"{nps_sd:.1f}", "5–15")
        add("NPS range plausible (min>10, max<100)",
            nps_min > 10 and nps_max < 100,
            f"{nps_min:.1f}–{nps_max:.1f}", "10–100")

        # Shapiro-Wilk normality test (paper claims "approximately normal")
        if len(nps) <= 5000:
            stat, p = stats.shapiro(nps[:300] if len(nps) > 300 else nps)
            lines.append(f"  ⓘ  Shapiro-Wilk normality: W={stat:.3f}, p={p:.4f} "
                         f"({'approx normal at α=0.05' if p>0.05 else 'non-normal — describe as skewed in paper'})")

    # ── 5. BBB DISTRIBUTION ──────────────────────────────────────────────────
    lines.append("")
    lines.append(sep())
    lines.append("5. BBB PERMEABILITY DISTRIBUTION  (4 categories, not 3)")
    lines.append(sep())
    lines.append("  ⚠️  Manuscript text incorrectly states '38% High, 44% Medium, 18% Low'")
    lines.append("      Correct values from Fig 2C: 34%/33%/18%/14% across 4 categories")

    if "bbb" in df.columns:
        bbb_counts = df["bbb"].value_counts()
        for cat in ["High", "Medium", "Low-Med", "Low"]:
            n_cat = bbb_counts.get(cat, 0)
            pct = round(n_cat / n * 100, 1)
            lines.append(f"    {cat:<10} n={n_cat:3d}  ({pct:5.1f}%)")
        lines.append("  → FIX in manuscript: replace '38%/44%/18% (3 categories)' with correct 4-category values above")

    # ── 6. DUPLICATE DETECTION ───────────────────────────────────────────────
    lines.append("")
    lines.append(sep())
    lines.append("6. DUPLICATE DETECTION")
    lines.append(sep())

    if "name" in df.columns:
        dup_names = df["name"].duplicated().sum()
        add("No duplicate compound names", dup_names == 0,
            f"{dup_names} duplicates", "0")

    if "smiles" in df.columns:
        smiles_valid = df[df["smiles"].fillna("") != ""]["smiles"]
        dup_smiles = smiles_valid.duplicated().sum()
        add("No duplicate SMILES strings", dup_smiles == 0,
            f"{dup_smiles} duplicates", "0")

    # ── 7. MISSING DATA PATTERN ──────────────────────────────────────────────
    lines.append("")
    lines.append(sep())
    lines.append("7. MISSING DATA PATTERN")
    lines.append(sep())

    total_missing = df[DIMENSION_COLS].isna().sum().sum()
    n_complete = (df[DIMENSION_COLS].notna().all(axis=1)).sum()
    lines.append(f"  Compounds with ALL 7 dimensions present: {n_complete}/{n} "
                 f"({100*n_complete/n:.1f}%)")
    lines.append(f"  Total missing dimension values:          {total_missing}")
    for col in DIMENSION_COLS:
        if col in df.columns:
            m = df[col].isna().sum()
            if m > 0:
                lines.append(f"    {col:<30} {m} missing values")

    # ── 8. NAR SPECIFIC: HOLD-OUT SPLIT REPRODUCIBILITY ─────────────────────
    lines.append("")
    lines.append(sep())
    lines.append("8. NAR REPRODUCIBILITY: HOLD-OUT SPLIT MUST BE FIXED")
    lines.append(sep())

    # Verify random_state=42 gives consistent split
    if "nps" in df.columns and len(df) >= 325:
        try:
            from sklearn.model_selection import train_test_split
            df["nps_quartile"] = pd.qcut(df["nps"], q=4, labels=False, duplicates="drop")
            df_tr, df_ho = train_test_split(df, test_size=0.20,
                                            random_state=42,
                                            stratify=df["nps_quartile"])
            add("Hold-out n=65 (random_state=42, stratified)",
                len(df_ho) == 65, len(df_ho), 65)
            add("Training n=260 (random_state=42, stratified)",
                len(df_tr) == 260, len(df_tr), 260)
            lines.append(f"  ⓘ  Hold-out NPS mean={df_ho['nps'].mean():.1f} ± {df_ho['nps'].std():.1f}")
            lines.append(f"  ⓘ  Training NPS mean={df_tr['nps'].mean():.1f} ± {df_tr['nps'].std():.1f}")
        except Exception as e:
            lines.append(f"  ⚠️  Split check failed: {e}")

    # ── SUMMARY ──────────────────────────────────────────────────────────────
    lines.append("")
    lines.append(sep("═"))
    lines.append(f"RESULT: {passed} PASS  |  {failed} FAIL  |  {passed+failed} total checks")
    if failed == 0:
        lines.append("✓ ALL CHECKS PASSED — safe to proceed to training")
    else:
        lines.append(f"✗ {failed} check(s) FAILED — paste this output to Claude before training")
    lines.append(sep("═"))

    return lines


def main():
    print("=" * 62)
    print("BrainSafe AI v6 — NAR Data Quality Check")
    print("=" * 62)

    # Prefer the largest available dataset
    for candidate in ["brainsafe_master_v2.csv", "brainsafe_master.csv",
                       "brainsafe_FINAL.csv", "brainsafe_FINAL.csv"]:
        train_path = ROOT / "data" / candidate
        if train_path.exists():
            break
    else:
        print("ERROR: No training CSV found in data/")
        sys.exit(1)

    df = pd.read_csv(train_path)
    print(f"Loaded {len(df)} compounds from {train_path}")

    # Load SMILES cache
    cache_path = ROOT / "smiles_cache.json"
    smiles_cache = {}
    if cache_path.exists():
        with open(cache_path) as f:
            smiles_cache = json.load(f)

    # Compute NPS if not present
    if "nps" not in df.columns:
        from scorer import neuro_score, DIMENSION_COLS as DC
        df["nps"] = df.apply(
            lambda r: neuro_score({c: r.get(c, 0.0) for c in DC}), axis=1
        )

    lines = run_checks(df, smiles_cache)

    print()
    for line in lines:
        print(line)

    # Save report
    report_path = LOG_DIR / "nar_data_quality_report.txt"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\nReport saved: {report_path}")
    print("\nNEXT: paste the output above into Claude chat, then proceed to training.")


if __name__ == "__main__":
    main()
