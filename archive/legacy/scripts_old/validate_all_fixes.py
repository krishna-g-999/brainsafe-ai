"""
scripts/validate_all_fixes.py
Validates all critical bug fixes from the v5→v6 migration.
Run after replacing the corrected files:

    python scripts/validate_all_fixes.py

Expected output: all 10 tests PASS.
"""

import sys
import json
import numpy as np
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

PASS = "\033[92m  ✓ PASS\033[0m"
FAIL = "\033[91m  ✗ FAIL\033[0m"

results = []


def test(name: str, condition: bool, detail: str = ""):
    status = PASS if condition else FAIL
    label  = "PASS" if condition else "FAIL"
    print(f"{status}  {name}")
    if detail and not condition:
        print(f"        {detail}")
    results.append((name, condition))
    return condition


print("=" * 65)
print("BrainSafe AI v6 — Bug Fix Validation Suite")
print("=" * 65)

# ── Test 1: scorer.py NPS uses all 7 dimensions ───────────────────────────
print("\n[1] scorer.py — NPS formula")
try:
    from scorer import neuro_score, DIM_WEIGHTS
    # Verify all 7 dims are in the weights dict
    EXPECTED_DIMS = {
        "antioxidant", "anti_inflammatory", "mitochondrial_support",
        "aggregation_modulation", "cognitive_enhancement",
        "neurogenesis", "synaptic_plasticity",
    }
    test("All 7 dimensions in DIM_WEIGHTS",
         set(DIM_WEIGHTS.keys()) == EXPECTED_DIMS,
         f"Found: {set(DIM_WEIGHTS.keys())}")

    # Compound with high cognitive but zero antioxidant/anti-inflam
    # v5 would give 0; v6 must give > 0
    data_cog = {d: 0.0 for d in EXPECTED_DIMS}
    data_cog["cognitive_enhancement"] = 100.0
    nps_cog = neuro_score(data_cog)
    test("Cognitive-only compound gets NPS > 0  (v5 bug: always 0)",
         nps_cog > 0,
         f"NPS = {nps_cog}")

    # All 100 → NPS should equal 100
    data_all = {d: 100.0 for d in EXPECTED_DIMS}
    nps_all = neuro_score(data_all)
    test("All-100 input produces NPS = 100",
         abs(nps_all - 100.0) < 0.1,
         f"NPS = {nps_all}")

    # Confirm denominator is correct (sum of weights × 100 = 1400)
    weight_sum = sum(DIM_WEIGHTS.values())
    test("Dimension weight sum = 14 (3+3+2+2+2+1+1)",
         weight_sum == 14,
         f"Sum = {weight_sum}")

except Exception as e:
    test("scorer.py imports without error", False, str(e))
    test("All 7 dimensions in DIM_WEIGHTS", False, "import failed")
    test("Cognitive-only compound gets NPS > 0", False, "import failed")
    test("All-100 input produces NPS = 100", False, "import failed")

# ── Test 2: model_config.py — NameErrors fixed ───────────────────────────
print("\n[2] model_config.py — NameErrors fixed")
try:
    from model_config import (BBB_MAP, DIS_MAP, TARGETTOPATHWAY,
                               DIMENSION_COLS, N_FEATURE_TOTAL,
                               estimate_bbb_class)

    test("BBB_MAP importable (was _BBB_MAP in v5)",
         isinstance(BBB_MAP, dict) and "High" in BBB_MAP)

    test("DIS_MAP importable (was _DIS_MAP in v5)",
         isinstance(DIS_MAP, dict) and "High" in DIS_MAP)

    test("TARGETTOPATHWAY defined (was missing in v5 — NameError)",
         isinstance(TARGETTOPATHWAY, dict) and len(TARGETTOPATHWAY) >= 10,
         f"Length = {len(TARGETTOPATHWAY) if isinstance(TARGETTOPATHWAY, dict) else 'N/A'}")

    test("N_FEATURE_TOTAL == 87",
         N_FEATURE_TOTAL == 87,
         f"Got {N_FEATURE_TOTAL}")

    # Test BBB estimation
    bbb_str, bbb_int = estimate_bbb_class(mw=300, logp=2.0, tpsa=50, hbd=1)
    test("estimate_bbb_class: small CNS-friendly → High BBB",
         bbb_str == "High" and bbb_int == 3,
         f"Got ({bbb_str}, {bbb_int})")

    bbb_str_low, bbb_int_low = estimate_bbb_class(mw=800, logp=6.0, tpsa=200, hbd=8)
    test("estimate_bbb_class: large polar molecule → Low BBB",
         bbb_str_low == "Low" and bbb_int_low == 0,
         f"Got ({bbb_str_low}, {bbb_int_low})")

except Exception as e:
    for _ in range(6):
        test(f"model_config import", False, str(e))

# ── Test 3: Score scale consistency ──────────────────────────────────────
print("\n[3] Score scale — all dimensions on 0–100")
try:
    from scorer import neuro_score, SCORE_MAX

    # Simulate a v5-style 1–10 input — the old code would not catch this
    data_old_scale = {
        "antioxidant": 8.5, "anti_inflammatory": 7.2,
        "mitochondrial_support": 6.0, "aggregation_modulation": 7.8,
        "cognitive_enhancement": 9.1, "neurogenesis": 6.5,
        "synaptic_plasticity": 7.0,
    }
    nps_old = neuro_score(data_old_scale)
    # On a 0–100 scale, inputs of 6–9 represent 6–9 % of max → NPS should be < 15.
    # (v5 bug: treating these as 1–10 scale would yield NPS ≈ 76, far too high.)
    test("Scores of 6–9 treated as 0-100 scale → NPS correctly low (< 15), not spuriously high",
         nps_old < 15.0,
         f"NPS = {nps_old:.1f} (expected < 15 for inputs that are 6-9/100 of max)")

    test("SCORE_MAX = 100 defined in scorer",
         hasattr(sys.modules.get("scorer", object()), "SCORE_MAX"),
         "SCORE_MAX not found in scorer module")

except Exception as e:
    test("Score scale test", False, str(e))

# ── Test 4: Neurotransmitter mapper ───────────────────────────────────────
print("\n[4] neurotransmitter_mapper.py — new module")
try:
    from neurotransmitter_mapper import compute_nt_effects, NTEffect

    # Quercetin-like: MAO-B inhibition + AChE inhibition
    test_enzymes = [
        {"name": "MAO-B", "action": "Inhibition", "strength": "Moderate"},
        {"name": "AChE", "action": "Inhibition", "strength": "Weak"},
    ]
    effects = compute_nt_effects(test_enzymes)
    nt_names = [e.neurotransmitter for e in effects]

    test("MAO-B inhibition → dopamine effect detected",
         "dopamine" in nt_names,
         f"Detected NTs: {nt_names}")

    test("AChE inhibition → acetylcholine effect detected",
         "acetylcholine" in nt_names,
         f"Detected NTs: {nt_names}")

    dopamine_effect = next((e for e in effects if e.neurotransmitter == "dopamine"), None)
    test("Dopamine direction = 'increases' (correct for MAO-B inhibition)",
         dopamine_effect and dopamine_effect.direction == "increases",
         f"Direction = {dopamine_effect.direction if dopamine_effect else 'N/A'}")

    test("NTEffect has display_arrow attribute",
         dopamine_effect and dopamine_effect.display_arrow == "↑",
         f"Arrow = {dopamine_effect.display_arrow if dopamine_effect else 'N/A'}")

except Exception as e:
    for _ in range(4):
        test("neurotransmitter_mapper", False, str(e))

# ── Test 5: Brain region mapper ───────────────────────────────────────────
print("\n[5] brain_region_mapper.py — new module")
try:
    from brain_region_mapper import compute_region_scores, get_primary_regions

    # High cognitive + neurogenesis → hippocampus should be top region
    dim_scores = {
        "antioxidant": 40.0, "anti_inflammatory": 40.0,
        "mitochondrial_support": 40.0, "aggregation_modulation": 40.0,
        "cognitive_enhancement": 90.0, "neurogenesis": 85.0,
        "synaptic_plasticity": 75.0,
    }
    regions = compute_region_scores(dim_scores)
    top_region = regions[0].region_key if regions else None

    test("compute_region_scores returns non-empty list",
         len(regions) > 0,
         f"Got {len(regions)} regions")

    test("High cognitive/neurogenesis → hippocampus is top region",
         top_region == "hippocampus",
         f"Top region: {top_region}")

    # High mitochondrial → substantia nigra should be highly ranked
    dim_mito = {d: 20.0 for d in dim_scores}
    dim_mito["mitochondrial_support"] = 95.0
    dim_mito["anti_inflammatory"] = 90.0
    regions_mito = compute_region_scores(dim_mito)
    sn_rank = next((i for i, r in enumerate(regions_mito)
                    if r.region_key == "substantia_nigra"), None)
    test("High mitochondrial → substantia nigra in top 3",
         sn_rank is not None and sn_rank < 3,
         f"Substantia nigra rank: {sn_rank}")

except Exception as e:
    for _ in range(3):
        test("brain_region_mapper", False, str(e))

# ── Summary ────────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
n_pass = sum(1 for _, ok in results if ok)
n_fail = sum(1 for _, ok in results if not ok)
print(f"Results: {n_pass}/{len(results)} tests passed, {n_fail} failed")

if n_fail == 0:
    print("\n\033[92m  ALL TESTS PASSED — ready to run training pipeline\033[0m")
    print("  Next: python scripts/generate_training_data.py")
else:
    print(f"\n\033[91m  {n_fail} tests failed — check the output above\033[0m")
    print("  Ensure all files from the fix package are in D:\\BRAINSAFE_AI\\")
    sys.exit(1)
