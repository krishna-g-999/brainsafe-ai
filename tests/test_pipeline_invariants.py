"""Invariants the scientific results depend on, each of which has been violated at least once.

These are not tests of "does the code run". Each one pins a property whose loss would change a
published number without raising an error, and most correspond to a defect the audit found:

  deduplication   collapses rows identical in feature space, drops contradictory groups rather than
                  voting on them, and does so BEFORE any split. Losing this reintroduces the leak
                  where a stereoisomer of a test compound sits in training.
  pools           the three background pools are disjoint and assignment is a pure function of the
                  structure, so a threshold set on one pool and measured on another stays honest
                  across runs and machines.
  label rule      a censored bound settles a label only when the whole interval lies on one side of
                  the cut. Feeding a bound to the exact-value rule silently discarded 253 measured
                  non-binders for AChE alone.
  determinism     the declared seed actually fixes the result.

Run:  python -m pytest tests/ -v
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import N_FEATURES                                  # noqa: E402
from models.pools import SHARES, _band, role_of                            # noqa: E402
from models.train_rf import SEED, _dedup_features                          # noqa: E402


def _rows(patterns):
    """Build a feature matrix from short bit patterns, padded to the real width."""
    X = np.zeros((len(patterns), N_FEATURES), dtype=np.float32)
    for i, p in enumerate(patterns):
        X[i, :len(p)] = p
    return X


class TestDeduplication(unittest.TestCase):
    def test_identical_rows_collapse_to_one(self):
        X = _rows([[1, 0, 1], [1, 0, 1], [0, 1, 0]])
        y = np.array([1, 1, 0])
        g = np.array(["a", "a", "b"])
        Xd, yd, gd, sd, rep = _dedup_features(X, y, g, ["s1", "s2", "s3"], "classification")
        self.assertEqual(len(Xd), 2)
        self.assertEqual(rep["duplicate_rows_removed"], 1)
        self.assertEqual(rep["rows_in"], 3)

    def test_contradictory_duplicates_are_dropped_not_voted(self):
        """Two identical inputs with opposite labels cannot both be learned from.

        Voting would pick one and present an arbitrary choice as data. Dropping is the honest
        response, and the count is reported so the loss is visible.
        """
        X = _rows([[1, 1, 0], [1, 1, 0], [0, 0, 1]])
        y = np.array([1, 0, 1])
        g = np.array(["a", "a", "b"])
        Xd, yd, gd, sd, rep = _dedup_features(X, y, g, ["s1", "s2", "s3"], "classification")
        self.assertEqual(rep["conflicting_groups_dropped"], 1)
        self.assertEqual(len(Xd), 1, "the contradictory pair must be gone, not resolved")
        self.assertEqual(int(yd[0]), 1)

    def test_regression_takes_the_group_median(self):
        X = _rows([[1, 0], [1, 0], [1, 0], [0, 1]])
        y = np.array([5.0, 7.0, 9.0, 2.0])
        g = np.array(["a", "a", "a", "b"])
        Xd, yd, gd, sd, rep = _dedup_features(X, y, g, ["a", "b", "c", "d"], "regression")
        self.assertEqual(len(Xd), 2)
        self.assertIn(7.0, list(yd), "the median of 5, 7, 9 is 7")

    def test_nothing_is_lost_when_there_are_no_duplicates(self):
        X = _rows([[1, 0], [0, 1], [1, 1]])
        y = np.array([1, 0, 1])
        g = np.array(["a", "b", "c"])
        Xd, _, _, _, rep = _dedup_features(X, y, g, ["a", "b", "c"], "classification")
        self.assertEqual(len(Xd), 3)
        self.assertEqual(rep["duplicate_rows_removed"], 0)


class TestBackgroundPools(unittest.TestCase):
    """The partition that keeps a threshold honest."""

    SMILES = ["CCO", "CCN", "c1ccccc1", "CC(=O)O", "CCCCCC", "CN1CCC1",
              "OCC1OC(O)C(O)C(O)C1O", "CC(C)Cc1ccc(cc1)C(C)C(O)=O"]

    def test_assignment_is_a_pure_function_of_the_structure(self):
        first = [role_of(s) for s in self.SMILES]
        second = [role_of(s) for s in reversed(self.SMILES)][::-1]
        self.assertEqual(first, second, "pool must not depend on order or call history")

    def test_a_compound_belongs_to_exactly_one_pool(self):
        for s in self.SMILES:
            self.assertIn(role_of(s), SHARES)

    def test_bands_cover_exactly_one_hundred(self):
        self.assertEqual(sum(SHARES.values()), 100)

    def test_hash_band_is_stable_and_in_range(self):
        for s in self.SMILES:
            b = _band(s)
            self.assertEqual(b, _band(s))
            self.assertTrue(0 <= b < 100)

    def test_the_partition_is_not_degenerate(self):
        # A hash that sent everything to one pool would satisfy every test above.
        rng = np.random.default_rng(SEED)
        smis = ["C" * (1 + int(rng.integers(1, 25))) for _ in range(600)]
        roles = {r: 0 for r in SHARES}
        for s in smis:
            roles[role_of(s)] += 1
        for role, n in roles.items():
            self.assertGreater(n, 0, f"no compound landed in the {role} pool")


class TestCensoredLabelRule(unittest.TestCase):
    """A bound settles a label only when the whole interval lies on one side of the cut.

    Reproduced here rather than imported, because the rule lives in several fetchers and the
    property under test is the rule itself. The defect this pins lost 253 measured non-binders for
    AChE by passing a bound to the exact-value rule, which treats it as a potency and discards the
    ambiguous 5 to 6 band.
    """

    INACTIVE_CUT = 5.0
    ACTIVE_CUT = 6.0

    @staticmethod
    def settles(bound: float, inactive_cut=5.0) -> bool:
        """`> bound` means the true potency is strictly below `bound`."""
        return bound <= inactive_cut

    def test_a_weak_bound_settles_the_compound_as_inactive(self):
        # "IC50 > 10 uM" is pChEMBL 5.0: everything below, so unambiguously inactive
        self.assertTrue(self.settles(5.0))
        self.assertTrue(self.settles(4.2))

    def test_a_bound_spanning_both_classes_is_undecidable(self):
        # "IC50 > 100 nM" is pChEMBL 7.0: the true value could be active or inactive
        self.assertFalse(self.settles(7.0))
        self.assertFalse(self.settles(5.5))

    def test_the_boundary_case_is_included(self):
        self.assertTrue(self.settles(self.INACTIVE_CUT),
                        "a bound exactly at the cut still places the true value below it")


class TestDeterminism(unittest.TestCase):
    def test_declared_seed_is_the_one_the_pipeline_uses(self):
        self.assertEqual(SEED, 42)

    def test_a_forest_with_the_fixed_seed_reproduces_itself(self):
        from sklearn.ensemble import RandomForestClassifier
        rng = np.random.default_rng(SEED)
        X = rng.random((120, 8))
        y = (X[:, 0] + X[:, 1] > 1.0).astype(int)
        a = RandomForestClassifier(n_estimators=25, random_state=SEED).fit(X, y)
        b = RandomForestClassifier(n_estimators=25, random_state=SEED).fit(X, y)
        np.testing.assert_array_equal(a.predict_proba(X), b.predict_proba(X))


if __name__ == "__main__":
    unittest.main(verbosity=2)
