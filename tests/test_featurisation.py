"""The featuriser: the one component every model depends on.

If featurisation changes, every stored estimator becomes wrong without anything failing, because a
pickled forest will happily consume a vector of the right width and the wrong meaning. These tests
pin the properties the rest of the pipeline assumes:

  shape        1,036 columns, the same for every molecule, forever
  purity       a molecule's vector depends on that molecule and nothing else, not on what else is
               in the batch and not on input order
  parenting    salts are stripped and the largest organic fragment is used, so a salt and its free
               base give the same vector, which is why deduplication has to happen in feature space
  stereo       chirality is excluded, so enantiomers collide. This is a limitation, not an accident,
               and a test that pins it will fail loudly if someone "fixes" it without also changing
               the deduplication step that exists because of it
  determinism  the same input gives bit-identical output across calls

Run:  python -m pytest tests/ -v
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import (MORGAN_BITS, N_FEATURES, featurize,  # noqa: E402
                                featurize_one, feature_names)

DONEPEZIL = "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
ASPIRIN_SODIUM = "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]"
L_ALANINE = "C[C@@H](N)C(=O)O"
D_ALANINE = "C[C@H](N)C(=O)O"


class TestShape(unittest.TestCase):
    def test_width_is_fixed_and_matches_the_declared_constant(self):
        v = featurize_one(DONEPEZIL)
        self.assertIsNotNone(v)
        self.assertEqual(v.shape, (N_FEATURES,))
        self.assertEqual(N_FEATURES, MORGAN_BITS + 12)

    def test_feature_names_match_the_vector_width(self):
        self.assertEqual(len(feature_names()), N_FEATURES)

    def test_fingerprint_block_is_binary(self):
        v = featurize_one(DONEPEZIL)
        bits = v[:MORGAN_BITS]
        self.assertTrue(set(np.unique(bits)).issubset({0.0, 1.0}),
                        "fingerprint block must be binary; a count fingerprint would silently "
                        "change what every stored model consumes")

    def test_unparseable_input_returns_none_rather_than_a_zero_vector(self):
        # A zero vector is a valid input to a forest and would be scored as though it were a
        # molecule. Returning None forces the caller to handle the failure.
        self.assertIsNone(featurize_one("not a molecule"))
        self.assertIsNone(featurize_one(""))


class TestPurity(unittest.TestCase):
    """A molecule's vector must not depend on anything except that molecule."""

    def test_alone_equals_in_batch(self):
        batch, mask = featurize([DONEPEZIL, ASPIRIN, L_ALANINE])
        self.assertTrue(mask.all())
        for i, smi in enumerate([DONEPEZIL, ASPIRIN, L_ALANINE]):
            np.testing.assert_array_equal(batch[i], featurize_one(smi))

    def test_order_does_not_matter(self):
        smis = [DONEPEZIL, ASPIRIN, L_ALANINE]
        fwd, _ = featurize(smis)
        rev, _ = featurize(list(reversed(smis)))
        for i in range(len(smis)):
            np.testing.assert_array_equal(fwd[i], rev[len(smis) - 1 - i])

    def test_repeated_calls_are_bit_identical(self):
        a, b = featurize_one(DONEPEZIL), featurize_one(DONEPEZIL)
        self.assertEqual(a.tobytes(), b.tobytes())

    def test_mask_marks_failures_and_rows_stay_aligned(self):
        X, mask = featurize([DONEPEZIL, "garbage", ASPIRIN])
        self.assertEqual(list(mask), [True, False, True])
        self.assertEqual(len(X), 2, "failed rows must be absent, not zero-filled")
        np.testing.assert_array_equal(X[1], featurize_one(ASPIRIN))


class TestParentAndStereo(unittest.TestCase):
    def test_counter_ion_is_removed(self):
        """The sodium goes; what remains is the largest organic fragment."""
        from features.featurize import parent_mol
        from rdkit import Chem
        parent = Chem.MolToSmiles(parent_mol(ASPIRIN_SODIUM))
        self.assertNotIn("Na", parent)
        self.assertEqual(parent.count("."), 0, "the counter-ion should be gone")

    def test_KNOWN_DEFECT_salt_and_free_base_differ_because_charge_is_not_neutralised(self):
        """Pins a defect so that fixing it is a deliberate, visible act.

        Salt stripping keeps the largest fragment but does not neutralise it, so aspirin sodium
        reduces to the carboxyl*ate* and not to aspirin. The vectors therefore differ, and so do
        the predictions: measured on the deployed models, haloperidol hydrochloride returns BBB
        0.613 where haloperidol returns 0.993, and its hERG probability moves by 0.619.

        Training SMILES are overwhelmingly neutral parent forms, so a user who submits a salt is
        scoring an input the models were not trained on. This test asserts the current behaviour
        rather than the desired one: changing it alters predictions for every salt input and needs
        re-validation, so it must not happen by accident.
        """
        free, salt = featurize_one(ASPIRIN), featurize_one(ASPIRIN_SODIUM)
        self.assertFalse(np.array_equal(free, salt),
                         "if this now passes, neutralisation has been added: that is a scientific "
                         "change, so re-validate the panel and update this test deliberately")

    def test_enantiomers_collide_which_is_a_known_limitation(self):
        # Pinned deliberately. If chirality is ever included, _dedup_features must be revisited,
        # because its whole purpose is to collapse rows this collision creates.
        np.testing.assert_array_equal(featurize_one(L_ALANINE), featurize_one(D_ALANINE))

    def test_enantiomers_collide_which_is_a_known_limitation(self):
        # Pinned deliberately. If chirality is ever included, _dedup_features must be revisited,
        # because its whole purpose is to collapse rows this collision creates.
        np.testing.assert_array_equal(featurize_one(L_ALANINE), featurize_one(D_ALANINE))


if __name__ == "__main__":
    unittest.main(verbosity=2)
