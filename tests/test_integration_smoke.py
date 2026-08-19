"""End-to-end: load the deployed models and score reference compounds through the real entry points.

The unit tests pin components. This pins the assembled system, through `app.predict_all` and
`app.disease_scores`, the same functions the web interface calls, so a change that leaves every
component correct and the composition wrong is still caught.

What is asserted is deliberately coarse-grained where the science is uncertain and exact where it is
not. Donepezil must return acetylcholinesterase as its driver, because that is pharmacology rather
than a tuning artefact and a panel that loses it is broken. The probability is asserted only to lie
in a wide band, because pinning it to three decimals would turn every legitimate retrain into a test
failure and train people to update the expected values without reading them.

Skipped rather than failed when models_rf/ is absent, so the suite still runs on a fresh clone or in
CI where the 0.84 GB archive has not been fetched. A skip is visible; a silent pass is not.

Run:  python -m pytest tests/test_integration_smoke.py -v
"""
from __future__ import annotations

import logging
import sys
import unittest
from pathlib import Path

logging.getLogger("streamlit").setLevel(logging.ERROR)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

MODELS = ROOT / "models_rf"
HAVE_MODELS = MODELS.exists() and any(MODELS.glob("*.joblib"))

# name -> (SMILES, expected driving target, expected condition)
REFERENCE = {
    "donepezil": ("COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2", "AChE", "Alzheimer's disease"),
    "haloperidol": ("O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1", "D2",
                    "Psychosis / schizophrenia"),
    "morphine": ("CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "OPRM1",
                 "Chronic pain"),
    "fluoxetine": ("CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1", "SERT", "Depression / anxiety"),
}
PERIPHERAL = {
    "atorvastatin": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
    "hydrochlorothiazide": "NS(=O)(=O)c1cc2c(cc1Cl)NCNS2(=O)=O",
}
REPORT_THRESHOLD = 0.30


@unittest.skipUnless(HAVE_MODELS, "models_rf/ absent; run python model_fetch.py first")
class TestDeployedPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import app
        cls.app = app
        cls.models = app.load_models()

    def test_every_reference_compound_featurises_and_scores(self):
        for name, (smi, _drv, _dz) in REFERENCE.items():
            with self.subTest(compound=name):
                r = self.app.predict_all(smi, self.models)
                self.assertIsNotNone(r, f"{name} failed to featurise")
                self.assertIn("targets", r)
                self.assertIn("BBB", r["targets"])

    def test_probabilities_are_probabilities(self):
        r = self.app.predict_all(REFERENCE["donepezil"][0], self.models)
        for ep, p in r["targets"].items():
            with self.subTest(endpoint=ep):
                self.assertGreaterEqual(p, 0.0)
                self.assertLessEqual(p, 1.0)

    def test_reference_drugs_return_the_pharmacologically_correct_driver(self):
        for name, (smi, driver, disease) in REFERENCE.items():
            with self.subTest(compound=name):
                r = self.app.predict_all(smi, self.models)
                bbb, _neuro, dz = self.app.disease_scores(r)
                self.assertTrue(dz, f"{name} produced no disease rows")
                top = dz[0]
                self.assertEqual(top["disease"], disease,
                                 f"{name} should surface {disease}, got {top['disease']}")
                self.assertIsNotNone(top.get("driver"), f"{name} has no named driver")
                self.assertEqual(top["driver"][0], driver,
                                 f"{name} should be driven by {driver}, "
                                 f"got {top['driver'][0]}")

    def test_central_drugs_are_predicted_to_reach_the_brain(self):
        for name, (smi, _d, _z) in REFERENCE.items():
            with self.subTest(compound=name):
                r = self.app.predict_all(smi, self.models)
                self.assertGreater(r["targets"]["BBB"], 0.5,
                                   f"{name} is a central drug and should cross")

    def test_peripheral_compounds_stay_silent(self):
        """The complementary behaviour: no disease call for a compound that does not arrive."""
        for name, smi in PERIPHERAL.items():
            with self.subTest(compound=name):
                r = self.app.predict_all(smi, self.models)
                _bbb, _neuro, dz = self.app.disease_scores(r)
                top = dz[0]["gated"] if dz else 0.0
                self.assertLess(top, REPORT_THRESHOLD,
                                f"{name} acts peripherally and should produce no call")

    def test_scoring_is_deterministic_to_floating_point_tolerance(self):
        """Repeated scoring agrees, but not bit-for-bit, and the distinction is worth stating.

        The forests are fitted with a fixed seed and are not re-fitted here, so the only source of
        variation is the order in which tree votes are summed: with n_jobs=-1 that order is not
        fixed, and floating-point addition is not associative. Measured difference is of order
        1e-16, the last bit or two of a float64.

        Asserting bit-equality would be asserting something the deployment does not provide, and
        would fail intermittently on a machine with a different core count. Asserting agreement to
        1e-12 says what is true: two identical queries give the same answer to far beyond the
        precision anyone reports or acts on.
        """
        smi = REFERENCE["haloperidol"][0]
        a = self.app.predict_all(smi, self.models)
        b = self.app.predict_all(smi, self.models)
        self.assertEqual(set(a["targets"]), set(b["targets"]))
        for ep in a["targets"]:
            with self.subTest(endpoint=ep):
                self.assertAlmostEqual(a["targets"][ep], b["targets"][ep], delta=1e-12)

    def test_withdrawn_endpoints_are_not_offered(self):
        """Five endpoints are withdrawn; a regression that re-deploys any of them must fail here.

        Nav1.1 and Cav3.2 were withdrawn by the deployed-specificity audit, for firing on trivial
        metabolites at every threshold. NRF2, NFKB1 and NR3C1 were added later to test whether the
        panel's natural-product gap could be closed by adding the targets natural products are
        actually assayed against; they carry some cross-validated signal (scaffold AUROC 0.719,
        0.711, 0.596) but cannot be given a threshold that recovers actives without firing on
        unrelated chemistry, at sensitivities of 0.250, 0.048 and 0.000. Re-deploying any of them
        would put an endpoint on the server that cannot answer the question it is asked.

        This set is deliberately pinned rather than derived: the point of the test is that a change
        to it must be argued for, not absorbed.
        """
        import json
        modes = json.loads((MODELS / "binder_modes.json").read_text(encoding="utf-8"))
        withdrawn = {k for k, v in modes.items() if not v.get("deployed", True)}
        self.assertEqual(withdrawn, {"Nav1_1", "Cav3_2", "NRF2", "NFKB1", "NR3C1"},
                         "the withdrawal set changed; if deliberate, update this test and say why")

    def test_every_withdrawal_records_why_and_on_what_evidence(self):
        """A withdrawal without a recorded reason is an undocumented scientific decision.

        The three natural-product withdrawals had no reason at all until this audit, and the reason
        was then hand-written into a generated file, where the next pipeline run would have erased
        it. Both halves are checked: the reason exists, and it points at evidence that is not the
        deployed-specificity audit for the three endpoints that audit never saw.
        """
        import json
        modes = json.loads((MODELS / "binder_modes.json").read_text(encoding="utf-8"))
        for ep, rec in modes.items():
            if rec.get("deployed", True):
                continue
            with self.subTest(endpoint=ep):
                self.assertTrue((rec.get("withdrawn_reason") or "").strip(),
                                f"{ep} is withdrawn with no recorded reason")
                self.assertTrue((rec.get("withdrawn_evidence") or "").strip(),
                                f"{ep} is withdrawn with no evidence pointer")
        for ep in ("NRF2", "NFKB1", "NR3C1"):
            with self.subTest(endpoint=ep, check="evidence is not the audit"):
                self.assertNotIn("deployed_specificity_audit", modes[ep]["withdrawn_evidence"],
                                 f"{ep} was never deployed, so that audit cannot be its evidence")

    def test_unparseable_input_is_rejected_rather_than_scored(self):
        for bad in ("", "   ", "not a molecule"):
            with self.subTest(value=repr(bad)):
                self.assertIsNone(self.app.predict_all(bad, self.models))


if __name__ == "__main__":
    unittest.main(verbosity=2)
