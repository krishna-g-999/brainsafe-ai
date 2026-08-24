"""The served panel must be the registered panel, and the graph must only name models that exist.

The panel registry reconciles three views: the endpoint tables, the fitted models, and the registry
itself. It does not see app.py, and that gap let a defect sit undetected. `BINDER_TARGETS` was a
hardcoded list and had drifted from the registry in both directions at once:

  GluA2  had been withdrawn from the panel and was still in the list, so it was still scored, still
         carried an Epilepsy edge in the knowledge graph at weight 0.6, and could therefore drive a
         reported condition from a model the project had judged unfit to deploy.
  Cav3_2 is deployed and validated and was absent from the list, so it was never scored at all.

Neither raised an error. Nothing in the freshness graph or the panel check could see either, because
both compare the registry against the models rather than against the code that serves them. These
tests close that gap from the other side.

Run:  python -m pytest tests/test_panel_app_consistency.py -v
"""
from __future__ import annotations

import sys
import unittest
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))

import panel  # noqa: E402


class TestServedPanelMatchesRegistry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import app
        cls.app = app
        cls.deployed = set(panel.names(deployed=True))
        cls.withdrawn = {e.name for e in panel.withdrawn()}

    def test_binder_targets_are_exactly_the_registry_less_the_regressors(self):
        expected = self.deployed - set(self.app.RECEPTOR_REGRESSORS)
        self.assertEqual(
            set(self.app.BINDER_TARGETS), expected,
            "the served binder list has drifted from the registry; extra entries are scored "
            "without being deployed, missing ones are deployed without being scored")

    def test_no_withdrawn_endpoint_is_scored(self):
        served = set(self.app.BINDER_TARGETS) | set(self.app.TARGET_CLASSIFIERS)
        leaked = sorted(self.withdrawn & served)
        self.assertEqual(leaked, [], f"withdrawn endpoints are still being scored: {leaked}")

    def test_no_withdrawn_endpoint_holds_a_disease_edge(self):
        leaked = sorted(self.withdrawn & set(self.app.KNOWLEDGE_GRAPH))
        self.assertEqual(
            leaked, [],
            f"withdrawn endpoints still claim a mechanism in the knowledge graph: {leaked}. A "
            f"model judged unfit to deploy must not be able to drive a reported condition")

    def test_every_graph_target_has_a_scoring_path(self):
        scorable = (set(self.app.BINDER_TARGETS) | set(self.app.TARGET_CLASSIFIERS)
                    | set(self.app.RECEPTOR_REGRESSORS) | {"NEURO"})
        orphans = sorted(set(self.app.KNOWLEDGE_GRAPH) - scorable)
        self.assertEqual(
            orphans, [],
            f"graph targets with no model behind them: {orphans}. Before the guard in "
            f"target_signal this raised a KeyError in the middle of a live query")

    def test_every_served_binder_has_a_binder_scoring_kind(self):
        wrong = sorted(t for t in self.app.BINDER_TARGETS
                       if self.app.TARGET_KIND.get(t) != "binder")
        self.assertEqual(
            wrong, [],
            f"served binders scored by the wrong branch: {wrong}. The enrichment branch reads a "
            f"base rate these endpoints do not have, and would fall through to a KeyError")


class TestOrphanTargetDegradesToSilence(unittest.TestCase):
    """The guard must return zero rather than raise, whatever the graph says."""

    def test_unknown_target_scores_zero(self):
        import app
        r = {"targets": {}, "receptor_binder": {}}
        self.assertEqual(app.target_signal(r, 0.0, "a_target_that_does_not_exist"), 0.0)


if __name__ == "__main__":
    unittest.main()
