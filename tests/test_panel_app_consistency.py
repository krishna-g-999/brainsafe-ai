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


class TestManuscriptCountsMatchThePanel(unittest.TestCase):
    """The manuscript states the panel's shape in prose, and prose does not regenerate.

    A prose pass found the abstract claiming 63 molecular targets where the server serves 54, and
    three places saying the pathway graph holds 52 targets when it holds 51, the graph having lost
    GluA2 when that endpoint was withdrawn from scoring. Neither number is produced by a training
    run, so nothing rebuilt them and nothing complained.

    check_manuscript_numbers.py cannot catch this. It asks whether a computed value appears
    somewhere in the documents, and a bare integer like 51 matches inside 0.951 or 151, so a small
    count is effectively unverifiable by substring search. This asserts the counts directly.
    """

    @classmethod
    def setUpClass(cls):
        import app
        cls.app = app
        # Both manuscripts are checked. The condensed draft is what would go to NAR at 4 to 5
        # printed pages; the full draft is the extended version the technical report draws on. Two
        # documents stating the same counts is two places for one of them to fall behind.
        cls.text = "\n".join(
            (ROOT / "manuscript" / f).read_text(encoding="utf-8")
            for f in ("NAR_WebServer_BrainSafe_draft.md", "NAR_condensed_draft.md")
            if (ROOT / "manuscript" / f).exists())

    def test_pathway_graph_target_count(self):
        import re
        n = len(self.app.KNOWLEDGE_GRAPH)
        stated = set(int(m) for m in re.findall(r"the (\d+) targets in the pathway graph", self.text))
        self.assertTrue(stated, "the manuscript no longer states a pathway-graph target count")
        self.assertEqual(
            stated, {n},
            f"the manuscript states {sorted(stated)} pathway-graph targets; the graph holds {n}")

    def test_molecular_target_count(self):
        import re
        served = (set(self.app.TARGET_CLASSIFIERS) | set(self.app.BINDER_TARGETS)
                  | set(self.app.RECEPTOR_REGRESSORS)) - {"BBB"}
        stated = re.search(r"engagement of (\d+) molecular targets", self.text)
        self.assertIsNotNone(stated, "the abstract no longer states a molecular-target count")
        self.assertEqual(
            int(stated.group(1)), len(served),
            f"the abstract claims {stated.group(1)} molecular targets; the server serves "
            f"{len(served)} excluding the barrier model")


class TestPanelCountsReconcile(unittest.TestCase):
    """The panel was being counted six ways, and the counts must add up.

    The interface quoted 47, 52, 55, 63, 70 and 75 in different sections. Each was correct for a
    different question and none was reconciled, which reads as carelessness whatever the arithmetic.
    panel_shape() is now the single source, and these assert the two identities that make its
    figures a partition rather than a list.
    """

    @classmethod
    def setUpClass(cls):
        import app
        cls.sh = app.panel_shape()
        cls.app = app

    def test_quantities_partition_into_targets_exposure_and_other(self):
        s = self.sh
        self.assertEqual(
            s["targets"] + s["exposure"] + s["other"], s["quantities"],
            "the panel breakdown does not sum to the number of distinct predicted quantities")

    def test_estimators_exceed_quantities_only_by_dual_modelled_proteins(self):
        s = self.sh
        self.assertEqual(
            s["quantities"] + len(s["dual_model"]), s["deployed"],
            f"{s['deployed']} deployed estimators cover {s['quantities']} quantities, which is "
            f"reconcilable only if exactly {s['deployed'] - s['quantities']} proteins carry two "
            f"models; {len(s['dual_model'])} do: {s['dual_model']}")

    def test_trained_is_deployed_plus_withdrawn(self):
        s = self.sh
        self.assertEqual(s["deployed"] + s["withdrawn"], s["trained"])

    def test_barrier_model_is_not_counted_as_a_target(self):
        """BBB is exposure. Counting it among the targets is how 54 became 55."""
        served = (set(self.app.TARGET_CLASSIFIERS) | set(self.app.BINDER_TARGETS)
                  | set(self.app.RECEPTOR_REGRESSORS))
        self.assertIn("BBB", served, "the barrier model is no longer served")
        self.assertEqual(
            self.sh["targets"], len(served) - 1,
            "the molecular-target count must exclude the barrier model, which is an exposure term")


class TestOrphanTargetDegradesToSilence(unittest.TestCase):
    """The guard must return zero rather than raise, whatever the graph says."""

    def test_unknown_target_scores_zero(self):
        import app
        r = {"targets": {}, "receptor_binder": {}}
        self.assertEqual(app.target_signal(r, 0.0, "a_target_that_does_not_exist"), 0.0)


if __name__ == "__main__":
    unittest.main()
