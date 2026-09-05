import tempfile
import unittest
from pathlib import Path

import pandas as pd

from atlas_reader import AtlasReader


class DashboardCandidatesTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.reader = AtlasReader(Path(self.directory.name))
        self.path = self.reader.root / "results/cmap/integrated_evidence/ATLAS_integrated_evidence_matrix.csv"
        self.path.parent.mkdir(parents=True)

    def test_preview_matches_final_order_and_evidence(self):
        # Earlier scores and ambiguous columns must not override final evidence.
        pd.DataFrame([
            {
                "priority_rank": 10 - i,
                "compound_classification": "identified_compound",
                "pubchem_status": "resolved",
                "active_target_gene_ids": "1234 | 5678",
                "integrated_target_network_score": i,
                "vina_mode_n": 9,
                "reference_best_affinity_kcal_mol": -12.0,
                "pert_iname": f"compound-{i}",
                "mean_tau": -60 - i,
                "validated_target_symbol": f"TARGET{i}",
                "best_affinity_kcal_mol": -9.327 if i == 0 else (-6.528 if i == 1 else None),
                "experimental_priority_score": 20 - i,
                "final_evidence_category": "EXPERIMENTAL_VALIDATION_WITH_CAUTION",
                "04u_rank": i + 1,
            }
            for i in range(10)
        ]).to_csv(self.path, index=False)

        candidates = self.reader.candidates()["rows"]
        dashboard = self.reader.dashboard()
        preview = dashboard["top_candidates"]

        self.assertEqual(len(preview), 8)
        self.assertEqual([r["name"] for r in preview], [r["pert_iname"] for r in candidates[:8]])
        for actual, source in zip(preview, candidates):
            self.assertEqual(actual["connectivity_score"], source["mean_tau"])
            self.assertEqual(actual["target"], source["validated_target_symbol"])
            self.assertEqual(actual["docking_score"], source["best_affinity_kcal_mol"])
            self.assertEqual(actual["final_score"], source["experimental_priority_score"])
            self.assertEqual(actual["status"], source["final_evidence_category"])
        self.assertEqual(next(m["value"] for m in dashboard["metrics"] if m["key"] == "docking"), 2)

    def test_missing_score_does_not_use_pose_count_or_reference_affinity(self):
        pd.DataFrame([{
            "pert_iname": "candidate",
            "active_target_gene_ids": "1234",
            "vina_mode_n": 9,
            "reference_best_affinity_kcal_mol": -12.0,
        }]).to_csv(self.path, index=False)

        preview = self.reader.dashboard()["top_candidates"]

        self.assertEqual(len(preview), 1)
        self.assertIsNone(preview[0]["target"])
        self.assertIsNone(preview[0]["docking_score"])

    def test_explicit_legacy_score_column_is_supported(self):
        pd.DataFrame([{"pert_iname": "candidate", "vina_mode_n": 9, "docking_score": 0.0}]).to_csv(self.path, index=False)
        self.assertEqual(self.reader.dashboard()["top_candidates"][0]["docking_score"], 0.0)

    def test_missing_candidates_produces_empty_preview(self):
        self.assertEqual(self.reader.candidates()["rows"], [])
        self.assertEqual(self.reader.dashboard()["top_candidates"], [])


if __name__ == "__main__":
    unittest.main()
