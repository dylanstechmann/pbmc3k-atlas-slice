from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import analyze  # noqa: E402


class ModuleTests(unittest.TestCase):
    def test_filter_module_splits_present_and_missing(self) -> None:
        available = {"CDKN1A", "GADD45A", "IL6"}
        present, missing = analyze.filter_module(
            ["CDKN1A", "GADD45A", "MDM2", "IL6"], available
        )
        self.assertEqual(present, ["CDKN1A", "GADD45A", "IL6"])
        self.assertEqual(missing, ["MDM2"])

    def test_pick_cell_type_prefers_highest_mean(self) -> None:
        means = {"B cell": 0.2, "T cell": 1.1, "NK cell": 0.4}
        self.assertEqual(analyze.pick_cell_type(means), "T cell")

    def test_pick_cell_type_empty(self) -> None:
        self.assertIsNone(analyze.pick_cell_type({}))

    def test_marker_sets_are_nonempty_pairs(self) -> None:
        for cell_type, markers in analyze.MARKER_SETS.items():
            self.assertTrue(1 <= len(markers) <= 4, cell_type)
            self.assertTrue(all(isinstance(m, str) and m for m in markers), cell_type)

    def test_modules_are_nonempty(self) -> None:
        for name, genes in analyze.MODULES.items():
            self.assertGreater(len(genes), 1, name)
            self.assertEqual(len(genes), len(set(genes)), name)

    def test_methods_report_selected_seed(self) -> None:
        receipt = {
            "started": "2026-09-26",
            "dataset_sha256": "abc",
            "versions": {"scanpy": "1", "anndata": "1", "numpy": "1"},
            "qc": {
                "min_genes": 200,
                "max_genes": 2500,
                "max_mito_fraction": 0.05,
                "min_cells_per_gene": 3,
            },
        }
        methods = analyze.render_methods(receipt, 17)
        self.assertIn("seed 17", methods)
        self.assertIn("--seed 17", methods)
        self.assertNotIn("seed 0", methods)

    def test_permutation_null_detects_a_planted_difference_and_not_noise(self) -> None:
        import numpy as np

        rng = np.random.default_rng(1)
        labels = np.repeat(["A", "B", "C"], 60)
        planted = rng.normal(0, 1, 180) + np.where(labels == "A", 3.0, 0.0)
        noise = rng.normal(0, 1, 180)
        hit = analyze.label_permutation_null(planted, labels, n_permutations=200, seed=3)
        null = analyze.label_permutation_null(noise, labels, n_permutations=200, seed=3)
        self.assertLess(hit["p_add_one"], 0.01)
        self.assertGreater(null["p_add_one"], 0.05)
        self.assertEqual(hit, analyze.label_permutation_null(planted, labels, n_permutations=200, seed=3))
        self.assertGreaterEqual(hit["p_add_one"], 1 / 201)

    def test_permutation_null_rejects_unusable_input(self) -> None:
        with self.assertRaises(ValueError):
            analyze.label_permutation_null([1.0, 2.0], ["A"])
        with self.assertRaises(ValueError):
            analyze.label_permutation_null([1.0, 2.0], ["A", "A"])

    def test_cell_type_counts(self) -> None:
        self.assertEqual(analyze.cell_type_counts(["B", "A", "B"]), {"A": 1, "B": 2})


if __name__ == "__main__":
    unittest.main()
