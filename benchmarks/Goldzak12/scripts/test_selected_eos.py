import copy
import unittest

import numpy as np

from fit_eos import birch_murnaghan
from fit_selected_eos import analyze, select_curves
from analyze_selected_eos import statistics


class SelectedEosTests(unittest.TestCase):
    def test_signed_statistics(self):
        result = statistics([
            {'method': 'test', 'delta_a_A': -2.0, 'delta_B0_GPa': 3.0},
            {'method': 'test', 'delta_a_A': 4.0, 'delta_B0_GPa': -3.0},
        ], 'method')[0]
        self.assertEqual(result['N'], 2)
        self.assertEqual(result['ME_a_A'], 1.0)
        self.assertEqual(result['MAE_a_A'], 3.0)
        self.assertEqual(result['ME_B0_GPa'], 0.0)
        self.assertEqual(result['RMSE_B0_GPa'], 3.0)

    def setUp(self):
        self.selection = {
            "selection_date": "2026-09-05", "solids": ["AlN", "BN"],
            "methods": ["gapwxc-gth", "gapw-ae", "hybrid-direct", "hybrid-one-center"],
            "default_points": [f"v{i:02d}" for i in range(1, 11)],
            "point_overrides": {"AlN": [f"v{i:02d}" for i in range(1, 10)]},
            "shared_ae_solids": ["BN"], "remaining_checks": ["grid convergence"],
        }
        self.rows = []
        volumes = np.linspace(90.0, 110.0, 10)
        energies = birch_murnaghan(volumes, -10.0, 100.0, 0.01, 4.0)
        for solid in self.selection["solids"]:
            methods = self.selection["methods"] if solid == "AlN" else self.selection["methods"][:2]
            for method in methods:
                for i, (v, e) in enumerate(zip(volumes, energies), 1):
                    self.rows.append({"method": method, "solid": solid, "point": f"v{i:02d}",
                                      "volume_A3": str(v), "energy_Ha": str(e), "execution_accepted": "true"})

    def test_common_window_and_shared_ae(self):
        curves = select_curves(self.rows, self.selection)
        self.assertEqual(len(curves), 8)
        for (method, solid), (source, rows) in curves.items():
            self.assertEqual(len(rows), 9 if solid == "AlN" else 10)
            if solid == "BN" and method.startswith("hybrid-"):
                self.assertEqual(source, "gapw-ae")

    def test_reject_duplicate_and_missing(self):
        with self.assertRaisesRegex(ValueError, "Duplicate source"):
            select_curves(self.rows + [self.rows[0]], self.selection)
        with self.assertRaisesRegex(ValueError, "Missing selected"):
            select_curves(self.rows[1:], self.selection)

    def test_reject_unaccepted_nonfinite_and_volume_mismatch(self):
        for field, value in [("execution_accepted", "false"), ("energy_Ha", "nan"), ("volume_A3", "91")]:
            rows = copy.deepcopy(self.rows)
            rows[0][field] = value
            with self.assertRaises(ValueError):
                select_curves(rows, self.selection)

    def test_known_eos_and_no_automatic_quality_release(self):
        fits, summary = analyze(self.rows, self.selection)
        self.assertEqual(summary["unique_selected_points"], 56)
        self.assertEqual(summary["independent_curves"], 6)
        self.assertTrue(summary["all_minima_bracketed"])
        self.assertFalse(summary["quality_release"])
        self.assertFalse(summary["cohesive_energies_computed"])
        for fit in fits:
            self.assertAlmostEqual(fit["V0_A3_cell"], 100.0, places=5)
            self.assertEqual(fit["quality_status"], "numerical_review_pending")


if __name__ == "__main__":
    unittest.main()
