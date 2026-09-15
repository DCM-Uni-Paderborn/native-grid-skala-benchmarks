"""Numerical and provenance regressions for the reported basis controls."""
import copy
import json
import unittest
from unittest.mock import patch
import basis_sensitivity as basis


class BasisSensitivityTests(unittest.TestCase):
    def test_paired_crystals_and_same_ice_population(self):
        report = basis.assess()
        self.assertAlmostEqual(report["crystals"]["CO2"]["basis_shift_kjmol"], 4.3094260471, places=8)
        self.assertAlmostEqual(report["crystals"]["NH3"]["qzvpp"]["lattice_energy_kjmol"], -38.1845506711, places=8)
        self.assertEqual(len(report["ice"]), 6)
        self.assertAlmostEqual(report["ice_mae"]["qzvpp"], 1.0950635539, places=8)
        self.assertAlmostEqual(report["ice_mae"]["tzvpp"], 22.3649283298, places=8)

    def test_outlier_is_reported_not_hidden_or_promoted(self):
        report = basis.assess()
        self.assertEqual(len(report["molecular_reactions"]), 5)
        row = next(r for r in report["molecular_reactions"] if r["subset"] == "BHROT27")
        self.assertEqual(row["scientific_status"], "unresolved_QZ_outlier")
        self.assertAlmostEqual(row["qz_kcal_mol"], 28.0093226301, places=8)
        self.assertIn(r"28.009$^{*}$", basis.tables(report)["molecular-basis-table-si.tex"])

    def test_wrong_input_hash_rejected(self):
        index = json.loads((basis.DATA / "index.json").read_text())
        case = copy.deepcopy(index["cases"][0])
        case["files_sha256"]["actual-input.inp"] = "0" * 64
        with self.assertRaises(AssertionError):
            basis.verified(case)

    def test_wrong_runtime_rejected(self):
        original = basis.json.loads
        def altered(text):
            record = original(text)
            if "runtime" in record:
                record["runtime"]["revision"] = "different-runtime"
            return record
        case = original((basis.DATA / "index.json").read_text())["cases"][0]
        with patch.object(basis.json, "loads", altered), self.assertRaises(AssertionError):
            basis.verified(case)


if __name__ == "__main__":
    unittest.main()
