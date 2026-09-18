"""Numerical and provenance regressions for the reported basis controls."""
import copy
import json
import re
import unittest
from unittest.mock import patch
import basis_sensitivity as basis


class BasisSensitivityTests(unittest.TestCase):
    def test_paired_crystals_and_same_ice_population(self):
        report = basis.assess()
        self.assertAlmostEqual(report["crystals"]["CO2"]["basis_shift_kjmol"], 4.3094260471, places=8)
        self.assertAlmostEqual(report["crystals"]["NH3"]["qzvpp"]["lattice_energy_kjmol"], -38.1845506711, places=8)
        self.assertEqual(len(report["ice"]), 13)
        self.assertIn("XIII", {r["phase"] for r in report["ice"]})
        self.assertAlmostEqual(report["ice_mae"]["qzvpp"], 1.1607398827, places=8)
        self.assertAlmostEqual(report["ice_mae"]["tzvpp"], 21.3491532798, places=8)
        self.assertAlmostEqual(report["ice_relative_mae"]["qzvpp"], .8214372870, places=8)
        self.assertAlmostEqual(report["ice_relative_mae"]["tzvpp"], 4.1725627527, places=8)

    def test_xiii_pair_and_absolute_relative_error_distinction(self):
        report = basis.assess()
        row = next(r for r in report["ice"] if r["phase"] == "XIII")
        relative = next(r for r in report["ice_relative"] if r["phase"] == "XIII")
        self.assertAlmostEqual(row["qzvpp"], -57.0356719856, places=8)
        self.assertAlmostEqual(row["tzvpp"], -79.9420202762, places=8)
        self.assertAlmostEqual(relative["qzvpp"] - relative["dmc_kjmol"], -1.4282092168, places=8)
        inputs = []
        for level in ("tzvpp", "qzvpp"):
            folder = basis.ROOT / f"benchmarks/DMC-ICE13/basis-controls/XIII/{level}"
            inp = basis.invariant((folder / "actual-input.inp").read_text())
            inputs.append(re.sub(r"(?m)^\s*(SYMMETRY_BACKEND|SYMMETRY_REDUCTION_METHOD|EPS_SYMMETRY)\s+[^\n]+\n", "", inp))
        self.assertEqual(*inputs)
        stats = report["ice_error_statistics"]
        self.assertEqual(stats["absolute"]["tzvpp"]["negative_error_count"], 13)
        self.assertEqual(stats["absolute"]["qzvpp"]["negative_error_count"], 1)
        self.assertEqual(stats["relative_ih"]["qzvpp"]["count"], 12)
        self.assertEqual(stats["relative_ih"]["qzvpp"]["maximum_error_phase"], "VII")

    def test_relative_energies_cancel_monomer_reference(self):
        report = basis.assess()
        index = json.loads((basis.DATA / "index.json").read_text())
        for level in ("tzvpp", "qzvpp"):
            solids = {c["system"]: basis.verified(c)[0] for c in index["cases"]
                      if (c["kind"], c["level"], c["phase"]) == ("ice", level, "solid")}
            ih = solids["Ih"]["energy_hartree"] / solids["Ih"]["molecules"]
            for row in report["ice_relative"]:
                run = solids[row["phase"]]
                direct = (run["energy_hartree"] / run["molecules"] - ih) * basis.HA_KJ
                self.assertAlmostEqual(row[level], direct, places=9)

    def test_outlier_is_reported_not_hidden_or_promoted(self):
        report = basis.assess()
        self.assertEqual(len(report["molecular_reactions"]), 6)
        row = next(r for r in report["molecular_reactions"] if r["subset"] == "BHROT27")
        self.assertEqual(row["scientific_status"], "initial_higher_energy_solution")
        self.assertAlmostEqual(row["qz_kcal_mol"], 28.0093226301, places=8)
        self.assertIn(r"28.009$^{*}$", basis.tables(report)["molecular-basis-table-si.tex"])
        self.assertAlmostEqual(report["ethane_state_controls"]["restart_barrier_kcal_mol"], 2.7973013720, places=8)
        self.assertAlmostEqual(report["ethane_state_controls"]["fine_grid_atomic_barrier_kcal_mol"], 27.5528354505, places=8)
        self.assertFalse(report["ethane_state_controls"]["included_in_uniform_production_statistics"])
        water = next(r for r in report["molecular_reactions"] if r["subset"] == "WATER27")
        self.assertAlmostEqual(water["qz_kcal_mol"], 29.1825487215, places=8)

    def test_three_solid_controls_preserve_volume_specific_settings(self):
        rows = basis.assess()["solid_basis_controls"]
        self.assertEqual([r["solid"] for r in rows], ["Si", "C", "MgO"])
        self.assertTrue(all(len(r["points"]) == 5 for r in rows))
        self.assertAlmostEqual(rows[2]["fits"]["qzvpp"]["a0_angstrom"], 4.0904453933, places=8)
        self.assertAlmostEqual(rows[2]["fits"]["qzvpp"]["B0_GPa"], 197.4532810081, places=6)
        self.assertNotEqual(basis.solid_invariant("SCHEME MONKHORST-PACK 5 5 5"),
                            basis.solid_invariant("SCHEME MONKHORST-PACK 6 6 6"))

    def test_generated_tables_match_manuscript_snapshot(self):
        for name, text in basis.tables(basis.assess()).items():
            self.assertEqual(text.strip(), (basis.ROOT / "paper/pccp" / name).read_text().strip())

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
