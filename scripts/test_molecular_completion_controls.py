"""Checks for the additional reaction controls and non-promoted precursors."""
import copy
import json
import unittest

import molecular_completion_controls as controls


class CompletionControlTests(unittest.TestCase):
    def test_reactions_and_scope(self):
        result = controls.assess()
        self.assertFalse(result["included_in_uniform_statistics"])
        self.assertFalse(result["precursor_energies_used"])
        self.assertEqual(len(result["final_species"]), 5)
        carb, ea = result["reactions"]
        self.assertAlmostEqual(carb["qz_kcal_mol"], 3.78373337248, places=8)
        self.assertGreater(abs(carb["qz_error_kcal_mol"]), abs(carb["tz_error_kcal_mol"]))
        self.assertAlmostEqual(ea["qz_kcal_mol"], 54.97633941006, places=8)
        self.assertLess(abs(ea["qz_error_kcal_mol"]), abs(ea["tz_error_kcal_mol"]))
        self.assertEqual(result["final_species"]["EA_25"]["electrons"], [18, 17])

    def test_only_newest_block_is_used(self):
        text = "PROGRAM STARTED\nSCF run converged\nPROGRAM ENDED\nPROGRAM STARTED\n[ABORT]"
        parsed = controls.parse_output(text)
        self.assertFalse(parsed["converged"])
        self.assertFalse(parsed["ended"])
        self.assertTrue(parsed["abort"])

    def test_physical_normalization_preserves_non_scf_settings(self):
        text = "&SCF\n &OT\n &END OT\n&END SCF\nCUTOFF 640\nCHARGE 0\n"
        self.assertNotEqual(controls.without_scf(text), controls.without_scf(text.replace("640", "600")))
        self.assertNotEqual(controls.without_scf(text), controls.without_scf(text.replace("CHARGE 0", "CHARGE 1")))

    def test_hash_mismatch_rejected(self):
        case = copy.deepcopy(json.loads((controls.DATA / "index.json").read_text())["cases"][0])
        case["hashes"]["actual-input.inp"] = "0" * 64
        with self.assertRaises(AssertionError):
            controls.verified(case)

    def test_generated_table_matches_snapshot(self):
        self.assertEqual(controls.table(controls.assess()),
                         (controls.ROOT / "paper/pccp/molecular-additional-basis-table-si.tex").read_text())


if __name__ == "__main__":
    unittest.main()
