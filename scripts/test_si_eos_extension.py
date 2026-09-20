"""Regression checks for the separate Si sampling control."""
import unittest
from unittest.mock import patch

from si_eos_extension import assess, table, verify_new, DATA
from molecular_completion_controls import parse_output
import json


class SiliconEOS(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = assess()

    def test_population_and_basis_shift(self):
        report = self.report
        self.assertEqual(len(report["points"]), 10)
        self.assertEqual(report["new_scf_points"], 5)
        self.assertEqual(report["reused_qz_points"], 5)
        self.assertFalse(report["included_in_uniform_statistics"])
        fits = report["fits"]["ten_all"]
        self.assertAlmostEqual(fits["qzvpp"]["a0_A"], 5.354665382752, places=7)
        self.assertAlmostEqual(fits["qzvpp"]["B0_GPa"], 106.891047773, places=4)
        for pair in report["fits"].values():
            self.assertLess(pair["qzvpp"]["a0_A"], pair["tzvpp"]["a0_A"])
            self.assertGreater(pair["qzvpp"]["B0_GPa"], pair["tzvpp"]["B0_GPa"])
            self.assertTrue(all(f["bracketed"] for f in pair.values()))

    def test_newest_block_and_diagonalizers(self):
        text = "PROGRAM STARTED\n SCF run converged\n PROGRAM ENDED\nPROGRAM STARTED\n"
        text += " 1 NoMix/Diag. 0.10 1.0 0.2 -1.0 -0.1\n 2 MBroy/Diag. 0.10 1.0 0.1 -1.1 -0.1\n[ABORT]\n"
        result = parse_output(text)
        self.assertEqual(len(result["steps"]), 2)
        self.assertFalse(result["converged"])
        self.assertFalse(result["ended"])
        self.assertTrue(result["abort"])

    def test_changed_archival_hash_rejected(self):
        case = json.loads((DATA / "index.json").read_text())["cases"][0]
        with patch("si_eos_extension.sha", return_value="mismatch"):
            with self.assertRaises(AssertionError):
                verify_new(case)

    def test_failed_scheduler_evidence_rejected(self):
        case = json.loads((DATA / "index.json").read_text())["cases"][0]
        case["slurm"]["job_and_steps"]["0"] = "FAILED|1:0"
        with self.assertRaises(AssertionError):
            verify_new(case)

    def test_table_population(self):
        content = table(self.report)
        self.assertEqual(content.count(" & TZVPP & "), 5)
        self.assertEqual(content.count(" & QZVPP & "), 5)
        self.assertIn("uniform ten-solid statistics are unchanged", content)


if __name__ == "__main__":
    unittest.main()
