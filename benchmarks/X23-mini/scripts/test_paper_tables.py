"""Check manuscript tables against the accepted analysis, without SCF runs."""

import copy
import json
import unittest

from build_paper_tables import ROOT, render


class PaperTableTests(unittest.TestCase):
    def setUp(self):
        self.data = json.loads((ROOT / "results/lattice-energies.json").read_text())

    def test_archived_tables_regenerate_exactly(self):
        for name, text in render(self.data).items():
            self.assertEqual((ROOT / "paper" / name).read_text(), text)

    def test_incomplete_or_duplicate_pairs_are_rejected(self):
        for change in ("count", "missing", "duplicate", "replacement"):
            with self.subTest(change=change):
                data = copy.deepcopy(self.data)
                if change == "count":
                    data["accepted_base_cases"] = 23
                elif change == "missing":
                    data["complete_pairs"].pop()
                elif change == "duplicate":
                    data["complete_pairs"].append(data["complete_pairs"][0])
                else:
                    data["complete_pairs"][-1] = data["complete_pairs"][0]
                with self.assertRaises(ValueError):
                    render(data)

    def test_missing_or_duplicate_paired_grid_control_is_rejected(self):
        for controls in ([], self.data["paired_numerical_controls"] * 2):
            data = dict(self.data, paired_numerical_controls=controls)
            with self.assertRaises(ValueError):
                render(data)


if __name__ == "__main__":
    unittest.main()
