"""Focused checks of provenance, normalization and publication selection."""
import copy
import json
import unittest
from unittest.mock import patch
import analyze_basis_controls as basis

class BasisTests(unittest.TestCase):
    def test_paired_shifts(self):
        r=basis.assessment()
        self.assertAlmostEqual(r["paired_basis_change_kjmol"],14.4154079916,places=8)
        self.assertAlmostEqual(r["paired_grid_change_kjmol"],-.2053177775,places=8)
        self.assertAlmostEqual(r["pairs"]["qzvpp-grid200-974"]["lattice_energy_kjmol"],-110.4274668659,places=8)

    def test_only_urea_ae_is_replaced(self):
        base=json.loads((basis.ROOT/"results/lattice-energies.json").read_text())["complete_pairs"]
        before=copy.deepcopy(base)
        selected=basis.selected_pairs(base)
        self.assertEqual(base,before)
        self.assertEqual(len(selected),12)
        for old,new in zip(base,selected):
            if (old["method"],old["system"])==("gapw-ae","urea"):
                self.assertNotEqual(old["execution_sha256"],new["execution_sha256"])
                self.assertEqual(new["superseded_base_lattice_energy_kjmol"],old["lattice_energy_kjmol"])
            else:
                self.assertEqual(old,new)

    def test_mismatched_provenance_rejected(self):
        real_sha=basis.sha
        def bad_sha(path):
            return "bad" if path.name=="actual-input.inp" else real_sha(path)
        with patch.object(basis,"sha",bad_sha):
            with self.assertRaises(AssertionError):
                basis.assessment()

if __name__=="__main__":
    unittest.main()

