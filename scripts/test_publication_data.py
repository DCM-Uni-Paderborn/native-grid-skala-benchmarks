"""Regression checks for the curated numerical populations and failure boundaries."""
import json
from pathlib import Path
import tempfile
import unittest

from data_checks import ROOT, close, compare_nested, metrics, sha, verify_files
from analyze_band_gaps import band_edges, gap_type, statistics
from additional_checks import ae_cutoff, adjoint, crystal_dft, paper_coverage

SAMPLE = '''# Point 1 Spin 1: 0 0 0
1 -10 2
2 -1 2
3 1 0
# Point 2 Spin 1: 0.5 0 0
1 -9 2
2 0 2
3 2 0
'''


class BandParsingTests(unittest.TestCase):
    def test_minimum_and_direct_are_distinct(self):
        result = band_edges(SAMPLE, 2, 2)
        self.assertEqual(result['sampled_gap_eV'], 1)
        self.assertEqual(result['sampled_direct_gap_eV'], 2)
        self.assertEqual(result['vbm_k'], [.5, 0, 0])

    def test_only_reviewed_core_indices_are_allowed(self):
        text = SAMPLE.replace('1 -10 2', '1 ************** 2')
        with self.assertRaises(ValueError):
            band_edges(text, 2, 2)
        self.assertFalse(band_edges(text, 2, 2, [1])['all_energies_readable'])

    def test_frontier_censoring_is_rejected(self):
        with self.assertRaises(ValueError):
            band_edges(SAMPLE.replace('2 -1 2', '2 ************** 2'), 2, 2, [2])

    def test_wrong_sampling_and_occupations_are_rejected(self):
        for text, points in [(SAMPLE, 3), (SAMPLE.replace('3 1 0', '3 1 1'), 2),
                             (SAMPLE.replace('3 1 0', '3 nan 0'), 2)]:
            with self.subTest(text=text, points=points), self.assertRaises(ValueError):
                band_edges(text, 2, points)

    def test_reported_classification_thresholds(self):
        self.assertEqual([gap_type(1, d) for d in (1, 1.00000642, 1.01)], ['D', 'N', 'I'])

    def test_missing_reference_is_never_zero_filled(self):
        with self.assertRaises(KeyError):
            statistics({'A': 1}, {}, ['A'])


class PublicationTests(unittest.TestCase):
    def test_current_float_coverage(self):
        self.assertEqual(paper_coverage(), {'tables': 33, 'figures': 4, 'citations': 60})

    def test_case_insensitive_paths_remain_unique(self):
        data = json.loads((ROOT / 'benchmarks/band-gaps/dataset.json').read_text())
        paths = [r['directory'].casefold() for r in data['executions'].values()]
        self.assertEqual(len(paths), len(set(paths)))
        for material in ('BAs', 'BaS'):
            self.assertIn(material, {r['material'] for r in data['primary']})

    def test_native_scope_is_exact(self):
        data = json.loads((ROOT / 'benchmarks/dietGMTKN55/dataset.json').read_text())
        self.assertEqual(len(data['reactions']), 62)
        self.assertEqual(len(data['executions']), 420)
        self.assertEqual(set(data['methods']), {'ae-tz', 'ae-qz', 'gth'})

    def test_energy_checks_are_not_relaxed_with_eos_tolerances(self):
        compare_nested({'B0_GPa': 197.4532732}, {'B0_GPa': 197.4532810})
        with self.assertRaises(ValueError):
            compare_nested({'energy_hartree': -76.00001}, {'energy_hartree': -76.0})
        with self.assertRaises(ValueError):
            compare_nested({'B0_GPa': 197.45}, {'B0_GPa': 197.46})

    def test_hash_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / 'output.out'
            path.write_text('original')
            record = {'directory': '.', 'files_sha256': {'output.out': sha(path)}}
            verify_files(record, root)
            path.write_text('changed')
            with self.assertRaises(ValueError):
                verify_files(record, root)

    def test_empty_or_nonfinite_statistics_are_rejected(self):
        for values in ([], [float('nan')], [float('inf')]):
            with self.assertRaises(ValueError):
                metrics(values)

    def test_ae_cutoff_reference(self):
        rows = ae_cutoff()
        close(rows[0]['difference_kJ_mol'], -.141062, 5e-7)
        self.assertEqual(rows[-1]['difference_kJ_mol'], 0)

    def test_adjoint_tolerance_is_not_a_measured_maximum(self):
        self.assertIsNone(adjoint()['dot_product_identity']['measured_maximum'])

    def test_crystal_literature_maes(self):
        rows = {r['method']: r for r in crystal_dft()}
        self.assertEqual(f"{rows['PBE+D3']['statistics']['MAE']:.2f}", '3.76')
        self.assertEqual(f"{rows['SCAN+rVV10']['statistics']['MAE']:.2f}", '3.54')


if __name__ == '__main__':
    unittest.main()
