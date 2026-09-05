import math
import unittest

from compare_selected_literature import error_statistics, matched_statistics
from fit_eos import ROOT, read_rows


class LiteratureComparisonTests(unittest.TestCase):
    def test_signed_statistics(self):
        r = error_statistics([1.0, 4.0], [2.0, 2.0])
        self.assertEqual(r['N'], 2)
        self.assertEqual(r['ME'], .5)
        self.assertEqual(r['MAE'], 1.5)
        self.assertAlmostEqual(r['RMSE'], math.sqrt(2.5))
        self.assertEqual(r['MARE_percent'], 75)

    def test_missing_property_is_not_zero(self):
        records = [dict(source=s, method=m, solid='C', a_A=3.5, B0_GPa=b)
                   for s, m, b in [('NativeSkala', 'native', 400), ('Schimka2011', 'PBE', '')]]
        rows = matched_statistics(records, {'C': dict(a_A=3.6, B0_GPa=450)}, ['native'], ['C'])
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(r['property'] == 'a_A' for r in rows))

    def test_duplicate_rejected(self):
        r = dict(source='NativeSkala', method='native', solid='C', a_A=3.5, B0_GPa=400)
        with self.assertRaises(ValueError):
            matched_statistics([r, r], {}, ['native'], ['C'])

    def test_curated_source_columns(self):
        zhang = {(r['method'], r['solid']): r for r in read_rows(ROOT / 'reference/zhang2018_selected.csv')}
        self.assertEqual(len(zhang), 54)
        self.assertEqual(float(zhang['PBE', 'C']['a_A']), 3.572)
        self.assertEqual(float(zhang['PBE', 'C']['a_ZPE_included_A']), 3.586)
        mo = {r['solid']: r for r in read_rows(ROOT / 'reference/mo2017_selected.csv')}
        self.assertEqual(float(mo['AlP']['B0_GPa']), 89.3)
        self.assertEqual(float(mo['BP']['B0_GPa']), 171.5)

    def test_identical_sets_for_each_ranking(self):
        rows = read_rows(ROOT / 'results/eos-literature-matched-statistics.csv')
        for group in ('Goldzak2022', 'Zhang2018', 'Schimka2011', 'Mo2017', 'PeriodicGFN2Manuscript'):
            for prop in ('a_A', 'B0_GPa'):
                subset = [r for r in rows if r['comparison_set'] == group and r['property'] == prop]
                self.assertLessEqual(len({r['solids'] for r in subset}), 1)
                self.assertLessEqual(len({r['N'] for r in subset}), 1)


if __name__ == '__main__':
    unittest.main()
