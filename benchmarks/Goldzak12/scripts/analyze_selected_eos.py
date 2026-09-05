"""Build structural comparisons and endpoint sensitivity without new SCF runs."""
import json
from collections import defaultdict

from fit_eos import ROOT, HA_PER_A3_TO_GPA, aggregate, fit_curve, read_rows, write_rows
from fit_selected_eos import analyze, select_curves


def statistics(rows, group):
    result = []
    groups = defaultdict(list)
    for row in rows:
        groups[row[group]].append(row)
    for label, values in sorted(groups.items()):
        record = {group: label, 'N': len(values)}
        for prop in ('a_A', 'B0_GPa'):
            for name, value in zip(('ME', 'MAE', 'RMSE', 'MaxAE'),
                                   aggregate([r['delta_' + prop] for r in values])):
                record[name + '_' + prop] = value
        result.append(record)
    return result


def main():
    selection = json.loads((ROOT / 'protocol/paper-eos-selection.json').read_text())
    rows = read_rows(ROOT / 'results/eos-selected-source.csv')
    fits, summary = analyze(rows, selection)
    lookup = {(r['method'], r['solid']): r for r in fits}
    refs = {r['solid']: r for r in read_rows(ROOT / 'reference/goldzak2022.csv')
            if r['method'] == 'experiment'}
    errors = []
    for row in fits:
        ref = refs[row['solid']]
        errors.append({'method': row['method'], 'solid': row['solid'],
                       'delta_a_A': row['a0_A'] - float(ref['a_A']),
                       'delta_B0_GPa': row['B0_GPa'] - float(ref['B0_GPa']),
                       'quality_status': 'numerical_review_pending'})
    pairs = []
    for left, right in [('gapwxc-gth', 'gapw-ae'), ('hybrid-direct', 'gapw-ae'),
                        ('hybrid-one-center', 'gapw-ae'),
                        ('hybrid-one-center', 'hybrid-direct')]:
        for solid in selection['solids']:
            a, b = lookup[left, solid], lookup[right, solid]
            pairs.append({'comparison': left + '_minus_' + right, 'solid': solid,
                          'delta_a_A': a['a0_A'] - b['a0_A'],
                          'delta_B0_GPa': a['B0_GPa'] - b['B0_GPa'],
                          'shared_ae_reuse': solid in selection['shared_ae_solids'],
                          'quality_status': 'numerical_review_pending'})
    sensitivity, seen = [], set()
    for (_, solid), (source, chosen) in select_curves(rows, selection).items():
        if (source, solid) in seen:
            continue
        seen.add((source, solid))
        base = lookup[source, solid]
        for label, subset in [('omit_smallest', chosen[1:]), ('omit_largest', chosen[:-1])]:
            p, _, _ = fit_curve(subset)
            sensitivity.append({'method': source, 'solid': solid, 'test': label,
                                'n_points': len(subset),
                                'delta_a_A': float(p[1] ** (1 / 3) - base['a0_A']),
                                'delta_B0_GPa': float(p[2] * HA_PER_A3_TO_GPA - base['B0_GPa']),
                                'minimum_bracketed': bool(float(subset[0]['volume_A3']) < p[1] < float(subset[-1]['volume_A3']))})
    tables = {'eos-selected-reference-errors.csv': errors,
              'eos-selected-reference-statistics.csv': statistics(errors, 'method'),
              'eos-selected-method-differences.csv': pairs,
              'eos-selected-method-statistics.csv': statistics(pairs, 'comparison'),
              'eos-selected-window-sensitivity.csv': sensitivity}
    tables['eos-selected-method-statistics-nonshared.csv'] = statistics(
        [r for r in pairs if not r['shared_ae_reuse']], 'comparison')
    literature = []
    for row in read_rows(ROOT / 'reference/goldzak2022.csv'):
        if row['solid'] in selection['solids']:
            literature.append({'source': 'Goldzak2022', 'method': row['method'],
                               'solid': row['solid'], 'a_A': row['a_A'], 'B0_GPa': row['B0_GPa']})
    for row in read_rows(ROOT / 'reference/periodic_gfn2_lc10.csv'):
        if row['solid'] in selection['solids'] and row['a_A']:
            literature.append({'source': 'PeriodicGFN2Manuscript', 'method': row['method'],
                               'solid': row['solid'], 'a_A': row['a_A'], 'B0_GPa': ''})
    tables['eos-selected-literature-values.csv'] = literature
    for name, values in tables.items():
        write_rows(ROOT / 'results' / name, list(values[0]), values)
    summary['maximum_endpoint_delta_a_A'] = max(abs(r['delta_a_A']) for r in sensitivity)
    summary['maximum_endpoint_delta_B0_GPa'] = max(abs(r['delta_B0_GPa']) for r in sensitivity)
    summary['endpoint_test_is_numerical_convergence_estimate'] = False
    (ROOT / 'results/eos-selected-analysis.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary, allow_nan=False))


if __name__ == '__main__':
    main()
