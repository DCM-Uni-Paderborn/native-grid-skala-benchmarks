"""Matched LC10 literature comparisons used in the current manuscript."""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

from fit_eos import ROOT, read_rows, write_rows

SOURCES = {'Goldzak2022': 'https://doi.org/10.1063/5.0119633',
           'Zhang2018': 'https://pure.mpg.de/rest/items/item_2599505_10/component/file_2618826/content'}


def error_statistics(values, refs):
    if not values or len(values) != len(refs):
        raise ValueError('Empty or mismatched comparison')
    if any(not math.isfinite(v) for v in values + refs) or any(r <= 0 for r in refs):
        raise ValueError('Invalid structural value')
    d = [v-r for v, r in zip(values, refs)]
    return dict(N=len(d), ME=sum(d)/len(d), MAE=sum(map(abs, d))/len(d),
                RMSE=math.sqrt(sum(x*x for x in d)/len(d)), MaxAE=max(map(abs, d)),
                MARE_percent=100*sum(abs(x/r) for x, r in zip(d, refs))/len(d))


def matched_statistics(records, references, native_methods, solids):
    available = defaultdict(dict)
    for row in records:
        key = (row['source'], row['method'])
        if row['solid'] in available[key]:
            raise ValueError('Duplicate method/solid: ' + str(key))
        available[key][row['solid']] = row
    groups = []
    for source in ('Goldzak2022', 'Zhang2018'):
        keys = [k for k in available if k[0] == source and k[1] != 'experiment']
        for prop in ('a_A', 'B0_GPa'):
            if not keys:
                continue
            subset = [s for s in solids if all(s in available[k] and available[k][s][prop] != ''
                                              for k in keys)]
            if not subset:
                continue
            for key in [('NativeSkala', m) for m in native_methods] + keys:
                vals = [float(available[key][s][prop]) for s in subset]
                refs = [float(references[s][prop]) for s in subset]
                groups.append(dict(comparison_set=source, source=key[0], method=key[1],
                                   property=prop, solids=';'.join(subset),
                                   reference='Goldzak2022 experiment, static-lattice convention',
                                   **error_statistics(vals, refs)))
    return groups


def analyze(fits=None):
    selection = json.loads((ROOT / 'protocol/paper-eos-selection.json').read_text())
    solids, methods = selection['solids'], selection['methods']
    if fits is None:
        fits = read_rows(ROOT / 'results/eos-selected-fits.csv')
    records = [dict(source='NativeSkala', method=r['method'], solid=r['solid'],
                    a_A=float(r['a0_A']), B0_GPa=float(r['B0_GPa'])) for r in fits]
    for filename, source in [('goldzak2022.csv', 'Goldzak2022'), ('zhang2018_selected.csv', 'Zhang2018')]:
        for r in read_rows(ROOT / 'reference' / filename):
            if r['solid'] in solids and r['a_A']:
                records.append(dict(source=source, method=r['method'], solid=r['solid'],
                                    a_A=float(r['a_A']), B0_GPa=float(r['B0_GPa']) if r.get('B0_GPa') else ''))
    refs = {r['solid']: r for r in records if r['source'] == 'Goldzak2022' and r['method'] == 'experiment'}
    stats = matched_statistics(records, refs, methods, solids)
    alternative = {s: dict(r) for s, r in refs.items()}
    for s, b in {'BN': 410.2, 'BP': 168.0, 'MgO': 169.8}.items():
        alternative[s]['B0_GPa'] = b
    alternate_stats = [r for r in matched_statistics(records, alternative, methods, solids)
                       if r['comparison_set'] == 'Goldzak2022' and r['property'] == 'B0_GPa']
    for row in alternate_stats:
        row['reference'] = 'Goldzak Table III alternate B0 for BN/BP/MgO; others primary'
    return {'values': records, 'statistics': stats, 'alternate_B0_statistics': alternate_stats}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for name, rows in analyze().items():
        write_rows(args.output / (name + '.csv'), list(rows[0]), rows)


if __name__ == '__main__':
    main()
