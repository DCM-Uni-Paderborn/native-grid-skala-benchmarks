"""Generate manuscript tables from the closed structural analysis, without SCF."""
import json
import re
from collections import defaultdict

from fit_eos import ROOT, read_rows, write_rows


LABELS = {
    'gapwxc-gth': r'GX', 'gapw-ae': 'AE',
    'hybrid-direct': 'HD', 'hybrid-one-center': 'HOC',
}
SOURCES = {
    'Goldzak2022': ('Goldzak2022Solids', 'ten-solid set'),
    'Zhang2018': ('Zhang2018Solids', 'nine solids, excluding AlN'),
    'Schimka2011': ('Schimka2011HSEsol', 'nine solids, excluding MgS'),
    'Mo2017': ('Mo2017Solids', 'eight solids, excluding AlN and BN'),
    'PeriodicGFN2Manuscript': ('Alizadeh2026PeriodicGFN2', 'nine solids, excluding MgO'),
}


def table(caption, label, headers, rows):
    return '\n'.join([
        r'\begin{table*}[t]', r'\centering', r'\small',
        r'\caption{' + caption + '}', r'\label{' + label + '}',
        r'\begin{tabular}{' + 'l' + 'r' * (len(headers)-1) + '}',
        r'\toprule', ' & '.join(headers) + r' \\', r'\midrule',
        *[' & '.join(str(x) for x in row) + r' \\' for row in rows],
        r'\bottomrule', r'\end{tabular}', r'\end{table*}', '',
    ])


def values(text, keyword):
    return tuple(sorted(set(re.findall(r'^\s*' + keyword + r'\s+([^!\n]+)', text, re.M))))


def main():
    output = ROOT / 'paper'
    output.mkdir(exist_ok=True)
    selection = json.loads((ROOT / 'protocol/paper-eos-selection.json').read_text())
    fits = read_rows(ROOT / 'results/eos-selected-fits.csv')
    refs = {r['solid']: r for r in read_rows(ROOT / 'reference/goldzak2022.csv')
            if r['method'] == 'experiment'}
    summary = read_rows(ROOT / 'results/eos-selected-reference-statistics.csv')
    main_rows = [[LABELS[r['method']], r['N']] +
                 [f"{float(r[k]):.4f}" for k in ('ME_a_A', 'MAE_a_A')] +
                 [f"{float(r[k]):.2f}" for k in ('ME_B0_GPa', 'MAE_B0_GPa')]
                 for r in summary]
    (output / 'periodic-summary.tex').write_text(table(
        r'Structural errors on the same ten solids relative to the static-lattice experimental references of Ref.~\citenum{Goldzak2022Solids}. '
        r'GX denotes \texttt{GAPW\_XC}--\GTH{}, AE denotes \GAPW{}-AE, and HD/HOC denote mixed \GAPW{}-AE/\GTH{} with direct/one-center GTH fields. '
        r'Errors are calculation minus reference; $a_0$ is in \AA{} and $B_0$ in GPa. These statistics describe the selected data, not a bound on discretization error.',
        'tab:goldzak12_plan', ['Representation', '$N$', 'ME($a_0$)', 'MAE($a_0$)', 'ME($B_0$)', 'MAE($B_0$)'], main_rows))
    lines = []
    for prop, unit, decimals in [('a0_A', r'\AA{}', 4), ('B0_GPa', 'GPa', 2)]:
        ref_prop = 'a_A' if prop == 'a0_A' else prop
        rows = []
        for s in selection['solids']:
            rows.append([s, f"{float(refs[s][ref_prop]):.{decimals}f}"] + [
                f"{float(next(r[prop] for r in fits if r['method'] == m and r['solid'] == s)):.{decimals}f}"
                for m in selection['methods']])
        symbol = '$a_0$' if prop == 'a0_A' else '$B_0$'
        lines.append(table('Selected ' + symbol + ' values (' + unit + r'). Experimental references follow Ref.~\citenum{Goldzak2022Solids}. '
                           r'BN and C reuse the AE curve for HD and HOC. These entries are not independent calculations.',
                           'tab:selected-' + prop, ['Solid', 'Experiment', 'GX', 'AE', 'HD', 'HOC'], rows))
    fit_rows = [[LABELS[r['method']], r['solid'], r['n_points'], f"{float(r['B0_prime']):.3f}",
                 f"{float(r['fit_rms_meV_atom']):.4f}"] for r in fits if r['method'] == r['source_method']]
    for i in range(0, len(fit_rows), 18):
        lines.append(table(r'EOS fit diagnostics, part ' + str(i//18+1) + r'. All full-window minima are bracketed. RMS residuals are in meV per atom.',
                           'tab:eos-fit-' + str(i//18+1), ['Method', 'Solid', '$N_V$', "$B'_0$", 'RMS'], fit_rows[i:i+18]))
    stats = read_rows(ROOT / 'results/eos-literature-matched-statistics.csv')
    for source, (citation, scope) in SOURCES.items():
        subset = [r for r in stats if r['comparison_set'] == source]
        indexed = defaultdict(dict)
        for row in subset:
            indexed[(row['source'], row['method'])][row['property']] = row
        rows = []
        for (_, method), props in indexed.items():
            a = props['a_A']
            b = props.get('B0_GPa')
            rows.append([LABELS.get(method, method), a['N'], f"{float(a['MAE']):.4f}",
                         f"{float(b['MAE']):.2f}" if b else '--'])
        lines.append(table(r'Matched comparison with Ref.~\citenum{' + citation + '}: ' + scope +
                           r'. Every row uses the identical intersection and the Goldzak static-lattice experimental reference. '
                           r'MAE($a_0$) is in \AA{} and MAE($B_0$) in GPa. A dash denotes unavailable data, not zero error. '
                           r'Native rows retain the numerical qualifications discussed in the text.',
                           'tab:literature-' + source, ['Method', '$N$', 'MAE($a_0$)', 'MAE($B_0$)'], rows))
    # Audit actual archived inputs; do not infer production settings from templates.
    records = json.loads((ROOT / 'results/eos-selected-input-provenance.json').read_text())
    settings = []
    meshes = defaultdict(set)
    for row in records:
        method, solid, point = row['key'].split('/')
        text = (ROOT / 'accepted-inputs' / row['key'] / 'actual-input.inp').read_text()
        entry = {'key': row['key']}
        for keyword in ('CUTOFF', 'REL_CUTOFF', 'RADIAL_GRID', 'LEBEDEV_GRID', 'SCHEME', 'EPS_SCF'):
            found = values(text, keyword)
            assert len(found) == 1, (row['key'], keyword, found)
            entry[keyword] = found[0].strip()
        for key, expected in [('CUTOFF', '800'), ('REL_CUTOFF', '60'), ('RADIAL_GRID', '150'), ('LEBEDEV_GRID', '770')]:
            assert entry[key] == expected, (row['key'], key)
        meshes[(solid, point)].add(entry['SCHEME'])
        settings.append(entry)
    assert all(len(x) == 1 for x in meshes.values()), 'Unmatched k meshes across methods'
    write_rows(ROOT / 'results/eos-selected-actual-settings.csv', list(settings[0]), sorted(settings, key=lambda x: x['key']))
    grid_rows = []
    for s in selection['solids']:
        groups = defaultdict(list)
        for (solid, point), mesh in sorted(meshes.items()):
            if solid == s:
                groups[next(iter(mesh)).split()[-1]].append(point[1:])
        description = '; '.join('$' + k + '^3$: ' + ','.join(points) for k, points in sorted(groups.items()))
        grid_rows.append([s, description])
    lines.append(table(r'Actual Gamma-centered k meshes by volume index. Indices 01--10 refer to the original evenly spaced volume ratios 0.90--1.10. '
                       r'Each selected volume uses the same mesh in all representations, but the mesh can change along an EOS.',
                       'tab:actual-k-meshes', ['Solid', 'Mesh and volume indices'], grid_rows))
    (output / 'periodic-tables-si.tex').write_text('\n'.join(lines))
    print(json.dumps({'tables_generated': 11, 'actual_inputs_audited': len(settings), 'same_mesh_per_paired_volume': True}))


if __name__ == '__main__':
    main()
