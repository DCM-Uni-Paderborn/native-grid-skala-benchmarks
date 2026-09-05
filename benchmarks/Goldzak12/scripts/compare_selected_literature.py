"""Compare existing EOS fits with literature on explicitly matched solid sets."""
import json
import math
from collections import defaultdict

from fit_eos import ROOT, read_rows, write_rows

LABELS = {'gapwxc-gth': 'Skala GAPW_XC-GTH', 'gapw-ae': 'Skala GAPW-AE',
          'hybrid-direct': 'Skala hybrid direct',
          'hybrid-one-center': 'Skala hybrid +1c'}
SOURCES = {
    'Goldzak2022': 'https://doi.org/10.1063/5.0119633',
    'Zhang2018': 'https://pure.mpg.de/rest/items/item_2599505_10/component/file_2618826/content',
    'Schimka2011': 'https://doi.org/10.1063/1.3524336',
    'Mo2017': 'https://doi.org/10.1103/PhysRevB.95.035118',
}


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
    for source in ('Goldzak2022', 'Zhang2018', 'Schimka2011', 'Mo2017',
                   'PeriodicGFN2Manuscript'):
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


def table(headers, rows):
    return ['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---']*len(headers)) + ' |'] + [
        '| ' + ' | '.join(str(c) for c in row) + ' |' for row in rows]


def main():
    selection = json.loads((ROOT / 'protocol/paper-eos-selection.json').read_text())
    solids, methods = selection['solids'], selection['methods']
    fits = read_rows(ROOT / 'results/eos-selected-fits.csv')
    records = [dict(source='NativeSkala', method=r['method'], solid=r['solid'],
                    a_A=r['a0_A'], B0_GPa=r['B0_GPa']) for r in fits]
    for filename, source in [('goldzak2022.csv', 'Goldzak2022'),
                             ('periodic_gfn2_lc10.csv', 'PeriodicGFN2Manuscript'),
                             ('zhang2018_selected.csv', 'Zhang2018'),
                             ('schimka2011_selected.csv', 'Schimka2011'),
                             ('mo2017_selected.csv', 'Mo2017')]:
        for r in read_rows(ROOT / 'reference' / filename):
            if r['solid'] in solids and r['a_A']:
                records.append(dict(source=source, method=r['method'], solid=r['solid'],
                                    a_A=r['a_A'], B0_GPa=r.get('B0_GPa', '')))
    refs = {r['solid']: r for r in records if r['source'] == 'Goldzak2022' and r['method'] == 'experiment'}
    stats = matched_statistics(records, refs, methods, solids)
    alternate_refs = {s: dict(r) for s, r in refs.items()}
    for s, b in {'BN': 410.2, 'BP': 168.0, 'MgO': 169.8}.items():
        alternate_refs[s]['B0_GPa'] = b
    alternate = [r for r in matched_statistics(records, alternate_refs, methods, solids)
                 if r['comparison_set'] == 'Goldzak2022' and r['property'] == 'B0_GPa']
    for r in alternate:
        r['reference'] = 'Goldzak Table III alternate B0 for BN/BP/MgO; others primary'
    write_rows(ROOT / 'results/eos-literature-alternate-reference-statistics.csv', list(alternate[0]), alternate)
    write_rows(ROOT / 'results/eos-literature-matched-statistics.csv', list(stats[0]), stats)
    write_rows(ROOT / 'results/eos-literature-comparison-values.csv', list(records[0]), records)
    residuals = []
    for r in records:
        if r['method'] == 'experiment':
            continue
        for prop in ('a_A', 'B0_GPa'):
            if r[prop] != '':
                residuals.append(dict(source=r['source'], method=r['method'], solid=r['solid'],
                                      property=prop, value=float(r[prop]), reference=float(refs[r['solid']][prop]),
                                      error=float(r[prop])-float(refs[r['solid']][prop])))
    write_rows(ROOT / 'results/eos-literature-signed-errors.csv', list(residuals[0]), residuals)

    lines = ['# Structural results and literature comparison', '',
             'Campaign closed; existing energies only. Numerical release remains pending.', '',
             '## Scope and conventions', '',
             'Ten cubic solids: ' + ', '.join(solids) + '. LiF and LiH are excluded consistently across methods.',
             '352 unique energies, 36 independent curves and 40 method/solid rows. AlN/AlP use v01-v09; others v01-v10.',
             'BN and C reuse the AE curve in both hybrid representations. The hybrid label denotes mixed AE/GTH treatment, not a hybrid exchange-correlation functional.',
             'All signed errors below use calculation minus the same Goldzak experimental reference, corrected to the static-lattice convention. No method-specific reference correction is fitted.',
             'Results with different source protocols are kept separate. Each comparison uses the exact common solid set, printed below.',
             'This is a comparison with the retrievable relevant benchmark literature, not an exhaustive inventory of every published value.', '',
             '## Native results', '']
    for prop, title in [('a_A', 'Lattice constant a0 (angstrom)'), ('B0_GPa', 'Bulk modulus B0 (GPa)')]:
        lines += ['### ' + title, '']
        matrix = []
        for s in solids:
            vals = [next(r[prop] for r in records if r['source'] == 'NativeSkala' and r['method'] == m and r['solid'] == s) for m in methods]
            matrix.append([s, f"{float(refs[s][prop]):.4f}"] + [f'{float(v):.4f}' for v in vals])
        lines += table(['Solid', 'Experiment', 'GX', 'AE', 'Hybrid direct', 'Hybrid +1c'], matrix) + ['']

    for source in ('Goldzak2022', 'Zhang2018', 'Schimka2011', 'Mo2017', 'PeriodicGFN2Manuscript'):
        lines += ['## ' + source, '']
        if source in SOURCES:
            lines += ['[Primary source](' + SOURCES[source] + ')', '']
        else:
            lines += ['Own periodic GFN manuscript, Table S3; not the original GFN2-xTB paper.', '']
        if source == 'Zhang2018':
            lines += ['SI Tables III/V: use the Uncorr. (static electronic) columns, not Corr. columns that include zero-point motion. SCAN was evaluated on PBE orbitals/densities. AlN is absent.', '']
        if source == 'Schimka2011':
            lines += ['Table III: lattice constants only; MgS absent. Bulk-modulus values were not extracted from its separate supplement.', '']
        if source == 'Mo2017':
            lines += ['Tables II/III: self-consistent Tao-Mo results; AlN and BN absent. Experimental references in that paper differ from Goldzak; errors here are recomputed against Goldzak.', '']
        if source == 'PeriodicGFN2Manuscript':
            lines += ['Nine-solid intersection: the manuscript LC10 excludes MgO/LiH, whereas this selected set excludes LiF/LiH. No bulk-modulus values are available in the curated table.', '']
        for prop in ('a_A', 'B0_GPa'):
            rows = [r for r in stats if r['comparison_set'] == source and r['property'] == prop]
            if not rows:
                continue
            lines += ['### ' + prop + ': matched errors', '', 'Solids: ' + rows[0]['solids'].replace(';', ', ') + '.', '']
            lines += table(['Method', 'N', 'ME', 'MAE', 'RMSE', 'MaxAE', 'MARE (%)'], [
                [LABELS.get(r['method'], r['method']), r['N']] + [f'{r[k]:.5f}' for k in ('ME', 'MAE', 'RMSE', 'MaxAE', 'MARE_percent')] for r in rows]) + ['']
            literature = [r for r in records if r['source'] == source and r['method'] != 'experiment' and r[prop] != '']
            names = list(dict.fromkeys(r['method'] for r in literature))
            lines += ['### ' + prop + ': literature values', '']
            lines += table(['Solid'] + names, [[s] + [f"{float(next(r[prop] for r in literature if r['method'] == m and r['solid'] == s)):.4f}" for m in names] for s in rows[0]['solids'].split(';')]) + ['']

    lines += ['## Alternate experimental B0 sensitivity', '',
              'Use the three alternate values listed in Goldzak Table III simultaneously, keeping all other primary references and all ten solids unchanged. This is a reference-choice sensitivity test, not another fitted result.', '']
    lines += table(['Method', 'N', 'MAE B0 primary (GPa)', 'MAE B0 alternate (GPa)'], [
        [LABELS.get(r['method'], r['method']), r['N'],
         f"{next(t['MAE'] for t in stats if t['comparison_set'] == 'Goldzak2022' and t['property'] == 'B0_GPa' and t['source'] == r['source'] and t['method'] == r['method']):.4f}",
         f"{r['MAE']:.4f}"] for r in alternate]) + ['']
    lines += ['## Interpretation and limits', '',
              '- Every native lattice constant is below the selected experimental reference. The contraction is systematic, not confined to one outlier.',
              '- Similarity among representations does not establish accuracy against experiment or validate the implementation independently. Functional, basis/pseudopotential and numerical contributions remain entangled.',
              '- The hybrid one-center correction changes a0 by at most 0.000542 angstrom; this does not remove the contraction. Its maximum B0 effect is 5.752 GPa (MgO).',
              '- N=10 one-center averages include two reused AE zeros. The N=8 independent comparison is in eos-selected-method-statistics-nonshared.csv.',
              '- Endpoint omission shifts a0 by up to 0.002444 angstrom and B0 by up to 9.997 GPa. These are fit-window sensitivity tests, not cutoff/k-point error bars.',
              '- Goldzak lists alternate experimental B0 values for BN (410.2 versus 388.5 GPa), BP (168.0 versus 176.5) and MgO (169.8 versus 173.0). Reference sensitivity must not be confused with numerical convergence.',
              '- Cohesive energies are not compared: no complete selected atomic-reference analysis is being claimed.', '',
              '## Aggregate-only literature is not a matched ranking', '',
              'Goldzak Table I lists PBE/PBEsol/SCAN MAEs of 0.061/0.030/0.030 angstrom and 12.2/7.8/7.4 GPa. These coincide with the 44-solid statistics in Tran, Stelzl and Blaha Table II; they must not be treated as newly calculated errors on this ten-solid set. The old aggregate metadata labeled them Goldzak12; this scope has been corrected to contextual/unspecified in that file.',
              '[Tran et al., JCP 144, 204120 (2016)](https://doi.org/10.1063/1.4948636). Its system-resolved supplement was not retrieved; no missing numbers were inferred from the aggregate.',
              'The system-resolved Zhang comparison above provides an independently sourced matched PBE/PBEsol/SCAN comparison instead.', '',
              '## Reproduction', '',
              '`python3 -B benchmarks/Goldzak12/scripts/compare_selected_literature.py`', '',
              'Inputs: selected fits, selection JSON and five curated reference CSVs. All comparison CSVs and this report are generated offline. No CP2K jobs, energy ledgers or acceptance criteria are changed.', '']
    (ROOT / 'results/LITERATURE-COMPARISON.md').write_text('\n'.join(lines))
    plot_residuals(records, refs, solids, methods)
    print(json.dumps(dict(values=len(records), statistics=len(stats), solids=len(solids), quality_release=False)))


def plot_residuals(records, refs, solids, methods):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), layout='constrained')
    palette = ['#087e8b', '#b23a48', '#425b2b', '#735d9d']
    for ax, prop, xlabel in zip(axes, ('a_A', 'B0_GPa'),
                               ('a0 - reference (angstrom)', 'B0 - reference (GPa)')):
        for i, m in enumerate(methods):
            values = {r['solid']: float(r[prop]) for r in records if r['source'] == 'NativeSkala' and r['method'] == m}
            ax.scatter([values[s]-float(refs[s][prop]) for s in solids],
                       [j+(i-1.5)*.12 for j in range(len(solids))],
                       label=LABELS[m], color=palette[i], s=28, marker=['o', 's', '^', 'x'][i], zorder=3)
        scs = {r['solid']: float(r[prop]) for r in records if r['source'] == 'Goldzak2022' and r['method'] == 'SCS-MP2'}
        ax.scatter([scs[s]-float(refs[s][prop]) for s in solids], range(len(solids)),
                   label='SCS-MP2 (Goldzak)', color='#292929', s=24, marker='D', zorder=3)
        ax.axvline(0, color='#555555', linewidth=.8)
        ax.set_yticks(range(len(solids)), solids)
        ax.invert_yaxis()
        ax.set_xlabel(xlabel)
        ax.grid(axis='x', color='#dddddd', linewidth=.6)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].legend(fontsize=8, loc='lower left')
    fig.suptitle('Selected ten-solid EOS results | numerical review pending', fontsize=13)
    fig.savefig(ROOT / 'results/eos-literature-residuals.png', dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
