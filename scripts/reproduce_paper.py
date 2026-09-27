#!/usr/bin/env python3
"""Recompute the published numerical comparisons without CP2K or network access."""
import argparse
import csv
import json
from pathlib import Path
import sys

from data_checks import ROOT
from analyze_molecular import analyze as molecular
from analyze_band_gaps import analyze as band_gaps
from basis_sensitivity import assess
from additional_checks import ae_cutoff, wavefunction_comparison, crystal_dft
from paper_figures import ice_data, ice_figure, molecular_figure, eos_figure, reconstruction_figure

sys.path.insert(0, str(ROOT / 'benchmarks/Goldzak12/scripts'))
from fit_selected_eos import analyze as eos
from compare_selected_literature import analyze as literature
sys.path.insert(0, str(ROOT / 'benchmarks/X23-mini/scripts'))
from analyze_basis_controls import assessment, selected_pairs, matched_basis_errors


def read(path):
    return json.loads((ROOT / path).read_text())


def csv_rows(path):
    with (ROOT / path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--figures', action='store_true')
    args = parser.parse_args()
    destination = args.output.resolve()
    if destination == ROOT or ROOT in destination.parents:
        parser.error('Choose an output directory outside the immutable data repository')
    destination.mkdir(parents=True, exist_ok=True)
    mol, bands, basis = molecular(), band_gaps(), assess()
    selection = read('benchmarks/Goldzak12/protocol/paper-eos-selection.json')
    fits, fit_summary = eos(csv_rows('benchmarks/Goldzak12/results/eos-selected-source.csv'), selection)
    crystal = read('benchmarks/X23-mini/results/lattice-energies.json')
    results = {'molecular': mol, 'band_gaps': bands, 'basis_controls': basis,
               'ae_cutoff': ae_cutoff(), 'LC10_fits': fits, 'LC10_fit_summary': fit_summary,
               'LC10_literature': literature(fits), 'crystals': selected_pairs(crystal['complete_pairs']),
               'crystal_basis': matched_basis_errors(assessment()),
               'crystal_wavefunction': wavefunction_comparison(), 'crystal_DFT': crystal_dft()}
    for name, data in results.items():
        (destination / (name + '.json')).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    write_csv(destination / 'molecular.csv', [{'reaction': r['reaction'], 'reference_kcal_mol': r['reference_kcal_mol'],
                                             **r['energies_kcal_mol']} for r in mol['rows']])
    write_csv(destination / 'LC10_fits.csv', fits)
    ice, ice_statistics = ice_data({'rows': basis['ice']})
    write_csv(destination / 'ice13.csv', ice)
    (destination / 'ice13-statistics.json').write_text(json.dumps(ice_statistics, indent=2)+'\n')
    improvements = []
    for row in mol['rows']:
        item = {'reaction': row['reaction']}
        for name, reference in [('reference', row['reference_kcal_mol']), ('gauxc', row['energies_kcal_mol']['gauxc_ae'])]:
            item[f'improvement_{name}_kj_mol'] = (abs(row['energies_kcal_mol']['ae-tz']-reference)
                                                 - abs(row['energies_kcal_mol']['ae-qz']-reference))*4.184
        improvements.append(item)
    write_csv(destination / 'molecular-basis-improvement.csv', improvements)
    if args.figures:
        import matplotlib.pyplot as plt
        counts = {name: sum(r[f'improvement_{name}_kj_mol'] > 0 for r in improvements)
                  for name in ('reference', 'gauxc')}
        figures = {'ice13-energy-errors': ice_figure(ice)[0],
                   'common62-basis-figure': molecular_figure(improvements, counts)[0],
                   'eos-literature-residuals': eos_figure(results['LC10_literature']['values'], selection),
                   'native-reconstruction': reconstruction_figure()}
        for name, figure in figures.items():
            figure.savefig(destination / (name + '.pdf'))
            figure.savefig(destination / (name + '.png'), dpi=200)
            plt.close(figure)
    print(json.dumps({'output': str(destination), 'reactions': len(mol['rows']),
                      'band_gaps': bands['primary_gaps'], 'ice_phases': len(ice), 'LC10_fits': len(fits)}))


if __name__ == '__main__':
    main()
