"""Offline checks for the additional numerical controls and publication sources."""
import json
from pathlib import Path
import re

from data_checks import ROOT, close, metrics, sha, total_energy


def ae_cutoff():
    data = json.loads((ROOT / 'convergence/ae-cutoff/dataset.json').read_text())
    energies = {r['cutoff_Ry']: total_energy(r) for r in data['cases']}
    assert set(energies) == {400, 600, 800, 1000, 1200}
    reference = energies[data['reference_cutoff_Ry']]
    return [{'cutoff_Ry': cutoff, 'total_energy_Ha': energy,
             'difference_kJ_mol': (energy-reference)*data['hartree_to_kjmol']}
            for cutoff, energy in sorted(energies.items())]


def adjoint():
    folder = ROOT / 'convergence/adjoint'
    evidence = json.loads((folder / 'evidence.json').read_text())
    for name, digest in evidence['source_reports_sha256'].items():
        assert sha(folder / name) == digest
    assert sha(folder / evidence['fixture']['file']) == evidence['fixture']['sha256']
    for name, digest in evidence['reports']['optimized']['generated_sources_sha256'].items():
        actual = evidence['fixture']['file'] if name == 'gapw_atom_adjoint_fixture.F90' else name
        assert sha(folder / actual) == digest
    for report in evidence['reports'].values():
        assert report['exit_code'] == 0 and not report['stderr']
        assert 'PASS: transpose' in report['output']
    assert evidence['dot_product_identity']['measured_maximum'] is None
    return {k: evidence[k] for k in ('dot_product_identity', 'projection_variants')}


def wavefunction_comparison():
    data = json.loads((ROOT / 'benchmarks/X23-mini/reference/wavefunction-comparison.json').read_text())
    rows = []
    for method in data['methods']:
        values = [v * method['sign_multiplier'] for v in method['source_energies']]
        assert values == method['lattice_energies'] and len(values) == 3
        error = metrics([a-b for a, b in zip(values, data['reference']['lattice_energies'])])
        assert f"{error['MAE']:.2f}" == method['mae_display']
        rows.append({'method': method['label'], 'energies_kJ_mol': values, 'statistics': error})
    return rows


def crystal_dft():
    data = json.loads((ROOT / 'benchmarks/X23-mini/reference/dft-comparison.json').read_text())
    return [{'method': name, 'energies_kJ_mol': energies,
             'statistics': metrics([v-r for v, r in zip(energies, data['DMC_reference'])])}
            for name, energies in data['methods'].items()]


def crystal_comparison():
    data = json.loads((ROOT / 'benchmarks/X23-mini/results/lattice-energies.json').read_text())
    reference = json.loads((ROOT / 'benchmarks/X23-mini/reference/wavefunction-comparison.json').read_text())
    population = reference['population']
    dmc = reference['reference']['lattice_energies']
    pairs = {(r['method'], r['system']): r for r in data['paper_pairs']}
    rows = []
    for method, label in [('gapw-ae', 'GAPW-AE'),
                          ('gapw-gth-one-center', 'GAPW-GTH one-centre'),
                          ('gapw-gth-direct', 'GAPW-GTH direct'),
                          ('gapwxc-gth', 'GAPW-XC/GTH')]:
        energies = [pairs[method, system]['lattice_energy_kjmol'] for system in population]
        for system, value in zip(population, dmc):
            close(pairs[method, system]['dmc_kjmol'], value)
        rows.append({'method': label, 'energies_kJ_mol': energies,
                     'statistics': metrics([v-r for v, r in zip(energies, dmc)])})
    rows.extend(crystal_dft())
    rows.extend(wavefunction_comparison())
    return rows


def check_crystal_table(text, rows):
    # Read the simple numeric tabular following this stable publication label.
    table = text.split(r'\label{tab:crystals}', 1)[1].split(r'\end{tabular}', 1)[0]
    cells = []
    for line in table.splitlines():
        if '&' not in line or line.strip().startswith('Method &'):
            continue
        values = [cell.strip().removesuffix(r'\\').strip() for cell in line.split('&')]
        values[0] = re.sub(r'\\cite\w*\{[^}]+\}', '', values[0])
        cells.append(values)
    expected_names = [r['method'] for r in rows] + ['DMC reference']
    if [r[0] for r in cells] != expected_names:
        raise ValueError('Crystal table and deposited method populations differ')
    for printed, row in zip(cells[:-1], rows):
        values = row['energies_kJ_mol'] + [row['statistics']['MAE']]
        if len(printed) != len(values) + 1:
            raise ValueError('Unexpected crystal table columns')
        for cell, value in zip(printed[1:], values):
            decimals = len(cell.split('.')[1])
            close(float(cell), value, 0.5 * 10**(-decimals))
    reference = json.loads((ROOT / 'benchmarks/X23-mini/reference/wavefunction-comparison.json').read_text())['reference']
    expected_dmc = [f'{energy:.1f}({round(uncertainty * 10)})' for energy, uncertainty in
                    zip(reference['lattice_energies'], reference['statistical_uncertainties'])]
    if cells[-1][1:] != expected_dmc + ['--']:
        raise ValueError('Crystal table and deposited DMC references differ')


def paper_coverage():
    folder = ROOT / 'paper/pccp'
    snapshot = json.loads((folder / 'source-snapshot.json').read_text())
    for name, digest in snapshot['files'].items():
        assert sha(folder / name) == digest, name
    texts = {name: (folder / name).read_text() for name in ('main.tex', 'supplementary_information.tex')}
    check_crystal_table(texts['main.tex'], crystal_comparison())
    assert len(list(folder.glob('*.tex'))) == 2
    labels = {name: re.findall(r'\\label\{([^}]+)\}', text) for name, text in texts.items()}
    coverage = json.loads((folder / 'data-coverage.json').read_text())
    floats = {label for group in labels.values() for label in group if label.startswith(('tab:', 'fig:'))}
    assert floats == coverage['floats'].keys(), (floats ^ coverage['floats'].keys())
    for entry in coverage['floats'].values():
        for path in entry['paths']:
            assert (ROOT / path).exists(), path
    bibliography = (folder / 'references.bib').read_text()
    bibkeys = re.findall(r'^@\w+\{([^,\s]+),', bibliography, re.M)
    assert len(bibkeys) == len(set(bibkeys))
    for name, text in texts.items():
        assert len(labels[name]) == len(set(labels[name]))
        assert not re.search(r'\\(?:input|include)\{', text)
        other = 'supplementary_information.tex' if name == 'main.tex' else 'main.tex'
        prefix = 'esi-' if name == 'main.tex' else 'paper-'
        for label in re.findall(r'\\(?:ref|eqref|pageref)\{([^}]+)\}', text):
            assert label in labels[name] or label.startswith(prefix) and label[len(prefix):] in labels[other], label
        for group in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}', text):
            assert set(group.split(',')) <= set(bibkeys), group
        for filename in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', text):
            assert (folder / 'figures' / filename).is_file(), filename
    return {'tables': sum(x.startswith('tab:') for x in floats),
            'figures': sum(x.startswith('fig:') for x in floats), 'citations': len(bibkeys)}
