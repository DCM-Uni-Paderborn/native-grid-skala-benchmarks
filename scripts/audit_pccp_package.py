#!/usr/bin/env python3
"""Audit only the published PCCP data, without networking or new calculations."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics

ROOT = Path(__file__).resolve().parents[1]
HA_KCAL = 627.5094740631
HA_KJ = 2625.4996394799


def read(path):
    return json.loads((ROOT / path).read_text())


def rows(path):
    with (ROOT / path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def close(a, b, tolerance=1e-9):
    assert math.isfinite(float(a)) and abs(float(a)-float(b)) <= tolerance, (a, b, tolerance)


def output(path, expected=None):
    text = (ROOT / path).read_text(errors='replace')
    start = text.rfind('PROGRAM STARTED')
    assert start >= 0, path
    block = text[start:]
    assert 'SCF run converged' in block and 'PROGRAM ENDED' in block, path
    assert '[ABORT]' not in block and 'SCF run NOT converged' not in block, path
    values = re.findall(r'ENERGY\| Total FORCE_EVAL.*?([-\d.]+)\s*$', block, re.M)
    assert values, path
    energy = float(values[-1])
    assert math.isfinite(energy), path
    timer = ROOT / path.parent / 'time.txt'
    if timer.exists():
        exits = re.findall(r'Exit status:\s*(-?\d+)', timer.read_text())
        assert exits and exits[-1] == '0', str(timer)
    execution = ROOT / path.parent / 'execution.txt'
    if execution.exists():
        exits = re.findall(r'^exit(?:_status|_code)?=(\d+)$', execution.read_text(), re.M)
        # Early launcher records omit an exit field; GNU time above records it.
        if exits:
            assert exits[-1] == '0', str(execution)
    if expected is not None:
        close(energy, expected)
    return energy


def inventory():
    manifest = read('paper/pccp/file-manifest.json')
    files = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
             if p.is_file() and '.git' not in p.relative_to(ROOT).parts
             and '__pycache__' not in p.parts and p.name != '.DS_Store'}
    files.discard('paper/pccp/file-manifest.json')
    expected = set(manifest['files'])
    assert files == expected, {'missing': sorted(expected-files), 'extra': sorted(files-expected)}
    for path, record in manifest['files'].items():
        assert sha(path) == record['sha256'], path
    return len(files)


def molecular():
    root = Path('benchmarks/dietGMTKN55/production-p25')
    data = root / 'paper-common-70'
    species = {r['route']+'/'+r['digest']: r for r in rows(data/'species-results.csv')}
    assert len(species) == 332
    actual = rows(root/'accepted-inputs/manifest.csv')
    assert len(actual) == 332
    for r in actual:
        key = r['route']+'/'+r['digest']
        assert r['availability'] == 'exact'
        assert sha(r['archived_input']) == r['actual_sha256'], key
        assert sha(r['canonical_input']) == r['canonical_sha256'], key
        output(root/'accepted-inputs'/key/'output.out', species[key]['total_energy_ha'])
    reference = read(data/'reference.json')
    index = read(data/'reaction-index.json')
    assert sha(data/'reference.json') == index['reference_subset_sha256']
    tabulated = {(r['subset'], int(r['reaction_id'])): r for r in rows(data/'native-protocol-comparison.csv')}
    assert len(tabulated) == len(index['reactions']) == 70
    for reaction in index['reactions']:
        key = reaction['subset'], reaction['reaction_id']
        row = tabulated[key]
        close(reference[key[0]][str(key[1])]['Energy'], row['reference_kcal_mol'])
        for protocol in ['gapwxc_gth','hybrid_ae_gth_direct','hybrid_ae_gth_one_center']:
            energy = sum(s['count'] * float(species[s['accepted_results'][protocol]['result_key']]['total_energy_ha'])
                         for s in reaction['species']) * HA_KCAL
            close(energy, row[protocol+'_kcal_mol'], 2e-6)
            close(energy-float(row['reference_kcal_mol']), row[protocol+'_error_kcal_mol'], 2e-6)
    metrics = {}
    for prefix in ['gapwxc_gth','hybrid_ae_gth_direct','hybrid_ae_gth_one_center']:
        errors = [float(r[prefix+'_error_kcal_mol']) for r in tabulated.values()]
        metrics[prefix] = {'MAE_kcal_mol':statistics.mean(map(abs,errors)),
                           'median_absolute_error_kcal_mol':statistics.median(map(abs,errors))}
    gauxc = [r for r in rows(data/'native-and-gauxc-paper-comparison.csv') if r['paper_gpw_gth_available']=='yes']
    assert len(gauxc) == 65
    difference = [float(r['gapwxc_gth_kcal_mol'])-float(r['paper_gpw_gth_with_current_d3_kcal_mol']) for r in gauxc]
    close(statistics.mean(map(abs,difference)), .4624, .00005)
    return {'reactions':70,'unique_species_routes':332,'GauXC_intersection':65,'metrics':metrics}


def solid_eos():
    root = Path('benchmarks/Goldzak12')
    selected = read(root/'protocol/paper-eos-selection.json')
    records = read(root/'results/eos-selected-input-provenance.json')
    data = rows(root/'results/eos-selected-source.csv')
    energies = {'/'.join(r[k] for k in ['method','solid','point']): r for r in data}
    assert len(energies) == len(records) == 352
    assert {r['solid'] for r in data} == set(selected['solids']) and len(selected['solids']) == 10
    for r in records:
        folder = root/'accepted-inputs'/r['key']
        for name, digest in r['files_sha256'].items():
            assert sha(folder/name) == digest, str(folder/name)
        assert sha(folder/'output.out') == r['output_sha256'], str(folder)
        output(folder/'output.out', energies[r['key']]['energy_Ha'])
    fits = rows(root/'results/eos-selected-fits.csv')
    assert len(fits) == 40
    return {'solids':selected['solids'],'unique_energies':352,'independent_EOS':36,'method_solid_comparisons':40}


def crystals():
    root = Path('benchmarks/X23-mini')
    data = read(root/'results/lattice-energies.json')
    records = {}
    for path in (ROOT/root/'results/validated').glob('*/*/*/validation.json'):
        r = json.loads(path.read_text())
        records[r['case_id']] = r
    assert len(records) == 24 and len(data['complete_pairs']) == 12
    for rel, digest in data['validation_file_sha256'].items():
        path = root/rel
        assert sha(path) == digest, str(path)
        r = read(path)
        for name, expected in r['final_hashes'].items():
            if (ROOT/path.parent/name).is_file():
                assert sha(path.parent/name) == expected, str(path.parent/name)
        assert sha(path.parent/'execution.json') == r['execution_sha256']
        output(path.parent/'output.out', r['energy_hartree'])
    for pair in data['complete_pairs']:
        key = pair['method']+'/'+pair['system']
        solid, molecule = records[key+'/solid'], records[key+'/molecule']
        energy = (solid['energy_hartree']/solid['molecules']-molecule['energy_hartree'])*HA_KJ
        close(energy,pair['lattice_energy_kjmol'])
        close(energy-pair['dmc_kjmol'],pair['signed_deviation_kjmol'])
    assert len(data['numerical_controls']) == 7
    for c in data['numerical_controls']:
        paths = [p for p in (ROOT/root/'results/controls').glob('*/validation.json')
                 if (r:=json.loads(p.read_text()))['case_id']==c['case_id'] and r['control_name']==c['control_name']]
        assert len(paths)==1
        r=json.loads(paths[0].read_text()); parent=records[c['case_id']]
        close((r['energy_hartree']-parent['energy_hartree'])/parent['molecules']*HA_KJ,
              c['delta_control_minus_parent_kjmol_per_molecule'])
    return {'base_executions':24,'pairs':12,'numerical_controls':7}


def numerical_controls():
    root = Path('convergence')
    for folder, table in [('cutoff/aconf8-b200','cutoff/aconf8-b200/results.csv'),
                           ('atom-grid/h2o','atom-grid/results/h2o-grid-results.csv')]:
        data = rows(root/table)
        assert all(r['converged']=='true' for r in data)
        for r in data: output(root/folder/r['path']/'output.out',r['total_energy_ha'])
    expected = {('gapwxc-gth','C'):.0200,('gapw-ae','C'):.0229,
                ('hybrid-direct','C'):.0229,('hybrid-one-center','C'):.0229,
                ('gapwxc-gth','MgO'):.0291,('hybrid-direct','MgO'):.0493,
                ('hybrid-one-center','MgO'):.0378}
    symmetry = {}
    for (method,solid), target in expected.items():
        base=root/'symmetry'/method/solid
        delta=abs(output(base/'full/output.out')-output(base/'reduced/output.out'))*HA_KJ
        close(delta,target,.00005)
        symmetry[method+'/'+solid]=delta
    derivative=root/'derivatives/periodic-water-20260909'
    report=read(derivative/'validated-results.json')
    force_errors=[]; stress_errors=[]; count=0
    for method,record in report['methods'].items():
        cases=record['cases']
        for name,r in cases.items():
            for filename,digest in r['sha256'].items():
                assert sha(derivative/'cases'/method/name/filename)==digest
            output(derivative/'cases'/method/name/'output.out',r['energy_Ha']);count+=1
        for kind,steps in [('force',[.001,.0005]),('strain',[.0001,.00005])]:
            for h in steps:
                fd=-(cases[f'{kind}-plus-{h}']['energy_Ha']-cases[f'{kind}-minus-{h}']['energy_Ha'])/(2*h)
                if kind=='force': force_errors.append(abs(fd-cases['central']['forces_Ha_bohr'][1][0]))
                else: stress_errors.append(abs(fd*4359.7447222071/512-cases['central']['stress_xx_GPa'])*1000)
    assert count==45
    assert max(force_errors)<2.857e-6 and max(stress_errors)<.473
    padding=Path('diagnostics/molecular-padding-convergence')
    manifest=read(padding/'manifest.json')
    for variant in manifest['variants']:
        for path in (ROOT/padding/variant).glob('*/output.out'): output(path.relative_to(ROOT))
    return {'cutoff_cases':36,'quadrature_cases':16,'padding_cases':9,'derivative_cases':45,
            'maximum_force_error_Ha_bohr':max(force_errors),'maximum_stress_error_MPa':max(stress_errors),
            'symmetry_pairs_kJ_mol_cell':symmetry}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-only',action='store_true',help='Skip the immutable file-inventory check')
    args=parser.parse_args()
    report={}
    if not args.data_only: report['hashed_files']=inventory()
    report.update(molecular=molecular(),LC10=solid_eos(),molecular_crystals=crystals(),numerical_checks=numerical_controls())
    for r in read('paper/pccp/execution-provenance.json'):
        for name,digest in r['files'].items():
            assert sha(Path(r['target'])/name)==digest, (r['target'],name)
    assert sum(1 for p in ROOT.rglob('output.out') if '.git' not in p.parts)==835
    for name,digest in read('paper/pccp/overleaf-source-snapshot.json')['files'].items():
        assert sha(Path('paper/pccp')/name)==digest
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
