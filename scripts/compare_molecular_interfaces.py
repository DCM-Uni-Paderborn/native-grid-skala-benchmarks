#!/usr/bin/env python3
"""Reproduce the reaction-matched native/GauXC/PySCF comparison without SCF runs.

The default calculation uses the archived, selected collaborator values. To
recreate that extract, supply --source-repo pointing to Molecular-Skala-in-CP2K
(requires numpy and openpyxl). No source workbook or production result is edited.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'benchmarks/dietGMTKN55/production-p25/paper-common-70'
NATIVE = ['gapwxc_gth', 'hybrid_ae_gth_direct', 'hybrid_ae_gth_one_center']
EXTERNAL = ['gauxc_ae', 'pyscf_unit', 'gauxc_gpw']
SHEETS = {'gauxc_ae': 'GAPW', 'gauxc_gpw': 'GPW',
          'pyscf_unit': 'PySCF unit'}
EH_KCAL_NATIVE = 627.5094740631
EH_KCAL_SOURCE = 627.509473777537


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def stats(values):
    return {'n': len(values), 'mean_signed': statistics.mean(values),
            'mean_absolute': statistics.mean(map(abs, values)),
            'median_absolute': statistics.median(map(abs, values)),
            'root_mean_square': statistics.mean(v*v for v in values)**0.5,
            'maximum_absolute': max(map(abs, values))}


def extract(source_repo):
    import numpy as np
    import openpyxl
    source = source_repo / 'raw/dietgmtkn55/validation'
    workbook = openpyxl.load_workbook(source/'energies.xlsx', data_only=True)
    index = json.loads((DATA/'reaction-index.json').read_text())
    reference = json.loads((DATA/'reference.json').read_text())

    def records(sheet):
        values = list(workbook[sheet].values)
        return [dict(zip(values[0], row)) for row in values[1:]
                if isinstance(row[0], str) and row[0].isdigit()]

    molecules = {}
    for key, sheet in SHEETS.items():
        molecules[key] = collections.defaultdict(list)
        for row in records(sheet+' molecular'):
            molecules[key][row['Reaction']].append(row)
    source_reactions = {r['Reaction']: r for r in records('Comparison')}
    geometries = {}
    for rows in molecules['pyscf_unit'].values():
        for row in rows:
            name = row['Input file']
            lines = (source/name).read_text().splitlines()
            atoms = [line.split() for line in lines[2:] if line.strip()]
            assert len(atoms) == int(lines[0])
            geometries[name] = ([a[0] for a in atoms],
                np.array([[float(x) for x in a[1:4]] for a in atoms]),
                *map(int, lines[1].split()[:2]))

    def difference(mol, ref):
        elements, xyz, charge, multiplicity = geometries[mol['Input file']]
        if (elements != ref['Elements'] or charge != ref['Charge']
                or multiplicity != ref['UHF']+1 or mol['Coefficient'] != ref['Count']):
            return float('inf')
        original = np.array(ref['Positions'])
        return float(np.max(abs(xyz-xyz.mean(0) - (original-original.mean(0)))))

    selected = []
    for reaction in index['reactions']:
        ref = reference[reaction['subset']][str(reaction['reaction_id'])]
        matches = []
        for source_id, mols in molecules['pyscf_unit'].items():
            if len(mols) != len(ref['Species']):
                continue
            assignments = [[i for i, mol in enumerate(mols) if difference(mol, spec) < 5e-5]
                           for spec in ref['Species'].values()]
            if all(len(a) == 1 for a in assignments) and len({a[0] for a in assignments}) == len(mols):
                matches.append((source_id, [a[0] for a in assignments]))
        assert len(matches) == 1, (reaction['subset'], reaction['reaction_id'], matches)
        source_id, assignment = matches[0]
        sr = source_reactions[source_id]
        assert abs(sr['Reference reaction energy (kcal/mol)']-reaction['reference_kcal_mol']) < 1e-8
        assert abs(sr['dietGMTKN55 weight']-reaction['weight']) < 1e-8
        record = {'subset': reaction['subset'], 'reaction_id': reaction['reaction_id'],
                  'source_reaction_id': source_id, 'source_subset': sr['Subset'],
                  'reference_kcal_mol': reaction['reference_kcal_mol'], 'weight': reaction['weight'],
                  'species': []}
        for (name, spec), i in zip(ref['Species'].items(), assignment):
            py = molecules['pyscf_unit'][source_id][i]
            item = {'native_species_name': name, 'coefficient': spec['Count'],
                    'maximum_centered_coordinate_difference_angstrom': difference(py, spec),
                    'methods': {}}
            for key in EXTERNAL:
                mol = molecules[key][source_id][i]
                assert mol['Coefficient'] == spec['Count'] and mol['Molecule'] == py['Molecule']
                item['methods'][key] = {
                    'converged': mol.get('Calculation status', 'converged') == 'converged',
                    'electronic_energy_ha': mol['Electronic energy (Eh)'],
                    'dispersion_energy_ha': mol['D3(BJ) correction (Eh)'],
                    'total_energy_ha': mol['Total energy (Eh)'],
                    'source_input': 'raw/dietgmtkn55/validation/'+mol['Input file'],
                    'source_output': 'raw/dietgmtkn55/validation/'+mol['Output file'],
                    'input_sha256': sha(source/mol['Input file']),
                    'output_sha256': sha(source/mol['Output file'])}
            record['species'].append(item)
        selected.append(record)
    assert len(selected) == 70
    return {'source_repository': 'https://github.com/DCM-Uni-Paderborn/Molecular-Skala-in-CP2K',
            'source_revision': subprocess.check_output(['git','rev-parse','HEAD'],cwd=source_repo,text=True).strip(),
            'source_workbook': 'raw/dietgmtkn55/validation/energies.xlsx',
            'source_workbook_sha256': sha(source/'energies.xlsx'),
            'reference_sha256': sha(DATA/'reference.json'),
            'coordinate_match_tolerance_angstrom': 5e-5,
            'matching': 'ordered atomic geometry after translation, charge, multiplicity, and signed stoichiometric coefficients; references and weights checked independently',
            'reactions': selected}


def calculate(source):
    assert sha(DATA/'reference.json') == source['reference_sha256']
    native = {(r['subset'],int(r['reaction_id'])): r for r in read_rows(DATA/'native-protocol-comparison.csv')}
    index = {(r['subset'],r['reaction_id']): r for r in json.loads((DATA/'reaction-index.json').read_text())['reactions']}
    reference = json.loads((DATA/'reference.json').read_text())
    species = {r['route']+'/'+r['digest']: r for r in read_rows(DATA/'species-results.csv')}
    result = []
    for src in source['reactions']:
        key = src['subset'],src['reaction_id']
        row = native[key]
        assert abs(float(row['reference_kcal_mol'])-src['reference_kcal_mol']) < 1e-8
        r = {k: src[k] for k in ['subset','reaction_id','source_reaction_id','source_subset','reference_kcal_mol','weight']}
        for method in NATIVE:
            r[method] = float(row[method+'_kcal_mol'])
            d3 = 0.
            for s in index[key]['species']:
                value = species[s['accepted_results'][method]['result_key']]['dispersion_energy_ha']
                if value == '':
                    assert reference[key[0]][str(key[1])]['Species'][s['name']]['Number'] == 1
                    value = '0'
                d3 += s['count'] * float(value)
            r[method+'_d3'] = d3 * EH_KCAL_NATIVE
        for method in EXTERNAL:
            complete = all(s['methods'][method]['converged'] for s in src['species'])
            for suffix, field in [('', 'total_energy_ha'), ('_d3', 'dispersion_energy_ha')]:
                r[method+suffix] = (sum(s['coefficient'] * s['methods'][method][field] for s in src['species'])
                                    * EH_KCAL_SOURCE if complete else None)
        result.append(r)
    assert len(result) == len({(r['subset'],r['reaction_id']) for r in result}) == 70
    reference_stats = {}
    for method in NATIVE+EXTERNAL:
        selected = [r for r in result if r[method] is not None]
        reference_stats[method] = stats([r[method]-r['reference_kcal_mol'] for r in selected])
    pairs = {}
    for a,b in [('hybrid_ae_gth_direct','gauxc_ae'),('hybrid_ae_gth_one_center','gauxc_ae'),
                ('hybrid_ae_gth_direct','pyscf_unit'),('gauxc_ae','pyscf_unit'),('gapwxc_gth','gauxc_gpw')]:
        selected = [r for r in result if r[a] is not None and r[b] is not None]
        pairs[a+'_minus_'+b] = {
            'as_reported': stats([r[a]-r[b] for r in selected]),
            'electronic_only': stats([r[a]-r[a+'_d3']-r[b]+r[b+'_d3'] for r in selected]),
            'dispersion_difference': stats([r[a+'_d3']-r[b+'_d3'] for r in selected])}
    assert reference_stats['gauxc_gpw']['n'] == 65
    return result, {'energy_unit': 'kcal/mol', 'error_convention':'calculated minus reference',
                    'source_revision':source['source_revision'],
                    'native_energies': 'unchanged recorded total energies including D3',
                    'external_energies':'source workbook totals including separately tabulated D3',
                    'reference_errors':reference_stats, 'pairwise_differences':pairs,
                    'native_65_reference_errors': {m:stats([r[m]-r['reference_kcal_mol'] for r in result if r['gauxc_gpw'] is not None]) for m in NATIVE},
                    'd3_discrepancies_gauxc_ae_minus_pyscf_unit':[
                        {'subset':r['subset'], 'reaction_id':r['reaction_id'],
                         'difference_kcal_mol':r['gauxc_ae_d3']-r['pyscf_unit_d3']}
                        for r in result if abs(r['gauxc_ae_d3']-r['pyscf_unit_d3'])>1e-7]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-repo',type=Path)
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    source = extract(args.source_repo) if args.source_repo else json.loads((DATA/'gauxc-source-common-70.json').read_text())
    result, summary = calculate(source)
    if args.write:
        if args.source_repo:
            (DATA/'gauxc-source-common-70.json').write_text(json.dumps(source,indent=2)+'\n')
        (DATA/'molecular-interface-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        with (DATA/'molecular-interface-comparison.csv').open('w',newline='') as handle:
            writer = csv.DictWriter(handle,fieldnames=list(result[0]),lineterminator='\n')
            writer.writeheader();writer.writerows(result)
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
