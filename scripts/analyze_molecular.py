"""Recompute the paper's 62-reaction comparison from native outputs."""
import json
import math

from data_checks import ROOT, close, metrics, total_energy

DATA = ROOT / 'benchmarks/dietGMTKN55'
EXTERNAL_FACTOR = 627.509473777537


def analyze():
    data = json.loads((DATA / 'dataset.json').read_text())
    external = json.loads((DATA / 'external-comparison.json').read_text())
    expected = json.loads((DATA / 'paper-results.json').read_text())
    refs = {r['reaction']: r for r in expected['rows']}
    sources = {f"{r['subset']}/{r['reaction_id']}": r for r in external['reactions']}
    energies = {key: total_energy(record) for key, record in data['executions'].items()}
    assert len(energies) == 420 and len(data['reactions']) == len(refs) == len(sources) == 62
    used = set()
    rows = []
    for reaction in data['reactions']:
        label = reaction['reaction']
        source = sources[label]
        close(source['reference_kcal_mol'], reaction['reference_kcal_mol'])
        values = {}
        for method in data['methods']:
            keys = [method + '/' + s['digest'] for s in reaction['species']]
            used.update(keys)
            values[method] = math.fsum(s['count'] * energies[key]
                                      for s, key in zip(reaction['species'], keys)) * data['hartree_to_kcal_mol']
        assert {s['name'] for s in reaction['species']} == {s['native_species_name'] for s in source['species']}
        for native in reaction['species']:
            item = next(s for s in source['species'] if s['native_species_name'] == native['name'])
            assert item['coefficient'] == native['count']
            assert item['maximum_centered_coordinate_difference_angstrom'] < external['coordinate_match_tolerance_angstrom']
        for method in ('gauxc_gpw', 'gauxc_ae', 'pyscf_unit'):
            for s in source['species']:
                result = s['methods'][method]
                assert result['converged']
                close(result['electronic_energy_ha'] + result['dispersion_energy_ha'], result['total_energy_ha'])
            values[method] = math.fsum(s['coefficient'] * s['methods'][method]['total_energy_ha']
                                      for s in source['species']) * EXTERNAL_FACTOR
        for method, value in values.items():
            close(value, refs[label]['energies'][method], 1e-6)
        rows.append({'reaction': label, 'reference_kcal_mol': reaction['reference_kcal_mol'],
                     'energies_kcal_mol': values})
    assert used == energies.keys(), 'Unreferenced native execution'
    statistics = {method: metrics([r['energies_kcal_mol'][method] - r['reference_kcal_mol'] for r in rows])
                  for method in rows[0]['energies_kcal_mol']}
    pairs = {a + '_minus_' + b: metrics([r['energies_kcal_mol'][a] - r['energies_kcal_mol'][b] for r in rows])
             for a, b in [('ae-tz', 'gauxc_ae'), ('ae-qz', 'gauxc_ae'), ('gth', 'gauxc_gpw'),
                          ('gauxc_ae', 'pyscf_unit'), ('ae-qz', 'ae-tz')]}
    return {'units': 'kcal/mol', 'error_convention': 'calculation minus reference',
            'weighting': 'unweighted', 'rows': rows, 'statistics': statistics,
            'pairwise_differences': pairs, 'unique_native_executions': len(energies)}


if __name__ == '__main__':
    result = analyze()
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2))
