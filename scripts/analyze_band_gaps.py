"""Recompute sampled band edges and all matched comparisons, entirely offline."""
from collections import Counter
import json
import math
import re

from data_checks import ROOT, close, compare_nested, total_energy, verify_files

DATA = ROOT / 'benchmarks/band-gaps'


def band_edges(text, occupied, expected_points, reviewed_indices=()):
    """Censoring is permitted only at the explicitly reviewed indices of this run."""
    if not isinstance(occupied, int) or occupied < 1:
        raise ValueError('Invalid occupied-band count')
    reviewed = set(reviewed_indices)
    if any(i < 1 or i >= occupied for i in reviewed):
        raise ValueError('Censoring cannot include a frontier or unoccupied band')
    points = []
    point = None
    overflow = set()
    for line in text.splitlines():
        match = re.match(r'#\s+Point\s+(\d+)\s+Spin\s+(\d+):\s+([^#]+)', line)
        if match:
            if point is not None:
                points.append(point)
            k = [float(x) for x in match[3].split()[:3]]
            if len(k) != 3 or not all(math.isfinite(x) for x in k) or int(match[2]) != 1:
                raise ValueError('Invalid k point or spin')
            point = {'k': k, 'bands': []}
        elif point is not None and line.strip() and not line.lstrip().startswith('#'):
            number, token, occupation = line.split()
            number, occupation = int(number), float(occupation)
            close(occupation, 2 if number <= occupied else 0, 1e-6)
            if token == '*' * 14:
                if number not in reviewed:
                    raise ValueError('Unreviewed censored band')
                energy = None
                overflow.add(number)
            else:
                energy = float(token)
                if not math.isfinite(energy):
                    raise ValueError('Nonfinite eigenvalue')
            point['bands'].append((number, energy))
    if point is not None:
        points.append(point)
    if not points or len(points) != expected_points or overflow != reviewed:
        raise ValueError('Sampling or reviewed censoring differs from the record')
    count = len(points[0]['bands'])
    for p in points:
        bands = p['bands']
        if count <= occupied or [b[0] for b in bands] != list(range(1, count + 1)):
            raise ValueError('Invalid band count or indices')
        finite = [energy for _, energy in bands if energy is not None]
        if any(a > b + 1e-7 for a, b in zip(finite, finite[1:])):
            raise ValueError('Unordered eigenvalues')
        p['vbm'], p['cbm'] = bands[occupied - 1][1], bands[occupied][1]
        if p['vbm'] is None or p['cbm'] is None:
            raise ValueError('Unreadable frontier')
    v = max(points, key=lambda p: p['vbm'])
    c = min(points, key=lambda p: p['cbm'])
    direct = min(points, key=lambda p: p['cbm'] - p['vbm'])
    return {'sampled_gap_eV': c['cbm'] - v['vbm'],
            'sampled_direct_gap_eV': direct['cbm'] - direct['vbm'],
            'vbm_eV': v['vbm'], 'cbm_eV': c['cbm'], 'vbm_k': v['k'], 'cbm_k': c['k'],
            'sampled_direct_k': direct['k'], 'points': len(points), 'occupied_bands': occupied,
            'censored_core_band_indices': sorted(overflow), 'all_energies_readable': not overflow}


def gap_type(minimum, direct):
    if not math.isfinite(minimum) or not math.isfinite(direct) or direct < minimum - 1e-8:
        raise ValueError('Invalid minimum/direct gaps')
    difference = direct - minimum
    return 'D' if difference < 1e-6 else 'N' if difference < 1e-3 else 'I'


def statistics(values, references, materials, override=None):
    errors = {m: values[m] - (override or {}).get(m, references[m]['experiment']) for m in materials}
    if not errors or not all(math.isfinite(x) for x in errors.values()):
        raise ValueError('Empty or nonfinite comparison')
    worst = max(errors, key=lambda m: abs(errors[m]))
    return {'n': len(errors), 'ME_eV': sum(errors.values()) / len(errors),
            'MAE_eV': sum(abs(e) for e in errors.values()) / len(errors),
            'RMSE_eV': math.sqrt(sum(e*e for e in errors.values()) / len(errors)),
            'max_abs_error_eV': abs(errors[worst]), 'max_error_material': worst,
            'signed_errors_eV': errors}


def comparisons(rows, references):
    values = {method: {r['material']: r['gap_eV'] for r in rows if r['method'] == method}
              for method in ('AE', 'GTH', 'mixed_1C')}
    assert not values['AE'].keys() & values['mixed_1C'].keys()
    values['AE_or_mixed_1C'] = {**values['AE'], **values['mixed_1C']}

    def compare(materials, source, methods=('AE', 'GTH'), override=None):
        materials = sorted(materials)
        ref = references[source]
        stats = {'Skala_' + m: statistics(values[m], ref['records'], materials, override) for m in methods}
        for f in ref['functionals']:
            stats[f] = statistics({m: ref['records'][m][f] for m in materials}, ref['records'], materials, override)
        return {'materials': materials, 'n': len(materials), 'statistics': stats}

    local, hybrid = 'Lee2021', 'Lee2022_truncated_Coulomb'
    pair = values['AE'].keys() & values['GTH'].keys()
    hpair = pair & references[hybrid]['records'].keys()
    mixed = values['mixed_1C'].keys() & values['GTH'].keys()
    combined = values['AE_or_mixed_1C'].keys() & values['GTH'].keys()
    hcombined = combined & references[hybrid]['records'].keys()
    for material in hcombined:
        close(references[local]['records'][material]['experiment'], references[hybrid]['records'][material]['experiment'])
    methods = ('AE_or_mixed_1C', 'GTH')
    result = {
        'paired_AE_GTH_Lee2021': compare(pair, local),
        'paired_AE_GTH_Lee2022': compare(hpair, hybrid),
        'paired_MgO_7_83_sensitivity': compare(pair, local, override={'MgO': 7.83}),
        'paired_mixed_GTH': compare(mixed, local, ('mixed_1C', 'GTH')),
        'paired_AE_or_mixed_GTH': compare(combined, local, methods),
        'paired_AE_or_mixed_GTH_Lee2022': compare(hcombined, hybrid, methods),
        'Lee2021_on_paired_Lee2022_population': compare(hcombined, local, methods),
        'paired_combined_MgO_7_83_sensitivity': compare(combined, local, methods, {'MgO': 7.83}),
        'paired_AE_GTH_fully_readable_bands': compare(pair & {r['material'] for r in rows
            if r['method'] == 'AE' and not r['core_output_censored']}, local),
    }
    for method in values:
        material_set = values[method].keys()
        common = material_set & references[hybrid]['records'].keys()
        if method != 'AE_or_mixed_1C':
            result['all_available_' + method] = compare(material_set, local, (method,))
        result[f'maximum_available_{method}_Lee2021'] = compare(material_set, local, (method,))
        result[f'maximum_available_{method}_Lee2022'] = compare(common, hybrid, (method,))
        result[f'Lee2021_on_{method}_Lee2022_population'] = compare(common, local, (method,))
    return result, values


def analyze():
    data = json.loads((DATA / 'dataset.json').read_text())
    expected = json.loads((DATA / 'paper-comparisons.json').read_text())
    edges = {}
    for case, record in data['executions'].items():
        folder = verify_files(record)
        total_energy(record)
        old = record['band_edges']
        parsed = band_edges((folder / 'bands.bs').read_text(), old['occupied_bands'], old['points'],
                            record['reviewed_censored_indices'])
        for key in parsed.keys() & old.keys():
            compare_nested(parsed[key], old[key], key)
        execution = json.loads((folder / 'execution.json').read_text())
        assert execution['execution_exit'] == execution['time_exit'] == 0
        close(execution['case']['expected_electrons'], 2 * parsed['occupied_bands'])
        edges[case] = parsed
    rows = data['primary']
    assert len(rows) == len({(r['material'], r['method']) for r in rows}) == 50
    assert len({r['material'] for r in rows}) == 28
    assert Counter(r['method'] for r in rows) == {'AE': 20, 'GTH': 23, 'mixed_1C': 7}
    assert sum(r['core_output_censored'] for r in rows) == 11
    for row in rows:
        close(edges[row['case']]['sampled_gap_eV'], row['gap_eV'], 1e-8)
        assert row['core_output_censored'] == bool(edges[row['case']]['censored_core_band_indices'])
    calculated, values = comparisons(rows, data['references'])
    compare_nested(calculated, expected)
    sensitivity = data['reference_sensitivities']
    alternative = sensitivity['alternative_reference']
    for key, comparison in sensitivity['comparisons'].items():
        refs = data['references'][comparison['source']]['records']
        assert comparison['materials'] == calculated[key]['materials']
        for method, saved in comparison['statistics'].items():
            vals = values[method.removeprefix('Skala_')] if method.startswith('Skala_') else {
                m: refs[m][method] for m in comparison['materials']}
            result = statistics(vals, refs, comparison['materials'], {alternative['material']: alternative['gap_eV']})
            for metric, value in saved['alternative_reference'].items():
                close(result[metric], value)
    pairs = []
    used = {r['case'] for r in rows}
    for pair in data['direct_one_center_pairs']:
        used.update((pair['direct'], pair['one_center']))
        difference = edges[pair['one_center']]['sampled_gap_eV'] - edges[pair['direct']]['sampled_gap_eV']
        pairs.append({**pair, 'one_center_minus_direct_eV': difference})
    assert used == edges.keys() and len(pairs) == 7
    classification = Counter(gap_type(edges[r['case']]['sampled_gap_eV'],
                                      edges[r['case']]['sampled_direct_gap_eV']) for r in rows)
    assert classification == {'D': 25, 'I': 24, 'N': 1}
    assert max(abs(p['one_center_minus_direct_eV']) for p in pairs) < 0.009
    return {'primary_gaps': 50, 'materials': 28, 'fully_readable_band_outputs': 39,
            'reviewed_censored_core_outputs': 11, 'classification': dict(classification),
            'edges': edges, 'comparisons': calculated, 'direct_one_center_pairs': pairs}


if __name__ == '__main__':
    result = analyze()
    print(json.dumps({k: v for k, v in result.items() if k not in ('edges', 'comparisons')}, indent=2))
