#!/usr/bin/env python3
"""Validate the small periodic derivative dataset and generate its SI table."""
import argparse
import csv
import json
from pathlib import Path
import re

from run_derivative_checks import HA_A3_GPA, make_input, parse, sha

METHODS = ['gpw-gth', 'gapwxc-gth', 'gapw-ae', 'gapw-gth-direct', 'gapw-gth-one-center']
LABELS = ['GPW/GTH', 'GAPW-XC/GTH', 'GAPW-AE', 'GAPW-GTH direct', 'GAPW-GTH one centre']
CASES = [('central', 0., 0.)]
for h in [0.001, 0.0005]:
    CASES += [(f'force-plus-{h}', h, 0.), (f'force-minus-{h}', -h, 0.)]
for eta in [0.0001, 0.00005]:
    CASES += [(f'strain-plus-{eta}', 0., eta), (f'strain-minus-{eta}', 0., -eta)]


def collect(root):
    result = {'scope': 'Periodic H2O in an 8 Angstrom cube, Gamma only, no D3',
              'collector_sha256': sha(__file__), 'methods': {}}
    for method in METHODS:
        folder = root / method
        records = [json.loads(p.read_text()) for p in folder.glob('*.json')]
        assert records, f'Missing execution provenance for {method}'
        keys = ['host', 'binary_sha256', 'model_sha256', 'basis_sha256', 'potential_sha256', 'source_revision']
        provenance = {k: records[0][k] for k in keys}
        assert all(all(r[k] == provenance[k] for k in keys) for r in records)
        wfns = list((folder / 'central').glob('*RESTART.wfn'))
        assert len(wfns) == 1
        provenance['central_wfn_sha256'] = sha(wfns[0])
        provenance['driver_sha256_values'] = sorted({r['script_sha256'] for r in records})
        values = {}
        for name, displacement, strain in CASES:
            case = folder / name
            assert (case / 'input.inp').read_text() == make_input(method, displacement, strain, name != 'central')
            output = (case / 'output.out').read_text(errors='replace')
            assert output.count('PROGRAM STARTED AT') == 1, case
            parsed = parse(case / 'output.out')
            assert parsed['converged'] and parsed['ended'] and not parsed['abort'], case
            assert len(parsed['forces_Ha_bohr']) == 3 and parsed['stress_xx_GPa'] is not None
            assert -100 < parsed['energy_Ha'] < -1
            block = output.rsplit('PROGRAM STARTED AT', 1)[1]
            electrons = [int(v) for v in re.findall(r'Number of electrons:\s+(\d+)', block)]
            assert electrons == ([10] if method == 'gapw-ae' else [8]), case
            assert re.search(r'DFT\| Charge\s+0\s*$', block, re.M)
            assert re.search(r'DFT\| Multiplicity\s+1\s*$', block, re.M)
            assert 'Spin restricted Kohn-Sham' in block
            steps = re.findall(r'^\s*(\d+)\s+OT\s+\S+\s+\S+\s+\S+\s+([-+\d.Ee]+)', block, re.M)
            assert steps and float(steps[-1][1]) <= 1e-9
            timer = (case / 'time.txt').read_text()
            assert re.search(r'Exit status:\s+0\s*$', timer)
            execution_records = [c for r in records for c in r['cases'] if c['name'] == name]
            assert execution_records and all(c['exit_code'] == 0 for c in execution_records)
            hashes = {f: sha(case / f) for f in ['input.inp', 'output.out', 'execution.log', 'time.txt']}
            for rec in execution_records:
                assert all(rec[k] == hashes[f] for k, f in [('input_sha256', 'input.inp'), ('output_sha256', 'output.out'), ('execution_sha256', 'execution.log'), ('time_sha256', 'time.txt')])
            if name != 'central':
                assert sha(case / 'seed.wfn') == provenance['central_wfn_sha256']
            grids = [int(v) for v in re.findall(r'PW_GRID\|\s+Bounds\s+\d+[^\n]*Points:\s+(\d+)', block)]
            assert grids and len(grids) % 3 == 0
            final_wfns = list(case.glob('*RESTART.wfn'))
            assert len(final_wfns) == 1
            parsed.pop('tail')
            parsed.update(sha256=hashes, final_wfn_sha256=sha(final_wfns[0]),
                          execution_exit=0, time_exit=0, electrons=electrons[0],
                          last_step=int(steps[-1][0]), last_residual=float(steps[-1][1]),
                          grid_dimensions=[grids[i:i+3] for i in range(0, len(grids), 3)],
                          wall_time=re.search(r'Elapsed \(wall clock\) time[^\n]*?:\s+([^\n]+)', timer)[1],
                          maximum_rss_KiB=int(re.search(r'Maximum resident set size \(kbytes\):\s+(\d+)', timer)[1]))
            values[name] = parsed
        assert all(v['grid_dimensions'] == values['central']['grid_dimensions'] for v in values.values())
        result['methods'][method] = {'provenance': provenance, 'cases': values}
    return result


def analyse(report):
    rows = []
    for method in METHODS:
        values = report['methods'][method]['cases']
        central = values['central']
        for kind, steps in [('force', [0.001, 0.0005]), ('strain', [0.0001, 0.00005])]:
            for step in steps:
                plus = values[f'{kind}-plus-{step}']['energy_Ha']
                minus = values[f'{kind}-minus-{step}']['energy_Ha']
                fd = -(plus-minus)/(2*step)
                if kind == 'strain':
                    fd *= HA_A3_GPA / 512.
                    analytic = central['stress_xx_GPa']
                else:
                    analytic = central['forces_Ha_bohr'][1][0]
                rows.append(dict(method=method, derivative=kind, step=step,
                                 analytic=analytic, finite_difference=fd,
                                 signed_error=fd-analytic,
                                 units='GPa' if kind == 'strain' else 'Ha/bohr'))
    return rows


def write_outputs(report, out, tex):
    out.mkdir(parents=True, exist_ok=True)
    if (out/'cases').exists():
        for method, entry in report['methods'].items():
            for name, rec in entry['cases'].items():
                for filename, digest in rec['sha256'].items():
                    assert sha(out/'cases'/method/name/filename) == digest
    (out/'validated-results.json').write_text(json.dumps(report, indent=2)+'\n')
    rows = analyse(report)
    with (out/'derivative-errors.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    if tex:
        body = []
        for method, label in zip(METHODS, LABELS):
            selected = [r for r in rows if r['method'] == method]
            errors = [r['signed_error'] * (1e6 if r['derivative'] == 'force' else 1e3) for r in selected]
            body.append(label + ' & ' + ' & '.join(f'{v:+.3f}' for v in errors) + r' \\')
        text = r'''Self-consistent checks use one water molecule in a periodically repeated 8~\AA{} cube at Gamma, without D3, to isolate the native electronic derivatives. The grid settings are 400/60~Ry, five multigrids, and 100 radial/434 angular points. The basis/core and reconstruction definitions match the corresponding production representations. The first hydrogen is displaced along $x$ by $h=10^{-3}$ and $5\times10^{-4}$ bohr. Homogeneous $xx$ strains are $\eta=10^{-4}$ and $5\times10^{-5}$. Every perturbed calculation starts from the exact central wavefunction and converges independently to an OT residual below $10^{-9}$. All 45 executions have normal termination and zero execution/time exits. The regular-grid dimensions remain unchanged under these perturbations.

\begin{table}[htbp]
\centering\small
\caption{Signed FD minus analytic derivatives. Force errors are in $\mu E_h$/bohr and stress errors in MPa. Columns give the two displacement/strain steps.}
\label{tab:derivatives}
\begin{tabular}{lrrrr}
\toprule
Representation & $\Delta F(10^{-3})$ & $\Delta F(5\times10^{-4})$ & $\Delta\sigma(10^{-4})$ & $\Delta\sigma(5\times10^{-5})$\\
\midrule
'''
        text += '\n'.join(body) + '\n'
        text += r'''\bottomrule
\end{tabular}
\end{table}

These checks test the self-consistent total-energy derivative, including the reconstruction and moving-grid terms, rather than only model differentiation at fixed fields. The smaller steps do not uniformly reduce the error: subtraction of converged energies amplifies finite SCF and numerical noise. Consequently, no ideal quadratic step-size scaling is claimed. The test covers one Cartesian force component and one diagonal strain component in this periodic Gamma-point system. It is not a validation of every force component, shear strain, finite-k-point derivative, or symmetry-breaking displacement under point-group reduction. Exact inputs, outputs, timing records, and file/model/binary hashes accompany the machine-readable analysis.
'''
        tex.write_text(text)
    if (out/'cases').exists():
        files = sorted(p for p in out.rglob('*') if p.is_file() and p.name != 'SHA256SUMS')
        (out/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(out)}\n' for p in files))
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--collect', type=Path, help='Remote case root with exact WFN files')
    source.add_argument('--report', type=Path, help='Previously validated, hash-addressed report')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--tex', type=Path)
    args = parser.parse_args()
    report = collect(args.collect) if args.collect else json.loads(args.report.read_text())
    write_outputs(report, args.out, args.tex)
