#!/usr/bin/env python3
"""Small self-consistent native-field derivative checks; no production inputs edited."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import time

BOHR_A = 0.529177210903
HA_A3_GPA = 4359.7447222071
ROOT = Path(__file__).resolve().parent
STACK = Path('/home/kuehne88/work/native-skala-terok')
SOURCE = STACK / 'x23-mini/runtime/source-37aedacbb3'
BINARY = STACK / 'x23-mini/runtime/build-37aedacbb3/bin/cp2k.psmp'
MODEL = STACK / 'pilot-source/models/skala-1.1-rev1.fun'
EXPECTED_BINARY = '4c4cdcf4ac51510d1d876310efa724ae7efa8bdab7fdd124d803e0220fb770fb'
EXPECTED_MODEL = '7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            digest.update(chunk)
    return digest.hexdigest()


def make_input(method, displacement=0.0, strain=0.0, restart=False):
    ae = method == 'gapw-ae'
    qs = 'GPW' if method == 'gpw-gth' else ('GAPW_XC' if method == 'gapwxc-gth' else 'GAPW')
    representation = 'DIRECT_VALENCE' if method in ('gpw-gth', 'gapw-gth-direct') else 'PAW_ONE_CENTER'
    coords = [['O', 4., 4., 4.], ['H', 4.76, 4., 4.59], ['H', 3.76, 4.74, 4.59]]
    coords[1][1] += displacement * BOHR_A
    for row in coords:
        row[1] *= 1 + strain
    coord_text = '\n'.join(f'{s} {x:.13f} {y:.13f} {z:.13f}' for s, x, y, z in coords)
    kinds = ''
    for element, charge, radius in [('O', 6, 1.4), ('H', 1, 1.0)]:
        basis = 'TZVPP-MOLOPT-PBE-ae' if ae else f'TZV2P-MOLOPT-PBE-GTH-q{charge}'
        potential = 'ALL' if ae else f'GTH-PBE-q{charge}'
        kinds += f'&KIND {element}\n BASIS_SET {basis}\n POTENTIAL {potential}\n RADIAL_GRID 100\n LEBEDEV_GRID 434\n HARD_EXP_RADIUS {radius}\n&END KIND\n'
    augmentation = '' if qs == 'GPW' else 'GAPW_1C_BASIS EXT_SMALL\n GAPW_ACCURATE_XCINT TRUE\n EPSFIT 1E-4\n EPSISO 1E-12\n EPSRHO0 1E-6\n LMAXN0 4'
    seed = 'WFN_RESTART_FILE_NAME seed.wfn' if restart else ''
    guess = 'RESTART' if restart else 'ATOMIC'
    return f'''&GLOBAL
 PROJECT fd_{method}
 RUN_TYPE ENERGY_FORCE
 PRINT_LEVEL MEDIUM
&END GLOBAL
&FORCE_EVAL
 METHOD QUICKSTEP
 STRESS_TENSOR ANALYTICAL
 &PRINT
  &FORCES ON
  &END FORCES
  &STRESS_TENSOR ON
  &END STRESS_TENSOR
 &END PRINT
 &SUBSYS
  &CELL
   ABC {8*(1+strain):.13f} 8.0 8.0
   PERIODIC XYZ
  &END CELL
  &COORD
{coord_text}
  &END COORD
{kinds}
 &END SUBSYS
 &DFT
  BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH_2026.2
  POTENTIAL_FILE_NAME POTENTIAL_UZH_2026.2
  CHARGE 0
  MULTIPLICITY 1
  UKS FALSE
  {seed}
  &QS
   METHOD {qs}
   EPS_DEFAULT 1E-12
   {augmentation}
  &END QS
  &MGRID
   CUTOFF 400
   REL_CUTOFF 60
   NGRIDS 5
  &END MGRID
  &POISSON
   PERIODIC XYZ
  &END POISSON
  &SCF
   SCF_GUESS {guess}
   EPS_SCF 1E-9
   MAX_SCF 200
   &OT
    MINIMIZER DIIS
    PRECONDITIONER FULL_ALL
    STEPSIZE 0.05
   &END OT
   &OUTER_SCF
    EPS_SCF 1E-9
    MAX_SCF 2
   &END OUTER_SCF
   &PRINT
    &RESTART ON
     BACKUP_COPIES 0
    &END RESTART
   &END PRINT
  &END SCF
  &XC
   &XC_FUNCTIONAL
    &GAUXC
     MODEL SKALA
     NATIVE_GRID TRUE
     NATIVE_GRID_LAYOUT ATOM_COMPOSITE
     NATIVE_GRID_ATOM_PARTITION SMOOTH
     NATIVE_GRID_DIAGNOSTICS TRUE
     NATIVE_GRID_USE_CUDA FALSE
     NATIVE_GRID_GAPW_DENSITY_PARTITION HARD_MINUS_SOFT
     PSEUDOPOTENTIAL_GAPW_REPRESENTATION {representation}
    &END GAUXC
   &END XC_FUNCTIONAL
  &END XC
 &END DFT
&END FORCE_EVAL
'''


def parse(output):
    text = output.read_text(errors='replace')
    start = text.rfind('PROGRAM STARTED AT')
    if start < 0:
        raise ValueError(f'Missing execution start in {output}')
    block = text[start:]
    energy = re.findall(r'ENERGY\|\s+Total FORCE_EVAL[^\n]*?\[hartree\]\s+([-+\d.Ee]+)', block)
    force_blocks = re.findall(r'ATOMIC FORCES in \[a\.u\.\](.*?)(?:SUM OF ATOMIC FORCES|ENERGY)', block, re.S)
    forces = []
    if force_blocks:
        for line in force_blocks[-1].splitlines():
            parts = line.split()
            if len(parts) == 6 and parts[0].isdigit():
                forces.append([float(x) for x in parts[-3:]])
    current_forces = re.findall(r'^\s*FORCES\|\s+\d+\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+[-+\d.Ee]+', block, re.M)
    if current_forces:
        forces = [[float(x) for x in row] for row in current_forces]
    stress_blocks = re.findall(r'STRESS\| Analytical stress tensor \[(bar|GPa)\](.*?)(?=STRESS\| Eigenvectors|\Z)', block, re.S)
    stress = []
    if stress_blocks:
        unit, tensor = stress_blocks[-1]
        stress = re.findall(r'STRESS\|\s+x\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)', tensor)
    return {'converged': bool(re.search(r'SCF run converged', block)),
            'ended': 'PROGRAM ENDED' in block, 'abort': '[ABORT]' in block,
            'energy_Ha': float(energy[-1]) if energy else None,
            'forces_Ha_bohr': forces,
            'stress_xx_GPa': float(stress[0][0]) * (1e-4 if unit == 'bar' else 1.) if stress else None,
            'tail': block[-1500:]}


def run(args):
    assert sha(BINARY) == EXPECTED_BINARY
    assert sha(MODEL) == EXPECTED_MODEL
    folder = ROOT / 'derivative-checks' / args.method
    folder.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(GAUXC_SKALA_MODEL=str(MODEL), CUDA_VISIBLE_DEVICES='-1',
               CP2K_DATA_DIR=str(SOURCE / 'data'), OMP_NUM_THREADS='8',
               OMP_DYNAMIC='FALSE', OMP_PROC_BIND='FALSE', OMP_STACKSIZE='256M',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env['PATH'] = '/home/kuehne88/micromamba/envs/cp2k-skala-gauxc/bin:' + env['PATH']
    cases = [('central', 0., 0.)]
    for step in [0.001, 0.0005]:
        cases.extend([(f'force-plus-{step}', step, 0.), (f'force-minus-{step}', -step, 0.)])
    for step in [0.0001, 0.00005]:
        cases.extend([(f'strain-plus-{step}', 0., step), (f'strain-minus-{step}', 0., -step)])
    if args.case:
        if args.case not in [case[0] for case in cases[1:]]:
            raise ValueError('Unknown perturbation case')
        cases = [case for case in cases if case[0] in ('central', args.case)]
        if not (folder / 'central' / 'time.txt').exists():
            raise RuntimeError('A completed central execution is required')
    report = {'host': socket.gethostname(), 'method': args.method, 'binary_sha256': sha(BINARY),
              'model_sha256': sha(MODEL), 'script_sha256': sha(__file__), 'cpus': args.cpus,
              'purpose': 'Derivative consistency of periodic H2O at Gamma, not a benchmark energy',
              'source_revision': '37aedacbb33677307818301a47a09ae435712293',
              'basis_sha256': sha(SOURCE / 'data/BASIS_MOLOPT_UZH_2026.2'),
              'potential_sha256': sha(SOURCE / 'data/POTENTIAL_UZH_2026.2'), 'cases': []}
    seed = None
    for name, displacement, strain in cases:
        case_dir = folder / name
        inp = case_dir / 'input.inp'
        expected_input = make_input(args.method, displacement, strain, seed is not None)
        reused = case_dir.exists()
        if reused:
            if not args.resume or inp.read_text() != expected_input:
                raise RuntimeError(f'Existing case requires matching input and --resume: {case_dir}')
            time_record = (case_dir / 'time.txt').read_text()
            if not re.search(r'Exit status:\s+0\s*$', time_record):
                raise RuntimeError(f'Existing execution not cleanly finished: {case_dir}')
            exit_code, wall_seconds = 0, None
        else:
            case_dir.mkdir(exist_ok=False)
            if seed:
                shutil.copyfile(seed, case_dir / 'seed.wfn')
            inp.write_text(expected_input)
            command = ['taskset', '-c', args.cpus, '/home/kuehne88/micromamba/envs/cp2k-skala-gauxc/bin/time', '-v', '-o', 'time.txt', str(BINARY), '-i', 'input.inp', '-o', 'output.out']
            started = time.time()
            with (case_dir / 'execution.log').open('w') as log:
                result = subprocess.run(command, cwd=case_dir, env=env, stdout=log, stderr=log, timeout=3600)
            exit_code, wall_seconds = result.returncode, time.time()-started
        rec = parse(case_dir / 'output.out')
        rec.update(name=name, displacement_bohr=displacement, strain_xx=strain, exit_code=exit_code,
                   wall_seconds=wall_seconds, reused_completed_execution=reused, input_sha256=sha(inp), output_sha256=sha(case_dir / 'output.out'),
                   execution_sha256=sha(case_dir/'execution.log'), time_sha256=sha(case_dir/'time.txt'))
        report['cases'].append(rec)
        report_path = folder / (f'report-{args.case}.json' if args.case else 'summary.json')
        report_path.write_text(json.dumps(report, indent=2)+'\n')
        print(args.method, name, rec['converged'], rec['energy_Ha'], flush=True)
        if not (rec['converged'] and rec['ended'] and not rec['abort'] and exit_code == 0 and rec['energy_Ha'] is not None):
            break
        if name == 'central':
            wfns = list(case_dir.glob('*RESTART.wfn'))
            if len(wfns) != 1 or len(rec['forces_Ha_bohr']) != 3 or rec['stress_xx_GPa'] is None:
                raise RuntimeError('Missing derivative or restart record')
            seed = wfns[0]
            report['seed_sha256'] = sha(seed)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('method', choices=['gpw-gth','gapwxc-gth','gapw-ae','gapw-gth-direct','gapw-gth-one-center'])
    parser.add_argument('--cpus', required=True)
    parser.add_argument('--resume', action='store_true', help='Reuse only exactly matching, cleanly ended executions')
    parser.add_argument('--case', help='Execute one unstarted perturbation using the exact completed central WFN')
    run(parser.parse_args())
