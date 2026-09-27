"""Shared offline validation for the archived numerical evidence."""
import hashlib
import fnmatch
import math
from pathlib import Path
import re
import statistics

ROOT = Path(__file__).resolve().parents[1]


def package_files(root=ROOT):
    """Honor the package's simple artifact globs, including on case-insensitive hosts."""
    patterns = [line.strip().rstrip('/').casefold() for line in (root / '.gitignore').read_text().splitlines()
                if line.strip() and not line.startswith('#')]
    for path in root.rglob('*'):
        if not path.is_file():
            continue
        parts = path.relative_to(root).parts
        if '.git' in parts or any(fnmatch.fnmatchcase(part.casefold(), pattern)
                                 for part in parts for pattern in patterns):
            continue
        yield path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(actual, expected, tolerance=1e-9):
    if not math.isfinite(float(actual)) or not math.isfinite(float(expected)):
        raise ValueError(f'Nonfinite comparison: {actual}, {expected}')
    if abs(float(actual) - float(expected)) > tolerance:
        raise ValueError(f'{actual} != {expected}, tolerance {tolerance}')


def verify_files(record, root=ROOT):
    folder = root / record['directory']
    for name, expected in record['files_sha256'].items():
        if sha(folder / name) != expected:
            raise ValueError(f'Hash mismatch: {folder / name}')
    return folder


def total_energy(record, root=ROOT):
    folder = verify_files(record, root)
    text = (folder / 'output.out').read_text(errors='replace')
    if 'PROGRAM STARTED' not in text:
        raise ValueError(f'Missing program start: {folder}')
    block = text.rsplit('PROGRAM STARTED', 1)[1]
    if ('SCF run converged' not in block or 'PROGRAM ENDED' not in block
            or '[ABORT]' in block or 'SCF run NOT converged' in block):
        raise ValueError(f'Incomplete SCF output: {folder}')
    energies = re.findall(r'ENERGY\| Total FORCE_EVAL.*?\[hartree\]\s+([-+\d.Ee]+)', block)
    if not energies:
        raise ValueError(f'Missing total energy: {folder}')
    value = float(energies[-1])
    close(value, record['energy_hartree'], 1e-10)
    return value


def metrics(values):
    if not values or not all(math.isfinite(v) for v in values):
        raise ValueError('Empty or nonfinite population')
    return {'n': len(values), 'ME': statistics.mean(values),
            'MAE': statistics.mean(map(abs, values)),
            'RMSE': math.sqrt(statistics.mean(v*v for v in values)),
            'median_absolute': statistics.median(map(abs, values)),
            'maximum_absolute': max(map(abs, values))}


def compare_nested(actual, expected, path=''):
    """Only fitted EOS parameters receive numerical-library tolerances."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or actual.keys() != expected.keys():
            raise ValueError(f'Keys differ: {path}')
        for key in expected:
            compare_nested(actual[key], expected[key], path + '.' + key)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise ValueError(f'Lengths differ: {path}')
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare_nested(a, b, path + f'[{i}]')
    elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
        tolerance = {'a0_angstrom': 2e-7, 'B0_GPa': 1e-4,
                     'Bprime': 1e-4, 'rms_meV_atom': 1e-7}.get(path.rsplit('.', 1)[-1], 1e-9)
        close(actual, expected, tolerance)
    elif actual != expected:
        raise ValueError(f'Value differs: {path}: {actual!r} != {expected!r}')
