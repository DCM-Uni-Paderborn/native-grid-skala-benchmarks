"""Hash the explicitly curated paper package, excluding campaign scratch data."""
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]


def paper_files():
    files = [ROOT / 'README.md']
    for directory in ('accepted-inputs', 'reference', 'results', 'paper'):
        files.extend(p for p in (ROOT / directory).rglob('*') if p.is_file())
    for name in ('paper-eos-selection.json', 'methods.csv', 'systems.csv', 'eos_volumes.csv'):
        files.append(ROOT / 'protocol' / name)
    for name in ('fit_eos.py', 'fit_selected_eos.py', 'analyze_selected_eos.py',
                 'compare_selected_literature.py', 'build_paper_tables.py',
                 'build_paper_manifest.py', 'import_zhang_reference.py',
                 'verify_selected_package.py', 'test_selected_eos.py',
                 'test_literature_comparison.py'):
        files.append(ROOT / 'scripts' / name)
    for path in files:
        if path.is_symlink() or not path.is_file():
            raise ValueError(path)
        if path.suffix.lower() in ('.wfn', '.pyc') or '__pycache__' in path.parts:
            raise ValueError(path)
    return sorted(set(files))


def main():
    files = paper_files()
    lines = [hashlib.sha256(p.read_bytes()).hexdigest() + '  ' +
             p.relative_to(REPO).as_posix() for p in files]
    (ROOT / 'SHA256SUMS').write_text('\n'.join(lines) + '\n')
    print(f'{len(files)} paper files; {sum(p.stat().st_size for p in files)} bytes')


if __name__ == '__main__':
    main()
