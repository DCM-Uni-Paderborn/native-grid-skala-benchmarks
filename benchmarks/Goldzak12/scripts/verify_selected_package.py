"""Check archived selected files against their manifest and execution hashes."""
import hashlib
import json
import re

from fit_eos import ROOT


def main():
    records = json.loads((ROOT / 'results/eos-selected-input-provenance.json').read_text())
    assert len(records) == len({r['key'] for r in records}) == 352
    checked = 0
    for row in records:
        root = ROOT / 'accepted-inputs' / row['key']
        for name, expected in row['files_sha256'].items():
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, str(root / name)
            checked += 1
        fields = dict(line.split('=', 1) for line in (root / 'execution.txt').read_text().splitlines() if '=' in line)
        assert fields.get('exit_status') == '0', row['key']
        assert row['files_sha256']['input.inp'] == fields['frozen_input_sha256'], row['key']
        actual = fields.get('actual_input', 'actual-input.inp')
        assert row['files_sha256'][actual] == fields['actual_input_sha256'], row['key']
        time = (root / 'time.txt').read_text()
        assert re.search(r'Exit status:\s*0\b', time), row['key']
    print(json.dumps({'selected_packages': len(records), 'verified_files': checked,
                      'execution_time_exits_zero': True, 'input_hashes_match_execution': True}))


if __name__ == '__main__':
    main()
