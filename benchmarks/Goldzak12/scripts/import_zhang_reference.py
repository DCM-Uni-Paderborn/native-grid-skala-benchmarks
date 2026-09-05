"""Extract the selected structural data from Zhang et al. SI, Tables III/V.

Optional reference refresh only; analysis uses the checked-in CSV offline.
Requires pypdf. The PDF is read in memory, not retained in the benchmark.
"""
import hashlib
import io
import json
import urllib.request
from pathlib import Path

from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://pure.mpg.de/rest/items/item_2599505_10/component/file_2618826/content'
METHODS = ('LDA', 'PBE', 'PBEsol', 'M06-L', 'SCAN', 'HSE06')
SOLIDS = ('AlP', 'BN', 'BP', 'C', 'LiCl', 'MgO', 'MgS', 'Si', 'SiC')


def extract_page(text):
    found = {}
    for line in text.splitlines():
        fields = line.split()
        if fields and fields[0] in SOLIDS and len(fields) == 14:
            if fields[0] in found:
                raise ValueError('Duplicate reference row: ' + fields[0])
            found[fields[0]] = [float(x) for x in fields[1:]]
    if set(found) != set(SOLIDS):
        raise ValueError('Incomplete PDF table: ' + str(sorted(found)))
    return found


def main():
    import csv
    data = urllib.request.urlopen(URL, timeout=30).read()
    pdf = PdfReader(io.BytesIO(data))
    a = extract_page(pdf.pages[3].extract_text())
    b = extract_page(pdf.pages[5].extract_text())
    assert a['C'][:4] == [3.532, 3.545, 3.572, 3.586]
    assert b['MgO'][:4] == [171.7, 166.9, 149.4, 145.2]
    rows = []
    for solid in SOLIDS:
        for i, method in enumerate(METHODS):
            rows.append(dict(source='Zhang2018', method=method, solid=solid,
                             a_A=a[solid][2*i], B0_GPa=b[solid][2*i],
                             a_ZPE_included_A=a[solid][2*i+1],
                             B0_ZPE_included_GPa=b[solid][2*i+1],
                             experimental_with_ZPE_a_A=a[solid][-1],
                             experimental_with_ZPE_B0_GPa=b[solid][-1]))
    with (ROOT / 'reference/zhang2018_selected.csv').open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    provenance = dict(url=URL, doi='10.1088/1367-2630/aac7f0',
                      pdf_sha256=hashlib.sha256(data).hexdigest(),
                      tables=['III', 'V'], rows=len(rows),
                      column_policy='Uncorr. is the static electronic result; Corr. includes ZPE. Compare Uncorr. with a consistently ZPE-subtracted experimental reference.',
                      scan_protocol='Non-self-consistent SCAN on PBE orbitals/density, as stated on SI page 2.')
    (ROOT / 'reference/zhang2018_selected-provenance.json').write_text(
        json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance))


if __name__ == '__main__':
    main()
