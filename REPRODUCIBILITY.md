# Reproducing the Paper

Run from the repository root with Python 3.10 or newer and NumPy.
Matplotlib is needed only for figures. The commands below do not run CP2K,
connect to clusters or change the archived data.

```sh
python3 -B scripts/audit_pccp_package.py
python3 -B -m unittest discover -s scripts -p 'test_*.py'
python3 -B -m unittest discover -s benchmarks/Goldzak12/scripts -p 'test_*.py'
python3 -B -m unittest discover -s benchmarks/X23-mini/scripts -p 'test_*.py'
python3 -B benchmarks/Goldzak12/scripts/verify_selected_package.py
python3 -B scripts/reproduce_paper.py --output /tmp/native-skala-analysis --figures
```

The audit verifies the complete file inventory and hashes, selected CP2K
completion records, total energies, molecular stoichiometry, lattice energies,
band occupations and sampled frontier energies. It recomputes the exact
matched comparison populations and checks the paper source graph. Do not use
Python's optimized mode (`-O`), which disables validation assertions.

The reproduction command writes JSON, CSV and four figure reconstructions
outside the repository. The deposited publication figures remain the accepted
assets. PDF metadata, fonts and insignificant fit digits can vary between
library versions. Fit-only tolerances are 2e-7 angstrom for lattice constants
and 1e-4 GPa for bulk moduli. They do not relax total-energy, file-hash or
band-edge checks.

## Numerical Sources

- **Molecules:** `benchmarks/dietGMTKN55/dataset.json` fixes the 62 reactions,
  stoichiometry, references and 420 selected native executions. The external
  comparison is a small pinned extract, not a copy of the GauXC dataset.
  Predecessors relevant to the selected molecular solutions are retained in
  `provenance`, outside numerical statistics.
- **Molecular crystals:** `results/lattice-energies.json` retains original and
  selected crystal/molecule energies. `analyze_basis_controls.py` reconstructs
  the QZVPP-ae selections and the basis/quadrature changes. Reference extracts
  identify the DFT and correlated methods and their different protocols.
- **Ice and basis controls:** `convergence/basis-sensitivity-20260915/index.json`
  links 56 executions covering CO2/NH3, all thirteen ice phases and Si/C/MgO.
  The four urea controls have their own X23 index. Relative ice energies cancel
  the monomer reference; the absolute-energy monomer geometry is documented.
- **LC10:** `protocol/paper-eos-selection.json` fixes materials, volumes and
  representation reuse. `results/eos-selected-input-provenance.json` identifies
  the 352 original outputs. The selected EOS and literature scripts reproduce
  structural values and same-population reference errors.
- **Band gaps:** `benchmarks/band-gaps/dataset.json` maps 50 primary gaps and
  seven direct/one-centre comparisons to 61 distinct executions. The parser
  reconstructs minimum and direct gaps from the sampled bands. Deep-core
  censoring is accepted only at explicitly reviewed indices, never at a band
  edge. Reference sensitivities remain separate from the baseline values.
- **Numerical controls:** the audit recomputes the AE and GTH cutoff results,
  quadrature and symmetry differences, and periodic-water finite differences.
  Adjoint kernel sources and archived checks are in `convergence/adjoint`.

## Publication Sources

`paper/pccp/main.tex` and `supplementary_information.tex` contain the entire
texts, including tables. They share `references.bib` and `figures/`.
`source-snapshot.json` identifies their Overleaf version and hashes.
`data-coverage.json` maps all tables and figures to their numerical evidence.

Compile in a separate copy with LaTeX, latexmk and the RSC bibliography style.
The supplied `latexmkrc` refreshes both documents' cross-references. Build
products are not part of the data manifest.

## Electronic-Structure Reruns

Original input paths must be mapped to an appropriate CP2K installation.
Obtain bases, pseudopotentials and Skala model weights from their cited
distributions and match the recorded hashes and scientific settings.
No restart WFN or binary is distributed, so exact restart trajectories are
not reproduced by the offline analysis. No new electronic-structure run is
needed for any command above.
