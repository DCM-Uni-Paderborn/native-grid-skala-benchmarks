# Native-grid Skala: PCCP Supporting Data

Inputs, outputs, reference extracts and offline analysis for **Native implementation
of the machine-learned Skala exchange-correlation functional in CP2K: Unified
one-centre reconstruction for molecular and condensed-phase calculations**.

| Reported comparison | Data |
| --- | --- |
| 62 dietGMTKN55 reactions, three native variants | [Molecular benchmark](benchmarks/dietGMTKN55), 420 selected executions for 140 species |
| CO2, NH3 and urea, four representations | [Molecular crystals](benchmarks/X23-mini), including basis/grid controls and DFT, MP2 and multi-level CCSD(T) literature comparisons |
| All thirteen DMC-ICE13 phases | [Ice](benchmarks/DMC-ICE13), TZVPP-ae and QZVPP-ae crystal/molecule energies and relative energies to Ih |
| LC10 lattice constants and bulk moduli | [Solids](benchmarks/Goldzak12), 352 energies, 36 independent EOS curves and 40 representation/material entries |
| 50 primary band gaps across 28 materials | [Band gaps](benchmarks/band-gaps), sampled edges, reference comparisons and seven direct/one-centre controls |
| Cutoff, atom quadrature, adjoint, derivatives and symmetry | [Numerical controls](convergence) |
| Isolated-cell padding | [Padding controls](diagnostics/molecular-padding-convergence) |
| Current manuscript and ESI | [Compact sources](paper/pccp) |

The package contains the data used by the current paper and the evidence needed
to reproduce its analyses. Superseded comparisons and unused publication assets
remain available in Git history rather than the current tree. Original execution
records are not rewritten to conceal failures or replace missing metadata.

## Reproduce

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md). The complete check requires Python
and NumPy, but no CP2K executable, model weights, network connection or new SCF:

```sh
python3 -B scripts/audit_pccp_package.py
python3 -B scripts/reproduce_paper.py --output /tmp/native-skala-analysis --figures
```

Figure generation additionally requires Matplotlib. Output must be outside this
immutable data tree. [Data coverage](paper/pccp/data-coverage.json) maps every
manuscript/ESI table and figure to its sources. The [file manifest](paper/pccp/file-manifest.json)
records hashes for all deposited files except itself.

## External Molecular Comparison

The earlier [molecular GauXC study](https://github.com/DCM-Uni-Paderborn/Molecular-Skala-in-CP2K)
is not duplicated. Only its species-energy extract, geometry correspondence,
source hashes and matched 62-reaction comparison are included here. Its raw
inputs, outputs and full benchmark remain in that repository.

Exact native inputs retain their original machine-specific paths and full basis
identifiers. Model weights, executables and restart wavefunctions are identified
by provenance records rather than redistributed. These external files are not
needed to reanalyze the archived results.
