# Basis and Grid Controls

`index.json` links **56 executions** to their original inputs and outputs:
eight CO2/NH3 basis controls, 30 ice/monomer/symmetry calculations, and 18
Si/diamond/MgO basis/runtime controls. The four urea controls are indexed in
X23-mini and are not duplicated here.

`scripts/basis_sensitivity.py` verifies these records and reconstructs
`assessment.json`. It covers all thirteen ice phases at both basis levels,
absolute lattice energies, relative energies to Ih, and matched five-volume
solid basis tests. The current molecular basis comparison is instead part of
the complete 62-reaction dataset under `benchmarks/dietGMTKN55`.

All inputs retain their original computational identifiers. No model weights,
wavefunctions or executable binaries are distributed.
