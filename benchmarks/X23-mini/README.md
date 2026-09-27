# CO2, NH3 and Urea Lattice Energies

Four native representations are compared with DMC. The archive contains
24 original crystal/molecule executions, seven numerical controls, and
twelve AE basis/grid controls. All three selected AE lattice energies use
QZVPP-ae with 200/974 quadrature.

`results/lattice-energies.json` retains the original `complete_pairs` and
the selected `paper_pairs`. The lattice-energy convention is
`(E_crystal/Z - E_molecule) * 2625.4996394799`, in kJ/mol.
Inputs and outputs preserve the exact geometry and computational settings.

`scripts/analyze_basis_controls.py` and the shared
`scripts/basis_sensitivity.py` verify the basis/grid controls and reconstruct
the selected results. The 1200 Ry numerical control changes the crystal only;
it is not a crystal/molecule lattice-energy cutoff series.

`reference/dft-comparison.json` contains the reported DFT literature values.
`reference/wavefunction-comparison.json` contains MP2/CBS, SCS(MI)-MP2/CBS
and the multi-level LNO-CCSD(T)/RPA+ph/HF comparison, with source locations,
energy definitions and geometry protocols. None is a native rerun of those
methods. The common DMC values are unchanged.
