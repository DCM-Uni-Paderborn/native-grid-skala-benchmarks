# CO2, NH3 and Urea: PCCP Molecular Crystals

Exactly 24 base executions (four representations, three crystal/molecule
pairs) and seven accepted ESI controls are included.
`results/lattice-energies.json` gives the 12 paired lattice energies,
DMC deviations and numerical sensitivities.
Every accepted record includes frozen/actual input, output and completion
provenance. Failed and unfinished controls are excluded.

`scripts/analyze_results.py` uses
`E_latt = (E_crystal/Z - E_molecule) * 2625.4996394799` in kJ/mol.
Only the paired 200/974 quadrature check changes both phases. The 1200 Ry
check is crystal-only and does not establish lattice-energy cutoff
convergence. The PCCP ESI states the numerical and basis/BSSE limitations.
