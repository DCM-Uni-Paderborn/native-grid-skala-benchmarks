# CO2, NH3 and Urea: PCCP Molecular Crystals

Exactly 24 original base executions (four representations, three
crystal/molecule pairs), seven original ESI controls and four paired Urea AE
basis/grid controls are included. `results/lattice-energies.json` retains the
12 original `complete_pairs` and the 12 selected `paper_pairs` separately.
Only Urea AE uses the larger QZVPP basis with 200/974 quadrature in the selected
comparison. CO2 and NH3 retain their original AE settings.
Every accepted record includes frozen/actual input, output and completion
provenance. Failed and unfinished controls are excluded.

`scripts/analyze_results.py` uses
`E_latt = (E_crystal/Z - E_molecule) * 2625.4996394799` in kJ/mol.
The paired 200/974 CO2 and Urea quadrature checks change both phases. The 1200 Ry
check is crystal-only and does not establish lattice-energy cutoff
convergence. The PCCP ESI states the numerical and basis/BSSE limitations.

`scripts/analyze_basis_controls.py` verifies the four Urea inputs, outputs and
execution hashes, then reconstructs the grid shift (-0.2053 kJ/mol) and the
matched TZVPP-to-QZVPP basis shift (+14.4154 kJ/mol). The selected Urea AE
lattice energy is -110.4275 kJ/mol. The original values remain in the ESI basis
table and archive. No complete-basis extrapolation or counterpoise correction
is implied. Pending follow-up calculations are not included in this package.
