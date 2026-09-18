# Paired Ice Basis Controls

This directory supports the complete matched thirteen-phase basis comparison
in the PCCP manuscript and SI. It contains Ih, II, III, IV, VI, VII, VIII, IX,
XI, XIII, XIV, XV and XVII at both basis levels, the two isolated-water
references, and the same-basis spglib controls for IX and VII: 30 executions.
XIII is included at both basis levels following completion and validation of
its QZVPP calculation on 18 September 2026. Relative energies to Ih cancel
the molecular reference. Complete phase coverage does not establish a
complete-basis limit.

The absolute lattice-energy MAE is 21.349153 kJ/mol at TZVPP and 1.160740
kJ/mol at QZVPP. For the twelve nonzero differences to Ih, the corresponding
MAEs are 4.172563 and 0.821437 kJ/mol. QZVPP retains a mean signed absolute
error of +1.058003 kJ/mol and a mean relative error of -0.719912 kJ/mol.
For XIII, the QZVPP absolute error is +0.294328 kJ/mol, but its relative
error is -1.428209 kJ/mol. These distinct measures must not be conflated.

The exact inputs include crystal coordinates, basis and quadrature settings,
and the executed symmetry choice. QZVPP uses full point-group reduction of
the Gamma-centred 3 x 3 x 3 mesh. All runs use 32 CPU threads and preserve
source, executable, model, elapsed time and peak-memory provenance.

Reference: F. Della Pia et al., J. Chem. Phys. 157, 134701 (2022),
doi:10.1063/5.0102645. The molecular reference is reconstructed equilibrium
Partridge-Schwenke water (OH 0.95865 Angstrom, HOH 104.348 degrees), not a
directly obtained DMC monomer coordinate file. No energy alignment is applied.
This molecular-reference convention remains a limitation of the absolute comparison.

Run `python3 -B scripts/basis_sensitivity.py` from the repository root to
reproduce all thirteen paired energies, absolute and relative same-population
error statistics, and the two
symmetry differences. Inputs and output evidence are sufficient for the
analysis. External model weights and wavefunctions are not redistributed.
