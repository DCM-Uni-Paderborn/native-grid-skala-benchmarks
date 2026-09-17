# Paired Ice Basis Controls

This directory supports the matched twelve-phase basis comparison in the PCCP
SI, not a completed QZVPP DMC-ICE13 benchmark. It contains Ih, II, III, IV, VI,
VII, VIII, IX, XI, XIV, XV and XVII at both basis levels, the two isolated-water
references, and the same-basis spglib controls for IX and VII: 28 executions.
XIII is excluded from both averages because the QZ calculation requires longer
SCF convergence. Relative energies to Ih cancel the molecular reference.

The exact inputs include crystal coordinates, basis and quadrature settings,
and the executed symmetry choice. QZVPP uses full point-group reduction of
the Gamma-centred 3 x 3 x 3 mesh. All runs use 32 CPU threads and preserve
source, executable, model, elapsed time and peak-memory provenance.

Reference: F. Della Pia et al., J. Chem. Phys. 157, 134701 (2022),
doi:10.1063/5.0102645. The molecular reference is reconstructed equilibrium
Partridge-Schwenke water (OH 0.95865 Angstrom, HOH 104.348 degrees), not a
directly obtained DMC monomer coordinate file. No energy alignment is applied.
This convention and the incomplete phase coverage limit the absolute comparison.

Run `python3 -B scripts/basis_sensitivity.py` from the repository root to
reproduce all twelve paired energies, absolute and relative same-population
MAEs, and the two
symmetry differences. Inputs and output evidence are sufficient for the
analysis. External model weights and wavefunctions are not redistributed.
