# Matched Solid Basis Controls

The PCCP SI reports five QZVPP volumes for each of Si, diamond (C) and MgO,
plus one same-basis TZVPP runtime check per solid. The original TZ energies
remain in `../accepted-inputs/gapw-ae` and are reused rather than duplicated.
This is not a complete ten-solid QZVPP benchmark.

Each directory contains the exact executed input, output, timing, launcher
and stable process provenance. The calculations use traditional diagonalization
with EPS_SCF 5e-6, the recorded full symmetry backend and the same
volume-specific k-mesh as the paired TZ input. MgO v07/v09 use the completed
executions with corrected monitoring counts, without altering their inputs.

From the repository root, `python3 -B scripts/basis_sensitivity.py` verifies
the evidence and reproduces the matched fits, window tests and TZ runtime
reproduction errors. Sparse five-point fits do not justify precise pressure
derivatives. No wavefunctions, binaries or unreported runs are included.
