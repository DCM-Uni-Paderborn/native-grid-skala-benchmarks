# Periodic Native-Skala Derivative Checks

This dataset supports the finite-difference table in the PCCP supporting
information. It is a derivative-consistency test, not an additional molecular
crystal benchmark or an ice calculation.

## Scope

Five representations are tested: GPW/GTH, GAPW-XC/GTH, GAPW-AE, and ordinary
GAPW/GTH with direct and one-centre reconstructed valence fields. Each uses one
periodically repeated water molecule in an 8 Angstrom cubic cell at Gamma.
D3 is disabled to isolate the native electronic contribution. The central
geometry and every actual input are retained in `cases/`.

Each representation has nine independent self-consistent executions: the
central configuration, two signed displacements at each of two step sizes,
and two signed homogeneous strains at each of two step sizes. Displacements
move atom 2 (the first H) along x by 0.001 or 0.0005 bohr. The xx strain steps
are 0.0001 and 0.00005, with fractional coordinates fixed. Every perturbation
starts from that method's exact central WFN, whose hash is recorded.

The protocol uses 400/60 Ry, five multigrids, 100 radial/434 Lebedev points,
and OT with a 1e-9 residual threshold. This tighter threshold is specific to
numerical differentiation and does not alter any production result or setting.
The executable, model, basis, and potential hashes are in
`validated-results.json`. All cases use the same production executable at
source revision `37aedacbb33677307818301a47a09ae435712293`.

## Reproduce the Analysis

From the repository root:

```bash
python3 paper/pccp/analyse_derivative_checks.py \
  --report convergence/derivatives/periodic-water-20260909/validated-results.json \
  --out convergence/derivatives/periodic-water-20260909 \
  --tex paper/pccp/derivative-results.tex
```

The validation collector additionally checks each raw output's single newest
execution, SCF convergence, normal termination, charge/spin/electron count,
input equality, zero execution and timing exits, and unchanged regular-grid
dimensions. It records timing, peak RSS, and SHA-256 hashes for each input,
output, execution log, and time record. `SHA256SUMS` permits a separate local
integrity check of the curated files.

The separately authorized runner is `paper/pccp/run_derivative_checks.py`.
Re-running the electronic calculations requires the recorded CP2K/model/data
stack. Its host-specific runtime paths must be adapted explicitly; the central
calculation then generates the restart for its eight perturbations. WFN files
are preserved on Terok, not duplicated into this lightweight article dataset.
The exact original restart hashes are retained for provenance.

Force errors are finite difference minus analytic force, in Ha/bohr; stress
errors are in GPa. CP2K prints the stress tensor in bar, which is converted
before comparison. Halving a step need not reduce the error when subtraction
amplifies finite SCF/numerical noise. Only one force and one diagonal stress
component are tested; this is not an exhaustive test of shear, finite-k-point
derivatives, or symmetry-breaking displacements under point-group reduction.

No production energies were recalculated for this analysis. Launcher failures
before a CP2K execution and temporary parsing/coordination records are excluded
from the article dataset; the complete scientific executions are retained.
