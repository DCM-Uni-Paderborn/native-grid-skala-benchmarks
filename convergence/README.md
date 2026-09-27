# Numerical Controls

The deposited controls correspond to the current manuscript and ESI:

- `ae-cutoff/`: five H2O GAPW-AE/QZVPP-ae calculations at relative cutoff
  60 Ry and 50/50 atom quadrature, referenced to 1200 Ry.
- `cutoff/aconf8-b200/`: GTH GPW/GAPW-XC conformational-energy cutoff series.
- `atom-grid/`: independent radial/angular H2O quadrature tests.
- `adjoint/`: extracted interpolation/backprojection kernels and validation.
- `derivatives/periodic-water-20260909/`: 45 component-specific water
  force/stress finite-difference executions.
- `symmetry/`: seven full/reduced k-point comparisons.
- `basis-sensitivity-20260915/`: the index for crystal, ice and solid
  basis/grid controls.

Cell-padding data are in `diagnostics/molecular-padding-convergence`.
Molecular-crystal device, cutoff and symmetry controls remain with X23-mini.
The 62-reaction molecular dataset includes its own TZVPP-ae/QZVPP-ae comparison.
