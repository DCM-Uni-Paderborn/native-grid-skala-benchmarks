# Numerical convergence protocol

## Plane-wave grid

Plane-wave cutoffs are converged separately for each smooth density representation. The direct GPW-like valence field and the softer `METHOD GAPW_XC` field are not assigned a common cutoff by assumption. Relative energies, particle number, forces, and stress are compared against the tightest calculation of the same representation.

## Atom-centered composite grid

`RADIAL_GRID=50` and `LEBEDEV_GRID=50` are starting values, not accepted production settings. Radial and angular quadratures are separated by the sequence

1. fixed radial grid with the implemented 50-, 110-, 194-, 302-, 434-, and 590-point Lebedev grids as required;
2. fixed converged angular grid with 50, 100, 150, and 200 radial points as required;
3. a simultaneous tightening check;
4. transfer checks for C/H, O, a heavy pseudopotential kind, and a periodic solid.

CP2K interprets `LEBEDEV_GRID` as a requested minimum size and selects the smallest implemented grid that is at least this large. The benchmark inputs use the implemented sizes explicitly so that, for example, a request of 150 is not ambiguously reported instead of the realized 194-point grid.

The production grid is fixed at `RADIAL_GRID=100` and `LEBEDEV_GRID=434`. The latter is the realized CP2K rule for a nominal request near 400 points. This setting is validated further through reaction or conformational energies, integrated particle number, force finite differences, and virial/stress finite differences. Total-energy convergence alone is insufficient.

The first isolated-H2O series at 800 Ry demonstrates why the inherited 50/50 one-center grid is only a starting point. Relative to a 200/974 reference, the 50/50 grid changes the energy by -369.201 microhartree and integrates 0.0004486 excess electrons. The 100/434 and 100/590 grids reduce these residuals to 8.379 and 6.084 microhartree and to -3.031e-7 and -2.206e-7 electrons, respectively.

The present native implementation uses an unpruned radial-times-Lebedev product. Nominal point counts therefore grow from 2,500 points per atom at 50/50 to 43,400 at 100/434 and 194,800 at 200/974. On the B200 H2O pilot, however, the complete wall-time factors were only 1.00, 1.20, and 1.83 because atom rows are processed in accelerator batches and the remaining CP2K work is unchanged. The 100/434 grid is the production setting, 100/590 its tightening control, and 200/974 a numerical reference.

The finer angular requirement is not specific to a one-center correction. Bonding, neighboring basis functions, and the smooth atom partition make the density anisotropic on each target-atom grid. Skala then uses density, density gradient, and positive kinetic-energy density in nonlinear nonlocal descriptors, so angular quadrature errors can propagate through descriptor couplings. Molecular GauXC hides much of the corresponding point-count burden through shell-dependent grid pruning. Equivalent pruning of the native grid is a performance optimization only if forward selection and analytical adjoints remain identical.

The test is repeated for all-electron GAPW and pseudopotential `GAPW_XC` with `PAW_ONE_CENTER`. A direct-valence route is retained as a same-density reference. `GAPW_ACCURATE_XCINT` is enabled in production and a paired on/off control verifies that its classical smooth/soft switching weights do not define the atom-composite Skala field.

## Production gate

The complete dietGMTKN55 and LC10 calculations remain blocked until the plane-wave and atom-grid protocols are frozen. Production inputs are generated once from the validated settings and pass the structural audit in `scripts/audit_diet_inputs.py` before submission.
