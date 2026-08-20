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

The production grid is the least expensive setting for which further tightening leaves the reaction or conformational energy, integrated particle number, force finite differences, and virial/stress finite differences unchanged within the reported tolerances. Total-energy convergence alone is insufficient.

The test is repeated for all-electron GAPW and pseudopotential `GAPW_XC` with `PAW_ONE_CENTER`. A direct-valence route is retained as a same-density reference. `GAPW_ACCURATE_XCINT` is enabled in production and a paired on/off control verifies that its classical smooth/soft switching weights do not define the atom-composite Skala field.

## Production gate

The complete dietGMTKN55 and LC10 calculations remain blocked until the plane-wave and atom-grid protocols are frozen. Production inputs are generated once from the validated settings and pass the structural audit in `scripts/audit_diet_inputs.py` before submission.
