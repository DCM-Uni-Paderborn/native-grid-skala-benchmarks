# Convergence data

## Plane-wave cutoff

`cutoff/aconf8-b200` contains the paired ACONF-8 `ttt` and `gtg`
calculations for direct GPW-GTH and the smoother GAPW_XC-GTH one-center
representation. The 800 Ry conformational energies are 0.747457 and
0.745395 kcal/mol, respectively. At 400 Ry, the same-representation errors are
0.002152 and -0.000181 kcal/mol. This establishes 400 Ry as the initial
GAPW_XC pilot cutoff but does not impose it on direct GPW or all-electron GAPW.

## Atom-centered integration

`atom-grid/h2o` separates radial and Lebedev convergence for an isolated water
molecule. The 50/50 starting grid differs from the 100/590 result by
`3.74715e-4` Ha and integrates 8.00044858 electrons. At 100/590 the integrated
particle number is 7.99999978. The corresponding wall times are 2:43.26 and
3:22.34 on the recorded reference setup. These total-energy data motivate, but
do not by themselves determine, the production grid.

`atom-grid/aconf8` contains the first ACONF-8 diagnostics. The production
decision is based on paired `ttt`/`gtg` calculations so that convergence of a
relative energy is tested directly. The complete paired sequence and its
force/virial checks will be added after validation.

`atom-grid/h2o-fixed-density` separates quadrature changes from self-consistent
response. Rows with `converged=false` in the extracted CSV files are retained
as diagnostic records and are never used as reference values.
