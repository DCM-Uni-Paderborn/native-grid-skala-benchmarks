# Selected PCCP X23-mini lattice energies

24/24 base cases accepted; 12/12 crystal/molecule pairs complete.

| System | Method | E_latt (kJ/mol) | DMC (kJ/mol) | Signed deviation (kJ/mol) |
| --- | --- | ---: | ---: | ---: |
| CO2 | gapw-ae | -36.570 | -29.4 +/- 0.2 | -7.170 |
| NH3 | gapw-ae | -38.744 | -38.2 +/- 0.1 | -0.544 |
| urea | gapw-ae | -110.427 | -108.5 +/- 0.3 | -1.927 |
| CO2 | gapw-gth-direct | -37.376 | -29.4 +/- 0.2 | -7.976 |
| NH3 | gapw-gth-direct | -40.094 | -38.2 +/- 0.1 | -1.894 |
| urea | gapw-gth-direct | -113.805 | -108.5 +/- 0.3 | -5.305 |
| CO2 | gapw-gth-one-center | -37.392 | -29.4 +/- 0.2 | -7.992 |
| NH3 | gapw-gth-one-center | -40.076 | -38.2 +/- 0.1 | -1.876 |
| urea | gapw-gth-one-center | -113.719 | -108.5 +/- 0.3 | -5.219 |
| CO2 | gapwxc-gth | -37.408 | -29.4 +/- 0.2 | -8.008 |
| NH3 | gapwxc-gth | -40.074 | -38.2 +/- 0.1 | -1.874 |
| urea | gapwxc-gth | -113.725 | -108.5 +/- 0.3 | -5.225 |

Urea AE uses the validated QZVPP/200-974 crystal and molecule. The original TZVPP pair remains in complete_pairs in the JSON and in the basis-convergence evidence. Other selections are unchanged. Four additional Urea controls separate grid and basis sensitivity.

Negative deviations indicate stronger binding. DMC error bars are statistical; total DMC accuracy is estimated at about 1-2 kJ/mol.

These are preliminary single-point results with bounded representative numerical checks, not a comprehensive convergence study. The differences must not yet be attributed solely to the functional or density representation. Numerical controls below quantify only their stated changes and cannot be generalized to every system or representation. Basis/BSSE and untested effects remain unresolved.

## Matched representation differences

| System | Method A | Method B | A minus B (kJ/mol) |
| --- | --- | --- | ---: |
| CO2 | gapwxc-gth | gapw-ae | -0.838 |
| CO2 | gapwxc-gth | gapw-gth-direct | -0.032 |
| CO2 | gapwxc-gth | gapw-gth-one-center | -0.016 |
| CO2 | gapw-gth-one-center | gapw-gth-direct | -0.016 |
| NH3 | gapwxc-gth | gapw-ae | -1.330 |
| NH3 | gapwxc-gth | gapw-gth-direct | +0.020 |
| NH3 | gapwxc-gth | gapw-gth-one-center | +0.002 |
| NH3 | gapw-gth-one-center | gapw-gth-direct | +0.018 |
| urea | gapwxc-gth | gapw-ae | -3.298 |
| urea | gapwxc-gth | gapw-gth-direct | +0.080 |
| urea | gapwxc-gth | gapw-gth-one-center | -0.006 |
| urea | gapw-gth-one-center | gapw-gth-direct | +0.085 |

Only completed pairs for the same system enter these differences. GAPW_XC versus direct GAPW changes more than the one-center treatment; isolate that effect only with the two GAPW-GTH variants.

## Accepted numerical controls

| Case | Control | Control minus parent (kJ/mol of molecules) |
| --- | --- | ---: |
| gapwxc-gth/CO2/solid | cpu-model | -0.000127 |
| gapwxc-gth/CO2/solid | cutoff-1200 | -0.015499 |
| gapwxc-gth/CO2/molecule | atom-200-974 | -0.155241 |
| gapwxc-gth/CO2/solid | atom-200-974 | -0.110687 |
| gapwxc-gth/CO2/solid | kmesh-5 | +0.000391 |
| gapwxc-gth/NH3/solid | sym-tol-1e4 | -0.000429 |
| gapw-ae/urea/solid | cpu-model | -0.000379 |

Controls are not additional base cases. The CO2 CPU/GPU comparison tests the crystal model evaluation on Spark, not other hosts, quadrature or basis completeness.

For crystal-only k-mesh and symmetry-tolerance changes, the isolated molecule is unchanged, so the tabulated shift is also the lattice-energy change. Cutoff and quadrature changes require matched controls of both phases.

## Paired lattice-energy sensitivity

| System | Method | Control | Change in E_latt (kJ/mol) |
| --- | --- | --- | ---: |
| CO2 | gapwxc-gth | atom-200-974 | +0.044554 |

Each difference uses both accepted phases relative to their exact accepted parents. It assesses only the stated numerical change for that system and representation.

Reference: [Della Pia et al., PRL 133, 046401 (2024)](https://doi.org/10.1103/PhysRevLett.133.046401).
