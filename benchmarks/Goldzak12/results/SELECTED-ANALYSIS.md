# Selected structural comparison

The expanded [literature comparison](LITERATURE-COMPARISON.md) contains all
selected native values, matched-source statistics and system-resolved
HF/MP2/SCS/SOS, LDA/PBE/PBEsol/M06-L/SCAN/HSE06/HSEsol/TM and own GFN
comparisons. Each source is scored only on its exact intersection with the
selected set; different sample sizes are never combined into one ranking.

## Scope

The closed campaign is analyzed on AlN, AlP, BN, BP, C, LiCl, MgO, MgS, Si
and SiC. All methods use v01-v09 for AlN/AlP and v01-v10 otherwise. This is
352 unique energy points, 36 independent EOS curves and 40 method/solid
comparisons. BN and C reuse the AE curve in both hybrid representations.
No additional SCF or isolated-atom calculation is required for this analysis.

LiF is excluded because the GAPW_XC-GTH EOS is inconsistent despite formal
SCF convergence; its cause remains unresolved. LiH is excluded because the
AE/GX minima are outside the sampled range. The exclusions apply to all
methods, not according to agreement with experiment. All original data remain
preserved. The missing AlN/AlP v10 points are not needed for the common window.

## Reference comparison

Signed errors are native minus the experimental reference transcribed in
`../reference/goldzak2022.csv`; its source and zero-point conventions are
documented in `../reference/SOURCES.md`. Values below are descriptive results
of the selected data, not certified discretization-error estimates.

| Method | N | Lattice MAE (angstrom) | Bulk-modulus MAE (GPa) |
| --- | ---: | ---: | ---: |
| GAPW_XC-GTH | 10 | 0.08313 | 16.61 |
| GAPW-AE | 10 | 0.07284 | 16.89 |
| Hybrid GAPW-AE/GTH, direct | 10 | 0.08443 | 21.77 |
| Hybrid GAPW-AE/GTH, one-center | 10 | 0.08445 | 21.07 |

All four methods underestimate every selected reference lattice constant.
This common bias must not be described as quantitative agreement without
qualification. The current data alone do not separate functional bias from
basis, quadrature and other numerical contributions.

## Representation differences

The mean absolute lattice difference between GAPW_XC-GTH and GAPW-AE is
0.01893 angstrom (maximum 0.04029 angstrom). The corresponding bulk-modulus
mean absolute difference is 10.55 GPa. These are paired representation
differences, not errors relative to an exact all-electron calculation.

The one-center-minus-direct hybrid lattice difference is small on the
selected set: mean absolute 0.000135 angstrom and maximum 0.000542 angstrom.
The corresponding bulk-modulus differences are 1.11 GPa mean absolute and
5.75 GPa maximum. The common-set averages include two exactly reused AE
curves; `eos-selected-method-statistics-nonshared.csv` additionally reports
the eight solids with independent hybrid calculations. Small lattice effects
do not establish that every one-center energy contribution is negligible.

## Fit sensitivity and numerical limits

Every selected full-window minimum is bracketed and the maximum EOS fit RMS
is 0.674 meV/atom. Omitting either endpoint independently produces lattice
shifts up to 0.002444 angstrom and bulk-modulus shifts up to 9.997 GPa
(GAPW_XC-GTH/AlN, removing v09 from the selected nine-point window).
LiCl loses minimum bracketing if its smallest-volume point is removed, so
that point is retained. No endpoint has been discarded to improve statistics.

These tests measure fit-window sensitivity, not cutoff or k-point errors.
Consequently a universal 1 GPa bulk-modulus accuracy claim is unsupported.
The existing grid, cutoff and k-point controls require protocol-matched
curation before numerical release; neither SCF acceptance nor a smooth EOS
alone certifies the final production settings. No broad rerun is planned.

The structural comparison is the primary closeout product. Cohesive energies
are not reported here, and missing atomic reference runs do not reopen the
campaign. Final scientific release remains explicitly pending numerical
review. This document is a repository analysis note; the accompanying
manuscript tables and SI discussion are collected in `../paper/`.
