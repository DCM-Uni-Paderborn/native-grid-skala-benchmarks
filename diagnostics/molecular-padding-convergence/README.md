# Isolated-cell convergence beyond 20 angstrom

This diagnostic contains the matched ACONF-5 and ACONF-8 cell-convergence inputs used to select the isolated-molecule protocol. The production electronic settings are held fixed while the total Cartesian padding is increased. A padding value `p` means

```text
L_i = molecular_extent_i + p,
```

with centered coordinates, `PERIODIC NONE` in `CELL` and `POISSON`, and `POISSON_SOLVER ANALYTIC`.

## Validated native-grid result

| Padding (angstrom) | Maximum species change (microhartree) | Maximum reaction change (kcal/mol) | Maximum electron-count change | Status |
|---:|---:|---:|---:|---|
| 10 | 254054.088 | 162.659672 | 1.37e-6 | fail |
| 12 | 141083.364 | 118.802591 | 1.96e-7 | fail |
| 14 | 25210.150 | 15.036497 | 2.57e-7 | fail |
| 16 | 8044.599 | 5.014691 | 1.08e-7 | fail |
| 18 | 543.313 | 0.340330 | 4.25e-8 | fail |
| 20 | 15.899 | 0.010016 | 5.39e-8 | fail |
| 22 | 0.348 | 0.000298 | 9.72e-8 | pass |
| 24 | 0.189 | 0.000163 | 1.43e-7 | pass |
| 25 | 0.000 | 0.000000 | 0.00 | reference |

All changes are relative to 25 angstrom. The acceptance limits are 5 microhartree per species, 0.005 kcal/mol per reaction, and 1e-5 electrons. Production uses 25 angstrom of total padding as a conservative common setting. The tab-separated source table is `native-series-vs-p25.tsv`; the lower-padding inputs and GauXC controls are included in the curated diagnostic package and are evaluated with `scripts/check_padding_convergence.py`.

Only outputs with both `SCF run converged` and `PROGRAM ENDED AT` are scientifically valid. Interrupted attempts are excluded from all reported differences.
