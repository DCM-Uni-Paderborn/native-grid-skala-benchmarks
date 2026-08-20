# Native-Grid Skala Benchmarks

Reproducibility data for **Native-Grid Skala in CP2K for Molecular and Periodic Electronic-Structure Calculations** by Johann Potoschnig, Juerg Hutter, and Thomas D. Kuehne.

The repository is private while calculations and the manuscript are in progress. It will contain the inputs, source and model metadata, compact raw outputs, derived tables, and analysis scripts needed to reproduce the paper. Large wavefunction restart files and transient scheduler data are deliberately excluded.

## Scope

- plane-wave cutoff convergence for direct and `GAPW_XC` smooth fields;
- radial/Lebedev and one-center integration convergence;
- particle-number, density-matrix-adjoint, force, virial, and stress checks;
- molecular convergence to the GauXC reference route;
- full and symmetry-reduced k-point validation;
- dietGMTKN55 molecular benchmarks;
- LC10 periodic benchmarks;
- CPU/GPU timing and peak-memory metadata where scientifically relevant.

## Layout

```text
convergence/    Numerical convergence inputs, outputs, and extracted tables
metadata/       CP2K revision, build, model checksum, and machine information
protocols/      Frozen scientific protocols and acceptance criteria
scripts/        Input preparation, auditing, and result extraction
benchmarks/     Final dietGMTKN55 and LC10 data
```

Numerical settings are frozen only after energy differences, particle number, and analytical derivatives have converged together. A tighter setting already used for a calculation is never replaced by a lower-accuracy setting.
