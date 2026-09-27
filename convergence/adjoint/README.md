# Interpolation and Adjoint Checks

The four `.F90` kernels are the actual extracted CP2K production routines and
caller loops used by the archived validation, not independent numerical
implementations. Their hashes match `evidence.json`. The fixture retained here
matches the optimized report. The earlier checked report records its own
fixture hash and is not described as an identical rerun.

The dot-product identity passed a 2e-13 acceptance threshold; its measured
maximum was not printed. The separately measured projection-variant maxima
are 6.1088e-16 and 5.7513e-16. Those quantities must not be conflated.

The normal repository audit verifies the archived source/report hashes without
compiling or running anything. An optional kernel-only rerun can use an existing
compatible GNU/OpenMP CP2K CMake/Ninja build:

```sh
python3 convergence/adjoint/run_check.py --build-dir /path/to/cp2k-build \
  --work-dir /tmp/skala-adjoint-check --check
```

The build must expose `orbital_transformation_matrices_unittest`, whose resolved
link dependencies supply the real CP2K grid/harmonic types. This command does
not rebuild CP2K or perform an SCF. It uses the retained optimized fixture in
checked mode; it does not claim to recreate the unavailable earlier fixture.

Extracted CP2K code is Copyright CP2K developers group and distributed under
GPL-2.0-or-later; see `LICENSE`. Baseline and patched source hashes and compiler
settings are in `evidence.json`.
