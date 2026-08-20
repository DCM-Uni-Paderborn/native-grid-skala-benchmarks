# Software and hardware metadata

## CP2K reference build

- CP2K version: 2026.2 development
- source revision: `a044f58aaebbb9979b3034dee6aab15fdfa96238`
- source worktree: clean when the convergence calculations were started
- compiler: GCC 14.3.0
- principal build flags: MPI, OpenMP, CUDA offload, DBCSR acceleration,
  Spglib, OpenBLAS, LibTorch, and LibTorch CUDA
- CP2K executable SHA-256:
  `41d006f21158ec2725d6f1555468118fff0f0f597a3bf98d3f897c88f48bbb91`

## Skala artifacts

- Skala-1.1 Rev1 CPU artifact SHA-256:
  `7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd`
- Skala-1.1 Rev1 CUDA artifact SHA-256:
  `f848eae769dca91741a518ae7275d10caac398ab21db649f91bc1f136872f223`

The model artifacts are not duplicated in this repository. Their checksums bind
the reported calculations to the exact artifacts distributed separately under
their applicable license.

## B200 convergence node

- node class: HZDR Rosi NVIDIA GB200 partition
- CPU: two AMD EPYC 9755 sockets, 128 physical cores per socket
- GPU: four NVIDIA B200 accelerators, 183359 MiB reported memory per device
- NVIDIA driver: 595.71.05
- CUDA toolkit used for the build: 12.8
- interactive allocation memory limit: 375 GiB

The all-electron ACONF calculations require approximately 107--113 GiB of
host memory per process. The dynamic runner therefore limits this allocation to
three simultaneous calculations even when four GPUs are visible. This avoids
memory-pressure termination and is recorded separately from GPU occupancy.

Every timing table additionally records the executable revision, process/thread
layout, device count, and maximum resident set size. Absolute timings from
different accelerator classes are not combined into a scaling series.
