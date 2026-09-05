# Goldzak12 periodic benchmark

## Analysis and manuscript assets

The campaign is closed: no new SCF or atomic-reference runs are authorized.
The paper reports structural properties, not cohesive energies. Start with
the [structural analysis](results/SELECTED-ANALYSIS.md),
[matched literature comparison](results/LITERATURE-COMPARISON.md), and
[manuscript/SI assets](paper/README.md).

The selected package contains 352 unique energies, 36 independent curves and
40 method/solid comparisons. All executed inputs use 800/60 Ry and 150/770
radial/Lebedev grids. Gamma-centered meshes range from 4^3 to 7^3, can change
along a volume sequence, and match between methods at each volume; see
[the actual-input audit](results/eos-selected-actual-settings.csv).

All native lattice constants underestimate experiment. The small mixed-core
one-center lattice shifts do not remove this bias. Numerical review remains
pending; a complete final-setting error budget is not claimed.

Offline reproduction (Python 3, NumPy and Matplotlib, from repository root):

```bash
python3 -B benchmarks/Goldzak12/scripts/verify_selected_package.py
python3 -B benchmarks/Goldzak12/scripts/fit_selected_eos.py
python3 -B benchmarks/Goldzak12/scripts/analyze_selected_eos.py
python3 -B benchmarks/Goldzak12/scripts/compare_selected_literature.py
python3 -B benchmarks/Goldzak12/scripts/build_paper_tables.py
```

Only selected input packages, compact results/provenance, source tables and
offline analysis tools are included in the paper update on GitHub. Local
campaign queues and failed recovery scripts are not part of that package.
Large restart WFNs remain in the source archive; their recorded hashes do not
redistribute the binaries needed for identical restart trajectories.

## Current paper selection (2026-09-05)

The user-approved structural comparison uses ten solids: AlN, AlP, BN, BP,
C, LiCl, MgO, MgS, Si, and SiC. LiF and LiH are excluded consistently across
methods; their original inputs and results remain preserved. LiF has an
inconsistent GAPW_XC-GTH EOS despite formal SCF convergence, with unresolved
cause. LiH has AE/GAPW_XC-GTH minima below the sampled volume range.

The common fit window is v01-v09 for all AlN and AlP methods and v01-v10 for
the other eight solids. The exact selection is recorded in
[`protocol/paper-eos-selection.json`](protocol/paper-eos-selection.json).
Shared AE results are reused for the two hybrid protocols of BN and C;
these are not independent calculations. Atomic-reference completion and
publication of cohesive energies are separate from the structural comparison.

Execution acceptance is not a numerical-accuracy certificate. The selected
fits remain under review for fit-window sensitivity and final-setting grid,
cutoff, and k-point convergence. The existing control inventory documents
available evidence, not a completed universal error bound.

## Curated contents

- `accepted-inputs/`: 352 exact frozen/executed input packages and execution/time records.
- `results/`: selected energies, fits, paired shifts, endpoint tests, settings and provenance.
- `reference/`: system-resolved literature values and source conventions.
- `protocol/`: fixed selection, systems, volume sampling and method definitions.
- `paper/`: manuscript and SI assets.
- `scripts/`: offline fitting, literature comparison, table generation and integrity tests.

The intended full twelve-solid design remains in the local campaign archive.
The GitHub paper package contains only material supporting the selected analysis.
Its ten-solid set differs from the earlier GFN LC10: the matched GFN comparison
contains nine solids, excluding MgO. Missing literature values are never treated
as zero errors, and full-set aggregate statistics are not used as matched results.
