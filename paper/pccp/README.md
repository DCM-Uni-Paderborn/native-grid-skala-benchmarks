# Native-grid Skala: PCCP Manuscript

Online working copy: https://www.overleaf.com/project/6aa123c081b754e6cf561ee5

This separate project was copied from the current native-grid manuscript before
editing. The earlier molecular GauXC manuscript is not modified.

## Article Data

- Main text: native GPW/GAPW reconstruction, 70 common molecular reactions,
  and the electronic lattice energies of CO2, NH3, and urea.
- Theory: separate GPW, GAPW-AE, direct/one-centre ECP, and XC-specific
  representations; classical accurate integration versus joint Skala field
  reconstruction; periodic partitions and discrete adjoints. Implementation
  details were checked against source revision
  `37aedacbb33677307818301a47a09ae435712293` used by the derivative tests.
  The current online GauXC Memo (Overleaf project
  `6a4609286a3fce49c980bf95`) was also consulted on 2026-09-09, without
  modifying it. Density representation, quadrature layout, and core treatment
  are kept distinct; cross-term counts and input-selection details were
  independently checked rather than copied verbatim.
- Supporting information: numerical controls, cross-code comparisons,
  ten-solid EOS and literature comparisons, pressure derivatives, and
  compressibilities.
- Existing benchmark energies are reused. Additional finite-difference checks
  are distinct from the production benchmarks.

The EOS analysis uses 352 unique energies and 36 independent curves, yielding
40 method/material entries after explicitly identified AE reuse. No EOS is
inferred from the three molecular-crystal single-point pairs.

## Rebuilding

From this directory, run:

```bash
python3 build_submission_assets.py
python3 analyse_derivative_checks.py \
  --report ../../convergence/derivatives/periodic-water-20260909/validated-results.json \
  --out ../../convergence/derivatives/periodic-water-20260909 \
  --tex derivative-results.tex
latexmk -pdf main.tex supplementary_information.tex
```

Asset generation requires NumPy, Matplotlib, and the existing Goldzak12 fitting
dependencies. It reads the curated repository datasets and does not launch CP2K.
`source-data-sha256.json` identifies numerical inputs to the table generator.
The handwritten analysis includes remain separate from generated tables.

The shared `latexmkrc` first refreshes the independent ESI in `esi-build/`.
The main text imports its labels with an `esi-` prefix, so section and table
numbers remain synchronized when the ESI changes. PDF links target
`supplementary_information.pdf`. Compile in a separate staging copy to keep
temporary build products out of the curated repository inventory.

Abbreviations are defined independently in the abstract, main text, and ESI.
Table and figure captions use definitions introduced in the preceding prose.
Established names such as PBE, D3, BJ, CP2K, and GauXC are not expanded.

The script `run_derivative_checks.py` is the separately authorized derivative
test driver. Its runtime paths identify the preserved Terok production stack;
they must be adapted, with matching source/model/data provenance, for a new host.
It must not be confused with the manuscript/table build.

The 45 periodic-water derivative executions are complete. The curated dataset
is in `convergence/derivatives/periodic-water-20260909` at the repository root.
The maximum absolute errors over both steps and all five representations are
2.856e-6 Ha/bohr for the tested force component and 0.4722 MPa for the tested
diagonal stress component. These are component-specific consistency checks,
not full phonon, elasticity, or finite-k-point derivative validation.

## Submission Review

The PCCP text was checked against the original article brief on 2026-09-09.
The main narrative retains the full GPW/GAPW theory and joint one-centre
reconstruction, followed by the molecular comparison and CO2/NH3/urea binding.
Quantitative cutoff and other numerical-control details are in the ESI.
The ten-solid analysis remains supporting evidence rather than a second
main-text benchmark. Molecular accuracy and condensed-phase transfer are
distinguished from basis/core differences and numerical validation.
The CC4S outlook includes periodic coupled-cluster machine-learning work;
the QMC discussion explicitly includes periodic FCIQMC and the TurboRVB
kagome application with Azadi.

No production energy, benchmark selection, or model was changed in this
editorial revision, and no additional electronic-structure calculation was
launched. The ESI describes the LiF/LiH exclusions concisely through the
availability of consistent EOS data across all representations. Detailed
diagnostics remain in the underlying records; no exclusively SCF-based
explanation is substituted.

The affiliations and funding are adapted from
the current GauXC manuscript. The PCCP author list was subsequently revised
as requested, with unused affiliations removed and the remaining ones renumbered.
The equal-contribution statement and its author markers were removed as requested.
No conflict-of-interest
statement was present in its main text or SI; author confirmation is required
before submission. The manuscript marks that issue explicitly.

The public, versioned dataset release must be completed before submission. No
journal submission, git commit, or remote push is performed by these build steps.
