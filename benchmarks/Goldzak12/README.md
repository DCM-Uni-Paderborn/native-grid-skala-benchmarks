# LC10 Structural Benchmark

AlN, AlP, BN, BP, C, LiCl, MgO, MgS, Si and SiC form the reported LC10 set.
The directory name preserves the source benchmark's original name.

The package contains 352 unique total energies and 36 independent EOS curves,
yielding 40 representation/material entries through the documented reuse of
the AE BN/C curves. The authoritative selection is
`protocol/paper-eos-selection.json`; input/output hashes are in
`results/eos-selected-input-provenance.json`.

The scripts fit equilibrium lattice constants and bulk moduli and compute
matched literature comparisons using the references in `reference/`.
`basis-controls/` contains the separate Si, diamond and MgO basis tests.
Pressure derivatives remain fit parameters, not an additional reported
benchmark. Atomic cohesive-energy and unused literature comparisons are
outside this package.

Use the repository's `scripts/reproduce_paper.py` to regenerate the current
comparisons without changing archived files.
