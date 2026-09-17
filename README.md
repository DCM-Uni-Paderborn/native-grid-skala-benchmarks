# Native-grid Skala: PCCP Supporting Data

This repository contains only the results and reproducibility material used in
**Native implementation of the machine-learned Skala exchange-correlation functional
in CP2K: Unified one-centre reconstruction for molecular and condensed-phase
calculations** and its electronic supplementary information (ESI). The TeX sources in
[paper/pccp](paper/pccp) were checked against the
[online PCCP project](https://www.overleaf.com/project/6aa123c081b754e6cf561ee5)
on 16 September 2026.

| Paper / ESI content | Supporting data |
| --- | --- |
| 70 common dietGMTKN55 reactions; separate 65-reaction GauXC GPW-GTH and 70-reaction GauXC AE/ECP--PySCF comparisons | [Molecular data](benchmarks/dietGMTKN55/production-p25/paper-common-70), 332 unique native species-route executions and a pinned collaborator-data extract |
| CO2, NH3 and urea, four representations | [Molecular crystals](benchmarks/X23-mini), 24 original base executions, seven numerical controls and twelve paired AE basis/grid controls, with all three selected AE pairs at QZVPP |
| Six completed targeted molecular basis comparisons and ethane state/grid controls | [Paired molecular controls](benchmarks/dietGMTKN55/basis-controls), 33 executions, separate from the unchanged 70-reaction statistics |
| Same-twelve-phase ice basis comparison, relative energies to Ih and two symmetry checks | [Ice controls](benchmarks/DMC-ICE13), 28 executions including both molecular references, excluding XIII at both basis levels |
| Matched five-volume Si, diamond and MgO basis tests | [Solid basis controls](benchmarks/Goldzak12/basis-controls), 15 QZVPP points and three TZVPP runtime checks, separate from the uniform ten-solid comparison |
| Ten-solid structural comparison, in the ESI | [LC10](benchmarks/Goldzak12), 352 unique energies, 36 independent EOS, 40 method-solid comparisons |
| ACONF cutoff and water radial/angular quadrature | [Cutoff](convergence/cutoff/aconf8-b200) and [quadrature](convergence/atom-grid) |
| Isolated-cell 22/24/25 angstrom padding | [Padding](diagnostics/molecular-padding-convergence) |
| Force/stress finite differences | [45 periodic water cases](convergence/derivatives/periodic-water-20260909) |
| Seven full/reduced k-point comparisons | [14 symmetry executions](convergence/symmetry) |

Shared species and AE curves are counted once. Unreported campaigns, recovery
queues, failed diagnostics, atomic cohesive-energy calculations, unselected
solids, old cell sizes and earlier paper drafts are excluded.

## Reproducibility

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md). Run
`python3 -B scripts/audit_pccp_package.py` to verify file inventory, hashes,
selected outputs, reaction/lattice energies and numerical checks, without
running CP2K or contacting any cluster.

[File manifest](paper/pccp/file-manifest.json) maps every included file to its
paper/ESI scope and SHA-256 hash; it excludes itself to avoid a circular hash.
[Execution provenance](paper/pccp/execution-provenance.json) records the
original locations and hashes of the retrieved selected executions.

Exact inputs retain machine-specific paths as historical provenance. Model
weights, binaries and restart wavefunctions are identified by hashes rather
than redistributed. They are unnecessary for reproducing the reported
analysis, but an identical restart trajectory requires the original WFN.

The earlier molecular GauXC paper is not duplicated or changed here. Its
100-reaction aggregates are quoted in their original population; only the
matched comparison used here is included. Numerical qualifications in the
paper/ESI remain applicable; SCF convergence is not a universal accuracy bound.
