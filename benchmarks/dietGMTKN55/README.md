# Matched 62-Reaction Molecular Comparison

`dataset.json` defines 62 reactions and 140 distinct native species at three
levels: GAPW-AE/TZVPP-ae, GAPW-AE/QZVPP-ae and GAPW-XC/TZV2P-GTH.
Their 420 selected native inputs and outputs are in `executions/`.
The native AE series contains no mixed AE/GTH species.

`paper-results.json` fixes the reported comparison values.
`external-comparison.json` is the minimal matched GauXC/PySCF extract, with
stoichiometric and geometry correspondence, source revision and hashes.
The complete external inputs, outputs and benchmark are maintained in the
linked molecular GauXC repository, not copied here.

`provenance/` retains selection records and the original predecessors needed
to trace the selected molecular results. These are not additional benchmark
entries. The QZVPP ethane selection is explicit rather than silently replacing
its original solution. Obsolete 70-/65-reaction aggregates are not mixed with
the present 62-reaction statistics.

Run `python3 -B scripts/analyze_molecular.py` from the repository root to
verify hashes and recompute every reaction energy and matched statistic.

Native inputs preserve their full original basis names and paths. Publication
shorthand such as QZVPP-ae does not rename the computational input identifiers.
