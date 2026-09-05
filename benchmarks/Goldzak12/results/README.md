# Validated results

## Campaign closeout

The recovery campaign was closed by user decision on 2026-09-05. The two
remaining HD AlN/v10 executions were cancelled, not classified as scientific
failures. Their remote checkpoints, outputs and terminal closeout records are
preserved. The recurring computation monitor is removed; no new recovery runs
are authorized. The original accepted union remains 440/448.

See `SELECTED-ANALYSIS.md` for the structural comparison and its limitations.
The 352 selected points have exact input/execution packages in
`../accepted-inputs/`, indexed by `eos-selected-input-provenance.json`.
No remote WFN or bulk output was deleted or copied into the paper package.

```bash
python3 -B benchmarks/Goldzak12/scripts/fit_selected_eos.py
python3 -B benchmarks/Goldzak12/scripts/analyze_selected_eos.py
python3 -B benchmarks/Goldzak12/scripts/verify_selected_package.py
```

These commands only analyze existing energies. The reference errors, paired
method differences, aggregate statistics, retained literature values and 72
endpoint-omission tests are written to `eos-selected-*.csv`. The nonshared
method-statistics table excludes the reused BN/C AE curves (N=8); the common
table includes all ten solids (N=10). Cohesive energies are outside this
closeout analysis. The analysis scripts neither trigger new calculations nor
require missing atomic runs.

## Selected structural analysis

`existing-convergence-inventory.json` consolidates 208 previously extracted
clean control records for the retained solids. It preserves settings, energies,
binary/model identifiers and original paths. It is an evidence inventory,
not a final convergence certificate; only matched protocols may be compared.

The current ten-solid selection is in `../protocol/paper-eos-selection.json`.
`eos-selected-source.csv` contains compact execution-accepted periodic energies
and original remote paths, not a claim of final numerical certification.
`eos-selected-fits.csv` and `eos-selected-summary.json` are explicitly marked
`numerical_review_pending` and do not contain cohesive energies. Generate them
without rerunning CP2K or requiring atomic energies:

```bash
python3 -B benchmarks/Goldzak12/scripts/fit_selected_eos.py
```

AlN and AlP use the same nine volumes for all methods. Other selected solids
use ten volumes. Shared AE curves for BN/C are reused explicitly in the hybrid
comparison. Original excluded or unused volume points are not deleted. The
historical remote accepted-union ledger retains its original twelve-solid scope.

The [matched literature report](LITERATURE-COMPARISON.md) and
`eos-literature-*.csv` preserve source-specific intersections and static-lattice
reference conventions. `eos-selected-actual-settings.csv` audits the 352
executed inputs; `../paper/` contains the manuscript and SI tables.

No cohesive energies or missing isolated-atom references enter this analysis.
Final numerical release remains pending; no more SCF jobs are required to
reproduce the selected fits and comparison tables.
