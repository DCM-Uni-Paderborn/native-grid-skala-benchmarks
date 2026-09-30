# Publication Coverage

[data-coverage.json](data-coverage.json) maps all **32 tables and four figures**
in the current manuscript and ESI to deposited data and analysis. Labels,
rather than printed numbers, make this map stable under renumbering.

The coverage comprises the 62-reaction molecular comparison; molecular
crystals and their DFT/correlated reference context; thirteen ice phases;
LC10 structural results and basis controls; 50 primary band gaps and seven
direct/one-centre comparisons; and the cutoff, quadrature, derivative,
adjoint, symmetry and cell-padding controls.

The manuscript text, source references and figure assets are fixed by
[source-snapshot.json](source-snapshot.json). The audit verifies all reference
targets, citations and assets and rejects missing or extra table/figure
mappings. It also reconstructs the numerical comparisons from the archived
outputs. No raw molecular GauXC data are replicated.

The crystal comparison in manuscript Table 3 includes the selected literature
methods formerly shown in a separate ESI table. The audit checks the table's
method order, all energies and MAEs, and the DMC statistical uncertainties
against the deposited records. Its obsolete ESI mapping has been removed.
