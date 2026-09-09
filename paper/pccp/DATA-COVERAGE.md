# PCCP Data Coverage

The online main text and ESI were synchronized and checked on 10 September
2026. The following selection is fixed by those documents, not by the former
production queues. Original selected inputs and outputs have not been edited.

| Manuscript / ESI claim | Supporting dataset | Selected executions |
| --- | --- | ---: |
| Main molecular table; ESI complete reaction table and error distribution | dietGMTKN55 common-70 reaction index, references, energies and accepted-input manifest | 332 |
| ESI 65-reaction GauXC comparison | Matched comparison CSV, external source revision | Reuses the molecular set |
| Quoted 100-reaction GauXC/PySCF context | ESI context table and cited author-manuscript version | No additional native executions |
| Main molecular-crystal table and figure; ESI total energies | X23-mini accepted base outputs and validation records | 24 |
| ESI paired grid, cutoff, k-mesh, symmetry tolerance and CPU/GPU checks | X23-mini accepted controls with exact parents | 7 |
| ESI LC10 structural tables, literature comparisons, EOS plots, B-prime and compressibility | 352 selected outputs, exact inputs, EOS selection and offline fitting scripts | 352 |
| ESI paired ACONF cutoff table | Nine cutoff values, two geometries, GPW and GAPW-XC | 36 |
| ESI GPW water radial/angular table | Self-consistent water quadrature series | 16 |
| ESI isolated-cell padding | Three species at 22, 24 and 25 angstrom | 9 |
| ESI force/stress finite differences | Five representations, nine executions each | 45 |
| ESI full/reduced k-point table | Seven tight-SCF pairs, including original inputs and outputs | 14 |

Total: **835 selected executions**. Reference transcriptions, generated
tables/figures, source TeX, analysis software and hash/provenance records are
supporting dependencies, not additional benchmark campaigns.

The file manifest is an exact allowlist: the audit rejects both missing files
and extra files. It also reconstructs reaction and lattice energies from the
selected final outputs. Offline regeneration reproduces all seven generated
ESI table files byte-for-byte, and all 40 EOS fits and 72 endpoint-sensitivity
records numerically exactly in the verification environment.

Excluded material was preserved outside this repository before pruning.
Removal from the current tree does not erase old Git history. No history
rewrite, new SCF calculation, or modification of the earlier GauXC paper was
performed for this reconciliation.
