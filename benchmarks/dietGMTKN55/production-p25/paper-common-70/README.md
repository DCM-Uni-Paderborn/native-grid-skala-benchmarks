# Fixed Common Set: 70 Reactions

`reference.json` and `reaction-index.json` contain official references,
geometries, stoichiometry and the mapping to each accepted species-route.
`species-results.csv` contains 332 unique routes; all executed inputs are
available and hash-verified in `../accepted-inputs/manifest.csv`.

`native-protocol-comparison.csv` contains all 70 reaction energies, including
W4-11/138. Shared AE species are reused, not counted twice.
`native-and-gauxc-paper-comparison.csv` contains the 65 matched GauXC values
and marks the five unavailable comparisons explicitly. These are not the
earlier paper's 100-reaction AE/ECP--PySCF aggregate population.

No unselected results, quarantines or recovery queues are included.
Run `scripts/audit_pccp_package.py` from the repository root to reconstruct
the reaction sums. Source revisions are in `../../provenance`.
