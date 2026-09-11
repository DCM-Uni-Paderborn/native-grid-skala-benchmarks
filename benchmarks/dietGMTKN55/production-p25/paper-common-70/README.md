# Fixed Common Set: 70 Reactions

`reference.json` and `reaction-index.json` contain official references,
geometries, stoichiometry and the mapping to each accepted species-route.
`species-results.csv` contains 332 unique routes; all executed inputs are
available and hash-verified in `../accepted-inputs/manifest.csv`.

`native-protocol-comparison.csv` contains all 70 reaction energies, including
W4-11/138. Shared AE species are reused, not counted twice.
The molecular interface comparison separates the core protocols:

- Native GAPW-XC/GTH versus GauXC GPW-GTH: 65 shared reactions, reference
  MAEs 2.219302 and 2.181098 kcal/mol, direct MAD 0.463144 kcal/mol.
- Native mixed AE/GTH versus GauXC AE/def2-ECP and PySCF: all 70
  reactions. Reference MAEs are 1.890482/1.890510, 1.237761 and 1.231410
  kcal/mol, respectively. Mixed-direct versus GauXC AE/ECP has MAD 1.888876.

PySCF is the only PySCF variant used here. It denotes the primary reference
with the atomic-radius adjustment in Becke partitioning disabled (unit radii).
The selected extract, comparison files and paper tables all use this variant.

`gauxc-source-common-70.json` contains only the selected collaborator species
values and their reaction mapping. The source is Stefano Battaglia's dataset
archived in `DCM-Uni-Paderborn/Molecular-Skala-in-CP2K` at revision
`25e80baf68bee449d1dc2fd11cc359fa503272b8`. Input/output paths and SHA-256
hashes refer to that pinned source, avoiding duplication of the raw archive.
Matching uses translated atomic geometries, charge, multiplicity and signed
stoichiometry. References and weights are independently checked. Source
reaction identifiers are preserved, including the BH76/BH76RC alias.

`molecular-interface-comparison.csv` contains reaction totals and separate D3
contributions. `molecular-interface-summary.json` contains reference errors,
direct differences, and dispersion sensitivity. As-reported totals are used
for the headline comparison. The retained HEAVY28/5 GauXC AE/PySCF D3 mismatch
is documented, not silently corrected. Electronic-only differences are
reported separately. Different Gaussian bases and core Hamiltonians prevent
interpreting these as same-protocol numerical identity tests.

`native-and-gauxc-paper-comparison.csv` retains the earlier 65-reaction
comparison with the same native D3 term added to both electronic energies.
Its MAD 0.462396 kcal/mol is the dispersion sensitivity quoted in the ESI,
not the new as-reported total-energy MAD. The five unavailable comparisons
are explicit. The prior paper's 100-reaction aggregates are context only.

Reproduce the new comparison without CP2K or network access:

```bash
python3 -B scripts/compare_molecular_interfaces.py --write
```

To recreate the selected extract from the collaborator repository, add
`--source-repo ../Molecular-Skala-in-CP2K` (NumPy and openpyxl required).

No unselected results, quarantines or recovery queues are included.
Run `scripts/audit_pccp_package.py` from the repository root to reconstruct
the reaction sums and both interface intersections. Source revisions are in
the selected extract and `../../provenance`.
