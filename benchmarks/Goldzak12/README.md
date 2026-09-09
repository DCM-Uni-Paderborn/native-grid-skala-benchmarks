# LC10: PCCP ESI Structural Data

Only AlN, AlP, BN, BP, C, LiCl, MgO, MgS, Si and SiC are included.
The directory name retains the source benchmark provenance.
There are 352 unique energies and 36 independent EOS curves; shared AE curves
for BN and C yield 40 method-solid comparisons without duplicate calculations.

The authoritative selection is `protocol/paper-eos-selection.json`.
`accepted-inputs` contains original inputs, outputs and completion records;
`results/eos-selected-input-provenance.json` fixes their hashes.
The fitting scripts reproduce lattice constants, bulk moduli, pressure
derivatives and endpoint sensitivity. No atomic cohesive energies are included.
The scientific discussion and numerical limits are in the PCCP ESI.
