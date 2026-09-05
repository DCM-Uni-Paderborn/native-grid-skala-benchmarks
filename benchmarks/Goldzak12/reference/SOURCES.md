# Reference provenance

## Goldzak12

T. Goldzak, X. Wang, H.-Z. Ye, and T. C. Berkelbach, "Accurate
thermochemistry of covalent and ionic solids from spin-component-scaled
MP2," *J. Chem. Phys.* **157**, 174112 (2022), DOI
[`10.1063/5.0119633`](https://doi.org/10.1063/5.0119633).

`goldzak2022.csv` transcribes Tables II-IV of the article.  The experimental
lattice constants and cohesive energies are the zero-point-corrected values
used by the authors.  The primary bulk-modulus reference follows their Table
III; alternate experimental values are retained in `protocol/systems.csv` for
BN, BP, and MgO.

The same article lists aggregate PBE, PBEsol, and SCAN errors. These numbers
coincide with the 44-solid statistics of Tran, Stelzl and Blaha, JCP 144,
204120 (2016), DOI `10.1063/1.4948636`, Table II. They are not verified as a
Goldzak12 recomputation and must not be labeled as matched ten-solid errors.
Those values are retained as contextual in `aggregate_literature.csv`; system-resolved DFT
numbers from calculations with other numerical protocols are not mixed into
the primary reference table.

## Periodic GFN2 manuscript

V. Alizadeh et al., "Periodic GFN2-xTB in CP2K: Multipolar Ewald
Electrostatics, k-Point Sampling, and Transferability Benchmarks for Solids,"
accepted in *J. Chem. Theory Comput.*, manuscript `ct-2026-010347`.

`periodic_gfn2_lc10.csv` reproduces the GFN1-xTB and GFN2-xTB lattice constants
and cohesive energies from Table S3 for the ten systems used there.  The SI
records that MgO and LiH did not bracket valid equation-of-state minima before
the SCC branch collapsed.  These two systems are therefore retained with an
explicit non-result status, and no numerical values are inferred for them.

## Additional aggregate context

Y. Mo, R. Car, V. N. Staroverov, G. E. Scuseria, and J. Tao, "Assessment of
the Tao-Mo nonempirical semilocal density functional in applications to solids
and surfaces," *Phys. Rev. B* **95**, 035118 (2017), DOI
[`10.1103/PhysRevB.95.035118`](https://doi.org/10.1103/PhysRevB.95.035118).

The Tao-Mo row is deliberately labeled as a different-set aggregate.  It is
useful methodological context in the SI, but it is not scored together with
Goldzak12 and does not enter any native-Skala error statistic.

## System-resolved literature added for the closed campaign

`mo2017_selected.csv` retains the self-consistent TM columns of Tables II/III
for the eight overlapping solids (no AlN or BN). Values were checked against
the column layout of the [published author-hosted PDF](https://repository.rice.edu/server/api/core/bitstreams/42110573-5f78-4625-8490-5a1a4ba64fb3/content).
The original article uses different experimental conventions; the new matched
statistics use the same Goldzak reference for both native and literature data.
For example, its TM B0 values for BP and AlP are 171.5 and 89.3 GPa, not the
adjacent HSE06 values 178.4 and 94.3 GPa. Missing columns must not be collapsed.

G.-X. Zhang, A. M. Reilly, A. Tkatchenko, and M. Scheffler, "Performance of
various density-functional approximations for cohesive properties of 64 bulk
solids," New J. Phys. 20, 063020 (2018), DOI
[`10.1088/1367-2630/aac7f0`](https://doi.org/10.1088/1367-2630/aac7f0).
`zhang2018_selected.csv` contains the nine overlapping solids (AlN absent),
six DFAs, and both static and ZPE-included values from SI Tables III/V.
Use `Uncorr.` for the electronic EOS comparison against the Goldzak reference;
`Corr.` already includes zero-point motion and must not be corrected twice.
The SI states that SCAN was evaluated non-self-consistently on PBE orbitals.
`import_zhang_reference.py` reproduces the extraction; the URL and PDF SHA256
are recorded in `zhang2018_selected-provenance.json`.

L. Schimka, J. Harl, and G. Kresse, "Improved hybrid functional for solids:
The HSEsol functional," JCP 134, 024116 (2011), DOI
[`10.1063/1.3524336`](https://doi.org/10.1063/1.3524336).
`schimka2011_selected.csv` transcribes Table III for the nine overlapping
solids (MgS absent), including PBE, PBEsol, HSE06 and HSEsol. All are static
theoretical lattice constants. The paper's experimental lattice constants
are not substituted for the common Goldzak references. Its separate B0
supplement was not extracted; these cells remain blank.

Scope: the comparison covers the retrievable, relevant benchmark sources
above, not every published result for these compounds. Missing source values
are never reconstructed from aggregate errors. Other phases, cohesive
energies without complete native atomic references, and unmatched-set
aggregate rankings are outside this structural comparison.
