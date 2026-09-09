# Reference Provenance for the PCCP LC10 Selection

Only the ten native solids and the structural quantities used in the ESI are
retained. Cohesive energies and unmatched-set aggregate tables are excluded.

* `goldzak2022.csv`: Tables II and III of T. Goldzak, X. Wang, H.-Z. Ye and
  T. C. Berkelbach, JCP 157, 174112 (2022),
  [DOI](https://doi.org/10.1063/5.0119633). The experimental lattice constants
  are the authors' zero-point-corrected values. The primary bulk-modulus
  reference follows Table III. Its alternate BN/BP/MgO values,
  410.2/168.0/169.8 GPa, are explicit in `compare_selected_literature.py`.
* `periodic_gfn2_lc10.csv`: structural values from Table S3 of the accepted
  author manuscript by V. Alizadeh et al., "Periodic GFN2-xTB in CP2K:
  Multipolar Ewald Electrostatics, k-Point Sampling, and Transferability
  Benchmarks for Solids", JCTC manuscript ct-2026-010347. Only the overlap is
  retained; blank/non-result values are not inferred.
* `mo2017_selected.csv`: self-consistent TM columns in Tables II/III of
  Y. Mo et al., PRB 95, 035118 (2017),
  [DOI](https://doi.org/10.1103/PhysRevB.95.035118). Eight solids overlap
  (AlN and BN absent). The matched errors use Goldzak's references, not the
  different experimental conventions in this source.
* `zhang2018_selected.csv`: nine overlapping solids (AlN absent), six DFAs,
  from SI Tables III/V of G.-X. Zhang et al., NJP 20, 063020 (2018),
  [DOI](https://doi.org/10.1088/1367-2630/aac7f0). Source URL and PDF SHA-256
  are in `zhang2018_selected-provenance.json`. Static "Uncorr." values are
  used for the electronic EOS comparison; "Corr." includes zero-point
  motion and must not be corrected twice. SCAN was evaluated on PBE orbitals.
* `schimka2011_selected.csv`: Table III of L. Schimka, J. Harl and G. Kresse,
  JCP 134, 024116 (2011),
  [DOI](https://doi.org/10.1063/1.3524336). Nine solids overlap (MgS absent);
  PBE, PBEsol, HSE06 and HSEsol static lattice constants are included.
  Bulk moduli from the separate supplement were not extracted.

Only matched system-resolved values enter native-versus-literature scores.
Other phases, cohesive energies without native atomic references and
different-set aggregate rankings are outside this structural comparison.
