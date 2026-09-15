# Reproducing the PCCP Analysis

All commands below analyze existing data; none launches a new SCF calculation.
The audit uses Python's standard library. Fitting and figures additionally
require NumPy and Matplotlib. LaTeX needs the packages listed in the source
files and TeX Live's RSC bibliography style.

From the repository root:

```sh
python3 -B scripts/audit_pccp_package.py
python3 -B benchmarks/Goldzak12/scripts/verify_selected_package.py
python3 -B -m unittest discover -s benchmarks/Goldzak12/scripts -p 'test_*.py'
python3 -B -m unittest discover -s benchmarks/X23-mini/scripts -p 'test_*.py'
```

Regenerate numerical analysis and figures in a separate checkout, since the
immutable file audit intentionally detects changed artifacts:

```sh
python3 -B benchmarks/Goldzak12/scripts/fit_selected_eos.py
python3 -B benchmarks/Goldzak12/scripts/analyze_selected_eos.py
python3 -B benchmarks/Goldzak12/scripts/compare_selected_literature.py
python3 -B benchmarks/Goldzak12/scripts/build_paper_tables.py
python3 -B benchmarks/X23-mini/scripts/analyze_results.py
python3 -B scripts/basis_sensitivity.py --write
python3 -B benchmarks/X23-mini/scripts/analyze_basis_controls.py
python3 -B benchmarks/X23-mini/scripts/build_paper_tables.py
python3 -B paper/pccp/analyse_derivative_checks.py --report convergence/derivatives/periodic-water-20260909/validated-results.json --out convergence/derivatives/periodic-water-20260909 --tex paper/pccp/derivative-results.tex
python3 -B paper/pccp/build_submission_assets.py
cd paper/pccp
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary_information.tex
```

Numerical-library versions may affect final fit digits and rendering metadata.
The archived tables and figures are the manuscript snapshot.

## Original Inputs

* Molecular reactions: `production-p25/accepted-inputs/manifest.csv` maps all
  332 selected species-routes to frozen and actually executed inputs.
  `reaction-index.json` supplies stoichiometry and shared AE reuse.
* LC10: `results/eos-selected-input-provenance.json` identifies all 352
  selected executions; `protocol/paper-eos-selection.json` fixes the common
  solid/volume selection.
* Molecular crystals: each accepted record includes its original output,
  frozen/actual input, timing and execution evidence. Controls name their
  exact parent. The 1200 Ry check is crystal-only, not a paired lattice-energy
  convergence test.
  `results/lattice-energies.json` retains the original `complete_pairs` and
  records the manuscript selection separately as `paper_pairs`. All three AE
  entries use validated QZVPP crystal/molecule pairs. Twelve additional
  inputs and outputs in `results/basis-controls` reproduce the matched tests.
  `results/basis-convergence.json` retains the Urea quadrature/basis analysis.
* The additional CO2/NH3, molecular-reaction, and ice controls are indexed in
  `convergence/basis-sensitivity-20260915/index.json`. The independent
  `scripts/basis_sensitivity.py` verifies completed-block markers, energy,
  electron count, input/runtime hashes, and paired reaction stoichiometry.
  It reproduces the new SI tables and `assessment.json`. The QZ ethane
  outlier remains explicitly unresolved and never enters the 70-reaction
  production statistics. The ice average uses the same six phases at both
  basis levels and retains the documented molecular-reference caveat.
* Derivative and symmetry checks retain their exact inputs and outputs.
  Their component-specific scope is defined in the ESI.

To rerun CP2K, map machine-specific paths to your installation while preserving
the scientific settings and recorded source, basis, potential and model
versions. External CP2K data and Skala models must be obtained from their
cited distributions. Dataset provenance records their hashes. No WFN or
compiled binaries are distributed. An identical WFN restart cannot be
promised without those external files; final energies and analysis are
independently checkable from the included outputs.
