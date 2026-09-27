# Sampled Electronic Band Gaps

`dataset.json` contains 50 primary gaps for 28 materials: 20 pure AE,
23 GAPW-XC/GTH and seven mixed AE/GTH one-centre results. Primary selection
is fixed by core/basis availability, not by agreement with experiment.
All primary results use the specified triple-zeta bases.

The seven supplementary direct/one-centre comparisons reuse primary and
legacy controls and require 61 distinct archived executions in total.
`scripts/analyze_band_gaps.py` reconstructs extrema from `bands.bs`, checks
occupations and SCF evidence, and reproduces the exact reference populations.
Material directory names separate element symbols, so BAs and BaS remain
distinct on case-insensitive filesystems.

Thirty-nine primary outputs have fully readable bands. Eleven have finite,
reviewed frontier energies with explicitly indexed censored deep-core values.
The parser accepts only those reviewed indices and rejects missing frontier
bands. D/N/I indicate direct-minus-minimum splittings below 1e-6 eV,
from 1e-6 to below 1e-3 eV, and at least 1e-3 eV. These are reporting
thresholds, not physical accuracy bounds.

`paper-comparisons.json` records matched and maximum-available statistics.
Experimental and functional reference extracts include source attribution.
BAs and MgO sensitivities are separate from the unchanged baseline references.
Reference-review and terminal-outcome records remain in `provenance/`, outside
the numerical populations. Original execution records, contradictory status
fields and missing timing/RSS data are preserved, not repaired retrospectively.

No QZ band-gap campaign, active launch queue, model weights or WFN files are
included. No new electronic-structure calculation is needed for analysis.
