# Molecular error distribution for the paper

`native-protocol-comparison.csv` contains the 70-reaction intersection used
unchanged in the manuscript and SI. Errors are calculated reaction energy minus
the official reference, in kcal/mol. The total reaction energies already include
the archived D3(BJ) contribution; no dispersion term or fitted offset is added
by this analysis.

From the repository root, reproduce the additional SI table and summary with:

```bash
python3 scripts/summarize_molecular_error_distribution.py
```

The script writes `error-distribution-summary.json` (including the source SHA-256)
and `molecular-distribution-si.tex`. It checks that the source contains exactly
70 unique reactions and reports strict threshold counts, the largest error's
share of the sum of squared errors, and the union of each protocol's three
largest absolute errors. It performs no electronic-structure calculation.

GX denotes GAPW_XC-GTH. HD and HOC denote mixed GAPW-AE/GTH with direct and
one-center fields, respectively. The set spans 36 parent subsets; these are
unweighted common-set statistics, not full dietGMTKN55 weighted scores. The
lower HD/HOC MAE is accompanied by fewer reactions below 1 kcal/mol absolute
error (42 versus 46 for GX), emphasizing the importance of the error tails.
Most HD/HOC energies share their AE contributions, so their near equality is
not 70 independent tests of heavy-element reconstruction.

The generated TeX is included by the SI of the
[Native-grid Skala manuscript](https://www.overleaf.com/project/6a860f5e29d899a1050644a1).
