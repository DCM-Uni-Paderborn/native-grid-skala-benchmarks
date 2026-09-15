# Targeted AE Basis Controls

The SI reports five completed reaction pairs from a seven-reaction study
selected before launch. These 24 species/basis executions cover CARBHB12/4,
G2RC/15, ADIM6/2, BHROT27/1, and RG18/1. All completed pairs are retained,
including the anomalous QZVPP ethane torsional barrier. Formal SCF convergence
does not resolve that outlier's electronic-state and quadrature questions.

These data do not replace, filter, or modify the fixed 70-reaction production
dataset. No representative QZVPP MAE is claimed. The incomplete W4-11 and
WATER27 pairs, queued follow-ups, and operational diagnostics are not archived
here because their energies are not reported in the paper or SI.

Run `python3 -B scripts/basis_sensitivity.py` from the repository root. It
uses the original reaction index, checks each digest and coefficient, verifies
the paired input invariants and completed-output evidence, and reproduces
the SI comparison against official references, GauXC AE, and PySCF.
