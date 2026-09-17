# Targeted AE Basis Controls

The SI reports six completed reaction pairs from a seven-reaction study
selected before launch. These 30 species/basis executions cover CARBHB12/4,
G2RC/15, ADIM6/2, BHROT27/1, RG18/1 and WATER27/20. All completed pairs are
retained, including the initial high-energy QZVPP ethane solution. Three
additional ethane state/grid controls document a lower-energy restart and
show that a finer quadrature with atomic initial guesses alone does not
resolve the large initial torsional barrier.

These data do not replace, filter, or modify the fixed 70-reaction production
dataset. No representative QZVPP MAE is claimed. The incomplete W4-11
pair, queued follow-ups, and unreported operational diagnostics are not archived
here because their energies are not reported in the paper or SI.

Run `python3 -B scripts/basis_sensitivity.py` from the repository root. It
uses the original reaction index, checks each digest and coefficient, verifies
the paired input invariants and completed-output evidence, and reproduces
the SI comparison against official references, GauXC AE, and PySCF.
