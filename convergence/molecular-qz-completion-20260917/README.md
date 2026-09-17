# Additional Molecular QZ Controls

This package supports the two additional reaction-resolved basis controls in
the PCCP ESI: CARBHB12/11 and G21EA/25. All five final species are converged,
unsmeared QZVPP calculations at the frozen production geometry, Hamiltonian,
quadrature, cutoffs, charge, spin and dispersion settings. The orbital basis
and SCF strategy are the only changes relative to the archived TZVPP inputs.
The uniform 70-reaction statistics are unchanged.

For 3CL and EA_25, `precursor/` preserves the initial, nonconverged
diagonalization calculation with modified Broyden mixing. Its final orbitals
seeded `final/`, which converged with OT. No precursor energy is used in an
energy difference. These two precursor records are included only to reproduce
the initialization of the reported final results, not as benchmark results.
The original `diagnostic_only` flags and execution records are preserved
without retroactive rewriting. The separate assessment establishes their use
as reported numerical controls, not replacements for uniform production data.

The exact inputs, outputs, timings and process/runtime provenance are archived.
Remote executable, library, source, model, data and final-file hashes were
rechecked before retrieval. The execution records retain the exact source and
final wavefunction hashes. Large wavefunctions and binaries are not duplicated
in Git. Reproduction of the recovery path requires the recorded precursor
calculation followed by the unsmeared final input using its checkpoint as
`source.wfn`. Absolute runtime paths must be adapted for another installation.

Run `python3 -B scripts/molecular_completion_controls.py --write` from the
repository root to verify the archived inputs and execution chain and
recalculate the two reactions. This analysis neither launches CP2K nor uses
the literature reference to select an electronic state. The GauXC and PySCF
comparison values come from the existing reaction-matched source dataset.
