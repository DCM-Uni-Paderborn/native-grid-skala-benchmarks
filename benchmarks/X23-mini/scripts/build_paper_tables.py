"""Generate manuscript tables from the accepted paired-energy analysis only."""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SYSTEMS = ("CO2", "NH3", "urea")
METHODS = ("gapwxc-gth", "gapw-ae", "gapw-gth-direct", "gapw-gth-one-center")
LABELS = {
    "gapwxc-gth": r"\texttt{GAPW\_XC}--\GTH{}",
    "gapw-ae": r"\GAPW{}-AE",
    "gapw-gth-direct": r"\GAPW{}--\GTH{}, direct",
    "gapw-gth-one-center": r"\GAPW{}--\GTH{}, one center",
    "CO2": r"CO$_2$",
    "NH3": r"NH$_3$",
    "urea": "Urea",
}


def render(data):
    pairs = {(p["method"], p["system"]): p for p in data["complete_pairs"]}
    expected = {(method, system) for method in METHODS for system in SYSTEMS}
    if (data["accepted_base_cases"] != 24
            or len(data["complete_pairs"]) != 12 or set(pairs) != expected):
        raise ValueError("The manuscript table requires all 24 accepted base cases")
    main = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Fixed-geometry electronic lattice energies including D3(BJ), in kJ\,mol$^{-1}$ per molecule. Parentheses contain signed deviations from DMC, not numerical uncertainties. Negative deviations denote stronger binding. The DMC uncertainties are statistical; the reference's estimated total uncertainty is about 2~kJ\,mol$^{-1}$.\cite{DellaPia2024X23}}",
        r"\label{tab:x23_lattice}",
        r"\begin{tabular}{lrrr}",
        r"\toprule",
        r"Representation & CO$_2$ & NH$_3$ & Urea \\",
        r"\midrule",
    ]
    for method in METHODS:
        values = []
        for system in SYSTEMS:
            pair = pairs[method, system]
            values.append(f"{pair['lattice_energy_kjmol']:.3f} ({pair['signed_deviation_kjmol']:+.3f})")
        main.append(LABELS[method] + " & " + " & ".join(values) + r" \\")
    refs = [pairs[METHODS[0], s] for s in SYSTEMS]
    main += [
        r"\midrule",
        r"DMC\cite{DellaPia2024X23} & "
        + " & ".join(f"${p['dmc_kjmol']:.1f}\\pm{p['dmc_statistical_uncertainty_kjmol']:.1f}$" for p in refs)
        + r" \\",
        r"\bottomrule", r"\end{tabular}", r"\end{table*}", "",
    ]
    si = [
        r"\begin{table*}[t]", r"\centering",
        r"\caption{Accepted total energies including D3(BJ), in hartree, and paired lattice energies in kJ\,mol$^{-1}$. Crystal energies refer to the complete cell. Extra digits support reproducibility and are not an accuracy claim.}",
        r"\label{tab:x23_total_si}", r"\begin{tabular}{llrrrr}", r"\toprule",
        r"System & Representation & $Z$ & $E_{\rm solid}$ & $E_{\rm molecule}$ & $E_{\rm latt}$ \\",
        r"\midrule",
    ]
    for system in SYSTEMS:
        for method in METHODS:
            p = pairs[method, system]
            si.append(f"{LABELS[system]} & {LABELS[method]} & {p['molecules_per_cell']} & {p['solid_energy_hartree']:.12f} & {p['molecule_energy_hartree']:.12f} & {p['lattice_energy_kjmol']:.6f}" + r" \\")
    si += [r"\bottomrule", r"\end{tabular}", r"\end{table*}", ""]
    controls = {(c["case_id"], c["control_name"]): c for c in data["numerical_controls"]}
    grid = [c for c in data["paired_numerical_controls"] if c["system"] == "CO2" and c["method"] == "gapwxc-gth" and c["control_name"] == "atom-200-974"]
    if len(grid) != 1:
        raise ValueError("Missing or duplicate accepted paired quadrature control")
    rows = [(r"CO$_2$, \texttt{GAPW\_XC}--\GTH{}", "150/770 to 200/974", r"$\Delta E_{\rm latt}$, both phases", grid[0]["delta_lattice_energy_kjmol"])]
    for case, name, change, observable in (
        ("gapwxc-gth/CO2/solid", "cutoff-1200", "800 to 1200 Ry", r"$\Delta(E_{\rm solid}/Z)$ only"),
        ("gapwxc-gth/CO2/solid", "kmesh-5", r"$3^3$ to $5^3$", r"$\Delta E_{\rm latt}$"),
        ("gapwxc-gth/NH3/solid", "sym-tol-1e4", r"Symmetry tolerance $10^{-6}$ to $10^{-4}$", r"$\Delta E_{\rm latt}$"),
        ("gapwxc-gth/CO2/solid", "cpu-model", "CPU minus CUDA", r"$\Delta E_{\rm latt}$"),
        ("gapw-ae/urea/solid", "cpu-model", "CPU minus CUDA", r"$\Delta E_{\rm latt}$"),
    ):
        method, system, _ = case.split("/")
        value = controls[case, name]["delta_control_minus_parent_kjmol_per_molecule"]
        rows.append((LABELS[system] + ", " + LABELS[method], change, observable, value))
    si += [
        r"\begin{table*}[t]", r"\centering",
        r"\caption{Bounded numerical sensitivities in kJ\,mol$^{-1}$ per molecule. All unlisted scientific settings are held fixed. Seven accepted control executions produce six comparisons because the atom-grid refinement requires both phases.}",
        r"\label{tab:x23_controls_si}",
        r"\begin{tabular}{@{}p{0.27\textwidth}p{0.29\textwidth}p{0.25\textwidth}r@{}}",
        r"\toprule", r"System and representation & Change & Observable & Shift \\", r"\midrule",
    ]
    for label, change, observable, value in rows:
        si.append(f"{label} & {change} & {observable} & ${value:+.6f}$" + r" \\")
    si += [r"\bottomrule", r"\end{tabular}", r"\end{table*}", ""]
    return {"x23-tables-si.tex": "\n".join(si)}


def main():
    data = json.loads((ROOT / "results/lattice-energies.json").read_text())
    for name, text in render(data).items():
        path = ROOT / "paper" / name
        if not path.exists() or path.read_text() != text:
            path.write_text(text)
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
