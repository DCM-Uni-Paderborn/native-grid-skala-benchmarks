#!/usr/bin/env python3
"""Validate paired Urea basis controls and select the reported finite-basis result."""
import copy
import hashlib
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
HA_KJ = 2625.4996394799

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def verified_run(level, phase):
    folder = ROOT / "results/basis-controls/urea" / level / phase
    record = json.loads((folder / "execution.json").read_text())
    validation = json.loads((folder / "validation.json").read_text())
    assert validation["status"] == "accepted_completed_control"
    assert validation["execution_sha256"] == sha(folder / "execution.json")
    assert record["case"]["id"] == level + "/" + phase
    assert record["clean_scf_markers"] and record["execution_exit"] == record["time_exit"] == 0
    assert record["runtime"]["revision"] == "37aedacbb33677307818301a47a09ae435712293"
    assert record["runtime"]["binary_sha256"] == "4c4cdcf4ac51510d1d876310efa724ae7efa8bdab7fdd124d803e0220fb770fb"
    assert record["model_sha256"] == "7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd"
    for name, digest in record["final_hashes"].items():
        path = folder / name
        if path.exists():
            assert sha(path) == digest, path
    assert sha(folder / "actual-input.inp") == record["case"]["input_sha256"]
    text = (folder / "actual-input.inp").read_text()
    expected_basis = level.split("-")[0].upper()
    assert text.count("BASIS_SET " + expected_basis + "-MOLOPT-PBE-ae") == 4
    assert text.count("RADIAL_GRID 200") == text.count("LEBEDEV_GRID 974") == 4
    assert "NATIVE_GRID_USE_CUDA .FALSE." in text and "EPS_SCF 5.0E-7" in text
    output = (folder / "output.out").read_text()
    start = output.rfind("PROGRAM STARTED AT")
    assert start >= 0
    block = output[start:]
    assert "SCF run converged" in block and "PROGRAM ENDED AT" in block
    assert "[ABORT]" not in block and "SCF run NOT converged" not in block
    energies = re.findall(r"ENERGY\| Total FORCE_EVAL.*?\[hartree\]\s+([-+\d.Ee]+)", block)
    energy = float(energies[-1])
    assert math.isfinite(energy) and energy == validation["energy_hartree"]
    electrons = 64 if phase == "solid" else 32
    assert [int(v) for v in re.findall(r"Number of electrons:\s+(\d+)", block)] == [electrons]
    assert validation["residual"] <= 5e-7
    assert re.search(r"Exit status:\s+0\s*$", (folder / "time.txt").read_text())
    def invariant(inp):
        inp = re.sub(r"(?m)^\s*(PROJECT_NAME|WFN_RESTART_FILE_NAME|SCF_GUESS)\s+[^\n]*\n", "", inp)
        return inp.replace("QZVPP-MOLOPT-PBE-ae", "TZVPP-MOLOPT-PBE-ae")
    return {
        "energy_hartree": energy, "electrons": electrons,
        "molecules": 2 if phase == "solid" else 1,
        "residual": validation["residual"], "scf_steps": validation["scf_steps"],
        "source_directory": str(folder.relative_to(ROOT)),
        "execution_sha256": sha(folder / "execution.json"),
        "validation_sha256": sha(folder / "validation.json"),
        "output_sha256": sha(folder / "output.out"),
        "input_sha256": sha(folder / "actual-input.inp"),
    }, invariant(text)

def assessment():
    result = {"system": "urea", "method": "gapw-ae", "pairs": {}, "hartree_to_kjmol": HA_KJ}
    invariant = {}
    for level in ("tzvpp-grid200-974", "qzvpp-grid200-974"):
        phases = {}
        for phase in ("solid", "molecule"):
            phases[phase], normalized = verified_run(level, phase)
            if phase in invariant:
                assert invariant[phase] == normalized, "Basis control changed another physical setting"
            invariant[phase] = normalized
        result["pairs"][level] = {
            "basis": level.split("-")[0].upper() + "-MOLOPT-PBE-ae",
            "radial_lebedev": [200, 974], "phases": phases,
            "lattice_energy_kjmol": (phases["solid"]["energy_hartree"]/2 - phases["molecule"]["energy_hartree"]) * HA_KJ,
        }
    baseline = (ROOT / "results/controls/urea-ae-cpu/validation.json")
    molecule = ROOT / "results/validated/gapw-ae/urea/molecule/validation.json"
    e0 = (json.loads(baseline.read_text())["energy_hartree"]/2
          - json.loads(molecule.read_text())["energy_hartree"]) * HA_KJ
    tz = result["pairs"]["tzvpp-grid200-974"]["lattice_energy_kjmol"]
    qz = result["pairs"]["qzvpp-grid200-974"]["lattice_energy_kjmol"]
    result.update(baseline_cpu_150_770_kjmol=e0, paired_grid_change_kjmol=tz-e0,
                  paired_basis_change_kjmol=qz-tz,
                  selected_paper_pair="qzvpp-grid200-974",
                  interpretation="Material finite-basis sensitivity, not a complete-basis or counterpoise extrapolation")
    return result

def selected_pairs(base_pairs):
    pairs = copy.deepcopy(base_pairs)
    report = assessment()
    selected = report["pairs"][report["selected_paper_pair"]]
    for pair in pairs:
        if (pair["method"], pair["system"]) == ("gapw-ae", "urea"):
            pair["superseded_base_lattice_energy_kjmol"] = pair["lattice_energy_kjmol"]
            pair.update(solid_energy_hartree=selected["phases"]["solid"]["energy_hartree"],
                        molecule_energy_hartree=selected["phases"]["molecule"]["energy_hartree"],
                        lattice_energy_kjmol=selected["lattice_energy_kjmol"],
                        signed_deviation_kjmol=selected["lattice_energy_kjmol"]-pair["dmc_kjmol"],
                        basis=selected["basis"], radial_lebedev=[200,974],
                        selection_reason="Larger AE basis following the paired grid/basis sensitivity test",
                        execution_sha256={p:r["execution_sha256"] for p,r in selected["phases"].items()},
                        selected_source_directories={p:r["source_directory"] for p,r in selected["phases"].items()})
    return pairs

def table(report):
    lines = [
        r"\begin{table}[htbp]", r"\centering", r"\small",
        r"\caption{Paired Urea AE basis and quadrature checks. Total energies are in hartree and lattice energies in kJ~mol$^{-1}$. Both phases use the same basis and grid. All four new executions satisfy the OT threshold $5\times10^{-7}$ and end normally. The QZVPP pair supplies the Urea AE entry in Table~\ref{tab:x23_total_si}.}",
        r"\label{tab:urea-basis}",
        r"\begin{tabular}{llrrr}", r"\toprule",
        r"Basis & Radial/Lebedev & $E_{\rm solid}$ & $E_{\rm molecule}$ & $E_{\rm latt}$\\", r"\midrule",
    ]
    lines.append(r"TZVPP & 150/770 & -450.589911695840 & -225.247483909549 & -124.637557\\")
    for level,pair in report["pairs"].items():
        a,b = pair["phases"]["solid"],pair["phases"]["molecule"]
        lines.append(f"{level.split('-')[0].upper()} & 200/974 & {a['energy_hartree']:.12f} & {b['energy_hartree']:.12f} & {pair['lattice_energy_kjmol']:.6f}" + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)

def main():
    result=assessment()
    (ROOT/"results/basis-convergence.json").write_text(json.dumps(result,indent=2)+"\n")
    (ROOT/"paper/x23-basis-si.tex").write_text(table(result))
    print(json.dumps(result,indent=2))

if __name__ == "__main__":
    main()

