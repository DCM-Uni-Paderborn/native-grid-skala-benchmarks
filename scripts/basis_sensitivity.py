#!/usr/bin/env python3
"""Reproduce the reported paired basis controls from exact completed executions."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import statistics

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "convergence/basis-sensitivity-20260915"
HA_KJ = 2625.4996394799
HA_KCAL = 627.5094740631
BINARIES = {"rosi": "3bd708ff4769f7fc98a15d9c692c44c4586a9d1d3351fe84e3842d1b6e5f20c9",
            "terok": "4c4cdcf4ac51510d1d876310efa724ae7efa8bdab7fdd124d803e0220fb770fb"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verified(case):
    folder = ROOT / case["target"]
    for name, digest in case["files_sha256"].items():
        assert sha(folder / name) == digest, (folder, name)
    record = json.loads((folder / "execution.json").read_text())
    validation = case["validation"]
    assert sha(folder / "execution.json") == validation["execution_sha256"]
    assert sha(folder / "output.out") == validation["output_sha256"]
    assert record["execution_exit"] == record["time_exit"] == 0
    assert record["runtime"]["revision"] == "37aedacbb33677307818301a47a09ae435712293"
    assert record["runtime"]["binary_sha256"] == BINARIES[case["host"]]
    assert record["model_sha256"] == "7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd"
    assert sha(folder / "actual-input.inp") == record["case"]["input_sha256"]
    inp = (folder / "actual-input.inp").read_text()
    assert set(re.findall(r"(?m)^\s*BASIS_SET\s+(\S+)", inp)) == {case["level"].upper() + "-MOLOPT-PBE-ae"}
    assert set(re.findall(r"(?m)^\s*POTENTIAL\s+(\S+)", inp)) == {"ALL"}
    assert set(float(x) for x in re.findall(r"(?m)^\s*EPS_SCF\s+(\S+)", inp)) == {5e-7}
    assert "GAPW_ACCURATE_XCINT .TRUE." in inp and "TYPE DFTD3(BJ)" in inp
    block = (folder / "output.out").read_text().rsplit("PROGRAM STARTED", 1)
    assert len(block) == 2
    block = block[1]
    assert "SCF run converged" in block and "PROGRAM ENDED" in block
    assert "[ABORT]" not in block and "SCF run NOT converged" not in block
    energy = float(re.findall(r"ENERGY\| Total FORCE_EVAL.*?\[hartree\]\s+([-+\d.Ee]+)", block)[-1])
    assert math.isfinite(energy) and energy < 0 and energy == validation["energy_hartree"]
    expected = record["case"].get("expected_spin_electrons", [record["case"]["expected_electrons"]])
    assert list(map(int, re.findall(r"Number of electrons:\s+(\d+)", block))) == expected
    assert validation["residual"] <= 5e-7 and record["peak_rss_kib"] > 0 and record["elapsed_wall_s"] > 0
    assert re.search(r"Exit status:\s+0\s*$", (folder / "time.txt").read_text())
    return {"energy_hartree": energy, "molecules": record["case"]["molecules"],
            "input_sha256": sha(folder / "actual-input.inp"),
            "execution_sha256": sha(folder / "execution.json"),
            "source_directory": case["target"], "residual": validation["residual"],
            "scf_steps": validation["scf_steps"]}, inp, record["case"]


def invariant(inp):
    inp = re.sub(r"(?m)^\s*(PROJECT_NAME|SCF_GUESS|WFN_RESTART_FILE_NAME)\s+[^\n]*\n", "", inp)
    return inp.replace("QZVPP-MOLOPT-PBE-ae", "TZVPP-MOLOPT-PBE-ae")


def assess():
    index = json.loads((DATA / "index.json").read_text())
    assert len(index["cases"]) == 48
    runs = {}
    for case in index["cases"]:
        key = (case["kind"], case["system"], case["level"], case["phase"])
        assert key not in runs
        runs[key] = verified(case)
    result = {"snapshot_utc": index["snapshot_utc"], "crystals": {}, "molecular_reactions": [], "ice": [], "ice_symmetry": []}
    for system in ("CO2", "NH3"):
        pairs = {}
        for level in ("tzvpp", "qzvpp"):
            phases = {phase: runs[("crystal", system, level, phase)][0] for phase in ("solid", "molecule")}
            pairs[level] = {"phases": phases, "lattice_energy_kjmol": (phases["solid"]["energy_hartree"] / phases["solid"]["molecules"] - phases["molecule"]["energy_hartree"]) * HA_KJ}
        for phase in ("solid", "molecule"):
            assert invariant(runs[("crystal", system, "tzvpp", phase)][1]) == invariant(runs[("crystal", system, "qzvpp", phase)][1])
        pairs["basis_shift_kjmol"] = pairs["qzvpp"]["lattice_energy_kjmol"] - pairs["tzvpp"]["lattice_energy_kjmol"]
        result["crystals"][system] = pairs
    reaction_index = json.loads((ROOT / "benchmarks/dietGMTKN55/production-p25/paper-common-70/reaction-index.json").read_text())
    for expected in index["molecular_reactions"]:
        reaction = next(r for r in reaction_index["reactions"] if (r["subset"], r["reaction_id"]) == (expected["subset"], expected["reaction_id"]))
        values = {}
        for level in ("tzvpp", "qzvpp"):
            energy = 0.0
            for species in reaction["species"]:
                system = f"{reaction['subset']}-{reaction['reaction_id']}/{species['name']}"
                run, inp, original = runs[("molecular", system, level, "molecule")]
                assert species["digest"] == original["digest"]
                assert invariant(inp) == invariant(runs[("molecular", system, "tzvpp", "molecule")][1])
                energy += species["count"] * run["energy_hartree"]
            values[level] = energy * HA_KCAL
            assert abs(values[level] - expected["tz_kcal_mol" if level == "tzvpp" else "qz_kcal_mol"]) < 1e-7
        result["molecular_reactions"].append(dict(expected, scientific_status="unresolved_QZ_outlier" if reaction["subset"] == "BHROT27" else "paired_basis_control"))
    for system, (reference, uncertainty) in index["ice_reference"]["values"].items():
        row = {"phase": system, "dmc_kjmol": reference, "dmc_statistical_uncertainty_kjmol": uncertainty}
        for level in ("tzvpp", "qzvpp"):
            solid = runs[("ice", system, level, "solid")][0]
            molecule = runs[("ice", "water-ps", level, "molecule")][0]
            row[level] = (solid["energy_hartree"] / solid["molecules"] - molecule["energy_hartree"]) * HA_KJ
        result["ice"].append(row)
    for system in ("IX", "VII"):
        old = runs[("ice", system, "tzvpp", "solid")][0]
        new = runs[("ice-symmetry", system, "tzvpp", "solid")][0]
        result["ice_symmetry"].append({"phase": system, "spglib_minus_original_kjmol": (new["energy_hartree"] - old["energy_hartree"]) / old["molecules"] * HA_KJ})
    result["ice_mae"] = {level: statistics.mean(abs(r[level] - r["dmc_kjmol"]) for r in result["ice"]) for level in ("tzvpp", "qzvpp")}
    result["ice_reference"] = index["ice_reference"]
    return result


def tables(report):
    lines = [r"\begin{table}[htbp]", r"\centering\small",
             r"\caption{Targeted paired molecular basis controls, in kcal~mol$^{-1}$. Both bases use the same geometry, Hamiltonian, quadrature, and dispersion. The starred QZVPP ethane barrier is an unresolved outlier despite SCF convergence, not a recommended replacement. These selected reactions do not define a revised 70-reaction MAE.}",
             r"\label{tab:molecular-basis-controls}", r"\begin{tabular}{lrrrrr}", r"\toprule",
             r"Reaction & Reference & TZVPP & QZVPP & GauXC AE & PySCF\\", r"\midrule"]
    for r in report["molecular_reactions"]:
        label = f"{r['subset']}/{r['reaction_id']}"
        star = r"$^{*}$" if r["scientific_status"] == "unresolved_QZ_outlier" else ""
        lines.append(f"{label} & {r['reference_kcal_mol']:.3f} & {r['tz_kcal_mol']:.3f} & {r['qz_kcal_mol']:.3f}{star} & {r['gauxc_ae_kcal_mol']:.3f} & {r['pyscf_kcal_mol']:.3f}" + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    ice = [r"\begin{table}[htbp]", r"\centering\small",
           r"\caption{Preliminary six-phase ice basis comparison. Electronic lattice energies and signed QZVPP--DMC residuals are in kJ~mol$^{-1}$ per water molecule. Parentheses give DMC statistical uncertainties. The averages concern these same six phases only, not the complete DMC-ICE13 benchmark.}",
           r"\label{tab:ice-basis}", r"\begin{tabular}{lrrrr}", r"\toprule",
           r"Phase & TZVPP & QZVPP & DMC & QZVPP--DMC\\", r"\midrule"]
    for r in report["ice"]:
        ice.append(f"{r['phase']} & {r['tzvpp']:.3f} & {r['qzvpp']:.3f} & {r['dmc_kjmol']:.2f}({round(r['dmc_statistical_uncertainty_kjmol']*100):d}) & {r['qzvpp']-r['dmc_kjmol']:+.3f}" + r"\\")
    ice += [r"\midrule", f"MAE & {report['ice_mae']['tzvpp']:.3f} & {report['ice_mae']['qzvpp']:.3f} & -- & --" + r"\\",
            r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return {"molecular-basis-table-si.tex": "\n".join(lines), "ice-basis-table-si.tex": "\n".join(ice)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = assess()
    if args.write:
        (DATA / "assessment.json").write_text(json.dumps(report, indent=2) + "\n")
        for name, content in tables(report).items():
            (ROOT / "paper/pccp" / name).write_text(content)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
