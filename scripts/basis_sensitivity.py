#!/usr/bin/env python3
"""Reproduce the reported paired basis controls from exact completed executions."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys

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
    threshold = 5e-6 if case["kind"] == "solid" else 5e-7
    assert set(float(x) for x in re.findall(r"(?m)^\s*EPS_SCF\s+(\S+)", inp)) == {threshold}
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
    assert validation["residual"] <= threshold and record["peak_rss_kib"] > 0 and record["elapsed_wall_s"] > 0
    assert re.search(r"Exit status:\s+0\s*$", (folder / "time.txt").read_text())
    return {"energy_hartree": energy, "molecules": record["case"].get("molecules"),
            "input_sha256": sha(folder / "actual-input.inp"),
            "execution_sha256": sha(folder / "execution.json"),
            "source_directory": case["target"], "residual": validation["residual"],
            "scf_steps": validation["scf_steps"]}, inp, record["case"]


def invariant(inp):
    inp = re.sub(r"(?m)^\s*(PROJECT_NAME|SCF_GUESS|WFN_RESTART_FILE_NAME)\s+[^\n]*\n", "", inp)
    return inp.replace("QZVPP-MOLOPT-PBE-ae", "TZVPP-MOLOPT-PBE-ae")


def assess():
    index = json.loads((DATA / "index.json").read_text())
    assert len(index["cases"]) == 56
    runs = {}
    for case in index["cases"]:
        key = (case["kind"], case["system"], case["level"], case["phase"])
        assert key not in runs
        runs[key] = verified(case)
    result = {"snapshot_utc": index["snapshot_utc"], "crystals": {}, "ice": [], "ice_symmetry": []}
    for system in ("CO2", "NH3"):
        pairs = {}
        for level in ("tzvpp", "qzvpp"):
            phases = {phase: runs[("crystal", system, level, phase)][0] for phase in ("solid", "molecule")}
            pairs[level] = {"phases": phases, "lattice_energy_kjmol": (phases["solid"]["energy_hartree"] / phases["solid"]["molecules"] - phases["molecule"]["energy_hartree"]) * HA_KJ}
        for phase in ("solid", "molecule"):
            assert invariant(runs[("crystal", system, "tzvpp", phase)][1]) == invariant(runs[("crystal", system, "qzvpp", phase)][1])
        pairs["basis_shift_kjmol"] = pairs["qzvpp"]["lattice_energy_kjmol"] - pairs["tzvpp"]["lattice_energy_kjmol"]
        result["crystals"][system] = pairs
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
    ih = next(r for r in result["ice"] if r["phase"] == "Ih")
    result["ice_relative"] = [{"phase": r["phase"], **{k: r[k] - ih[k] for k in ("tzvpp", "qzvpp", "dmc_kjmol")}}
                              for r in result["ice"]]
    result["ice_relative_mae"] = {level: statistics.mean(abs(r[level] - r["dmc_kjmol"])
        for r in result["ice_relative"] if r["phase"] != "Ih") for level in ("tzvpp", "qzvpp")}
    result["ice_error_statistics"] = {
        "absolute": ice_error_statistics(result["ice"]),
        "relative_ih": ice_error_statistics([r for r in result["ice_relative"] if r["phase"] != "Ih"]),
    }
    result["ice_reference"] = index["ice_reference"]
    result["solid_basis_controls"] = solid_assessment(runs)
    return result


def ice_error_statistics(rows):
    result = {}
    for level in ("tzvpp", "qzvpp"):
        errors = [r[level] - r["dmc_kjmol"] for r in rows]
        worst = max(range(len(rows)), key=lambda i: abs(errors[i]))
        result[level] = {
            "count": len(rows), "mae_kjmol": statistics.mean(map(abs, errors)),
            "mse_kjmol": statistics.mean(errors),
            "rmse_kjmol": math.sqrt(statistics.mean(e * e for e in errors)),
            "maximum_absolute_error_kjmol": abs(errors[worst]),
            "maximum_error_phase": rows[worst]["phase"],
            "negative_error_count": sum(e < 0 for e in errors),
        }
    return result


def solid_assessment(runs):
    sys.path.insert(0, str(ROOT / "benchmarks/Goldzak12/scripts"))
    from fit_eos import fit_curve, HA_PER_A3_TO_GPA

    def fit(rows):
        parameters, _, rms = fit_curve(rows)
        return {"a0_angstrom": float(parameters[1] ** (1 / 3)),
                "B0_GPa": float(parameters[2] * HA_PER_A3_TO_GPA),
                "Bprime": float(parameters[3]), "rms_meV_atom": float(rms * 27211.386245988 / 8),
                "bracketed": bool(min(float(r["volume_A3"]) for r in rows) < parameters[1]
                                  < max(float(r["volume_A3"]) for r in rows))}

    root = ROOT / "benchmarks/Goldzak12"
    with (root / "results/eos-selected-source.csv").open() as f:
        old = list(csv.DictReader(f))
    results = []
    for solid in ("Si", "C", "MgO"):
        points = []
        for point in ("v01", "v03", "v05", "v07", "v09"):
            run, inp, case = runs[("solid", solid, "qzvpp", point)]
            original = next(r for r in old if (r["solid"], r["method"], r["point"]) == (solid, "gapw-ae", point))
            assert float(original["volume_A3"]) == case["volume_A3"]
            assert float(original["energy_Ha"]) == case["old_tz_energy_hartree"]
            original_input = root / "accepted-inputs/gapw-ae" / solid / point / "actual-input.inp"
            assert sha(original_input) == case["source_input_sha256"]
            assert solid_invariant(inp) == solid_invariant(original_input.read_text())
            points.append({"point": point, "volume_A3": case["volume_A3"],
                           "tzvpp": float(original["energy_Ha"]), "qzvpp": run["energy_hartree"]})
        fits = {}
        for level in ("tzvpp", "qzvpp"):
            rows = [{"volume_A3": p["volume_A3"], "energy_Ha": p[level]} for p in points]
            fits[level] = fit(rows)
            assert fits[level]["bracketed"]
        controls = [v for k, v in runs.items() if k[:3] == ("solid", solid, "tzvpp")]
        assert len(controls) == 1
        control, _, case = controls[0]
        results.append({"solid": solid, "points": points, "fits": fits,
                        "tz_runtime_reproduction_error_hartree": control["energy_hartree"] - case["old_tz_energy_hartree"],
                        "scope": "Matched five-volume control, not a complete QZVPP ten-solid benchmark"})
    return results


def solid_invariant(text):
    text = re.sub(r"(?m)^(\s*PROJECT_NAME)\s+\S+", r"\1 PROJECT", text)
    text = re.sub(r"(?m)^\s*(SYMMETRY_BACKEND|SYMMETRY_REDUCTION_METHOD|INVERSION_SYMMETRY_ONLY)\s+[^\n]+\n", "", text)
    text = re.sub(r"BASIS_SET (?:TZVPP|QZVPP)-MOLOPT-PBE-ae", "BASIS_SET AE_BASIS", text)
    text = re.sub(r"(?m)^\s*&OUTER_SCF\s*\n.*?^\s*&END OUTER_SCF\s*\n", "", text, flags=re.S)
    return re.sub(r"(?m)^(\s*MAX_SCF)\s+\d+", r"\1 MAXIMUM", text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Update the derived assessment JSON')
    args = parser.parse_args()
    report = assess()
    if args.write:
        (DATA / 'assessment.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
