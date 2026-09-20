"""Reproduce the matched ten-volume Si basis control from archived executions."""
import argparse
import csv
from datetime import datetime
import json
import math
from pathlib import Path
import re
import sys

from basis_sensitivity import ROOT, BINARIES, sha, solid_invariant, verified
from molecular_completion_controls import parse_output

sys.path.insert(0, str(ROOT / "benchmarks/Goldzak12/scripts"))
from fit_eos import fit_curve, HA_PER_A3_TO_GPA

DATA = ROOT / "convergence/si-eos-extension-20260920"
SOLIDS = ROOT / "benchmarks/Goldzak12"


def rosi_path(path):
    # ROSI resolves this verified home symlink in /proc executable/library paths.
    return path.replace("/home/kuehne88/", "/data/home2/kuehne88/", 1)


def verify_new(case):
    directory = ROOT / case["target"]
    record = json.loads((directory / "execution.json").read_text())
    for name, digest in case["files_sha256"].items():
        assert sha(directory / name) == digest, (directory, name)
        if name != "execution.json":
            assert record["final_hashes"][name] == digest
    assert set(case["slurm"]["job_and_steps"].values()) == {"COMPLETED|0:0"}
    assert record["preflight_exit"] == record["execution_exit"] == record["time_exit"] == 0
    assert record["peak_rss_kib"] > 0 and record["elapsed_wall_s"] > 0
    assert re.search(r"Exit status:\s+0\s*$", (directory / "time.txt").read_text())
    runtime = record["runtime"]
    assert runtime["revision"] == "37aedacbb33677307818301a47a09ae435712293"
    assert runtime["binary_sha256"] == BINARIES["rosi"]
    assert set(runtime["cp2k_shared_libraries"].values()) == {"f9001e8375239e34dd68fe8a928a295af2ffebc5f8e5fa3c8aa437fb43889d88"}
    assert record["model_sha256"] == "7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd"
    original = record["case"]
    assert original["charge"] == 0 and original["multiplicity"] == 1
    inp = (directory / "actual-input.inp").read_text()
    assert sha(directory / "actual-input.inp") == original["input_sha256"]
    assert set(re.findall(r"(?m)^\s*EPS_SCF\s+(\S+)", inp)) == {"5e-6"}
    assert "&OUTER_SCF" not in inp and "&SMEAR" not in inp
    assert "SCHEME MONKHORST-PACK 4 4 4" in inp
    assert "BASIS_SET QZVPP-MOLOPT-PBE-ae" in inp
    assert "UKS .FALSE." in inp
    parsed = parse_output((directory / "output.out").read_text())
    assert parsed["converged"] and parsed["ended"] and not parsed["abort"]
    assert parsed["electrons"] == [112] and parsed["steps"][-1]["residual"] <= 5e-6
    assert math.isfinite(parsed["energy_hartree"])
    block = (directory / "output.out").read_text().rsplit("PROGRAM STARTED", 1)[1]
    assert re.findall(r"BRILLOUIN\| Symmetry backend\s+(\S+)", block)[-1] == "SPGLIB"
    assert int(re.findall(r"BRILLOUIN\| List of Kpoints.*?\s+(\d+)\s*$", block, re.M)[-1]) == 10
    assert not any("inversion" in line.lower() and "fallback" in line.lower() for line in block.splitlines())
    snapshots = json.loads((directory / "stable-launch.json").read_text())["snapshots"]
    assert len(snapshots) == 2 and len(record["affinity"]) == 32
    assert (datetime.fromisoformat(snapshots[1]["observed"]) - datetime.fromisoformat(snapshots[0]["observed"])).total_seconds() >= 2
    for snap in snapshots:
        for key in ("pid", "start_ticks", "cwd", "executable"):
            assert snap[key] == snapshots[0][key]
        assert snap["cwd"] == record["cwd"] and snap["executable"] == rosi_path(runtime["binary"])
        assert set(snap["loaded_cp2k_libraries"]) == {rosi_path(p) for p in runtime["cp2k_shared_libraries"]}
        assert snap["affinity"] == record["affinity"]
        assert all(set(mask) <= set(record["affinity"]) for mask in snap["thread_affinities"].values())
    return parsed, inp, original


def fit(points, level):
    rows = [{"volume_A3": p["volume_A3"], "energy_Ha": p[level]} for p in points]
    parameters, _, rms = fit_curve(rows)
    return dict(n=len(rows), a0_A=float(parameters[1] ** (1 / 3)),
                B0_GPa=float(parameters[2] * HA_PER_A3_TO_GPA), Bprime=float(parameters[3]),
                rms_meV_atom=float(rms * 27211.386245988 / 8),
                bracketed=bool(min(r["volume_A3"] for r in rows) < parameters[1] < max(r["volume_A3"] for r in rows)))


def assess():
    index = json.loads((DATA / "index.json").read_text())
    old_index = json.loads((ROOT / "convergence/basis-sensitivity-20260915/index.json").read_text())
    old = {c["phase"]: c for c in old_index["cases"] if (c["kind"], c["system"], c["level"]) == ("solid", "Si", "qzvpp")}
    new = {c["point"]: c for c in index["cases"]}
    assert set(old) == {"v01", "v03", "v05", "v07", "v09"}
    assert set(new) == {"v02", "v04", "v06", "v08", "v10"}
    with (SOLIDS / "results/eos-selected-source.csv").open() as stream:
        tz = {r["point"]: r for r in csv.DictReader(stream) if (r["method"], r["solid"]) == ("gapw-ae", "Si")}
    points = []
    for point in sorted(tz):
        parsed, inp, original = verified(old[point]) if point in old else verify_new(new[point])
        source = SOLIDS / "accepted-inputs/gapw-ae/Si" / point / "actual-input.inp"
        assert sha(source) == original["source_input_sha256"]
        assert solid_invariant(inp) == solid_invariant(source.read_text())
        assert float(tz[point]["volume_A3"]) == original["volume_A3"]
        assert float(tz[point]["energy_Ha"]) == original["old_tz_energy_hartree"]
        points.append(dict(point=point, volume_A3=original["volume_A3"], tzvpp=float(tz[point]["energy_Ha"]), qzvpp=parsed["energy_hartree"]))
    assert len(points) == 10
    windows = {"five_odd": points[::2], "ten_all": points, "nine_omit_smallest": points[1:],
               "nine_omit_largest": points[:-1], "eight_omit_both": points[1:-1]}
    fits = {name: {level: fit(rows, level) for level in ("tzvpp", "qzvpp")} for name, rows in windows.items()}
    assert all(f["bracketed"] for window in fits.values() for f in window.values())
    return dict(scope=index["scope"], points=points, fits=fits, new_scf_points=5, reused_qz_points=5,
                included_in_uniform_statistics=False, Bprime_status="exploratory; fit-window dependent")


def table(report):
    labels = {"five_odd": "Original five", "ten_all": "All ten", "nine_omit_smallest": "Omit smallest", "nine_omit_largest": "Omit largest", "eight_omit_both": "Omit both endpoints"}
    lines = [r"\begin{table}[htbp]", r"\centering\small",
             r"\caption{Matched Si EOS fits with different sampling and fit windows. Lattice constants are in \AA{}, bulk moduli in GPa, and fit RMS residuals in meV per atom. All minima are bracketed. Five new QZVPP points complete the ten-volume curve; the original five QZVPP points and all TZVPP energies are reused. The uniform ten-solid statistics are unchanged.}",
             r"\label{tab:si-eos-extension}", r"\begin{tabular}{lcrrrr}"]
    lines += [r"\toprule", r"Window & Basis & $N$ & $a_0$ & $B_0$ & RMS\\", r"\midrule"]
    for window, levels in report["fits"].items():
        for level, result in levels.items():
            lines.append(f"{labels[window]} & {level.upper()} & {result['n']} & {result['a0_A']:.5f} & {result['B0_GPa']:.2f} & {result['rms_meV_atom']:.3f}" + r"\\")
    return "\n".join(lines + [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = assess()
    if args.write:
        (DATA / "assessment.json").write_text(json.dumps(report, indent=2) + "\n")
        (ROOT / "paper/pccp/si-eos-extension-table-si.tex").write_text(table(report))
    print(json.dumps(report, indent=2))
