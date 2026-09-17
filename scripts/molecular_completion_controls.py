"""Reproduce two additional QZ reaction controls without changing benchmark data."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "convergence/molecular-qz-completion-20260917"
PRODUCTION = ROOT / "benchmarks/dietGMTKN55/production-p25"
HA_KCAL = 627.5094740631


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def without_scf(text):
    """Remove only SCF controls, retaining every physical input section."""
    lines = text.splitlines(keepends=True)
    result, depth, count = [], 0, 0
    for line in lines:
        tokens = line.strip().upper().split()
        token = tokens[0] if tokens else ""
        if depth == 0 and token == "&SCF":
            depth, count = 1, count + 1
            continue
        if depth:
            if token == "&END":
                depth -= 1
            elif token.startswith("&"):
                depth += 1
        else:
            result.append(line)
    assert depth == 0 and count == 1
    text = "".join(result)
    text = re.sub(r"(?m)^\s*(PROJECT_NAME|WFN_RESTART_FILE_NAME)\s+[^\n]*\n", "", text)
    return text.replace("QZVPP-MOLOPT-PBE-ae", "TZVPP-MOLOPT-PBE-ae")


def parse_output(text):
    assert "PROGRAM STARTED" in text
    block = text.rsplit("PROGRAM STARTED", 1)[1]
    steps = []
    for line in block.splitlines():
        fields = line.split()
        if len(fields) >= 7 and fields[0].isdigit() and (fields[1] == "OT" or "/Diag." in fields[1]):
            steps.append({"step": int(fields[0]), "residual": float(fields[-3]), "energy": float(fields[-2])})
    energies = re.findall(r"ENERGY\| Total FORCE_EVAL.*?\[hartree\]\s+([-+\d.Ee]+)", block)
    return dict(converged="SCF run converged" in block and "SCF run NOT converged" not in block,
                ended="PROGRAM ENDED" in block, abort="[ABORT]" in block,
                energy_hartree=float(energies[-1]) if energies else None, steps=steps,
                electrons=list(map(int, re.findall(r"Number of electrons:\s+(\d+)", block))))


def verified(case):
    directory = DATA / case["target"]
    for name, digest in case["hashes"].items():
        assert sha(directory / name) == digest, (directory, name)
    record = json.loads((directory / "execution.json").read_text())
    for name in case["hashes"]:
        if name != "execution.json":
            assert sha(directory / name) == record["final_hashes"][name]
    inp = (directory / "actual-input.inp").read_text()
    original = record["case"]
    assert sha(directory / "actual-input.inp") == original["input_sha256"]
    assert original["digest"] == case["digest"]
    source = ROOT / original["source_input"]
    assert sha(source) == original["source_input_sha256"]
    assert without_scf(inp) == without_scf(source.read_text()), case["name"]
    assert "&KPOINTS" not in inp and "&SMEAR" not in inp and "ADDED_MOS" not in inp
    assert set(re.findall(r"(?m)^\s*BASIS_SET\s+(\S+)", inp)) == {"QZVPP-MOLOPT-PBE-ae"}
    assert set(re.findall(r"(?m)^\s*POTENTIAL\s+(\S+)", inp)) == {"ALL"}
    assert re.findall(r"(?m)^\s*PERIODIC\s+(\S+)", inp) == ["NONE", "NONE"]
    runtime = record["runtime"]
    assert runtime["revision"] == "37aedacbb33677307818301a47a09ae435712293"
    assert runtime["binary_sha256"] == "4c4cdcf4ac51510d1d876310efa724ae7efa8bdab7fdd124d803e0220fb770fb"
    assert record["model_sha256"] == "7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd"
    launch = json.loads((directory / "launch-environment.json").read_text())
    assert launch["affinity"] == record["affinity"]
    assert set(launch["loaded_cp2k_libraries"]) == set(runtime["cp2k_shared_libraries"])
    assert int(launch["runtime_environment"]["OMP_NUM_THREADS"]) == len(record["affinity"])
    parsed = parse_output((directory / "output.out").read_text())
    assert parsed["electrons"] == original["expected_spin_electrons"]
    assert record["peak_rss_kib"] > 0 and record["elapsed_wall_s"] > 0
    if case["stage"] == "final":
        assert parsed["converged"] and parsed["ended"] and not parsed["abort"]
        assert parsed["steps"][-1]["residual"] <= 5e-7
        assert set(map(float, re.findall(r"(?m)^\s*EPS_SCF\s+(\S+)", inp))) == {5e-7}
        assert math.isfinite(parsed["energy_hartree"]) and parsed["energy_hartree"] < 0
        assert record["execution_exit"] == record["time_exit"] == 0
        assert re.search(r"Exit status:\s+0\s*$", (directory / "time.txt").read_text())
    else:
        assert not parsed["converged"] and record["execution_exit"] != 0
    return record, parsed


def assess():
    index = json.loads((DATA / "index.json").read_text())
    assert len(index["cases"]) == 7
    cases = {(c["name"], c["stage"]): c for c in index["cases"]}
    records = {key: verified(c) for key, c in cases.items()}
    for name in ("3CL", "EA_25"):
        final, precursor = records[name, "final"][0], records[name, "precursor"][0]
        seed = final["case"]["seed"]
        parent_hashes = cases[name, "precursor"]["hashes"]
        assert seed["parent_execution_sha256"] == parent_hashes["execution.json"]
        assert seed["parent_output_sha256"] == parent_hashes["output.out"]
        assert seed["sha256"] == precursor["final_hashes"][Path(seed["source"]).name]
        assert seed["sha256"] == final["final_hashes"]["source.wfn"]
        assert final["diagnostic_only"] is True
    reaction_index = json.loads((PRODUCTION / "paper-common-70/reaction-index.json").read_text())
    with (PRODUCTION / "paper-common-70/molecular-interface-comparison.csv").open() as stream:
        comparison = {(r["subset"], int(r["reaction_id"])): r for r in csv.DictReader(stream)}
    with (PRODUCTION / "paper-common-70/species-results.csv").open() as stream:
        species = {r["digest"]: r for r in csv.DictReader(stream) if r["route"] == "shared-gapw-ae"}
    rows = []
    for key in index["reactions"]:
        reaction = next(r for r in reaction_index["reactions"] if (r["subset"], r["reaction_id"]) == tuple(key))
        tz, qz = 0.0, 0.0
        for item in reaction["species"]:
            record, parsed = records[item["name"], "final"]
            assert record["case"]["charge"] == item["charge"]
            assert record["case"]["multiplicity"] == item["multiplicity"]
            assert record["case"]["digest"] == item["digest"]
            tz += item["count"] * float(species[item["digest"]]["total_energy_ha"])
            qz += item["count"] * parsed["energy_hartree"]
        old = comparison[tuple(key)]
        assert abs(tz * HA_KCAL - float(old["hybrid_ae_gth_direct"])) < 1e-7
        row = dict(subset=key[0], reaction_id=key[1], reference_kcal_mol=reaction["reference_kcal_mol"],
                   tz_kcal_mol=tz * HA_KCAL, qz_kcal_mol=qz * HA_KCAL,
                   gauxc_ae_kcal_mol=float(old["gauxc_ae"]), pyscf_kcal_mol=float(old["pyscf_unit"]))
        for level in ("tz", "qz"):
            row[level + "_error_kcal_mol"] = row[level + "_kcal_mol"] - row["reference_kcal_mol"]
        rows.append(row)
    return dict(scope=index["scope"], reactions=rows, included_in_uniform_statistics=False,
                precursor_energies_used=False,
                final_species={name: dict(energy_hartree=p["energy_hartree"], residual=p["steps"][-1]["residual"],
                                         scf_steps=len(p["steps"]), electrons=p["electrons"],
                                         input_sha256=r["case"]["input_sha256"])
                               for (name, stage), (r, p) in records.items() if stage == "final"})


def table(report):
    lines = [r"\begin{table}[htbp]", r"\centering\small",
             r"\caption{Additional reaction-resolved molecular basis controls, in kcal~mol$^{-1}$. QZVPP energies use the same physical input as the archived TZVPP production calculations, with the enlarged orbital basis and a converged unsmeared SCF solution. The GauXC and PySCF values retain their respective protocols. These controls do not replace entries in the uniform 70-reaction comparison.}",
             r"\label{tab:molecular-additional-basis}", r"\begin{tabular}{lrrrrr}", r"\toprule",
             r"Reaction & Reference & TZVPP & QZVPP & GauXC AE & PySCF\\", r"\midrule"]
    for r in report["reactions"]:
        lines.append(f"{r['subset']}/{r['reaction_id']} & {r['reference_kcal_mol']:.3f} & {r['tz_kcal_mol']:.3f} & {r['qz_kcal_mol']:.3f} & {r['gauxc_ae_kcal_mol']:.3f} & {r['pyscf_kcal_mol']:.3f}" + r"\\")
    return "\n".join(lines + [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = assess()
    if args.write:
        (DATA / "assessment.json").write_text(json.dumps(report, indent=2) + "\n")
        (ROOT / "paper/pccp/molecular-additional-basis-table-si.tex").write_text(table(report))
    print(json.dumps(report, indent=2))
