"""Checks specific to the paired molecular AE basis controls."""
import re


def value(text, keyword):
    found = re.findall(r"(?m)^\s*" + keyword + r"\s+([^\n]+)", text)
    if len(found) != 1:
        raise ValueError((keyword, found))
    return found[0].strip()


def normalize(text):
    text = re.sub(r"(?m)^(\s*PROJECT_NAME)\s+\S+", r"\1 PROJECT", text)
    text = re.sub(r"(?m)^(\s*EPS_SCF)\s+\S+", r"\1 EPS", text)
    return re.sub(r"BASIS_SET (?:TZVPP|QZVPP)-MOLOPT-PBE-ae", "BASIS_SET AE_BASIS", text)


def validate_input(text, case):
    assert int(value(text, "CHARGE")) == case["charge"]
    assert int(value(text, "MULTIPLICITY")) == case["multiplicity"]
    assert value(text, "UKS") == (".TRUE." if case["multiplicity"] != 1 else ".FALSE.")
    assert value(text, "SCF_GUESS") == "ATOMIC"
    assert value(text, "NATIVE_GRID_USE_CUDA") == ".FALSE."
    assert value(text, "NATIVE_GRID_LAYOUT") == "ATOM_COMPOSITE"
    assert value(text, "POISSON_SOLVER") == "ANALYTIC"
    assert re.findall(r"(?m)^\s*PERIODIC\s+(\S+)", text) == ["NONE", "NONE"]
    assert [float(x) for x in re.findall(r"(?m)^\s*EPS_SCF\s+(\S+)", text)] == [5e-7, 5e-7]
    assert value(text, "CUTOFF") == "640.0" and value(text, "REL_CUTOFF") == "60.0"
    assert value(text, "MINIMIZER") == "DIIS" and value(text, "PRECONDITIONER") == "FULL_ALL"
    assert value(text, "GAPW_ACCURATE_XCINT") == ".TRUE."
    assert value(text, "REFERENCE_FUNCTIONAL") == "B3LYP"
    assert "&KPOINTS" not in text and "WFN_RESTART_FILE_NAME" not in text
    assert re.findall(r"(?m)^\s*BASIS_SET\s+(\S+)", text) == [case["basis"]] * case["kinds"]
    assert re.findall(r"(?m)^\s*POTENTIAL\s+(\S+)", text) == ["ALL"] * case["kinds"]
    assert re.findall(r"(?m)^\s*RADIAL_GRID\s+(\S+)", text) == ["100"] * case["kinds"]
    assert re.findall(r"(?m)^\s*LEBEDEV_GRID\s+(\S+)", text) == ["434"] * case["kinds"]
    assert "&CENTER_COORDINATES" in text
    coord = re.search(r"&COORD\s*\n(.*?)&END COORD", text, re.S).group(1)
    atoms = [line.split() for line in coord.splitlines() if line.strip()]
    xyz = [[float(x) for x in row[1:]] for row in atoms]
    assert len(atoms) == case["atoms"] and all(len(r) == 3 for r in xyz)
    z = {"H": 1, "C": 6, "N": 7, "O": 8, "Ne": 10}
    assert sum(z[row[0]] for row in atoms) - case["charge"] == case["expected_electrons"]
    abc = [float(x) for x in value(text, "ABC").split()]
    for axis, cell in enumerate(abc):
        extent = max(p[axis] for p in xyz) - min(p[axis] for p in xyz)
        assert 25 - 1e-7 <= cell - extent < 26.000001
        assert abs(sum(p[axis] for p in xyz) / len(xyz)) < 1e-5
    return True


def electron_counts(case):
    n, spin = case["expected_electrons"], case["multiplicity"] - 1
    return [(n + spin) // 2, (n - spin) // 2] if spin else [n]
