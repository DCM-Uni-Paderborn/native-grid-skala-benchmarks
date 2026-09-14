#!/usr/bin/env python3
"""Form provisional lattice energies only from matched, accepted pilot results."""
import json
import math
import hashlib
from pathlib import Path
from analyze_basis_controls import selected_pairs

ROOT = Path(__file__).resolve().parents[1]
REVISION = "37aedacbb33677307818301a47a09ae435712293"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(path, expected):
    if sha(path) != expected:
        raise ValueError(f"Frozen input mismatch: {path}")

HARTREE_TO_KJMOL = 2625.4996394799


def pair_result(solid, molecule, reference):
    prefix, phase = solid["case_id"].rsplit("/", 1)
    if phase != "solid" or molecule["case_id"] != prefix + "/molecule":
        raise ValueError("Cannot mix systems or representations")
    if any(r["status"] != "accepted_scf_numerical_accuracy_unassessed" for r in (solid, molecule)):
        raise ValueError("Unaccepted result")
    z = solid["molecules"]
    if z < 1 or molecule["molecules"] != 1 or solid["electrons"] != z * molecule["electrons"]:
        raise ValueError("Inconsistent molecular or electron normalization")
    if any(r["source_revision"] != REVISION or r["charge"] != 0 or r["multiplicity"] != 1
           or not r["restricted"] or not math.isfinite(r["energy_hartree"]) for r in (solid, molecule)):
        raise ValueError("Incompatible physical provenance")
    energy = (solid["energy_hartree"] / z - molecule["energy_hartree"]) * HARTREE_TO_KJMOL
    return {"method": prefix.split("/")[0], "system": prefix.split("/")[1],
            "molecules_per_cell": z, "solid_energy_hartree": solid["energy_hartree"],
            "molecule_energy_hartree": molecule["energy_hartree"], "lattice_energy_kjmol": energy,
            "dmc_kjmol": reference[0], "dmc_statistical_uncertainty_kjmol": reference[1],
            "signed_deviation_kjmol": energy - reference[0],
            "status": "provisional_numerical_convergence_unassessed",
            "execution_sha256": {"solid": solid["execution_sha256"], "molecule": molecule["execution_sha256"]}}


def representation_comparisons(pairs):
    by_key = {(p["method"], p["system"]): p for p in pairs}
    if len(by_key) != len(pairs):
        raise ValueError("Duplicate lattice-energy pair")
    comparisons = []
    for system in sorted({p["system"] for p in pairs}):
        for a, b in (("gapwxc-gth", "gapw-ae"), ("gapwxc-gth", "gapw-gth-direct"),
                     ("gapwxc-gth", "gapw-gth-one-center"), ("gapw-gth-one-center", "gapw-gth-direct")):
            if (a, system) in by_key and (b, system) in by_key:
                comparisons.append({"system": system, "method_a": a, "method_b": b,
                                    "delta_a_minus_b_kjmol": by_key[a, system]["lattice_energy_kjmol"] - by_key[b, system]["lattice_energy_kjmol"],
                                    "status": "provisional_numerical_convergence_unassessed"})
    return comparisons


def control_difference(control, parent):
    if (control["status"] != "accepted_scf_numerical_accuracy_unassessed"
            or control["case_id"] != parent["case_id"]
            or control["parent_execution_sha256"] != parent["execution_sha256"]):
        raise ValueError("Unaccepted or mismatched numerical control")
    for key in ("source_revision", "binary_sha256", "electrons", "molecules", "charge", "multiplicity", "restricted"):
        if control[key] != parent[key]:
            raise ValueError(f"Control changes physical/runtime provenance: {key}")
    delta = (control["energy_hartree"] - parent["energy_hartree"]) / parent["molecules"] * HARTREE_TO_KJMOL
    if not math.isfinite(delta):
        raise ValueError("Nonfinite control difference")
    return {"case_id": control["case_id"], "control_name": control["control_name"],
            "delta_control_minus_parent_kjmol_per_molecule": delta,
            "control_execution_sha256": control["execution_sha256"],
            "parent_execution_sha256": parent["execution_sha256"]}


def main():
    base = json.loads((ROOT / "manifest.json").read_text())
    cases, accepted = {}, {}
    for name in ("manifest.json", "manifest-gapw-gth.json"):
        manifest_path = ROOT / name
        for case in json.loads(manifest_path.read_text())["cases"]:
            cases[case["id"]] = (case, sha(manifest_path))
    evidence = {}
    for path in sorted((ROOT / "results/validated").glob("*/*/*/validation.json")):
        record = json.loads(path.read_text())
        case, manifest_hash = cases[record["case_id"]]
        verify(ROOT / case["input"], case["input_sha256"])
        if record["frozen_manifest_sha256"] != manifest_hash or record["final_hashes"]["frozen-input.inp"] != case["input_sha256"]:
            raise ValueError("Frozen evidence mismatch")
        if record["status"] != "accepted_scf_numerical_accuracy_unassessed":
            raise ValueError("Unaccepted validation file")
        if record["case_id"] in accepted:
            raise ValueError("Duplicate accepted case")
        accepted[record["case_id"]] = record
        evidence[str(path.relative_to(ROOT))] = sha(path)
    pairs = []
    for case_id, solid in sorted(accepted.items()):
        prefix, phase = case_id.rsplit("/", 1)
        molecule = accepted.get(prefix + "/molecule")
        if phase == "solid" and molecule:
            pairs.append(pair_result(solid, molecule, base["references"][prefix.split("/")[1]]))
    contrasts = representation_comparisons(pairs)
    controls = []
    for path in sorted((ROOT / "results/controls").glob("*/validation.json")):
        record = json.loads(path.read_text())
        parent = accepted[record["case_id"]]
        controls.append(control_difference(record, parent))
        evidence[str(path.relative_to(ROOT))] = sha(path)
    paired_controls = paired_control_differences(controls)
    report = {"scope": base["scope"], "accepted_base_cases": len(accepted), "total_base_cases": len(cases),
              "missing_base_cases": sorted(set(cases) - set(accepted)), "complete_pairs": pairs,
              "representation_comparisons": contrasts,
              "paper_pairs": selected_pairs(pairs),
              "numerical_controls": controls,
              "paired_numerical_controls": paired_controls,
              "reference_doi": base["reference_doi"], "reference_uncertainty": base["reference_uncertainty"],
              "conversion_hartree_to_kjmol": HARTREE_TO_KJMOL, "validation_file_sha256": evidence,
              "caveat": "SCF/provenance accepted. Listed numerical controls quantify only their stated changes for the tested systems and representations, not a universal error bound. Basis/BSSE and untested numerical effects remain unresolved."}
    text = "# Selected PCCP X23-mini lattice energies\n\n"
    text += f"{len(accepted)}/{len(cases)} base cases accepted; {len(pairs)}/12 crystal/molecule pairs complete.\n\n"
    text += "| System | Method | E_latt (kJ/mol) | DMC (kJ/mol) | Signed deviation (kJ/mol) |\n| --- | --- | ---: | ---: | ---: |\n"
    for p in report["paper_pairs"]:
        text += f"| {p['system']} | {p['method']} | {p['lattice_energy_kjmol']:.3f} | {p['dmc_kjmol']:.1f} +/- {p['dmc_statistical_uncertainty_kjmol']:.1f} | {p['signed_deviation_kjmol']:+.3f} |\n"
    text += "\nUrea AE uses the validated QZVPP/200-974 crystal and molecule. The original TZVPP pair remains in complete_pairs in the JSON and in the basis-convergence evidence. Other selections are unchanged. Four additional Urea controls separate grid and basis sensitivity.\n"
    text += "\nNegative deviations indicate stronger binding. DMC error bars are statistical; total DMC accuracy is estimated at about 1-2 kJ/mol.\n\n"
    text += "These are preliminary single-point results with bounded representative numerical checks, not a comprehensive convergence study. "
    text += "The differences must not yet be attributed solely to the functional or density representation. "
    text += "Numerical controls below quantify only their stated changes and cannot be generalized to every system or representation. Basis/BSSE and untested effects remain unresolved.\n\n"
    if contrasts:
        text += "## Matched representation differences\n\n"
        text += "| System | Method A | Method B | A minus B (kJ/mol) |\n| --- | --- | --- | ---: |\n"
        for c in representation_comparisons(report["paper_pairs"]):
            text += f"| {c['system']} | {c['method_a']} | {c['method_b']} | {c['delta_a_minus_b_kjmol']:+.3f} |\n"
        text += "\nOnly completed pairs for the same system enter these differences. "
        text += "GAPW_XC versus direct GAPW changes more than the one-center treatment; isolate that effect only with the two GAPW-GTH variants.\n\n"
    if controls:
        text += "## Accepted numerical controls\n\n"
        text += "| Case | Control | Control minus parent (kJ/mol of molecules) |\n| --- | --- | ---: |\n"
        for c in controls:
            text += f"| {c['case_id']} | {c['control_name']} | {c['delta_control_minus_parent_kjmol_per_molecule']:+.6f} |\n"
        text += "\nControls are not additional base cases. The CO2 CPU/GPU comparison tests the crystal model evaluation on Spark, not other hosts, quadrature or basis completeness.\n\n"
        text += "For crystal-only k-mesh and symmetry-tolerance changes, the isolated molecule is unchanged, so the tabulated shift is also the lattice-energy change. Cutoff and quadrature changes require matched controls of both phases.\n\n"
    if paired_controls:
        text += "## Paired lattice-energy sensitivity\n\n"
        text += "| System | Method | Control | Change in E_latt (kJ/mol) |\n| --- | --- | --- | ---: |\n"
        for c in paired_controls:
            text += f"| {c['system']} | {c['method']} | {c['control_name']} | {c['delta_lattice_energy_kjmol']:+.6f} |\n"
        text += "\nEach difference uses both accepted phases relative to their exact accepted parents. It assesses only the stated numerical change for that system and representation.\n\n"
    text += f"Reference: [Della Pia et al., PRL 133, 046401 (2024)](https://doi.org/{base['reference_doi']}).\n"
    for path, content in ((ROOT / "results/lattice-energies.json", json.dumps(report, indent=2) + "\n"),
                          (ROOT / "results/lattice-energies.md", text)):
        if not path.exists() or path.read_text() != content:
            path.write_text(content)
    print(text)


def paired_control_differences(controls):
    grouped = {}
    for control in controls:
        prefix, phase = control["case_id"].rsplit("/", 1)
        name = "cutoff-1200" if control["control_name"] == "cutoff-1200-mem200" else control["control_name"]
        phases = grouped.setdefault((prefix, name), {})
        if phase not in ("solid", "molecule") or phase in phases:
            raise ValueError("Invalid or duplicate control phase")
        phases[phase] = control
    pairs = []
    for (prefix, name), phases in sorted(grouped.items()):
        if set(phases) == {"solid", "molecule"}:
            method, system = prefix.split("/")
            delta = (phases["solid"]["delta_control_minus_parent_kjmol_per_molecule"]
                     - phases["molecule"]["delta_control_minus_parent_kjmol_per_molecule"])
            pairs.append({"method": method, "system": system, "control_name": name,
                          "delta_lattice_energy_kjmol": delta,
                          "control_execution_sha256": {p: c["control_execution_sha256"] for p, c in phases.items()},
                          "parent_execution_sha256": {p: c["parent_execution_sha256"] for p, c in phases.items()}})
    return pairs


if __name__ == "__main__":
    main()
