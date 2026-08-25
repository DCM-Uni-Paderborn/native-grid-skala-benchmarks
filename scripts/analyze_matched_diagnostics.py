#!/usr/bin/env python3
"""Extract matched-route molecular and reaction energies from CP2K outputs."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


HARTREE_TO_KCAL_MOL = 627.5094740631
WAVELET_EDGE_WARNING = "Density non-zero on the"


def last_float(pattern: str, text: str) -> float | None:
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    return float(matches[-1]) if matches else None


def cell_lengths(text: str) -> list[float] | None:
    vectors = re.findall(
        r"CELL\| Vector [abc] \[angstrom\]:\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
        text,
    )
    if len(vectors) < 3:
        return None
    return [sum(float(value) ** 2 for value in vector) ** 0.5 for vector in vectors[:3]]


def parse_output(run_dir: Path) -> dict[str, object]:
    output_file = run_dir / "output.out"
    timing_file = run_dir / "time.txt"
    if not output_file.exists():
        return {"status": "missing"}

    text = output_file.read_text(errors="replace")
    clean = "PROGRAM ENDED AT" in text
    converged = "SCF run converged" in text
    result: dict[str, object] = {
        "status": (
            "invalid-wavelet-boundary"
            if WAVELET_EDGE_WARNING in text
            else "complete" if clean and converged else "incomplete"
        ),
        "wavelet_edge_warning": WAVELET_EDGE_WARNING in text,
        "scf_steps": len(re.findall(r"^\s*\d+\s+(?:OT|Broy\.|Mix\/Diag\.)", text, re.MULTILINE)),
        "energy_hartree": last_float(
            r"^\s*ENERGY\| Total FORCE_EVAL .*?\s+([-+0-9.Ee]+)\s*$", text
        ),
        "dispersion_hartree": last_float(
            r"^\s*Dispersion energy:\s+([-+0-9.Ee]+)\s*$", text
        ),
        "xc_hartree": last_float(
            r"^\s*SKALA_GPW\| Active atom-composite XC energy\s+([-+0-9.Ee]+)\s*$",
            text,
        ),
        "integrated_electrons": last_float(
            r"^\s*SKALA_GPW\| Active atom-composite electrons\s+([-+0-9.Ee]+)\s*$",
            text,
        ),
        "cell_angstrom": cell_lengths(text),
    }
    if result["xc_hartree"] is None:
        result["xc_hartree"] = last_float(
            r"^\s*Exchange-correlation energy:\s+([-+0-9.Ee]+)\s*$", text
        )
    if result["integrated_electrons"] is None:
        result["integrated_electrons"] = last_float(
            r"^\s*Number of electrons:\s+([-+0-9.Ee]+)\s*$", text
        )

    if timing_file.exists():
        timing = timing_file.read_text(errors="replace")
        memory = re.findall(r"Maximum resident set size \(kbytes\):\s+(\d+)", timing)
        wall = re.findall(r"Elapsed \(wall clock\) time[^:]*:\s+(.+)", timing)
        result["peak_memory_kib"] = int(memory[-1]) if memory else None
        result["wall_time"] = wall[-1].strip() if wall else None
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--json", type=Path, help="Optional path for structured results")
    args = parser.parse_args()

    manifest = json.loads((args.root / "manifest.json").read_text())
    results: dict[str, object] = {"species": {}, "reactions": {}}

    print("variant\tcomplete\ttotal\treaction\tenergy_kcal_mol")
    for variant in manifest["variants"]:
        species_results = {
            digest: parse_output(args.root / variant / digest)
            for digest in manifest["species"]
        }
        results["species"][variant] = species_results
        complete = sum(row["status"] == "complete" for row in species_results.values())

        reaction_results = {}
        for reaction, stoichiometry in manifest["reactions"].items():
            rows = [species_results[digest] for digest in stoichiometry]
            if all(row["status"] == "complete" for row in rows):
                energy = sum(
                    coefficient * float(species_results[digest]["energy_hartree"])
                    for digest, coefficient in stoichiometry.items()
                )
                energy *= HARTREE_TO_KCAL_MOL
                reaction_results[reaction] = energy
                print(
                    f"{variant}\t{complete}\t{len(species_results)}\t"
                    f"{reaction}\t{energy:.9f}"
                )
        results["reactions"][variant] = reaction_results
        if not reaction_results:
            print(f"{variant}\t{complete}\t{len(species_results)}\t-\t-")

    if args.json:
        args.json.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
