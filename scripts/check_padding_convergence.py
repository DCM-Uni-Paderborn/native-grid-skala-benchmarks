#!/usr/bin/env python3
"""Check molecular-cell convergence for the matched ACONF diagnostics."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


HARTREE_TO_KCAL_MOL = 627.5094740631
ENERGY_TOLERANCE_HARTREE = 5.0e-6
REACTION_TOLERANCE_KCAL_MOL = 5.0e-3
ELECTRON_TOLERANCE = 1.0e-5
WAVELET_EDGE_WARNING = "Density non-zero on the"
REACTIONS = {
    "ACONF-5": {
        "053d252e4457117cde6b791ea1b292d848da0bc6": -1,
        "18523213f9133ced8440d2c7b30dcb631edfce6f": 1,
    },
    "ACONF-8": {
        "053d252e4457117cde6b791ea1b292d848da0bc6": -1,
        "1e98e5428039d926bcd2f3279b449c25eb96ffbf": 1,
    },
}
FAMILIES = ("gth-gauxc", "gth-native-one-center-gapwxc")


def last_float(pattern: str, text: str) -> float | None:
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    return float(matches[-1]) if matches else None


def parse_output(path: Path) -> dict[str, float | bool | None]:
    if not path.exists():
        return {"valid": False, "energy": None, "electrons": None}
    text = path.read_text(errors="replace")
    energy = last_float(
        r"^\s*ENERGY\| Total FORCE_EVAL .*?\s+([-+0-9.Ee]+)\s*$", text
    )
    electrons = last_float(
        r"^\s*SKALA_GPW\| Active atom-composite electrons\s+([-+0-9.Ee]+)\s*$",
        text,
    )
    if electrons is None:
        electrons = last_float(r"^\s*Number of electrons:\s+([-+0-9.Ee]+)\s*$", text)
    return {
        "valid": (
            energy is not None
            and "SCF run converged" in text
            and "PROGRAM ENDED AT" in text
            and WAVELET_EDGE_WARNING not in text
        ),
        "energy": energy,
        "electrons": electrons,
    }


def variant_name(family: str, padding: int, solver: str) -> str:
    suffix = "-wavelet" if solver == "wavelet" else ""
    return f"{family}-p{padding}{suffix}"


def compare(
    root: Path,
    padding: int,
    reference: int,
    families: tuple[str, ...],
    solvers: tuple[str, ...],
) -> bool:
    passed = True
    print(
        "family\tsolver\tmax_species_dE_uEh\t"
        "max_reaction_dE_kcal_mol\tmax_dN\tstatus"
    )
    for family in families:
        for solver in solvers:
            candidate = variant_name(family, padding, solver)
            baseline = variant_name(family, reference, solver)
            candidate_rows = {
                digest: parse_output(root / candidate / digest / "output.out")
                for reaction in REACTIONS.values()
                for digest in reaction
            }
            baseline_rows = {
                digest: parse_output(root / baseline / digest / "output.out")
                for reaction in REACTIONS.values()
                for digest in reaction
            }
            valid = all(
                row["valid"]
                for row in (*candidate_rows.values(), *baseline_rows.values())
            )
            if not valid:
                print(f"{family}\t{solver}\t-\t-\t-\tincomplete")
                passed = False
                continue

            species_delta = max(
                abs(float(candidate_rows[digest]["energy"]) - float(baseline_rows[digest]["energy"]))
                for digest in candidate_rows
            )
            electron_delta = max(
                abs(
                    float(candidate_rows[digest]["electrons"])
                    - float(baseline_rows[digest]["electrons"])
                )
                for digest in candidate_rows
                if candidate_rows[digest]["electrons"] is not None
                and baseline_rows[digest]["electrons"] is not None
            )
            reaction_delta = 0.0
            for reaction in REACTIONS.values():
                delta = sum(
                    coefficient
                    * (
                        float(candidate_rows[digest]["energy"])
                        - float(baseline_rows[digest]["energy"])
                    )
                    for digest, coefficient in reaction.items()
                )
                reaction_delta = max(reaction_delta, abs(delta) * HARTREE_TO_KCAL_MOL)

            row_passed = (
                species_delta <= ENERGY_TOLERANCE_HARTREE
                and reaction_delta <= REACTION_TOLERANCE_KCAL_MOL
                and electron_delta <= ELECTRON_TOLERANCE
            )
            passed = passed and row_passed
            print(
                f"{family}\t{solver}\t"
                f"{species_delta * 1.0e6:.3f}\t{reaction_delta:.6f}\t"
                f"{electron_delta:.3e}\t{'pass' if row_passed else 'fail'}"
            )
    return passed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("padding", type=int)
    parser.add_argument("reference", type=int)
    parser.add_argument(
        "--family",
        action="append",
        choices=FAMILIES,
        help="Family to compare; repeat as needed (default: all families)",
    )
    parser.add_argument(
        "--solver",
        action="append",
        choices=("analytic", "wavelet"),
        help="Poisson solver to compare; repeat as needed (default: both)",
    )
    args = parser.parse_args()
    families = tuple(args.family) if args.family else FAMILIES
    solvers = tuple(args.solver) if args.solver else ("analytic", "wavelet")
    sys.exit(
        0
        if compare(args.root, args.padding, args.reference, families, solvers)
        else 1
    )


if __name__ == "__main__":
    main()
