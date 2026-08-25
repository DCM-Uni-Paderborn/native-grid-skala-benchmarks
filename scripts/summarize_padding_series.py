#!/usr/bin/env python3
"""Summarize the validated ACONF molecular-padding series."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from check_padding_convergence import HARTREE_TO_KCAL_MOL, REACTIONS, parse_output


def variant_root(low_root: Path, high_root: Path, family: str, padding: int) -> Path:
    root = low_root if padding <= 20 else high_root
    return root / f"{family}-p{padding}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("low_root", type=Path)
    parser.add_argument("high_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--family", default="gth-native-one-center-gapwxc"
    )
    parser.add_argument("--reference", type=int, default=25)
    parser.add_argument(
        "--paddings", type=int, nargs="+", default=[10, 12, 14, 16, 18, 20, 22, 24, 25]
    )
    args = parser.parse_args()

    digests = sorted({digest for reaction in REACTIONS.values() for digest in reaction})
    reference = {
        digest: parse_output(
            variant_root(
                args.low_root, args.high_root, args.family, args.reference
            )
            / digest
            / "output.out"
        )
        for digest in digests
    }
    if not all(row["valid"] for row in reference.values()):
        raise SystemExit("Reference padding is incomplete")

    rows = []
    for padding in args.paddings:
        candidate = {
            digest: parse_output(
                variant_root(args.low_root, args.high_root, args.family, padding)
                / digest
                / "output.out"
            )
            for digest in digests
        }
        if not all(row["valid"] for row in candidate.values()):
            rows.append([padding, "incomplete", "", ""])
            continue
        species = max(
            abs(float(candidate[digest]["energy"]) - float(reference[digest]["energy"]))
            for digest in digests
        )
        reaction = max(
            abs(
                sum(
                    coefficient
                    * (
                        float(candidate[digest]["energy"])
                        - float(reference[digest]["energy"])
                    )
                    for digest, coefficient in stoichiometry.items()
                )
            )
            * HARTREE_TO_KCAL_MOL
            for stoichiometry in REACTIONS.values()
        )
        electron = max(
            abs(
                float(candidate[digest]["electrons"])
                - float(reference[digest]["electrons"])
            )
            for digest in digests
        )
        rows.append([padding, f"{species * 1e6:.6f}", f"{reaction:.9f}", f"{electron:.9e}"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "padding_angstrom",
                "max_species_delta_microhartree",
                "max_reaction_delta_kcal_mol",
                "max_electron_delta",
            ]
        )
        writer.writerows(rows)
    for row in rows:
        print("\t".join(map(str, row)))


if __name__ == "__main__":
    main()
