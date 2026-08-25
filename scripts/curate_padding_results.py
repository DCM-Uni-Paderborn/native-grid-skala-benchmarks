#!/usr/bin/env python3
"""Collect only validated molecular-padding inputs and compact results."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from check_padding_convergence import REACTIONS


LOW_VARIANTS = [
    f"{family}-p{padding}"
    for family in ("gth-native-one-center-gapwxc", "gth-gauxc")
    for padding in (10, 12, 14, 16, 18, 20)
]
HIGH_VARIANTS = [
    "gth-native-one-center-gapwxc-p22",
    "gth-native-one-center-gapwxc-p24",
    "gth-native-one-center-gapwxc-p25",
    "gth-gauxc-p24",
]
FILES = ("input.inp", "output.out", "time.txt")


def copy_variant(source: Path, target: Path, variant: str) -> int:
    copied = 0
    for run_dir in sorted((source / variant).iterdir()):
        if not run_dir.is_dir():
            continue
        output = run_dir / "output.out"
        if not output.exists():
            continue
        text = output.read_text(errors="replace")
        if "SCF run converged" not in text or "PROGRAM ENDED AT" not in text:
            raise SystemExit(f"Refusing incomplete output: {output}")
        destination = target / variant / run_dir.name
        destination.mkdir(parents=True, exist_ok=True)
        for name in FILES:
            path = run_dir / name
            if path.exists():
                shutil.copy2(path, destination / name)
        copied += 1
    if copied != 3:
        raise SystemExit(f"Expected three validated cases for {variant}, found {copied}")
    return copied


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("low_root", type=Path)
    parser.add_argument("high_root", type=Path)
    parser.add_argument("target", type=Path)
    args = parser.parse_args()

    args.target.mkdir(parents=True, exist_ok=True)
    variants = {}
    for variant in LOW_VARIANTS:
        variants[variant] = copy_variant(args.low_root, args.target, variant)
    for variant in HIGH_VARIANTS:
        variants[variant] = copy_variant(args.high_root, args.target, variant)

    for name in ("README.md", "native-series-vs-p25.tsv", "gauxc-series-vs-p24.tsv"):
        shutil.copy2(args.high_root / name, args.target / name)
    manifest = {
        "status": "validated",
        "variants": variants,
        "accepted_total_padding_angstrom": 25,
        "accepted_route": "GAPW_XC-GTH/PAW_ONE_CENTER",
        "cell_definition": "L_i = molecular_extent_i + total_padding_angstrom",
        "boundary_conditions": {
            "cell_periodic": "NONE",
            "poisson_periodic": "NONE",
            "poisson_solver": "ANALYTIC",
            "coordinates": "centered",
        },
        "acceptance_limits": {
            "species_hartree": 5.0e-6,
            "reaction_kcal_mol": 5.0e-3,
            "electrons": 1.0e-5,
        },
        "reactions": REACTIONS,
        "files_per_case": list(FILES),
    }
    (args.target / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
