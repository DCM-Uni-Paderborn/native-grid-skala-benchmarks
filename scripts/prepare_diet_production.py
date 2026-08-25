#!/usr/bin/env python3
"""Build the three native-Skala dietGMTKN55 production protocols."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


HEAVY_GTH_ELEMENTS = {"Bi", "I", "Te"}

SCF_BLOCK = """      &SCF
         EPS_SCF 1e-07
         MAX_SCF 200
         SCF_GUESS ATOMIC
         &OUTER_SCF
            EPS_SCF 1e-07
            MAX_SCF 20
         &END OUTER_SCF
         &OT
            ENERGY_GAP 0.001
            MINIMIZER DIIS
            PRECONDITIONER FULL_ALL
            STEPSIZE 0.10
         &END OT
      &END SCF
"""


def replace_once(
    text: str, pattern: str, replacement: str, *, dotall: bool = False
) -> str:
    flags = re.MULTILINE | (re.DOTALL if dotall else 0)
    updated, count = re.subn(pattern, replacement, text, count=1, flags=flags)
    if count != 1:
        raise ValueError(f"Expected one match for {pattern!r}, found {count}")
    return updated


def set_total_padding(text: str, padding: int) -> str:
    coord_match = re.search(r"(?ms)^\s*&COORD\s*$\n(.*?)^\s*&END COORD\s*$", text)
    if coord_match is None:
        raise ValueError("Coordinate block not found")

    coordinates = []
    for line in coord_match.group(1).splitlines():
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            coordinates.append(tuple(float(value) for value in fields[-3:]))
        except ValueError:
            continue
    if not coordinates:
        raise ValueError("No Cartesian coordinates found")

    lengths = []
    for axis in range(3):
        values = [coordinate[axis] for coordinate in coordinates]
        lengths.append(float(math.ceil(max(values) - min(values) + padding)))
    abc = " ".join(f"{length:.1f}" for length in lengths)
    return replace_once(text, r"(?m)^(\s*ABC\s+).*$", rf"\g<1>{abc}")


def common_transform(
    text: str, cutoff: int, radial: int, lebedev: int, padding: int, project: str
) -> str:
    text = replace_once(
        text, r"^      &SCF\n.*?^      &END SCF\n", SCF_BLOCK, dotall=True
    )
    text = replace_once(text, r"(?m)^(\s*CUTOFF\s+)\S+", rf"\g<1>{cutoff}.0")
    text = replace_once(text, r"(?m)^(\s*PROJECT_NAME\s+)\S+", rf"\g<1>{project}")
    text, radial_count = re.subn(
        r"(?m)^(\s*RADIAL_GRID\s+)\S+", rf"\g<1>{radial}", text
    )
    text, lebedev_count = re.subn(
        r"(?m)^(\s*LEBEDEV_GRID\s+)\S+", rf"\g<1>{lebedev}", text
    )
    if radial_count == 0 or lebedev_count == 0:
        raise ValueError("Native atom-grid settings are missing")
    return set_total_padding(text, padding)


def set_qs_method(text: str, method: str) -> str:
    return replace_once(text, r"(?m)^(\s*METHOD\s+)(?:GAPW|GAPW_XC)$", rf"\g<1>{method}")


def convert_covered_kinds_to_ae(text: str) -> str:
    kind_pattern = re.compile(
        r"(?ms)^(\s*)&KIND\s+(\S+)\n(.*?)^\1&END KIND\n"
    )

    def convert(match: re.Match[str]) -> str:
        indent, element, body = match.groups()
        if element in HEAVY_GTH_ELEMENTS:
            return match.group(0)
        body, basis_count = re.subn(
            r"(?m)^(\s*BASIS_SET\s+)\S+", r"\g<1>TZVPP-MOLOPT-PBE-ae", body, count=1
        )
        body, potential_count = re.subn(
            r"(?m)^(\s*POTENTIAL\s+)\S+", r"\g<1>ALL", body, count=1
        )
        if basis_count != 1 or potential_count != 1:
            raise ValueError(f"Could not convert KIND {element} to all electron")
        return f"{indent}&KIND {element}\n{body}{indent}&END KIND\n"

    converted, count = kind_pattern.subn(convert, text)
    if count == 0:
        raise ValueError("No KIND sections found")
    return converted


def write_case(destination: Path, digest: str, route: str, text: str) -> str:
    relative = Path(route) / digest / "input.inp"
    (destination / relative).parent.mkdir(parents=True, exist_ok=True)
    (destination / relative).write_text(text)
    return str(relative)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, help="Directory containing diet-inputs/manifest.json")
    parser.add_argument("destination", type=Path)
    parser.add_argument("--ae-cutoff", type=int, required=True)
    parser.add_argument("--gapwxc-cutoff", type=int, default=400)
    parser.add_argument("--radial-grid", type=int, default=100)
    parser.add_argument("--lebedev-grid", type=int, default=434)
    parser.add_argument("--padding", type=int, default=12)
    args = parser.parse_args()

    manifest = json.loads((args.source / "manifest.json").read_text())
    output_manifest = {
        "protocol": {
            "shared_ae_cutoff_ry": args.ae_cutoff,
            "hybrid_cutoff_ry": args.ae_cutoff,
            "gapwxc_gth_cutoff_ry": args.gapwxc_cutoff,
            "radial_grid": args.radial_grid,
            "lebedev_grid": args.lebedev_grid,
            "molecular_padding_angstrom": args.padding,
            "eps_scf": 1.0e-7,
            "ot_stepsize": 0.10,
            "heavy_gth_elements": sorted(HEAVY_GTH_ELEMENTS),
        },
        "species": {},
    }

    for digest, record in manifest["species"].items():
        elements = set(record["elements"])
        heavy = bool(elements & HEAVY_GTH_ELEMENTS)
        routes: dict[str, str] = {}

        if not heavy:
            source = args.source / record["routes"]["gapw-ae"]
            text = common_transform(
                source.read_text(), args.ae_cutoff, args.radial_grid,
                args.lebedev_grid, args.padding, f"{digest[:12]}_gapw_ae"
            )
            text = text.replace("QZVPP-MOLOPT-PBE-ae", "TZVPP-MOLOPT-PBE-ae")
            routes["shared-gapw-ae"] = write_case(
                args.destination, digest, "shared-gapw-ae", text
            )
        else:
            for route, source_route, representation in (
                ("hybrid-direct-heavy", "gapw-gth-valence", "DIRECT_VALENCE"),
                ("hybrid-one-center-heavy", "gapw-gth-paw", "PAW_ONE_CENTER"),
            ):
                source = args.source / record["routes"][source_route]
                text = convert_covered_kinds_to_ae(source.read_text())
                text = common_transform(
                    text, args.ae_cutoff, args.radial_grid, args.lebedev_grid,
                    args.padding, f"{digest[:12]}_{route.replace('-', '_')}"
                )
                text = set_qs_method(text, "GAPW")
                if f"PSEUDOPOTENTIAL_GAPW_REPRESENTATION {representation}" not in text:
                    raise ValueError(f"Incorrect representation in {digest} {route}")
                routes[route] = write_case(args.destination, digest, route, text)

        source = args.source / record["routes"]["gapw-gth-paw"]
        text = common_transform(
            source.read_text(), args.gapwxc_cutoff, args.radial_grid,
            args.lebedev_grid, args.padding, f"{digest[:12]}_gapwxc_gth_paw"
        )
        text = set_qs_method(text, "GAPW_XC")
        routes["gapwxc-gth-one-center"] = write_case(
            args.destination, digest, "gapwxc-gth-one-center", text
        )

        output_manifest["species"][digest] = {
            "charge": record["charge"],
            "multiplicity": record["multiplicity"],
            "elements": record["elements"],
            "heavy_gth": heavy,
            "routes": routes,
        }

    args.destination.mkdir(parents=True, exist_ok=True)
    (args.destination / "manifest.json").write_text(
        json.dumps(output_manifest, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
