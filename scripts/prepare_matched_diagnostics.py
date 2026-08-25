#!/usr/bin/env python3
"""Prepare matched GauXC and native-grid Skala diagnostic calculations."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path


REACTIONS = {
    "ACONF-5": {
        "053d252e4457117cde6b791ea1b292d848da0bc6": -1,
        "18523213f9133ced8440d2c7b30dcb631edfce6f": 1,
    },
    "ACONF-8": {
        "053d252e4457117cde6b791ea1b292d848da0bc6": -1,
        "1e98e5428039d926bcd2f3279b449c25eb96ffbf": 1,
    },
    "CARBHB12-11": {
        "cb06d2a1d55acb75f8b039f41acbb3165684a9ef": -1,
        "435dcfec2baff1096f952206ad475f6f9b4242b2": 1,
        "6007ed0531318684dce37c642e2f2aacbb6d349d": 1,
    },
    "G2RC-15": {
        "5f7357571186a3960d758f7969a9dfa60688788e": -3,
        "f539e6c9a52a881fdd9cc6018b7b03df2f6d56e3": 2,
        "00333493a55a481e863241b5357e055f4054bb24": -1,
    },
}


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


GAUXC_BLOCK = """            &GAUXC
               GRID FINE
               MODEL SKALA
               NATIVE_GRID_USE_CUDA .TRUE.
               PRUNING_SCHEME ROBUST
               RADIAL_QUADRATURE MURAKNOWLES
            &END GAUXC"""


@dataclass(frozen=True)
class Variant:
    source_route: str
    cutoff: int
    padding: int
    radial: int
    lebedev: int
    method: str
    backend: str
    poisson_solver: str = "ANALYTIC"


VARIANTS = {
    "ae-native-p12-r100-l434": Variant(
        "gapw-ae", 640, 12, 100, 434, "GAPW", "native"
    ),
    "ae-gauxc-p12": Variant(
        "gapw-ae", 640, 12, 100, 434, "GAPW", "gauxc"
    ),
    "ae-gauxc-p30": Variant(
        "gapw-ae", 640, 30, 100, 434, "GAPW", "gauxc"
    ),
    "ae-gauxc-p12-wavelet": Variant(
        "gapw-ae", 640, 12, 100, 434, "GAPW", "gauxc", "WAVELET"
    ),
    "ae-native-p30-r100-l434": Variant(
        "gapw-ae", 640, 30, 100, 434, "GAPW", "native"
    ),
    "ae-native-p12-r100-l590": Variant(
        "gapw-ae", 640, 12, 100, 590, "GAPW", "native"
    ),
    "ae-native-p12-r100-l434-wavelet": Variant(
        "gapw-ae", 640, 12, 100, 434, "GAPW", "native", "WAVELET"
    ),
    "gth-native-direct-p12-r100-l434": Variant(
        "gapw-gth-valence", 640, 12, 100, 434, "GAPW", "native"
    ),
    "gth-gauxc-p12": Variant(
        "gapw-gth-valence", 640, 12, 100, 434, "GAPW", "gauxc"
    ),
    "gth-gauxc-p16": Variant(
        "gapw-gth-valence", 640, 16, 100, 434, "GAPW", "gauxc"
    ),
    "gth-gauxc-p20": Variant(
        "gapw-gth-valence", 640, 20, 100, 434, "GAPW", "gauxc"
    ),
    "gth-gauxc-p30": Variant(
        "gapw-gth-valence", 640, 30, 100, 434, "GAPW", "gauxc"
    ),
    "gth-gauxc-p12-wavelet": Variant(
        "gapw-gth-valence", 640, 12, 100, 434, "GAPW", "gauxc", "WAVELET"
    ),
    "gth-gauxc-p16-wavelet": Variant(
        "gapw-gth-valence", 640, 16, 100, 434, "GAPW", "gauxc", "WAVELET"
    ),
    "gth-gauxc-p20-wavelet": Variant(
        "gapw-gth-valence", 640, 20, 100, 434, "GAPW", "gauxc", "WAVELET"
    ),
    "gth-gauxc-p30-wavelet": Variant(
        "gapw-gth-valence", 640, 30, 100, 434, "GAPW", "gauxc", "WAVELET"
    ),
    "gth-native-one-center-gapw-p12": Variant(
        "gapw-gth-paw", 640, 12, 100, 434, "GAPW", "native"
    ),
    "gth-native-one-center-gapwxc-p12": Variant(
        "gapw-gth-paw", 400, 12, 100, 434, "GAPW_XC", "native"
    ),
    "gth-native-one-center-gapwxc-p16": Variant(
        "gapw-gth-paw", 400, 16, 100, 434, "GAPW_XC", "native"
    ),
    "gth-native-one-center-gapwxc-p20": Variant(
        "gapw-gth-paw", 400, 20, 100, 434, "GAPW_XC", "native"
    ),
    "gth-native-one-center-gapwxc-p30": Variant(
        "gapw-gth-paw", 400, 30, 100, 434, "GAPW_XC", "native"
    ),
    "gth-native-one-center-gapwxc-p12-wavelet": Variant(
        "gapw-gth-paw", 400, 12, 100, 434, "GAPW_XC", "native", "WAVELET"
    ),
    "gth-native-one-center-gapwxc-p16-wavelet": Variant(
        "gapw-gth-paw", 400, 16, 100, 434, "GAPW_XC", "native", "WAVELET"
    ),
    "gth-native-one-center-gapwxc-p20-wavelet": Variant(
        "gapw-gth-paw", 400, 20, 100, 434, "GAPW_XC", "native", "WAVELET"
    ),
    "gth-native-one-center-gapwxc-p30-wavelet": Variant(
        "gapw-gth-paw", 400, 30, 100, 434, "GAPW_XC", "native", "WAVELET"
    ),
}

for padding in (24, 25, 26):
    VARIANTS.setdefault(
        f"ae-native-p{padding}-r100-l434",
        Variant("gapw-ae", 640, padding, 100, 434, "GAPW", "native"),
    )
    VARIANTS.setdefault(
        f"ae-gauxc-p{padding}",
        Variant("gapw-ae", 640, padding, 100, 434, "GAPW", "gauxc"),
    )


PADDING_SCAN = (4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 25, 26, 28, 30)
PADDING_FAMILIES = {
    "gth-gauxc": ("gapw-gth-valence", 640, "GAPW", "gauxc"),
    "gth-native-one-center-gapwxc": (
        "gapw-gth-paw",
        400,
        "GAPW_XC",
        "native",
    ),
}

for padding in PADDING_SCAN:
    for prefix, (source_route, cutoff, method, backend) in PADDING_FAMILIES.items():
        VARIANTS.setdefault(
            f"{prefix}-p{padding}",
            Variant(source_route, cutoff, padding, 100, 434, method, backend),
        )
        VARIANTS.setdefault(
            f"{prefix}-p{padding}-wavelet",
            Variant(
                source_route,
                cutoff,
                padding,
                100,
                434,
                method,
                backend,
                "WAVELET",
            ),
        )



def replace_once(
    text: str, pattern: str, replacement: str, *, dotall: bool = False
) -> str:
    flags = re.MULTILINE | (re.DOTALL if dotall else 0)
    updated, count = re.subn(
        pattern, replacement, text, count=1, flags=flags
    )
    if count != 1:
        raise ValueError(f"Expected one match for {pattern!r}, found {count}")
    return updated


def set_total_padding(text: str, padding: int, *, cubic: bool) -> str:
    coord_match = re.search(r"(?ms)^\s*&COORD\s*$\n(.*?)^\s*&END COORD\s*$", text)
    if coord_match is None:
        raise ValueError("Coordinate block not found")

    coordinates = []
    for line in coord_match.group(1).splitlines():
        fields = line.split()
        if len(fields) >= 4:
            try:
                coordinates.append(tuple(float(value) for value in fields[-3:]))
            except ValueError:
                continue
    if not coordinates:
        raise ValueError("No Cartesian coordinates found")

    extents = [
        max(coordinate[axis] for coordinate in coordinates)
        - min(coordinate[axis] for coordinate in coordinates)
        for axis in range(3)
    ]
    requested_lengths = [extent + padding for extent in extents]
    if cubic:
        side = max(requested_lengths)
        lengths = [side, side, side]
    else:
        lengths = requested_lengths
    abc = " ".join(f"{length:.8f}" for length in lengths)
    return replace_once(text, r"(?m)^(\s*ABC\s+).*$", rf"\g<1>{abc}")


def ensure_center_coordinates(text: str) -> str:
    if re.search(r"(?m)^\s*&CENTER_COORDINATES(?:\s+.*)?$", text):
        return text

    topology = re.search(
        r"(?ms)^(?P<indent>[ \t]*)&TOPOLOGY\s*$.*?^"
        r"(?P=indent)&END TOPOLOGY\s*$",
        text,
    )
    if topology is not None:
        indent = topology.group("indent")
        end_start = text.rfind("&END TOPOLOGY", topology.start(), topology.end())
        center = (
            f"{indent}   &CENTER_COORDINATES\n"
            f"{indent}   &END CENTER_COORDINATES\n"
        )
        return text[:end_start] + center + text[end_start:]

    coord = re.search(r"(?m)^(?P<indent>[ \t]*)&COORD\s*$", text)
    if coord is None:
        raise ValueError("Coordinate block not found for molecular centering")
    indent = coord.group("indent")
    center = (
        f"{indent}&TOPOLOGY\n"
        f"{indent}   &CENTER_COORDINATES\n"
        f"{indent}   &END CENTER_COORDINATES\n"
        f"{indent}&END TOPOLOGY\n"
    )
    return text[:coord.start()] + center + text[coord.start():]


def extract_cell_metadata(text: str) -> dict[str, list[float]]:
    coord_match = re.search(r"(?ms)^\s*&COORD\s*$\n(.*?)^\s*&END COORD\s*$", text)
    cell_match = re.search(r"(?m)^\s*ABC\s+(\S+)\s+(\S+)\s+(\S+)\s*$", text)
    if coord_match is None or cell_match is None:
        raise ValueError("Coordinates or cell not found")
    coordinates = []
    for line in coord_match.group(1).splitlines():
        fields = line.split()
        if len(fields) >= 4:
            try:
                coordinates.append(tuple(float(value) for value in fields[-3:]))
            except ValueError:
                continue
    extents = [
        max(coordinate[axis] for coordinate in coordinates)
        - min(coordinate[axis] for coordinate in coordinates)
        for axis in range(3)
    ]
    cell = [float(value) for value in cell_match.groups()]
    return {
        "molecular_extent_angstrom": extents,
        "cell_angstrom": cell,
        "total_padding_angstrom": [
            cell_length - extent for cell_length, extent in zip(cell, extents)
        ],
        "vacuum_per_side_after_centering_angstrom": [
            0.5 * (cell_length - extent)
            for cell_length, extent in zip(cell, extents)
        ],
    }


def prepare_input(text: str, variant: Variant, project: str) -> str:
    text = ensure_center_coordinates(text)
    text = replace_once(
        text, r"^      &SCF\n.*?^      &END SCF\n", SCF_BLOCK, dotall=True
    )
    text = replace_once(text, r"(?m)^(\s*CUTOFF\s+)\S+", rf"\g<1>{variant.cutoff}.0")
    text = replace_once(text, r"(?m)^(\s*PROJECT_NAME\s+)\S+", rf"\g<1>{project}")
    text = replace_once(
        text, r"(?m)^(\s*METHOD\s+)(?:GAPW|GAPW_XC)$", rf"\g<1>{variant.method}"
    )
    text = replace_once(
        text,
        r"(?m)^(\s*POISSON_SOLVER\s+)\S+$",
        rf"\g<1>{variant.poisson_solver}",
    )
    text = re.sub(
        r"(?m)^(\s*RADIAL_GRID\s+)\S+", rf"\g<1>{variant.radial}", text
    )
    text = re.sub(
        r"(?m)^(\s*LEBEDEV_GRID\s+)\S+", rf"\g<1>{variant.lebedev}", text
    )
    text = set_total_padding(
        text,
        variant.padding,
        cubic=variant.poisson_solver == "WAVELET",
    )

    if variant.source_route == "gapw-ae":
        text = text.replace("QZVPP-MOLOPT-PBE-ae", "TZVPP-MOLOPT-PBE-ae")

    if variant.backend == "gauxc":
        text = replace_once(
            text,
            r"^([ \t]*)&GAUXC[^\n]*\n.*?^\1&END GAUXC[^\n]*$",
            GAUXC_BLOCK,
            dotall=True,
        )

    if "GAPW_ACCURATE_XCINT .TRUE." not in text:
        raise ValueError("GAPW_ACCURATE_XCINT must remain enabled")
    return text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, help="Directory containing source route folders")
    parser.add_argument("destination", type=Path)
    parser.add_argument(
        "--reaction",
        action="append",
        choices=sorted(REACTIONS),
        help="Reaction to include; repeat as needed (default: all diagnostic reactions)",
    )
    parser.add_argument(
        "--variant",
        action="append",
        choices=sorted(VARIANTS),
        help="Variant to include; repeat as needed (default: all variants)",
    )
    args = parser.parse_args()

    reaction_names = args.reaction or sorted(REACTIONS)
    variant_names = args.variant or sorted(VARIANTS)
    species = sorted(
        {
            digest
            for reaction in reaction_names
            for digest in REACTIONS[reaction]
        }
    )

    manifest = {
        "reactions": {name: REACTIONS[name] for name in reaction_names},
        "variants": {name: asdict(VARIANTS[name]) for name in variant_names},
        "species": species,
        "cells": {},
    }
    args.destination.mkdir(parents=True, exist_ok=True)

    for variant_name in variant_names:
        variant = VARIANTS[variant_name]
        manifest["cells"][variant_name] = {}
        for digest in species:
            source = args.source / variant.source_route / f"{digest}.inp"
            if not source.exists():
                raise FileNotFoundError(source)
            run_dir = args.destination / variant_name / digest
            run_dir.mkdir(parents=True, exist_ok=True)
            project = f"diag_{variant_name.replace('-', '_')}_{digest[:8]}"
            prepared = prepare_input(source.read_text(), variant, project)
            (run_dir / "input.inp").write_text(prepared)
            manifest["cells"][variant_name][digest] = extract_cell_metadata(prepared)

    (args.destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
