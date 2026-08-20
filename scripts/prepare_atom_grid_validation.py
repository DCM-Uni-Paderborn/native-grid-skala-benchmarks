#!/usr/bin/env python3
"""Prepare paired ACONF-8 atom-grid convergence calculations."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


ACONF = {
    "ttt": "053d252e4457117cde6b791ea1b292d848da0bc6",
    "gtg": "1e98e5428039d926bcd2f3279b449c25eb96ffbf",
}

GRIDS = (
    (50, 110),
    (50, 194),
    (50, 302),
    (50, 434),
    (100, 434),
    (150, 434),
    (150, 590),
)

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


def replace_once(text: str, pattern: str, replacement: str) -> str:
    updated, count = re.subn(
        pattern, replacement, text, count=1, flags=re.MULTILINE | re.DOTALL
    )
    if count != 1:
        raise ValueError(f"Expected one match for {pattern!r}, found {count}")
    return updated


def transform(
    text: str,
    *,
    route: str,
    cutoff: int,
    radial: int,
    lebedev: int,
    project: str,
) -> str:
    text = replace_once(text, r"^      &SCF\n.*?^      &END SCF\n", SCF_BLOCK)
    text = replace_once(text, r"(?m)^(\s*CUTOFF\s+)\S+", rf"\g<1>{cutoff}.0")
    text = replace_once(text, r"(?m)^(\s*PROJECT_NAME\s+)\S+", rf"\g<1>{project}")
    text = re.sub(r"(?m)^(\s*RADIAL_GRID\s+)\S+", rf"\g<1>{radial}", text)
    text = re.sub(r"(?m)^(\s*LEBEDEV_GRID\s+)\S+", rf"\g<1>{lebedev}", text)

    if route == "ae":
        text = text.replace("QZVPP-MOLOPT-PBE-ae", "TZVPP-MOLOPT-PBE-ae")
    elif route == "gapwxc-paw":
        text = replace_once(
            text, r"(?m)^(\s*METHOD\s+)(?:GAPW|GAPW_XC)$", r"\g<1>GAPW_XC"
        )
    else:
        raise ValueError(f"Unknown route: {route}")

    if "GAPW_ACCURATE_XCINT .TRUE." not in text:
        raise ValueError("GAPW_ACCURATE_XCINT must be enabled")
    return text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, help="Directory containing diet-inputs")
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()

    routes = {
        "ae": ("gapw-ae", 800),
        "gapwxc-paw": ("gapw-gth-paw", 400),
    }
    for route, (source_route, cutoff) in routes.items():
        for conformer, digest in ACONF.items():
            source = args.source / source_route / f"{digest}.inp"
            for radial, lebedev in GRIDS:
                label = f"aconf8-{route}-r{radial}-l{lebedev}-{conformer}"
                run_dir = args.destination / label
                run_dir.mkdir(parents=True, exist_ok=True)
                output = transform(
                    source.read_text(),
                    route=route,
                    cutoff=cutoff,
                    radial=radial,
                    lebedev=lebedev,
                    project=label.replace("-", "_"),
                )
                (run_dir / "input.inp").write_text(output)


if __name__ == "__main__":
    main()
