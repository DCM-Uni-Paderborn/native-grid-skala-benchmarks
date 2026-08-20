#!/usr/bin/env python3
"""Prepare the small production-protocol cutoff validation set."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


ACONF = {
    "ttt": "053d252e4457117cde6b791ea1b292d848da0bc6",
    "gtg": "1e98e5428039d926bcd2f3279b449c25eb96ffbf",
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


def replace_once(text: str, pattern: str, replacement: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE | re.DOTALL)
    if count != 1:
        raise ValueError(f"Expected one match for {pattern!r}, found {count}")
    return updated


def transform(text: str, cutoff: int, project: str, *, ae: bool) -> str:
    text = replace_once(text, r"^      &SCF\n.*?^      &END SCF\n", SCF_BLOCK)
    text = replace_once(text, r"(?m)^(\s*CUTOFF\s+)\S+", rf"\g<1>{cutoff}.0")
    text = replace_once(text, r"(?m)^(\s*PROJECT_NAME\s+)\S+", rf"\g<1>{project}")
    if ae:
        text = text.replace("QZVPP-MOLOPT-PBE-ae", "TZVPP-MOLOPT-PBE-ae")
    return text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, help="Prepared diet-inputs directory")
    parser.add_argument("destination", type=Path)
    parser.add_argument(
        "--routes",
        default="ae,direct,gapw-paw",
        help="Comma-separated subset of ae,direct,gapw-paw",
    )
    args = parser.parse_args()
    routes = set(args.routes.split(","))
    unknown = routes - {"ae", "direct", "gapw-paw"}
    if unknown:
        raise ValueError(f"Unknown routes: {', '.join(sorted(unknown))}")

    args.destination.mkdir(parents=True, exist_ok=True)
    for conformer, digest in ACONF.items():
        if "ae" in routes:
            ae_source = args.source / "gapw-ae" / f"{digest}.inp"
            for cutoff in (400, 600, 800):
                label = f"aconf8-gapw-ae-tzvpp-{cutoff}-{conformer}"
                run_dir = args.destination / label
                run_dir.mkdir(exist_ok=True)
                output = transform(ae_source.read_text(), cutoff, label.replace("-", "_"), ae=True)
                (run_dir / "input.inp").write_text(output)

        if "direct" in routes:
            direct_source = args.source / "gapw-gth-valence" / f"{digest}.inp"
            label = f"aconf8-gapw-gth-direct-800-{conformer}"
            run_dir = args.destination / label
            run_dir.mkdir(exist_ok=True)
            output = transform(direct_source.read_text(), 800, label.replace("-", "_"), ae=False)
            (run_dir / "input.inp").write_text(output)

        if "gapw-paw" in routes:
            paw_source = args.source / "gapw-gth-paw" / f"{digest}.inp"
            for cutoff in (400, 600, 800):
                label = f"aconf8-gapw-gth-paw-{cutoff}-{conformer}"
                run_dir = args.destination / label
                run_dir.mkdir(exist_ok=True)
                output = transform(paw_source.read_text(), cutoff, label.replace("-", "_"), ae=False)
                (run_dir / "input.inp").write_text(output)


if __name__ == "__main__":
    main()
