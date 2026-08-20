#!/usr/bin/env python3
"""Audit generated native-Skala dietGMTKN55 production inputs."""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path


HEAVY_GTH_ELEMENTS = {"Bi", "I", "Te"}


def kind_blocks(text: str) -> list[tuple[str, str]]:
    return re.findall(r"(?ms)^\s*&KIND\s+(\S+)\n(.*?)^\s*&END KIND\s*$", text)


def one_match(text: str, pattern: str, label: str) -> str:
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if len(matches) != 1:
        raise ValueError(f"expected one {label}, found {len(matches)}")
    return matches[0]


def audit_file(path: Path, route: str, protocol: dict[str, object]) -> None:
    text = path.read_text()
    method = one_match(text, r"^\s*METHOD\s+(GAPW(?:_XC)?)$", "QS method")
    cutoff = float(one_match(text, r"^\s*CUTOFF\s+([0-9.]+)$", "cutoff"))
    representations = re.findall(
        r"^\s*PSEUDOPOTENTIAL_GAPW_REPRESENTATION\s+(\S+)$", text, flags=re.MULTILINE
    )

    expected_method = "GAPW_XC" if route == "gapwxc-gth-one-center" else "GAPW"
    expected_cutoff = float(
        protocol["gapwxc_gth_cutoff_ry"]
        if route == "gapwxc-gth-one-center"
        else protocol["shared_ae_cutoff_ry"]
    )
    if method != expected_method:
        raise ValueError(f"method {method}, expected {expected_method}")
    if cutoff != expected_cutoff:
        raise ValueError(f"cutoff {cutoff}, expected {expected_cutoff}")
    if "GAPW_ACCURATE_XCINT .TRUE." not in text:
        raise ValueError("GAPW_ACCURATE_XCINT is not enabled")
    if len(re.findall(r"^\s*EPS_SCF\s+1e-07$", text, flags=re.MULTILINE)) != 2:
        raise ValueError("inner and outer EPS_SCF are not both 1e-07")
    if "MINIMIZER DIIS" not in text or "PRECONDITIONER FULL_ALL" not in text:
        raise ValueError("OT DIIS/FULL_ALL protocol is incomplete")
    if "STEPSIZE 0.10" not in text:
        raise ValueError("OT step size is not 0.10")

    blocks = kind_blocks(text)
    if not blocks:
        raise ValueError("no KIND blocks")

    if route == "gapwxc-gth-one-center":
        if representations != ["PAW_ONE_CENTER"]:
            raise ValueError(f"representations {representations}, expected PAW_ONE_CENTER")
        for element, body in blocks:
            if re.search(r"^\s*POTENTIAL\s+ALL$", body, flags=re.MULTILINE):
                raise ValueError(f"all-electron potential in all-GTH route for {element}")
            if "TZV2P-MOLOPT-PBE-GTH" not in body:
                raise ValueError(f"non-UZH-GTH basis for {element}")
        return

    expected_representation = {
        "hybrid-direct-heavy": "DIRECT_VALENCE",
        "hybrid-one-center-heavy": "PAW_ONE_CENTER",
    }.get(route)
    if expected_representation is not None and representations != [expected_representation]:
        raise ValueError(
            f"representations {representations}, expected {expected_representation}"
        )

    for element, body in blocks:
        potential_all = bool(re.search(r"^\s*POTENTIAL\s+ALL$", body, flags=re.MULTILINE))
        if element in HEAVY_GTH_ELEMENTS:
            if potential_all or "TZV2P-MOLOPT-PBE-GTH" not in body:
                raise ValueError(f"heavy kind {element} is not UZH GTH")
        else:
            if not potential_all or "TZVPP-MOLOPT-PBE-ae" not in body:
                raise ValueError(f"covered kind {element} is not UZH all electron")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()

    manifest = json.loads((args.root / "manifest.json").read_text())
    counts: Counter[str] = Counter()
    errors: list[str] = []
    for digest, record in manifest["species"].items():
        for route, relative in record["routes"].items():
            counts[route] += 1
            path = args.root / relative
            try:
                audit_file(path, route, manifest["protocol"])
            except (OSError, ValueError) as exc:
                errors.append(f"{digest} {route}: {exc}")

    print(json.dumps(dict(sorted(counts.items())), indent=2))
    if errors:
        print("\n".join(errors))
        raise SystemExit(1)
    print(f"Validated {sum(counts.values())} generated inputs.")


if __name__ == "__main__":
    main()
