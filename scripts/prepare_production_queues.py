#!/usr/bin/env python3
"""Create disjoint host queues for the frozen dietGMTKN55 production."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


COORD_RE = re.compile(r"^\s*[A-Z][a-z]?\s+[-+0-9]")


def cost_key(root: Path, relative: str) -> tuple[int, int, str]:
    text = (root / relative / "input.inp").read_text(encoding="utf-8")
    atoms = sum(bool(COORD_RE.match(line)) for line in text.splitlines())
    return atoms, len(text), relative


def write_list(path: Path, entries: list[str]) -> None:
    path.write_text("".join(f"{entry}\n" for entry in entries), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("production_root", type=Path)
    parser.add_argument("--spark", type=int, default=32)
    parser.add_argument("--terok", type=int, default=144)
    args = parser.parse_args()

    root = args.production_root.resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    routes = sorted(
        {
            relative.removesuffix("/input.inp")
            for species in manifest["species"].values()
            for relative in species["routes"].values()
        }
    )
    if len(routes) != 480:
        raise SystemExit(f"Expected 480 routes, found {len(routes)}")

    gth = sorted(
        (route for route in routes if route.startswith("gapwxc-gth-one-center/")),
        key=lambda route: cost_key(root, route),
    )
    shared_ae = sorted(
        (route for route in routes if route.startswith("shared-gapw-ae/")),
        key=lambda route: cost_key(root, route),
    )
    hybrid = sorted(
        (route for route in routes if route.startswith("hybrid-")),
        key=lambda route: cost_key(root, route),
    )
    if (len(gth), len(shared_ae), len(hybrid)) != (236, 228, 16):
        raise SystemExit(
            "Unexpected route counts: "
            f"GTH={len(gth)}, shared-AE={len(shared_ae)}, hybrid={len(hybrid)}"
        )

    spark = gth[: args.spark]
    terok = gth[args.spark : args.spark + args.terok]
    rosi_gth = gth[args.spark + args.terok :]
    queues = {
        "spark-gpu.txt": spark,
        "terok-cpu.txt": terok,
        "rosi-gpu-gth.txt": rosi_gth,
        "rosi-gpu-hybrid.txt": hybrid,
        "rosi-cpu-ae.txt": shared_ae,
    }

    queue_dir = root / "queues"
    queue_dir.mkdir(exist_ok=True)
    for name, entries in queues.items():
        write_list(queue_dir / name, entries)
    write_list(queue_dir / "terok-cpu-early.txt", terok[:6])
    write_list(queue_dir / "terok-cpu-main.txt", terok[6:])

    assigned = [entry for entries in queues.values() for entry in entries]
    if len(assigned) != len(set(assigned)) or set(assigned) != set(routes):
        raise SystemExit("Queues are not a disjoint and complete partition")

    summary = {
        "production_root": str(root),
        "total": len(assigned),
        "queues": {name: len(entries) for name, entries in queues.items()},
        "terok_operational_split": {"early": 6, "main": len(terok) - 6},
    }
    (queue_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
