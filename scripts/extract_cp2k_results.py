#!/usr/bin/env python3
"""Extract compact convergence metadata from CP2K run directories."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def last_match(pattern: str, text: str, flags: int = 0) -> str:
    matches = re.findall(pattern, text, flags)
    if not matches:
        return ""
    value = matches[-1]
    return value if isinstance(value, str) else " ".join(value)


def unique_matches(pattern: str, text: str, flags: int = 0) -> str:
    return ";".join(dict.fromkeys(re.findall(pattern, text, flags)))


def keyword_enabled(keyword: str, text: str) -> str:
    value = last_match(
        rf"^\s*{re.escape(keyword)}\s+(\S+)$", text, re.MULTILINE
    ).upper()
    return str(value in {"T", ".TRUE.", "TRUE", "YES", "ON"}).lower()


def parse_run(output_path: Path, root: Path) -> dict[str, str]:
    run_dir = output_path.parent
    input_path = run_dir / "input.inp"
    time_path = run_dir / "time.txt"
    output = output_path.read_text(errors="replace")
    input_text = input_path.read_text(errors="replace") if input_path.exists() else ""
    time_text = time_path.read_text(errors="replace") if time_path.exists() else ""

    return {
        "path": str(run_dir.relative_to(root)),
        "method": last_match(r"^\s*METHOD\s+(GAPW(?:_XC)?|GPW)$", input_text, re.MULTILINE),
        "representation": last_match(
            r"^\s*PSEUDOPOTENTIAL_GAPW_REPRESENTATION\s+(\S+)$",
            input_text,
            re.MULTILINE,
        ),
        "cutoff_ry": last_match(r"^\s*CUTOFF\s+([0-9.]+)$", input_text, re.MULTILINE),
        "relative_cutoff_ry": last_match(
            r"^\s*REL_CUTOFF\s+([0-9.]+)$", input_text, re.MULTILINE
        ),
        "radial_grid": unique_matches(r"^\s*RADIAL_GRID\s+(\d+)$", input_text, re.MULTILINE),
        "lebedev_grid": unique_matches(r"^\s*LEBEDEV_GRID\s+(\d+)$", input_text, re.MULTILINE),
        "hard_exp_radius": unique_matches(
            r"^\s*HARD_EXP_RADIUS\s+([0-9.]+)$", input_text, re.MULTILINE
        ),
        "accurate_xcint": keyword_enabled("GAPW_ACCURATE_XCINT", input_text),
        "eps_scf": last_match(r"^\s*EPS_SCF\s+(\S+)$", input_text, re.MULTILINE),
        "converged": str("SCF run converged" in output).lower(),
        "scf_steps": last_match(r"SCF run converged in\s+(\d+) steps", output),
        "total_energy_ha": last_match(
            r"ENERGY\| Total FORCE_EVAL.*?(-?\d+\.\d+)", output
        ),
        "skala_xc_energy_ha": last_match(
            r"Active atom-composite XC energy\s+([+-]?\d+\.\d+(?:E[+-]?\d+)?)",
            output,
        ),
        "composite_electrons": last_match(
            r"Active atom-composite electrons\s+([+-]?\d+\.\d+(?:E[+-]?\d+)?)",
            output,
        ),
        "wall_time": last_match(
            r"^.*Elapsed \(wall clock\) time \([^)]*\):\s*([0-9:.]+)\s*$",
            time_text,
            re.MULTILINE,
        ),
        "max_rss_kb": last_match(r"Maximum resident set size \(kbytes\):\s*(\d+)", time_text),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    rows = [parse_run(path, args.root) for path in sorted(args.root.rglob("output.out"))]
    if not rows:
        raise SystemExit(f"No output.out files found below {args.root}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
