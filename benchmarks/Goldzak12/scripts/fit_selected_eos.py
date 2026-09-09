#!/usr/bin/env python3
"""Fit the selected structural comparison without requiring atomic energies."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from fit_eos import HA_PER_A3_TO_GPA, HA_TO_EV, ROOT, fit_curve, read_rows, write_rows


def select_curves(rows, selection):
    indexed = {}
    for row in rows:
        key = (row["method"], row["solid"], row["point"])
        if key in indexed:
            raise ValueError(f"Duplicate source point: {key}")
        indexed[key] = row
    curves = {}
    for solid in selection["solids"]:
        points = selection["point_overrides"].get(solid, selection["default_points"])
        if len(points) < 5 or len(set(points)) != len(points):
            raise ValueError(f"Invalid point selection: {solid}")
        for method in selection["methods"]:
            source_method = method
            if solid in selection["shared_ae_solids"] and method.startswith("hybrid-"):
                source_method = "gapw-ae"
            chosen = []
            for point in points:
                key = (source_method, solid, point)
                if key not in indexed:
                    raise ValueError(f"Missing selected source point: {key}")
                row = indexed[key]
                if row["execution_accepted"].lower() != "true":
                    raise ValueError(f"Source not execution-accepted: {key}")
                values = [float(row["volume_A3"]), float(row["energy_Ha"])]
                if not np.all(np.isfinite(values)) or values[0] <= 0:
                    raise ValueError(f"Invalid energy or volume: {key}")
                chosen.append(row)
            chosen.sort(key=lambda row: float(row["volume_A3"]))
            if len({float(row["volume_A3"]) for row in chosen}) != len(chosen):
                raise ValueError(f"Duplicate volume: {source_method}/{solid}")
            curves[(method, solid)] = (source_method, chosen)
    # Paired representation errors require exactly the same physical volumes.
    volumes_by_solid = defaultdict(set)
    for (_, solid), (_, chosen) in curves.items():
        volumes_by_solid[solid].add(tuple(float(row["volume_A3"]) for row in chosen))
    if any(len(values) != 1 for values in volumes_by_solid.values()):
        raise ValueError("Methods do not share identical selected volumes")
    return curves


def analyze(rows, selection):
    curves = select_curves(rows, selection)
    results = []
    cache = {}
    for (method, solid), (source_method, chosen) in sorted(curves.items()):
        source_key = (source_method, solid)
        if source_key not in cache:
            parameters, _, rms = fit_curve(chosen)
            energy, volume, bulk, derivative = parameters
            minimum = min(float(row["volume_A3"]) for row in chosen)
            maximum = max(float(row["volume_A3"]) for row in chosen)
            cache[source_key] = {
                "n_points": len(chosen),
                "E0_Ha_cell": float(energy),
                "V0_A3_cell": float(volume),
                "a0_A": float(volume ** (1.0 / 3.0)),
                "B0_GPa": float(bulk * HA_PER_A3_TO_GPA),
                "B0_prime": float(derivative),
                "fit_rms_meV_atom": float(rms * HA_TO_EV * 1000.0 / 8.0),
                "minimum_bracketed": bool(minimum < volume < maximum),
            }
        results.append({
            "method": method, "solid": solid, "source_method": source_method,
            **cache[source_key], "quality_status": "numerical_review_pending",
        })
    unique_points = {
        (source, solid, row["point"])
        for (_, solid), (source, chosen) in curves.items() for row in chosen
    }
    summary = {
        "selection_date": selection.get("selection_date", "2026-09-05"),
        "solids": selection["solids"],
        "unique_selected_points": len(unique_points),
        "independent_curves": len(cache),
        "method_solid_comparisons": len(results),
        "all_minima_bracketed": all(row["minimum_bracketed"] for row in results),
        "maximum_fit_rms_meV_atom": max(row["fit_rms_meV_atom"] for row in results),
        "quality_release": False,
        "cohesive_energies_computed": False,
        "remaining_checks": selection.get("remaining_checks", ["Numerical qualifications are described in the PCCP ESI."]),
    }
    return results, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bulk", type=Path, default=ROOT / "results/eos-selected-source.csv")
    parser.add_argument("--selection", type=Path, default=ROOT / "protocol/paper-eos-selection.json")
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    results, summary = analyze(read_rows(args.bulk), json.loads(args.selection.read_text()))
    write_rows(args.output / "eos-selected-fits.csv", list(results[0]), results)
    (args.output / "eos-selected-summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="ascii"
    )
    print(json.dumps(summary, allow_nan=False))


if __name__ == "__main__":
    main()
