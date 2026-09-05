#!/usr/bin/env python3
"""Fit Goldzak12 equations of state and build paper-facing error tables."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
HA_TO_EV = 27.211386245988
HA_PER_A3_TO_GPA = 4359.7447222071


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def birch_murnaghan(
    volume: np.ndarray,
    energy0: float,
    volume0: float,
    bulk0: float,
    bulk_prime: float,
) -> np.ndarray:
    x = (volume0 / volume) ** (2.0 / 3.0)
    eta = x - 1.0
    return energy0 + (9.0 * volume0 * bulk0 / 16.0) * (
        bulk_prime * eta**3 + (6.0 - 4.0 * x) * eta**2
    )


def fit_curve(rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray, float]:
    volumes = np.array([float(row["volume_A3"]) for row in rows])
    energies = np.array([float(row["energy_Ha"]) for row in rows])

    def linear_fit(volume0: float) -> tuple[float, np.ndarray]:
        x = (volume0 / volumes) ** (2.0 / 3.0)
        eta = x - 1.0
        prefactor = 9.0 * volume0 / 16.0
        design = np.column_stack(
            (
                np.ones_like(volumes),
                prefactor * eta**3,
                prefactor * eta**2 * (6.0 - 4.0 * x),
            )
        )
        coefficients, _, _, _ = np.linalg.lstsq(design, energies, rcond=None)
        energy0, bulk_times_prime, bulk0 = coefficients
        if bulk0 <= 0.0:
            return math.inf, np.array((energy0, volume0, bulk0, math.nan))
        bulk_prime = bulk_times_prime / bulk0
        if not 0.0 < bulk_prime < 15.0:
            return math.inf, np.array((energy0, volume0, bulk0, bulk_prime))
        parameters = np.array((energy0, volume0, bulk0, bulk_prime))
        residual = energies - birch_murnaghan(volumes, *parameters)
        return float(np.dot(residual, residual)), parameters

    lower = volumes.min() * 0.8
    upper = volumes.max() * 1.2
    scan = np.linspace(lower, upper, 2001)
    objectives = np.array([linear_fit(value)[0] for value in scan])
    best = int(np.argmin(objectives))
    left = scan[max(best - 1, 0)]
    right = scan[min(best + 1, scan.size - 1)]
    golden = (math.sqrt(5.0) - 1.0) / 2.0
    c = right - golden * (right - left)
    d = left + golden * (right - left)
    fc = linear_fit(c)[0]
    fd = linear_fit(d)[0]
    for _ in range(100):
        if abs(right - left) < 1.0e-12 * max(1.0, abs(left), abs(right)):
            break
        if fc < fd:
            right, d, fd = d, c, fc
            c = right - golden * (right - left)
            fc = linear_fit(c)[0]
        else:
            left, c, fc = c, d, fd
            d = left + golden * (right - left)
            fd = linear_fit(d)[0]
    _, parameters = linear_fit(0.5 * (left + right))
    if not np.all(np.isfinite(parameters)):
        raise ValueError("Birch-Murnaghan fit did not yield physical parameters")

    residual = energies - birch_murnaghan(volumes, *parameters)
    rms = float(np.sqrt(np.mean(residual**2)))
    jacobian = np.empty((volumes.size, 4))
    for column, parameter in enumerate(parameters):
        step = 1.0e-6 * max(abs(parameter), 1.0)
        plus = parameters.copy()
        minus = parameters.copy()
        plus[column] += step
        minus[column] -= step
        jacobian[:, column] = (
            birch_murnaghan(volumes, *plus)
            - birch_murnaghan(volumes, *minus)
        ) / (2.0 * step)
    degrees = max(volumes.size - parameters.size, 1)
    variance = float(np.dot(residual, residual) / degrees)
    covariance = variance * np.linalg.pinv(jacobian.T @ jacobian)
    return parameters, covariance, rms


def composition(system: dict[str, str]) -> dict[str, int]:
    if system["sublattice_a"] == system["sublattice_b"]:
        return {system["sublattice_a"]: 8}
    return {system["sublattice_a"]: 4, system["sublattice_b"]: 4}


def aggregate(values: list[float]) -> tuple[float, float, float, float]:
    array = np.asarray(values, dtype=float)
    return (
        float(array.mean()),
        float(np.abs(array).mean()),
        float(np.sqrt(np.mean(array**2))),
        float(np.max(np.abs(array))),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bulk",
        type=Path,
        default=ROOT / "results" / "validated-bulk-energies.csv",
    )
    parser.add_argument(
        "--atoms",
        type=Path,
        default=ROOT / "results" / "validated-atomic-energies.csv",
    )
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    args = parser.parse_args()

    systems = {
        row["solid"]: row for row in read_rows(ROOT / "protocol" / "systems.csv")
    }
    references = {
        row["solid"]: row
        for row in read_rows(ROOT / "reference" / "goldzak2022.csv")
        if row["method"] == "experiment"
    }
    atomic_energies = {
        (row["method"], row["element"]): float(row["energy_Ha"])
        for row in read_rows(args.atoms)
        if row["accepted"].lower() == "true"
    }
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in read_rows(args.bulk):
        if row["accepted"].lower() == "true":
            grouped[(row["method"], row["solid"])].append(row)

    fit_rows: list[dict[str, object]] = []
    error_rows: list[dict[str, object]] = []
    errors_by_method: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for (method, solid), rows in sorted(grouped.items()):
        if len(rows) != 10:
            raise ValueError(f"{method}/{solid}: expected 10 accepted volumes, got {len(rows)}")
        parameters, covariance, fit_rms = fit_curve(rows)
        energy0, volume0, bulk0, bulk_prime = parameters
        lattice0 = volume0 ** (1.0 / 3.0)
        bulk_gpa = bulk0 * HA_PER_A3_TO_GPA
        atom_sum = 0.0
        for element, count in composition(systems[solid]).items():
            try:
                atom_sum += count * atomic_energies[(method, element)]
            except KeyError as exc:
                raise ValueError(f"Missing atom energy for {method}/{element}") from exc
        cohesive = (atom_sum - energy0) * HA_TO_EV / 8.0
        reference = references[solid]
        delta_a = lattice0 - float(reference["a_A"])
        delta_bulk = bulk_gpa - float(reference["B0_GPa"])
        delta_cohesive = cohesive - float(reference["E_coh_eV_atom"])
        fit_rows.append(
            {
                "method": method,
                "solid": solid,
                "E0_Ha_cell": f"{energy0:.15f}",
                "V0_A3_cell": f"{volume0:.10f}",
                "a0_A": f"{lattice0:.8f}",
                "B0_GPa": f"{bulk_gpa:.6f}",
                "B0_prime": f"{bulk_prime:.6f}",
                "E_coh_eV_atom": f"{cohesive:.8f}",
                "fit_rms_Ha": f"{fit_rms:.8e}",
                "sigma_V0_A3": f"{math.sqrt(max(covariance[1, 1], 0.0)):.8e}",
            }
        )
        error_rows.append(
            {
                "method": method,
                "solid": solid,
                "delta_a_A": f"{delta_a:.8f}",
                "delta_B0_GPa": f"{delta_bulk:.6f}",
                "delta_E_coh_eV_atom": f"{delta_cohesive:.8f}",
            }
        )
        errors_by_method[method]["a"].append(delta_a)
        errors_by_method[method]["B0"].append(delta_bulk)
        errors_by_method[method]["E_coh"].append(delta_cohesive)

    aggregate_rows: list[dict[str, object]] = []
    for method, properties in sorted(errors_by_method.items()):
        row: dict[str, object] = {"method": method, "N": len(properties["a"])}
        for label, values in properties.items():
            me, mae, rmse, maxae = aggregate(values)
            row[f"ME_{label}"] = f"{me:.8f}"
            row[f"MAE_{label}"] = f"{mae:.8f}"
            row[f"RMSE_{label}"] = f"{rmse:.8f}"
            row[f"MaxAE_{label}"] = f"{maxae:.8f}"
        aggregate_rows.append(row)

    write_rows(
        args.output / "eos-fits.csv",
        [
            "method",
            "solid",
            "E0_Ha_cell",
            "V0_A3_cell",
            "a0_A",
            "B0_GPa",
            "B0_prime",
            "E_coh_eV_atom",
            "fit_rms_Ha",
            "sigma_V0_A3",
        ],
        fit_rows,
    )
    write_rows(
        args.output / "reference-errors.csv",
        ["method", "solid", "delta_a_A", "delta_B0_GPa", "delta_E_coh_eV_atom"],
        error_rows,
    )
    write_rows(
        args.output / "aggregate-statistics.csv",
        [
            "method",
            "N",
            "ME_a",
            "MAE_a",
            "RMSE_a",
            "MaxAE_a",
            "ME_B0",
            "MAE_B0",
            "RMSE_B0",
            "MaxAE_B0",
            "ME_E_coh",
            "MAE_E_coh",
            "RMSE_E_coh",
            "MaxAE_E_coh",
        ],
        aggregate_rows,
    )

    comparison_rows: list[dict[str, object]] = []
    for row in read_rows(ROOT / "reference" / "goldzak2022.csv"):
        comparison_rows.append(
            {
                "source": "Goldzak2022",
                "method": row["method"],
                "solid": row["solid"],
                "a_A": row["a_A"],
                "B0_GPa": row["B0_GPa"],
                "E_coh_eV_atom": row["E_coh_eV_atom"],
                "status": "published",
                "note": row["reference_note"],
            }
        )
    for row in read_rows(ROOT / "reference" / "periodic_gfn2_lc10.csv"):
        comparison_rows.append(
            {
                "source": "PeriodicGFN2Manuscript",
                "method": row["method"],
                "solid": row["solid"],
                "a_A": row["a_A"],
                "B0_GPa": "",
                "E_coh_eV_atom": row["E_coh_eV_atom"],
                "status": row["status"],
                "note": row["source"],
            }
        )
    for row in fit_rows:
        comparison_rows.append(
            {
                "source": "this work",
                "method": row["method"],
                "solid": row["solid"],
                "a_A": row["a0_A"],
                "B0_GPa": row["B0_GPa"],
                "E_coh_eV_atom": row["E_coh_eV_atom"],
                "status": "strictly validated",
                "note": "native-grid Skala Goldzak12 protocol",
            }
        )
    write_rows(
        args.output / "all-comparison-values.csv",
        [
            "source",
            "method",
            "solid",
            "a_A",
            "B0_GPa",
            "E_coh_eV_atom",
            "status",
            "note",
        ],
        comparison_rows,
    )

    fit_lookup = {(row["method"], row["solid"]): row for row in fit_rows}
    pairs = (
        ("gapwxc-gth_minus_gapw-ae", "gapwxc-gth", "gapw-ae"),
        ("hybrid-direct_minus_gapw-ae", "hybrid-direct", "gapw-ae"),
        ("hybrid-one-center_minus_gapw-ae", "hybrid-one-center", "gapw-ae"),
        (
            "hybrid-one-center_minus_hybrid-direct",
            "hybrid-one-center",
            "hybrid-direct",
        ),
    )
    pair_rows: list[dict[str, object]] = []
    pair_errors: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for label, left, right in pairs:
        for solid in sorted(systems):
            if (left, solid) not in fit_lookup or (right, solid) not in fit_lookup:
                continue
            left_row = fit_lookup[(left, solid)]
            right_row = fit_lookup[(right, solid)]
            delta_a = float(left_row["a0_A"]) - float(right_row["a0_A"])
            delta_bulk = float(left_row["B0_GPa"]) - float(right_row["B0_GPa"])
            delta_cohesive = float(left_row["E_coh_eV_atom"]) - float(
                right_row["E_coh_eV_atom"]
            )
            pair_rows.append(
                {
                    "comparison": label,
                    "solid": solid,
                    "delta_a_A": f"{delta_a:.8f}",
                    "delta_B0_GPa": f"{delta_bulk:.6f}",
                    "delta_E_coh_eV_atom": f"{delta_cohesive:.8f}",
                }
            )
            pair_errors[label]["a"].append(delta_a)
            pair_errors[label]["B0"].append(delta_bulk)
            pair_errors[label]["E_coh"].append(delta_cohesive)
    write_rows(
        args.output / "native-method-differences.csv",
        [
            "comparison",
            "solid",
            "delta_a_A",
            "delta_B0_GPa",
            "delta_E_coh_eV_atom",
        ],
        pair_rows,
    )

    pair_stat_rows: list[dict[str, object]] = []
    for label, properties in sorted(pair_errors.items()):
        record: dict[str, object] = {"comparison": label, "N": len(properties["a"])}
        for property_name, values in properties.items():
            me, mae, rmse, maxae = aggregate(values)
            record[f"ME_{property_name}"] = f"{me:.8f}"
            record[f"MAE_{property_name}"] = f"{mae:.8f}"
            record[f"RMSE_{property_name}"] = f"{rmse:.8f}"
            record[f"MaxAE_{property_name}"] = f"{maxae:.8f}"
        pair_stat_rows.append(record)
    write_rows(
        args.output / "native-method-difference-statistics.csv",
        [
            "comparison",
            "N",
            "ME_a",
            "MAE_a",
            "RMSE_a",
            "MaxAE_a",
            "ME_B0",
            "MAE_B0",
            "RMSE_B0",
            "MaxAE_B0",
            "ME_E_coh",
            "MAE_E_coh",
            "RMSE_E_coh",
            "MaxAE_E_coh",
        ],
        pair_stat_rows,
    )


if __name__ == "__main__":
    main()
