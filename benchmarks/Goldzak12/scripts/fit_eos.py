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
