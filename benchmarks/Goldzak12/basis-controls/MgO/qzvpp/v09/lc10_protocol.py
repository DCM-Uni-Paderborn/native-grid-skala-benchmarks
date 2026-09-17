"""Strict protocol and newest-block checks for the LC10 paired basis study."""
import hashlib
import math
import os
from pathlib import Path
import re

from execution_helpers import dump, live_launch_snapshot, now


def values(text, key):
    return re.findall(r"(?m)^\s*" + re.escape(key) + r"\s+([^\n]+)", text)


def structural_text(text):
    text = re.sub(r"(?m)^(\s*PROJECT_NAME)\s+\S+", r"\1 PROJECT", text)
    text = re.sub(r"(?m)^\s*(SYMMETRY_BACKEND|SYMMETRY_REDUCTION_METHOD|INVERSION_SYMMETRY_ONLY)\s+[^\n]+\n", "", text)
    text = re.sub(r"BASIS_SET (?:TZVPP|QZVPP)-MOLOPT-PBE-ae", "BASIS_SET AE_BASIS", text)
    text = re.sub(r"(?m)^\s*&OUTER_SCF\s*\n.*?^\s*&END OUTER_SCF\s*\n", "", text, flags=re.S)
    text = re.sub(r"(?m)^(\s*MAX_SCF)\s+\d+", r"\1 MAXIMUM", text)
    return text


def validate_input(text, case):
    def expect(key, wanted):
        actual = [v.strip() for v in values(text, key)]
        if actual != wanted:
            raise ValueError((key, actual, wanted))
    expect("BASIS_SET", [case["basis"]] * case["kinds"])
    expect("POTENTIAL", ["ALL"] * case["kinds"])
    expect("RADIAL_GRID", ["150"] * case["kinds"])
    expect("LEBEDEV_GRID", ["770"] * case["kinds"])
    for key, value in {"CHARGE": "0", "MULTIPLICITY": "1", "UKS": ".FALSE.",
                       "EPS_SCF": "5e-6", "MAX_SCF": "600", "SCF_GUESS": "ATOMIC",
                       "CUTOFF": "800", "REL_CUTOFF": "60", "NGRIDS": "5",
                       "SYMMETRY": ".TRUE.", "FULL_GRID": ".FALSE.",
                       "INVERSION_SYMMETRY_ONLY": ".FALSE.", "EPS_SYMMETRY": "1.0E-8",
                       "SYMMETRY_BACKEND": case["backend"],
                       "SYMMETRY_REDUCTION_METHOD": case["backend"],
                       "SCHEME": case["kmesh"], "POISSON_SOLVER": "PERIODIC",
                       "NATIVE_GRID_USE_CUDA": ".FALSE.", "NATIVE_GRID_LAYOUT": "ATOM_COMPOSITE",
                       "NATIVE_GRID_ATOM_PARTITION": "SMOOTH", "GAPW_ACCURATE_XCINT": ".TRUE.",
                       "REFERENCE_FUNCTIONAL": "B3LYP", "ALPHA": "0.10",
                       "ALGORITHM": "STANDARD"}.items():
        expect(key, [value])
    expect("PERIODIC", ["XYZ", "XYZ"])
    if "&OUTER_SCF" in text or "&OT" in text or "WFN_RESTART_FILE_NAME" in text:
        raise ValueError("Unexpected solver or restart")
    if hashlib.sha256(structural_text(text).encode()).hexdigest() != case["structural_sha256"]:
        raise ValueError("Geometry or frozen numerical settings changed")
    return True


def trajectory(text):
    start = text.rfind("PROGRAM STARTED AT")
    block = text[start:] if start >= 0 else ""
    rows = []
    for line in block.splitlines():
        fields = line.split()
        if len(fields) < 7 or not fields[0].isdigit() or not re.match(r"(?:P_Mix|Broy|Pulay|Diag|OT)", fields[1]):
            continue
        try:
            rows.append(dict(step=int(fields[0]), residual=float(fields[-3]), energy=float(fields[-2])))
        except ValueError:
            continue
    finite = all(math.isfinite(r["residual"]) and math.isfinite(r["energy"]) for r in rows)
    best = min(rows, key=lambda r: r["residual"]) if rows and finite else None
    old = min((r["residual"] for r in rows[:-10]), default=None)
    factor = old / best["residual"] if old is not None and best and best["residual"] > 0 else None
    conv = "SCF run converged" in block and "SCF run NOT converged" not in block
    classification = "improving" if factor is not None and factor >= 1.1 else "early_trajectory"
    if factor is not None and factor < 1.1:
        tail = [r["residual"] for r in rows[-10:]]
        classification = "stationary_floor_oscillatory" if max(tail) > 2 * min(tail) else "stationary_floor"
    if not rows:
        classification = "initializing"
    if conv:
        classification = "converged_pending_validation"
    if not finite:
        classification = "divergent_nonfinite"
    energies = re.findall(r"ENERGY\| Total FORCE_EVAL.*?\[hartree\]\s+([-+\d.Ee]+)", block)
    nk = re.findall(r"BRILLOUIN\| List of Kpoints.*?\s+(\d+)\s*$", block, re.M)
    backend = re.findall(r"BRILLOUIN\| Symmetry backend\s+(\S+)", block)
    reduction = re.findall(r"BRILLOUIN\| Symmetry reduction method\s+(\S+)", block)
    fallback = any("inversion" in line.lower() and any(w in line.lower() for w in ("fallback", "falling back", "only")) for line in block.splitlines())
    return dict(emitted_steps=len(rows), best=best, last12=rows[-12:], ten_step_factor=factor,
                classification=classification, converged=conv, ended="PROGRAM ENDED AT" in block,
                abort="[ABORT]" in block, energy_hartree=float(energies[-1]) if energies else None,
                electrons=[int(n) for n in re.findall(r"Number of electrons:\s+(\d+)", block)],
                irreducible_kpoints=int(nk[-1]) if nk else None, backend=backend[-1] if backend else None,
                reduction=reduction[-1] if reduction else None, inversion_fallback=fallback,
                outer_messages=[l.strip() for l in block.splitlines() if "outer SCF" in l][-8:])


def capture_stable(directory, group, runtime, affinity):
    destination = directory / "stable-launch.json"
    snapshots = []
    if destination.exists():
        snapshots = __import__("json").loads(destination.read_text())["snapshots"]
        if len(snapshots) >= 2:
            return True
    for proc in Path("/proc").glob("[0-9]*"):
        try:
            snap = live_launch_snapshot(proc, directory, group, runtime)
            if snap is None:
                continue
            threads = {p.name: sorted(os.sched_getaffinity(int(p.name))) for p in (proc / "task").iterdir()}
            if snap["affinity"] != affinity or not all(set(c) <= set(affinity) for c in threads.values()):
                continue
            snap.update(thread_affinities=threads, start_ticks=(proc / "stat").read_text().split(") ", 1)[1].split()[19])
            if snapshots and snap["start_ticks"] != snapshots[0]["start_ticks"]:
                raise RuntimeError("PID identity changed")
            snapshots.append(snap)
            dump(destination, dict(snapshots=snapshots, observed=now()))
            return len(snapshots) >= 2
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return False
