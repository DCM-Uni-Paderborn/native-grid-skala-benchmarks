#!/usr/bin/env python3
"""Run the archived extracted CP2K kernels against a compatible GNU/OpenMP build.

The build must expose the orbital_transformation_matrices_unittest Ninja target.
No CP2K rebuild, model inference or SCF calculation is performed.
"""
# Adapted from CP2K's kernel-validation driver; GPL-2.0-or-later, see LICENSE.
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument(
        "--check", action="store_true", help="Bounds/FPE checks, no timing"
    )
    parser.add_argument("--rows", type=int, default=4000)
    args = parser.parse_args()
    source_dir = Path(__file__).resolve().parent
    build, work = args.build_dir.resolve(), args.work_dir.resolve()
    root = source_dir.parents[1]
    if work == root or root in work.parents:
        parser.error('Choose a work directory outside the immutable data repository')
    work.mkdir(parents=True, exist_ok=True)
    evidence = json.loads((source_dir / "evidence.json").read_text())
    sources = [source_dir / (name + ".F90") for name in ("original", "values_only", "screened", "parallel")]
    # Reuse the build's resolved compiler and link dependencies without guessing
    # site-specific library paths. Invoke argv directly, never through a shell.
    commands = subprocess.check_output(
        [
            "ninja",
            "-C",
            str(build),
            "-t",
            "commands",
            "orbital_transformation_matrices_unittest",
        ],
        text=True,
    )
    link = shlex.split(commands.strip().splitlines()[-1])
    if link[:2] != [":", "&&"] or link[-2:] != ["&&", ":"]:
        raise ValueError("Unsupported CMake link command layout")
    link = link[2:-2]
    objects = [x for x in link if x.endswith(".F.o")]
    if len(objects) != 1:
        raise ValueError("Expected one unit-test driver object")
    compiler = link[0]
    executable = work / ("adjoint_checked" if args.check else "adjoint_release")
    flags = [
        "-cpp",
        "-fopenmp",
        "-ffree-line-length-none",
        "-I",
        str(build / "src/mod_files"),
    ]
    flags += (
        ["-O1", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow"]
        if args.check
        else ["-O3", "-funroll-loops"]
    )
    sources.append(source_dir / "adjoint-validation-fixture.F90")
    compiled = []
    for source in sources:
        target = work / (source.stem + ".o")
        subprocess.run(
            [compiler, *flags, "-c", str(source), "-o", str(target)],
            cwd=work,
            check=True,
        )
        compiled.append(str(target))
    pos = link.index(objects[0])
    link[pos : pos + 1] = compiled
    link[link.index("-o") + 1] = str(executable)
    subprocess.run(link, cwd=build, check=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    run_command = [
        str(executable),
        "check" if args.check else "benchmark",
        str(args.rows),
    ]
    if Path("/usr/bin/time").exists():
        run_command = [
            "/usr/bin/time",
            "-l" if sys.platform == "darwin" else "-v",
            *run_command,
        ]
    completed = subprocess.run(
        run_command, cwd=work, env=env, text=True, capture_output=True
    )
    print(completed.stdout, end="")
    print(completed.stderr, end="")
    report = {
        "source_evidence": evidence["reports"]["optimized"]["baseline_revision"],
        "fixture_note": "Retained optimized fixture; checked mode uses this same fixture.",
        "compiler": subprocess.check_output(
            [compiler, "--version"], text=True
        ).splitlines()[0],
        "flags": flags,
        "exit_code": completed.returncode,
        "output": completed.stdout,
        "stderr": completed.stderr,
        "generated_sources_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources
        },
        "scope": "Extracted production kernels and exact caller loop, real CP2K types/harmonics; no SCF/model test",
    }
    (
        work / ("checked-results.json" if args.check else "benchmark-results.json")
    ).write_text(json.dumps(report, indent=2) + "\n")
    completed.check_returncode()


if __name__ == "__main__":
    main()
