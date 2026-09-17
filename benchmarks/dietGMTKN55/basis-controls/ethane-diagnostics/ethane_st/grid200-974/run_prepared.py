#!/usr/bin/env python3
"""Execute an immutable prepared case with the parallel native-grid runtime."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import socket
import subprocess
import time
from molecular_protocol import validate_input

from execution_helpers import (DATA_HASHES, MODELS, REVISION, capture_launch_provenance,
                               dump, memory, now, require_no_case_process, sha, trajectory, verify)

BINARIES = {
    "terok": "4c4cdcf4ac51510d1d876310efa724ae7efa8bdab7fdd124d803e0220fb770fb",
    "rosi": "3bd708ff4769f7fc98a15d9c692c44c4586a9d1d3351fe84e3842d1b6e5f20c9",
}
SOURCE = "11e24230b22a531652f05a5dc222c92903612086614446cc37c247833120e0e6"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifest", type=Path)
    p.add_argument("case")
    p.add_argument("--runtime", type=Path, required=True)
    p.add_argument("--host", choices=BINARIES, required=True)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--cpus")
    p.add_argument("--check-only", action="store_true")
    args = p.parse_args()
    root = args.manifest.resolve().parent
    manifest = json.loads(args.manifest.read_text())
    case = next(c for c in manifest["cases"] if c["id"] == args.case)
    if case.get("requires"):
        assessment = json.loads((root / "paired-grid-assessment.json").read_text())
        if not assessment.get("basis_stage_authorized"):
            raise RuntimeError("Paired grid assessment required before the basis stage")
    runtime = json.loads(args.runtime.read_text())
    if runtime["revision"] != REVISION or runtime["binary_sha256"] != BINARIES[args.host]:
        raise RuntimeError("Unexpected runtime")
    verify(runtime["binary"], BINARIES[args.host])
    verify(Path(runtime["source"]) / "src/qs_vxc_atom.F", SOURCE)
    for path, digest in runtime["cp2k_shared_libraries"].items():
        verify(path, digest)
    inp = root / case["input"]
    verify(inp, case["input_sha256"])
    validate_input(inp.read_text(), case)
    env = os.environ.copy()
    for name, digest in DATA_HASHES.items():
        verify(Path(env["CP2K_DATA_DIR"]) / name, digest)
    verify(env["GAUXC_SKALA_MODEL"], MODELS["cpu"])
    if args.cpus:
        cpus = set()
        for span in args.cpus.split(","):
            ends = list(map(int, span.split("-")))
            cpus.update(range(ends[0], ends[-1] + 1))
        if not cpus <= os.sched_getaffinity(0):
            raise RuntimeError("Requested CPUs outside allocation")
        os.sched_setaffinity(0, cpus)
    affinity = sorted(os.sched_getaffinity(0))
    if len(affinity) != args.threads:
        raise RuntimeError(f"Affinity must contain exactly {args.threads} CPUs: {affinity}")
    if args.host == "rosi" and (args.threads != 32 or env.get("SLURM_CPUS_PER_TASK") != "32"):
        raise RuntimeError("ROSI timing series requires exactly 32 allocated cores")
    if memory() < 80:
        raise RuntimeError("Insufficient free memory")
    env.update(OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC="FALSE", OMP_PROC_BIND="FALSE",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_STACKSIZE="256M",
               CUDA_VISIBLE_DEVICES="-1", PYTHONDONTWRITEBYTECODE="1")
    env.pop("OMP_PLACES", None)
    directory = root / ("preflight" if args.check_only else "runs") / args.case
    require_no_case_process(directory)
    directory.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(inp, directory / "actual-input.inp")
    shutil.copyfile(args.manifest, directory / "manifest-source.json")
    for name in ("run_prepared.py", "execution_helpers.py", "molecular_protocol.py", "baseline_protocol.py"):
        shutil.copyfile(Path(__file__).parent / name, directory / name)
    if case.get("seed"):
        parent = case["parent"]
        origin = Path(parent["remote_directory"])
        for name, digest in (("execution.json", parent["execution_sha256"]),
                             ("actual-input.inp", parent["actual_input_sha256"]),
                             ("output.out", parent["output_sha256"])):
            verify(origin / name, digest)
        previous = trajectory((origin / "output.out").read_text())
        old_record = json.loads((origin / "execution.json").read_text())
        if old_record.get("execution_exit") != 0 or not previous["converged"] or not previous["ended"] or previous["abort"]:
            raise RuntimeError("Seed parent is not process-clean")
        seed = case["seed"]
        verify(seed["source"], seed["sha256"])
        shutil.copyfile(seed["source"], directory / seed["name"])
        verify(directory / seed["name"], seed["sha256"])
    record = {"case": case, "started": now(), "host": socket.gethostname(), "runtime": runtime,
              "runtime_record_sha256": sha(args.runtime), "runner_sha256": sha(__file__),
              "helper_sha256": sha(Path(__file__).with_name("execution_helpers.py")),
              "manifest_sha256": sha(args.manifest), "actual_input_sha256": sha(inp),
              "model_sha256": MODELS["cpu"], "data_sha256": DATA_HASHES,
              "threads": args.threads, "mpi_ranks": 1, "affinity": affinity,
              "cwd": str(directory), "available_gib_start": memory(),
              "scheduler": {k: v for k, v in env.items() if k.startswith("SLURM_")},
              "lscpu": subprocess.check_output(["lscpu"], text=True), "status": "preflight"}
    dump(directory / "execution.json", record)
    with (directory / "preflight.txt").open("w") as log:
        check = subprocess.run([runtime["binary"], "--check", "-i", "actual-input.inp"],
                               cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
    record["preflight_exit"] = check.returncode
    dump(directory / "execution.json", record)
    if check.returncode or args.check_only:
        raise SystemExit(check.returncode)
    cmd = [env.get("TIME_BINARY", "/usr/bin/time"), "-v", "-o", "time.txt",
           runtime["binary"], "-i", "actual-input.inp", "-o", "output.out"]
    begin = time.monotonic()
    with (directory / "console.txt").open("w") as log:
        child = subprocess.Popen(cmd, cwd=directory, env=env, stdout=log,
                                 stderr=subprocess.STDOUT, start_new_session=True)
        record.update(process_group=child.pid, status="running", execution_started=now())
        dump(directory / "execution.json", record)
        def stop(reason):
            record.setdefault("stop_reason", reason)
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
        signal.signal(signal.SIGTERM, lambda *_: stop("external SIGTERM"))
        signal.signal(signal.SIGINT, lambda *_: stop("external SIGINT"))
        captured = False
        while child.poll() is None:
            try:
                child.wait(timeout=60 if captured else 1)
            except subprocess.TimeoutExpired:
                captured = captured or capture_launch_provenance(directory, child.pid, runtime)
                output = directory / "output.out"
                report = trajectory(output.read_text(errors="replace")) if output.exists() else {}
                report.update(observed=now(), available_gib=memory())
                dump(directory / "trajectory.json", report)
                classification = report.get("classification", "")
                if classification == "divergent_nonfinite":
                    stop(classification)
                if memory() < 16:
                    stop("low memory")
                if not captured and time.monotonic() - begin > 180:
                    stop("runtime process/library provenance not verified")
                if record.get("stop_reason"):
                    try:
                        child.wait(timeout=60)
                    except subprocess.TimeoutExpired:
                        os.killpg(child.pid, signal.SIGKILL)
                        child.wait()
    output = directory / "output.out"
    report = trajectory(output.read_text(errors="replace")) if output.exists() else {}
    timing = (directory / "time.txt").read_text()
    time_exit = re.findall(r"Exit status:\s+(\d+)", timing)
    rss = re.findall(r"Maximum resident set size \(kbytes\):\s+(\d+)", timing)
    record.update(ended=now(), elapsed_wall_s=time.monotonic() - begin, execution_exit=child.returncode,
                  time_exit=int(time_exit[-1]) if time_exit else None,
                  peak_rss_kib=int(rss[-1]) if rss else None, markers=report)
    verify(directory / "actual-input.inp", case["input_sha256"])
    energy = report.get("energy_hartree")
    clean = (child.returncode == 0 and record["time_exit"] == 0 and report.get("converged")
             and report.get("ended") and not report.get("abort") and captured
             and report.get("electrons") == case["expected_spin_electrons"]
             and energy is not None and math.isfinite(energy) and energy < 0 and bool(rss))
    require_no_case_process(directory)
    record["status"] = "converged_pending_scientific_validation" if clean else "stopped_or_failed"
    record["clean_scf_markers"] = bool(clean)
    dump(directory / "trajectory.json", report)
    record["final_hashes"] = {f.name: sha(f) for f in directory.iterdir() if f.is_file() and f.name != "execution.json"}
    dump(directory / "execution.json", record)
    raise SystemExit(0 if clean else 1)


if __name__ == "__main__":
    main()
