#!/usr/bin/env python3
"""Run an immutable X23 case, monitor its newest execution and retain provenance."""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import signal
import shutil
import socket
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REVISION = "37aedacbb33677307818301a47a09ae435712293"
DATA_HASHES = {
    "BASIS_MOLOPT_UZH_2026.2": "dd2a387f4f1ed8ac84696bea7472a5ffb63a7dc8f5d948872d137089e78c38ed",
    "POTENTIAL_UZH_2026.2": "8ab49093391f67f2705679f0b3bfc5b84bbeb4003955b556ded322e9d1d6093d",
    "dftd3.dat": "4aa1e246a1661bb875590a67229cbc0747e9ad5e11ec1b8e47bce30db866dd1d",
}
MODELS = {"cpu": "7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd",
          "gpu": "f848eae769dca91741a518ae7275d10caac398ab21db649f91bc1f136872f223"}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def verify(path, expected):
    actual = sha(path)
    if actual != expected:
        raise RuntimeError(f"Hash mismatch: {path}: {actual}")
    return actual


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def dump(path, data):
    temp = path.with_suffix(path.suffix + ".new")
    temp.write_text(json.dumps(data, indent=2) + "\n")
    temp.replace(path)


def memory():
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) / 1048576
    raise RuntimeError("No MemAvailable")


def actual_input(frozen, gpu=False, recovery=None, control=None):
    text = frozen
    if recovery is not None:
        supported = (recovery["name"] == "full-all-step005" and recovery["case_id"] in (
            "gapw-ae/CO2/molecule", "gapw-ae/NH3/molecule", "gapw-ae/urea/molecule"))
        cg = recovery["name"] == "full-all-cg-step005" and recovery["case_id"] == "gapw-ae/urea/molecule"
        diis_from_cg = recovery["name"] == "full-all-diis-from-cg" and recovery["case_id"] == "gapw-ae/urea/molecule"
        crystal_cg = recovery["name"] == "irac-cg-step005" and recovery["case_id"] == "gapw-ae/CO2/solid" and not gpu
        crystal_diis = recovery["name"] == "irac-diis-from-cg" and recovery["case_id"] == "gapw-ae/CO2/solid" and not gpu
        if not (supported or cg or diis_from_cg or crystal_cg or crystal_diis):
            raise ValueError("Unsupported recovery plan")
        seed_name = recovery_seed_name(recovery)
        replacements = (("SCF_GUESS ATOMIC", "SCF_GUESS RESTART"),
                        ("STEPSIZE 0.10", "STEPSIZE 0.05"),
                        ("   &DFT\n", f"   &DFT\n      WFN_RESTART_FILE_NAME {seed_name}\n"))
        if crystal_cg or crystal_diis:
            if text.count("PRECONDITIONER FULL_ALL") != 1 or text.count("ALGORITHM IRAC") != 1:
                raise ValueError("Crystal recovery requires the frozen IRAC/FULL_ALL input")
        else:
            replacements += (("PRECONDITIONER FULL_KINETIC", "PRECONDITIONER FULL_ALL"),)
        for old, new in replacements:
            if text.count(old) != 1:
                raise ValueError(f"Unexpected input for recovery: {old}")
            text = text.replace(old, new)
        if cg or crystal_cg:
            if text.count("MINIMIZER DIIS") != 1:
                raise ValueError("Unexpected minimizer for CG recovery")
            text = text.replace("MINIMIZER DIIS", "MINIMIZER CG")
    if gpu:
        text = text.replace("NATIVE_GRID_USE_CUDA .FALSE.", "NATIVE_GRID_USE_CUDA .TRUE.")
    if control is not None:
        cpu = (control["name"] == "cpu-model" and control["case_id"] in (
            "gapwxc-gth/CO2/solid", "gapw-ae/urea/solid") and not gpu)
        symmetry = control["name"] == "sym-tol-1e4" and control["case_id"] == "gapwxc-gth/NH3/solid" and gpu
        kmesh = control["name"] == "kmesh-5" and control["case_id"] == "gapwxc-gth/CO2/solid" and gpu
        quadrature = (control["name"] == "atom-200-974" and control["case_id"] in (
            "gapwxc-gth/CO2/solid", "gapwxc-gth/CO2/molecule")
            and gpu == control["case_id"].endswith("/solid"))
        cutoff = (control["name"] in ("cutoff-1200", "cutoff-1200-mem200") and control["case_id"] in (
            "gapwxc-gth/CO2/solid", "gapwxc-gth/CO2/molecule")
            and gpu == control["case_id"].endswith("/solid"))
        if control["name"] == "cutoff-1200-mem200" and control["case_id"] != "gapwxc-gth/CO2/molecule":
            raise ValueError("Memory retry is only authorized for the CO2 molecule")
        if recovery or not (cpu or quadrature or cutoff or symmetry or kmesh):
            raise ValueError("Unsupported numerical control")
        seed_name = control_seed_name(control)
        for old, new in (("SCF_GUESS ATOMIC", "SCF_GUESS RESTART"),
                         ("   &DFT\n", f"   &DFT\n      WFN_RESTART_FILE_NAME {seed_name}\n")):
            if text.count(old) != 1:
                raise ValueError(f"Unexpected control input: {old}")
            text = text.replace(old, new)
        if symmetry:
            old = "         SYMMETRY .TRUE.\n"
            if text.count(old) != 1 or "EPS_SYMMETRY" in text:
                raise ValueError("Unexpected NH3 symmetry input")
            text = text.replace(old, old + "         EPS_SYMMETRY 1.0E-4\n")
        if kmesh:
            old = "         SCHEME MONKHORST-PACK 3 3 3\n"
            if text.count(old) != 1:
                raise ValueError("Unexpected CO2 k-point input")
            text = text.replace(old, "         SCHEME MONKHORST-PACK 5 5 5\n")
        if quadrature:
            for old, new in (("RADIAL_GRID 150", "RADIAL_GRID 200"),
                             ("LEBEDEV_GRID 770", "LEBEDEV_GRID 974")):
                if text.count(old) != 2:
                    raise ValueError(f"Unexpected CO2 quadrature input: {old}")
                text = text.replace(old, new)
        if cutoff:
            old = "         CUTOFF 800\n"
            if text.count(old) != 1:
                raise ValueError("Unexpected CO2 cutoff input")
            text = text.replace(old, "         CUTOFF 1200\n")
    return text


def control_seed_name(plan):
    return "seed.kp" if plan["case_id"].endswith("/solid") else "seed.wfn"


def recovery_seed_name(plan):
    return "seed.kp" if plan["case_id"].endswith("/solid") else "seed.wfn"


def recovery_source(plan, case, runtime):
    if plan["case_id"] != case["id"]:
        raise ValueError("Recovery case mismatch")
    parent = (ROOT / plan["source_directory"]).resolve()
    expected_parent = ROOT / "runs" / case["id"]
    if plan["name"] == "full-all-cg-step005":
        expected_parent = expected_parent / "recoveries/full-all-step005"
    if plan["name"] == "full-all-diis-from-cg":
        expected_parent = expected_parent / "recoveries/full-all-cg-step005"
    if plan["name"] == "irac-diis-from-cg":
        expected_parent = expected_parent / "recoveries/irac-cg-step005"
    if parent != expected_parent.resolve():
        raise ValueError("Recovery source is not the designated parent")
    verify(parent / "execution.json", plan["source_execution_sha256"])
    record = json.loads((parent / "execution.json").read_text())
    if record["case"] != case or not record.get("ended") or record.get("execution_exit") in (None, 0):
        raise ValueError("Recovery requires a finished, failed matching case")
    if record["runtime"]["binary_sha256"] != runtime["binary_sha256"]:
        raise ValueError("Recovery would mix binaries")
    if plan["name"] == "irac-cg-step005" and not (
            case["id"] == "gapw-ae/CO2/solid" and record.get("stop_reason") == "stationary_floor"
            and record.get("execution_exit") == -15 and record["markers"]["emitted_steps"] == 25):
        raise ValueError("Crystal CG recovery requires the preserved CO2 stationary stop")
    if plan["name"] in ("full-all-diis-from-cg", "irac-diis-from-cg"):
        verify(parent / "output.out", record["final_hashes"]["output.out"])
        output = (parent / "output.out").read_text()
        start = output.rfind("PROGRAM STARTED AT")
        if not (record.get("execution_exit") == 1 and start >= 0
                and trajectory(output)["converged"]
                and "mp_world_finalize: assert failed: leaking communicators" in output[start:]):
            raise ValueError("DIIS cleanup requires the preserved converged CG/finalization failure")
    name = plan["source_wfn_name"]
    if Path(name).name != name or record["final_hashes"].get(name) != plan["source_wfn_sha256"]:
        raise ValueError("WFN not preserved in source execution")
    verify(parent / name, plan["source_wfn_sha256"])
    return parent / name


def require_no_case_process(case_root):
    for proc in Path("/proc").glob("[0-9]*"):
        try:
            if proc.stat().st_uid != os.getuid() or "cp2k" not in (proc / "comm").read_text().lower():
                continue
            state = (proc / "stat").read_text().split(") ")[1][0]
            if state != "Z" and (proc / "cwd").resolve().is_relative_to(case_root.resolve()):
                raise RuntimeError(f"Case still has a live CP2K process: {proc.name}")
        except FileNotFoundError:
            continue


def live_launch_snapshot(proc, directory, group_id, runtime):
    pid = int(proc.name)
    # Open MPI creates child process groups inside the launcher's owned session.
    if (proc.stat().st_uid != os.getuid() or os.getsid(pid) != group_id
            or "cp2k" not in (proc / "comm").read_text().lower()
            or (proc / "stat").read_text().split(") ")[1][0] == "Z"
            or (proc / "cwd").resolve() != directory.resolve()
            or (proc / "exe").resolve() != Path(runtime["binary"]).resolve()):
        return None
    libraries = sorted({line.split()[-1] for line in (proc / "maps").read_text().splitlines()
                        if "libcp2k.so" in line})
    if ({str(Path(p).resolve()) for p in libraries}
            != {str(Path(p).resolve()) for p in runtime["cp2k_shared_libraries"]}):
        return None
    wanted = {"CUDA_VISIBLE_DEVICES", "GAUXC_SKALA_MODEL", "GAUXC_SKALA_CUDA_MODEL",
              "CP2K_DATA_DIR", "OMP_NUM_THREADS", "OMP_PROC_BIND", "OPENBLAS_NUM_THREADS",
              "MKL_NUM_THREADS", "LD_LIBRARY_PATH", "LD_PRELOAD"}
    environment = {k: v for item in (proc / "environ").read_text().split("\0") if "=" in item
                   for k, v in [item.split("=", 1)] if k in wanted}
    return {"observed": now(), "host": socket.gethostname(), "pid": pid,
            "process_group": os.getpgid(pid), "launcher_session": group_id,
            "cwd": str(directory.resolve()),
            "executable": str((proc / "exe").resolve()),
            "affinity": sorted(os.sched_getaffinity(pid)),
            "resources": [line for line in (proc / "status").read_text().splitlines()
                          if line.startswith(("VmRSS:", "Threads:"))],
            "loaded_cp2k_libraries": libraries, "runtime_environment": environment}


def capture_launch_provenance(directory, group_id, runtime):
    destination = directory / "launch-environment.json"
    if destination.exists():
        return True
    for proc in Path("/proc").glob("[0-9]*"):
        try:
            snapshot = live_launch_snapshot(proc, directory, group_id, runtime)
        except (OSError, ProcessLookupError):
            continue
        if snapshot is not None:
            dump(destination, snapshot)
            return True
    return False


def control_source(plan, case, runtime, validation_path=None):
    if plan["case_id"] != case["id"] or plan["source_directory"] != "runs/" + case["id"]:
        raise ValueError("Control source case mismatch")
    parent = ROOT / plan["source_directory"]
    verify(parent / "execution.json", plan["source_execution_sha256"])
    validation_path = validation_path or parent / "validation.json"
    verify(validation_path, plan["source_validation_sha256"])
    accepted = json.loads(validation_path.read_text())
    record = json.loads((parent / "execution.json").read_text())
    if (accepted["status"] != "accepted_scf_numerical_accuracy_unassessed"
            or accepted["execution_sha256"] != plan["source_execution_sha256"]
            or record["case"] != case or record.get("execution_exit") != 0
            or record["runtime"]["binary_sha256"] != runtime["binary_sha256"]):
        raise ValueError("Control needs an accepted parent with the same binary")
    name = plan["source_wfn_name"]
    if Path(name).name != name or accepted["final_hashes"].get(name) != plan["source_wfn_sha256"]:
        raise ValueError("Control WFN does not match accepted parent")
    verify(parent / name, plan["source_wfn_sha256"])
    if plan["name"] == "cutoff-1200-mem200":
        failed = parent / "controls/cutoff-1200"
        verify(failed / "execution.json", plan["failed_control_execution_sha256"])
        previous = json.loads((failed / "execution.json").read_text())
        if (previous["execution_exit"] != 137 or previous["markers"]["emitted_steps"] != 0
                or previous["case"] != case
                or previous["runtime"]["binary_sha256"] != runtime["binary_sha256"]
                or previous["control_plan"]["source_wfn_sha256"] != plan["source_wfn_sha256"]):
            raise ValueError("Memory retry requires the recorded pre-SCF memory failure")
        verify(failed / "actual-input.inp", previous["actual_input_sha256"])
        expected = actual_input((parent / "frozen-input.inp").read_text(), control=plan)
        if expected != (failed / "actual-input.inp").read_text():
            raise ValueError("Memory retry must preserve the failed control input exactly")
    return parent / name


def trajectory(text):
    start = text.rfind("PROGRAM STARTED AT")
    block = text[start:] if start >= 0 else ""
    rows, outer = [], 1
    for line in block.splitlines():
        fields = line.split()
        if len(fields) >= 8 and fields[0].isdigit() and fields[1] == "OT" and fields[2] != "INIT":
            try:
                row = {"step": int(fields[0]), "residual": float(fields[-3]), "energy": float(fields[-2])}
            except ValueError:
                continue
            if rows and row["step"] <= rows[-1]["step"]:
                outer += 1
            row["outer"] = outer
            rows.append(row)
    finite = all(math.isfinite(r["residual"]) and math.isfinite(r["energy"]) for r in rows)
    best = min(rows, key=lambda r: r["residual"]) if rows and finite else None
    old = min((r["residual"] for r in rows[:-10]), default=None)
    factor = old / best["residual"] if old is not None and best and best["residual"] > 0 else None
    conv = "SCF run converged" in block and "SCF run NOT converged" not in block
    if not finite:
        classification = "divergent_nonfinite"
    elif conv and "[ABORT]" in block:
        classification = "scf_converged_but_aborted"
    elif conv:
        classification = "converged_pending_validation"
    elif len(rows) < 20:
        classification = "initializing" if not rows else "early_trajectory"
    elif factor is not None and factor < 1.1:
        tail = [r["residual"] for r in rows[-10:]]
        classification = "stationary_floor_oscillatory" if max(tail) > 2 * min(tail) else "stationary_floor"
    else:
        classification = "improving"
    energies = re.findall(r"ENERGY\| Total FORCE_EVAL.*?\[hartree\]\s+([-+\d.Ee]+)", block)
    electrons = re.findall(r"Number of electrons:\s+(\d+)", block)
    return {"emitted_steps": len(rows), "best": best, "last12": rows[-12:], "ten_step_factor": factor,
            "classification": classification, "converged": conv, "ended": "PROGRAM ENDED AT" in block,
            "abort": "[ABORT]" in block, "energy_hartree": float(energies[-1]) if energies else None,
            "electrons": [int(n) for n in electrons],
            "outer_messages": [l.strip() for l in block.splitlines() if "outer SCF" in l][-8:],
            "native_electron_diagnostics": [l.strip() for l in block.splitlines() if "electron" in l.lower() and ("atom" in l.lower() or "skala" in l.lower())][-8:]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case")
    p.add_argument("--manifest", default="manifest.json")
    p.add_argument("--threads", type=int, default=12)
    p.add_argument("--cpus", help="Comma separated, inclusive ranges allowed")
    p.add_argument("--gpu", action="store_true")
    p.add_argument("--check-only", action="store_true")
    plans = p.add_mutually_exclusive_group()
    plans.add_argument("--recovery-plan", help="Immutable in-campaign recovery plan")
    plans.add_argument("--control-plan", help="Immutable numerical comparison plan")
    args = p.parse_args()
    manifest_path = (ROOT / args.manifest).resolve()
    if not manifest_path.is_relative_to(ROOT.resolve()):
        raise RuntimeError("Manifest must be inside the campaign")
    manifest = json.loads(manifest_path.read_text())
    if "parent_manifest" in manifest:
        verify(ROOT / manifest["parent_manifest"], manifest["parent_manifest_sha256"])
    case = next(c for c in manifest["cases"] if c["id"] == args.case)
    runtime = json.loads((ROOT / "runtime/runtime-build.json").read_text())
    if runtime["revision"] != REVISION or runtime["status"] != "built_not_scientifically_validated":
        raise RuntimeError("Pinned runtime build has not succeeded")
    binary = runtime["binary"]
    verify(binary, runtime["binary_sha256"])
    for library, expected in runtime["cp2k_shared_libraries"].items():
        verify(library, expected)
    verify(Path(runtime["source"]) / "src/qs_vxc_atom.F", runtime["source_sha256"])
    source_input = ROOT / case["input"]
    verify(source_input, case["input_sha256"])
    env = os.environ.copy()
    data = Path(env["CP2K_DATA_DIR"])
    for name, expected in DATA_HASHES.items():
        verify(data / name, expected)
    mode = "gpu" if args.gpu else "cpu"
    model = env["GAUXC_SKALA_CUDA_MODEL" if args.gpu else "GAUXC_SKALA_MODEL"]
    verify(model, MODELS[mode])
    # Verify both models when using CUDA; the adapter may load the host model too.
    verify(env["GAUXC_SKALA_MODEL"], MODELS["cpu"])
    if args.cpus:
        cpus = set()
        for field in args.cpus.split(","):
            span = list(map(int, field.split("-")))
            cpus.update(range(span[0], span[-1] + 1))
        if not cpus <= os.sched_getaffinity(0):
            raise RuntimeError("Requested CPUs outside allowed affinity")
        os.sched_setaffinity(0, cpus)
    if args.threads < 1 or args.threads > len(os.sched_getaffinity(0)):
        raise RuntimeError("Invalid thread/affinity allocation")
    if memory() < (45 if args.gpu else 60):
        raise RuntimeError("Insufficient available memory for a new case")
    env.update(OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC="FALSE", OMP_PROC_BIND="FALSE",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_STACKSIZE="256M")
    env.pop("OMP_PLACES", None)
    directory = ROOT / ("preflight" if args.check_only else "runs") / args.case
    recovery, seed, control = None, None, None
    if args.recovery_plan:
        plan_path = (ROOT / args.recovery_plan).resolve()
        if not plan_path.is_relative_to(ROOT.resolve()):
            raise RuntimeError("Recovery plan outside campaign")
        recovery = json.loads(plan_path.read_text())
        seed = recovery_source(recovery, case, runtime)
        require_no_case_process(ROOT / "runs" / args.case)
        actual_input(source_input.read_text(), args.gpu, recovery)
        directory = directory / "recoveries" / recovery["name"]
    if args.control_plan:
        plan_path = (ROOT / args.control_plan).resolve()
        if not plan_path.is_relative_to(ROOT.resolve()):
            raise RuntimeError("Control plan outside campaign")
        control = json.loads(plan_path.read_text())
        if control["name"] == "cutoff-1200-mem200" and int(env.get("SLURM_MEM_PER_NODE", "0")) < 200 * 1024:
            raise RuntimeError("Memory retry requires a Slurm allocation of at least 200 GiB")
        actual_input(source_input.read_text(), args.gpu, control=control)
        seed = control_source(control, case, runtime)
        require_no_case_process(ROOT / "runs" / args.case)
        directory = directory / "controls" / control["name"]
    directory.mkdir(parents=True, exist_ok=False)
    actual = actual_input(source_input.read_text(), args.gpu, recovery, control)
    if recovery:
        seed_name = recovery_seed_name(recovery)
        shutil.copyfile(seed, directory / seed_name)
        verify(directory / seed_name, recovery["source_wfn_sha256"])
        (directory / "recovery-plan.json").write_bytes(plan_path.read_bytes())
    if control:
        seed_name = control_seed_name(control)
        shutil.copyfile(seed, directory / seed_name)
        verify(directory / seed_name, control["source_wfn_sha256"])
        (directory / "control-plan.json").write_bytes(plan_path.read_bytes())
        shutil.copyfile(ROOT / control["source_directory"] / "validation.json", directory / "parent-validation.json")
        verify(directory / "parent-validation.json", control["source_validation_sha256"])
    (directory / "frozen-input.inp").write_bytes(source_input.read_bytes())
    (directory / "actual-input.inp").write_text(actual)
    (directory / "manifest-source.json").write_bytes(manifest_path.read_bytes())
    (directory / "runner-source.py").write_bytes(Path(__file__).read_bytes())
    command = shlex.split(env.get("CP2K_LAUNCHER", "")) + [binary]
    record = {"case": case, "host": socket.gethostname(), "started": now(), "runner_pid": os.getpid(),
              "runtime": runtime, "manifest": str(manifest_path), "manifest_sha256": sha(manifest_path), "runner_sha256": sha(__file__),
              "actual_input_sha256": sha(directory / "actual-input.inp"), "data_sha256": DATA_HASHES,
              "model_sha256": MODELS[mode], "mode": mode, "threads": args.threads,
              "affinity": sorted(os.sched_getaffinity(0)), "cwd": str(directory), "available_gib_start": memory(),
              "slurm_job": env.get("SLURM_JOB_ID"), "changes_from_frozen": ["CUDA model evaluation"] if args.gpu else []}
    if recovery:
        record.update(attempt=recovery["name"], recovery_plan=recovery,
                      recovery_plan_sha256=sha(plan_path))
        record["changes_from_frozen"] += ["Exact saved WFN restart", "STEPSIZE 0.05"]
        if recovery["name"] not in ("irac-cg-step005", "irac-diis-from-cg"):
            record["changes_from_frozen"] += ["FULL_ALL preconditioner"]
        if recovery["name"] in ("full-all-cg-step005", "irac-cg-step005"):
            record["changes_from_frozen"] += ["MINIMIZER CG"]
    if control:
        record.update(control_name=control["name"], control_plan=control, control_plan_sha256=sha(plan_path))
        record["changes_from_frozen"] += ["Exact accepted WFN restart"]
        if control["name"] == "atom-200-974":
            record["changes_from_frozen"] += ["RADIAL_GRID 200", "LEBEDEV_GRID 974"]
        if control["name"] in ("cutoff-1200", "cutoff-1200-mem200"):
            record["changes_from_frozen"] += ["CUTOFF 1200 Ry"]
        if control["name"] == "sym-tol-1e4":
            record["changes_from_frozen"] += ["EPS_SYMMETRY 1.0E-4"]
        if control["name"] == "kmesh-5":
            record["changes_from_frozen"] += ["Gamma-centered 5x5x5 k-point mesh"]
    dump(directory / "execution.json", record)
    with (directory / "preflight.txt").open("w") as log:
        preflight = subprocess.run(command + ["--check", "-i", "actual-input.inp"], cwd=directory, env=env,
                                   stdout=log, stderr=subprocess.STDOUT, timeout=180)
    record["preflight_exit"] = preflight.returncode
    dump(directory / "execution.json", record)
    if preflight.returncode or args.check_only:
        raise SystemExit(preflight.returncode)
    command = [env.get("TIME_BINARY", "/usr/bin/time"), "-v", "-o", "time.txt"] + command + ["-i", "actual-input.inp", "-o", "output.out"]
    with (directory / "console.txt").open("w") as log:
        child = subprocess.Popen(command, cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        record["process_group"] = child.pid
        dump(directory / "execution.json", record)
        def stop(signum, frame):
            record["stop_reason"] = f"signal {signum}"
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
        signal.signal(signal.SIGTERM, stop)
        signal.signal(signal.SIGINT, stop)
        launch_captured = False
        while child.poll() is None:
            try:
                child.wait(timeout=60 if launch_captured else 1)
            except subprocess.TimeoutExpired:
                if not launch_captured:
                    launch_captured = capture_launch_provenance(directory, child.pid, runtime)
                out = directory / "output.out"
                report = trajectory(out.read_text(errors="replace")) if out.exists() else {}
                report.update(observed=now(), available_gib=memory())
                dump(directory / "trajectory.json", report)
                if memory() < 16 or report.get("classification", "").startswith(("stationary_floor", "divergent")):
                    record["stop_reason"] = "low memory" if memory() < 16 else report["classification"]
                    if child.poll() is None:
                        os.killpg(child.pid, signal.SIGTERM)
                if record.get("stop_reason"):
                    try:
                        child.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        os.killpg(child.pid, signal.SIGKILL)
                        child.wait()
        record["execution_exit"] = child.returncode
    record["ended"] = now()
    out = directory / "output.out"
    report = trajectory(out.read_text(errors="replace")) if out.exists() else {}
    dump(directory / "trajectory.json", report)
    record["markers"] = report
    record["clean_scf_markers"] = child.returncode == 0 and report.get("converged") and report.get("ended") and not report.get("abort")
    record["final_hashes"] = {p.name: sha(p) for p in directory.iterdir() if p.is_file() and p.name != "execution.json"}
    record["status"] = "converged_pending_scientific_validation" if record["clean_scf_markers"] else "stopped_or_failed"
    dump(directory / "execution.json", record)
    raise SystemExit(0 if record["clean_scf_markers"] else 1)


if __name__ == "__main__":
    main()
