#!/usr/bin/env bash
set -euo pipefail

root=${1:?Usage: run_b200_pool.sh RUN_DIRECTORY}
source "${SKALA_ENV:-/home/kuehne88/cp2k-native-skala/env.sh}"

cp2k_root=${CP2K_NATIVE_ROOT:-/home/kuehne88/cp2k-native-skala}
cp2k_bin=${CP2K_BIN:-${cp2k_root}/install/cp2k-b200/bin/cp2k.psmp}
data_dir=${CP2K_DATA_DIR:-${cp2k_root}/install/cp2k-b200/share/cp2k/data}
cpu_model=${GAUXC_SKALA_MODEL:-${cp2k_root}/models/skala-1.1-rev1.fun}
cuda_model=${GAUXC_SKALA_CUDA_MODEL:-${cp2k_root}/models/skala-1.1-rev1-cuda.fun}
max_concurrent=${MAX_B200_JOBS:-3}
min_headroom_bytes=${MIN_JOB_MEMORY_HEADROOM_BYTES:-130000000000}

mapfile -t allowed_cpus < <(python3 - <<'PY'
from pathlib import Path

entry = next(
    line for line in Path("/proc/self/status").read_text().splitlines()
    if line.startswith("Cpus_allowed_list:")
).split(":", 1)[1].strip()
for part in entry.split(","):
    bounds = [int(value) for value in part.split("-")]
    start, stop = (bounds[0], bounds[-1])
    for cpu in range(start, stop + 1):
        print(cpu)
PY
)
threads_per_case=$(( ${#allowed_cpus[@]} / 4 ))
if ((threads_per_case < 1)); then
  echo "Fewer than four CPUs are available" >&2
  exit 1
fi
declare -a cpu_sets
for slot in 0 1 2 3; do
  offset=$((slot * threads_per_case))
  cpu_sets[$slot]=$(IFS=,; echo "${allowed_cpus[*]:offset:threads_per_case}")
done

mapfile -t gpu_uuids < <(
  nvidia-smi --query-gpu=uuid --format=csv,noheader,nounits
)
if ((${#gpu_uuids[@]} != 4)); then
  echo "Expected four visible GPUs, found ${#gpu_uuids[@]}" >&2
  exit 1
fi

run_is_active() {
  local run_dir=$1
  local pid
  for pid in $(pgrep -u "${USER}" -f "${cp2k_bin}" || true); do
    if [[ "$(readlink -f "/proc/${pid}/cwd" 2>/dev/null || true)" == \
      "$(readlink -f "${run_dir}")" ]]; then
      return 0
    fi
  done
  return 1
}

mapfile -t queue < <(
  find "${root}" -type f -name input.inp -printf '%h\n' | sort -V | while read -r run_dir; do
    if { [[ ! -f "${run_dir}/output.out" ]] || \
      ! grep -q "SCF run converged" "${run_dir}/output.out"; } && \
      ! run_is_active "${run_dir}"; then
      printf '%s\n' "${run_dir}"
    fi
  done
)

active_compute_count() {
  nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits 2>/dev/null \
    | sed '/^$/d' | wc -l
}

cgroup_memory_headroom() {
  local path max current
  path=$(cut -d: -f3 /proc/self/cgroup)
  while [[ "${path}" != "/" ]]; do
    max=$(cat "/sys/fs/cgroup${path}/memory.max" 2>/dev/null || echo max)
    if [[ "${max}" != "max" ]]; then
      current=$(cat "/sys/fs/cgroup${path}/memory.current")
      echo $((max - current))
      return
    fi
    path=$(dirname "${path}")
  done
  echo 9223372036854775807
}

gpu_is_free() {
  local slot=$1
  local uuid=$2
  local used
  used=$(nvidia-smi --id="${slot}" --query-gpu=memory.used \
    --format=csv,noheader,nounits 2>/dev/null)
  [[ "${used}" =~ ^[0-9]+$ ]] || return 1
  ((used < 100)) || return 1
  ! nvidia-smi --query-compute-apps=gpu_uuid --format=csv,noheader,nounits 2>/dev/null \
    | grep -Fxq "${uuid}"
}

run_one() {
  local run_dir=$1
  local slot=$2
  (
    cd "${run_dir}"
    export CP2K_DATA_DIR="${data_dir}"
    export GAUXC_SKALA_MODEL="${cpu_model}"
    export GAUXC_SKALA_CUDA_MODEL="${cuda_model}"
    export CUDA_VISIBLE_DEVICES="${slot}"
    export OMP_NUM_THREADS="${threads_per_case}"
    export OMP_PROC_BIND=false
    unset OMP_PLACES
    /usr/bin/time -v taskset -c "${cpu_sets[$slot]}" \
      "${cp2k_bin}" -i input.inp -o output.out 2>time.txt
  )
}

declare -a pids=(0 0 0 0)
declare -a labels=("" "" "" "")
next=0
status=0

while true; do
  active=0
  for slot in 0 1 2 3; do
    pid=${pids[$slot]}
    if ((pid > 0)); then
      if kill -0 "${pid}" 2>/dev/null; then
        active=$((active + 1))
        continue
      fi
      if ! wait "${pid}"; then
        status=1
      fi
      printf '%s slot=%d finished case=%s\n' \
        "$(date -Is)" "${slot}" "${labels[$slot]}"
      pids[$slot]=0
      labels[$slot]=""
    fi

    if ((next < ${#queue[@]})) && \
      (( $(active_compute_count) < max_concurrent )) && \
      (( $(cgroup_memory_headroom) >= min_headroom_bytes )) && \
      gpu_is_free "${slot}" "${gpu_uuids[$slot]}"; then
      run_dir=${queue[$next]}
      next=$((next + 1))
      rm -f "${run_dir}/output.out" "${run_dir}/time.txt"
      run_one "${run_dir}" "${slot}" &
      pids[$slot]=$!
      labels[$slot]=$(basename "${run_dir}")
      active=$((active + 1))
      printf '%s slot=%d pid=%d case=%s\n' \
        "$(date -Is)" "${slot}" "${pids[$slot]}" "${labels[$slot]}"
    fi
  done

  if ((next >= ${#queue[@]} && active == 0)); then
    break
  fi
  sleep 5
done

exit "${status}"
