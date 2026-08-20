#!/usr/bin/env bash
set -euo pipefail

root=${1:?Usage: run_b200_pool.sh RUN_DIRECTORY}
source "${SKALA_ENV:-/home/kuehne88/cp2k-native-skala/env.sh}"

cp2k_root=${CP2K_NATIVE_ROOT:-/home/kuehne88/cp2k-native-skala}
cp2k_bin=${CP2K_BIN:-${cp2k_root}/install/cp2k-b200/bin/cp2k.psmp}
data_dir=${CP2K_DATA_DIR:-${cp2k_root}/install/cp2k-b200/share/cp2k/data}
cpu_model=${GAUXC_SKALA_MODEL:-${cp2k_root}/models/skala-1.1-rev1.fun}
cuda_model=${GAUXC_SKALA_CUDA_MODEL:-${cp2k_root}/models/skala-1.1-rev1-cuda.fun}
cpu_sets=(32-63 64-95 128-159 160-191)

mapfile -t gpu_uuids < <(
  nvidia-smi --query-gpu=uuid --format=csv,noheader,nounits
)
if ((${#gpu_uuids[@]} != 4)); then
  echo "Expected four visible GPUs, found ${#gpu_uuids[@]}" >&2
  exit 1
fi

mapfile -t queue < <(
  find "${root}" -mindepth 1 -maxdepth 1 -type d | sort -V | while read -r run_dir; do
    if [[ ! -f "${run_dir}/output.out" ]] || ! grep -q "SCF run converged" "${run_dir}/output.out"; then
      printf '%s\n' "${run_dir}"
    fi
  done
)

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
    export OMP_NUM_THREADS=32
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

    if ((next < ${#queue[@]})) && gpu_is_free "${slot}" "${gpu_uuids[$slot]}"; then
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
