#!/usr/bin/env bash
set -uo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 DIAGNOSTIC_ROOT CASE_LIST" >&2
  exit 2
fi

root=$1
case_list=$2
work_root=${TEROK_SKALA_ROOT:-"$HOME/work/native-skala-terok"}
binary=${CP2K_BINARY:-"$work_root/install-cpu/bin/cp2k.psmp"}
env_dir=${CP2K_SKALA_ENV:-"$HOME/micromamba/envs/cp2k-skala-gauxc"}
model=${GAUXC_SKALA_MODEL:-"$work_root/pilot-source/models/skala-1.1-rev1.fun"}
cuda_root=${CUDA_TOOLKIT_ROOT:-"$HOME/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4"}
cuda_math="$HOME/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64"
threads=${MATCHED_THREADS:-12}
slots=${MATCHED_SLOTS:-8}
cpu_offset=${MATCHED_CPU_OFFSET:-0}
seed_dir=${MATCHED_SEED_DIR:-}

export PATH="$env_dir/bin:$PATH"
export LD_LIBRARY_PATH="$work_root/install-cpu/lib:$env_dir/lib:$cuda_root/lib64:$cuda_math${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export LD_PRELOAD="$env_dir/lib/python3.12/site-packages/torch/lib/libtorch_cpu.so${LD_PRELOAD:+:$LD_PRELOAD}"
export CP2K_DATA_DIR="$work_root/cp2k/data"
export GAUXC_SKALA_MODEL="$model"
export CUDA_VISIBLE_DEVICES=-1
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=FALSE
export OMP_STACKSIZE=512M
export OMPI_MCA_hwloc_base_binding_policy=none
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
unset OMP_PLACES

mapfile -t cases < <(grep -Ev '^\s*(#|$)' "$case_list")

run_case() {
  local relative=$1 slot=$2
  local run_dir=$root/$relative
  local digest=${relative##*/}
  local source_input=input.inp
  local first=$((cpu_offset + slot * threads))
  local last=$((first + threads - 1))
  local stamp

  [[ -f $run_dir/input.restart.inp ]] && source_input=input.restart.inp
  sed 's/NATIVE_GRID_USE_CUDA[[:space:]]*\.TRUE\./NATIVE_GRID_USE_CUDA .FALSE./' \
    "$run_dir/$source_input" >"$run_dir/input.cpu.inp"
  if [[ -n $seed_dir && -f $seed_dir/$digest.wfn ]]; then
    cp "$seed_dir/$digest.wfn" "$run_dir/seed.wfn"
    awk '
      /^[[:space:]]*WFN_RESTART_FILE_NAME[[:space:]]/ { next }
      /^[[:space:]]*&DFT[[:space:]]*$/ && !inserted_wfn {
        print
        match($0, /^[[:space:]]*/)
        indent = substr($0, 1, RLENGTH)
        print indent "   WFN_RESTART_FILE_NAME seed.wfn"
        inserted_wfn = 1
        next
      }
      /^[[:space:]]*SCF_GUESS[[:space:]]/ {
        match($0, /^[[:space:]]*/)
        indent = substr($0, 1, RLENGTH)
        print indent "SCF_GUESS RESTART"
        next
      }
      { print }
    ' "$run_dir/input.cpu.inp" >"$run_dir/input.cpu.restart.inp"
    mv "$run_dir/input.cpu.restart.inp" "$run_dir/input.cpu.inp"
  fi

  if [[ -f $run_dir/output.out ]] \
    && grep -q "Density non-zero on the" "$run_dir/output.out"; then
    echo "Invalid WAVELET boundary density: $relative" >&2
    return 1
  fi

  if [[ -f $run_dir/output.out ]] \
    && grep -q "SCF run converged" "$run_dir/output.out" \
    && grep -q "PROGRAM ENDED AT" "$run_dir/output.out"; then
    echo "Already completed: $relative"
    return 0
  fi

  if [[ -f $run_dir/output.out || -f $run_dir/time.txt ]]; then
    stamp=$(date +%Y%m%dT%H%M%S)
    mkdir -p "$run_dir/incomplete-attempts/$stamp"
    [[ -f $run_dir/output.out ]] && mv "$run_dir/output.out" "$run_dir/incomplete-attempts/$stamp/"
    [[ -f $run_dir/time.txt ]] && mv "$run_dir/time.txt" "$run_dir/incomplete-attempts/$stamp/"
  fi

  echo "$(date --iso-8601=seconds) slot=$slot starting $relative CPUs=$first-$last"
  (
    set -e
    cd "$run_dir" || exit 1
    export OMP_NUM_THREADS=$threads
    setsid "$env_dir/bin/time" -v taskset -c "$first-$last" \
      "$binary" -i input.cpu.inp -o output.out 2>time.txt &
    run_pid=$!
    while kill -0 "$run_pid" 2>/dev/null; do
      if [[ -f output.out ]] && grep -q "Density non-zero on the" output.out; then
        kill -TERM -- "-$run_pid" 2>/dev/null || true
        wait "$run_pid" 2>/dev/null || true
        echo "Invalid WAVELET boundary density: $relative" >&2
        exit 1
      fi
      sleep 5
    done
    wait "$run_pid"
    grep -q "SCF run converged" output.out
    grep -q "PROGRAM ENDED AT" output.out
    if grep -q "Density non-zero on the" output.out; then
      echo "Invalid WAVELET boundary density: $relative" >&2
      exit 1
    fi
  )
}

declare -a pids labels
for ((slot=0; slot<slots; slot++)); do
  pids[$slot]=0
  labels[$slot]=""
done

next=0
active=0
failures=0
while ((next < ${#cases[@]} || active > 0)); do
  for ((slot=0; slot<slots; slot++)); do
    pid=${pids[$slot]}
    if ((pid > 0)) && ! kill -0 "$pid" 2>/dev/null; then
      if wait "$pid"; then
        echo "$(date --iso-8601=seconds) slot=$slot completed ${labels[$slot]}"
      else
        echo "$(date --iso-8601=seconds) slot=$slot failed ${labels[$slot]}" >&2
        failures=$((failures + 1))
      fi
      pids[$slot]=0
      labels[$slot]=""
      active=$((active - 1))
    fi

    if ((pids[$slot] == 0 && next < ${#cases[@]})); then
      relative=${cases[$next]}
      next=$((next + 1))
      run_case "$relative" "$slot" &
      pids[$slot]=$!
      labels[$slot]="$relative"
      active=$((active + 1))
    fi
  done
  sleep 10
done

echo "Finished matched validation with $failures failures at $(date --iso-8601=seconds)"
exit "$failures"
