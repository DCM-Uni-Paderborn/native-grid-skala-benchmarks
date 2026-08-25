#!/usr/bin/env bash
set -uo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 PRODUCTION_ROOT CASE_LIST" >&2
  exit 2
fi

root=$1
case_list=$2
binary=${CP2K_BINARY:?Set CP2K_BINARY}
model=${GAUXC_SKALA_MODEL:?Set GAUXC_SKALA_MODEL}
threads=${PRODUCTION_THREADS:-20}
log=${PRODUCTION_LOG:-"$root/gpu-queue.log"}

mkdir -p "$(dirname "$log")"
exec > >(tee -a "$log") 2>&1

export OMP_NUM_THREADS=$threads
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=FALSE
export OMP_STACKSIZE=512M
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none
export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}

failures=0
while IFS= read -r relative || [[ -n $relative ]]; do
  [[ -z $relative || $relative == \#* ]] && continue
  run_dir=$root/$relative
  output=$run_dir/output.out

  if [[ -f $output ]] \
    && grep -q "SCF run converged" "$output" \
    && grep -q "PROGRAM ENDED AT" "$output"; then
    echo "Already completed: $relative"
    continue
  fi

  if [[ -f $output || -f $run_dir/time.txt ]]; then
    stamp=$(date +%Y%m%dT%H%M%S)
    mkdir -p "$run_dir/incomplete-attempts/$stamp"
    [[ -f $output ]] && mv "$output" "$run_dir/incomplete-attempts/$stamp/"
    [[ -f $run_dir/time.txt ]] && mv "$run_dir/time.txt" "$run_dir/incomplete-attempts/$stamp/"
  fi

  printf 'host=%s\nstarted=%s\nroute=%s\nthreads=%s\n' \
    "$(hostname)" "$(date --iso-8601=seconds)" "$relative" "$threads" \
    >"$run_dir/execution.txt"
  echo "$(date --iso-8601=seconds) starting $relative"
  (
    cd "$run_dir" || exit 1
    /usr/bin/time -v "$binary" -i input.inp -o output.out 2>time.txt
    grep -q "SCF run converged" output.out
    grep -q "PROGRAM ENDED AT" output.out
  ) || {
    echo "$(date --iso-8601=seconds) failed $relative" >&2
    failures=$((failures + 1))
  }
done <"$case_list"

echo "Finished GPU queue with $failures failures at $(date --iso-8601=seconds)"
exit "$failures"
