#!/usr/bin/env bash
#SBATCH --job-name=skala-matched-validation
#SBATCH --partition=gpu-h100-casus
#SBATCH --account=casus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --gres=gpu:4
#SBATCH --mem=900G
#SBATCH --time=24:00:00
#SBATCH --output=/home/kuehne88/cp2k-native-skala/logs/matched-validation-%j.out
#SBATCH --error=/home/kuehne88/cp2k-native-skala/logs/matched-validation-%j.err

set -euo pipefail

root=/home/kuehne88/cp2k-native-skala
run_root=${MATCHED_ROOT:-${root}/validation-uzh/matched-diagnostics}
case_list=${MATCHED_CASE_LIST:-${run_root}/aconf-native-cell-cases.txt}
cp2k_prefix=${root}/install/cp2k-h100
dbcsr_prefix=${root}/install/dbcsr-h100
spglib_prefix=${root}/install/spglib
libtorch_prefix=${root}/install/libtorch
binary=${cp2k_prefix}/bin/cp2k.psmp

source "${root}/env.sh"
export LD_LIBRARY_PATH="${cp2k_prefix}/lib:${dbcsr_prefix}/lib:${spglib_prefix}/lib:${libtorch_prefix}/lib:${LD_LIBRARY_PATH}"
export CP2K_DATA_DIR="${cp2k_prefix}/share/cp2k/data"
export GAUXC_SKALA_MODEL="${root}/models/skala-1.1-rev1.fun"
export GAUXC_SKALA_CUDA_MODEL="${root}/models/skala-1.1-rev1-cuda.fun"
export OMP_NUM_THREADS=32
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_STACKSIZE=512M
export OPENBLAS_NUM_THREADS=1

mapfile -t cases < <(grep -Ev '^\s*(#|$)' "$case_list")

run_one() {
  local relative=$1
  local run_dir=${run_root}/${relative}
  local input=input.inp
  local stamp

  [[ -f ${run_dir}/input.restart.inp ]] && input=input.restart.inp

  if [[ -f ${run_dir}/output.out ]] \
    && grep -q "SCF run converged" "${run_dir}/output.out" \
    && grep -q "PROGRAM ENDED AT" "${run_dir}/output.out"; then
    echo "Already completed: ${relative}"
    return 0
  fi

  if [[ -f ${run_dir}/output.out || -f ${run_dir}/time.txt ]]; then
    stamp=$(date +%Y%m%dT%H%M%S)
    mkdir -p "${run_dir}/incomplete-attempts/${stamp}"
    [[ -f ${run_dir}/output.out ]] \
      && mv "${run_dir}/output.out" "${run_dir}/incomplete-attempts/${stamp}/"
    [[ -f ${run_dir}/time.txt ]] \
      && mv "${run_dir}/time.txt" "${run_dir}/incomplete-attempts/${stamp}/"
  fi

  cd "$run_dir"
  echo "$(date --iso-8601=seconds) starting ${relative}"
  /usr/bin/time -v "$binary" -i "$input" -o output.out >launcher.out 2>time.txt
  grep -q "SCF run converged" output.out
  grep -q "PROGRAM ENDED AT" output.out
  echo "$(date --iso-8601=seconds) completed ${relative}"
}
export -f run_one
export run_root binary

status=0
for ((offset=0; offset<${#cases[@]}; offset+=4)); do
  pids=()
  for ((slot=0; slot<4 && offset+slot<${#cases[@]}; slot++)); do
    relative=${cases[offset+slot]}
    srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=32 --gres=gpu:1 \
      --mem=215G --cpu-bind=cores \
      bash -lc "run_one '$relative'" &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    wait "$pid" || status=1
  done
done

exit "$status"
