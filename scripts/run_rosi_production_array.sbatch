#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/home/kuehne88/cp2k-native-skala/logs/p25-%A_%a.out
#SBATCH --error=/home/kuehne88/cp2k-native-skala/logs/p25-%A_%a.err

set -euo pipefail

run_root=${RUN_ROOT:?Set RUN_ROOT}
case_list=${CASE_LIST:?Set CASE_LIST}
mode=${PRODUCTION_MODE:?Set PRODUCTION_MODE to gpu or cpu}
native_root=/home/kuehne88/cp2k-native-skala
source "$native_root/env.sh"

mapfile -t cases < <(grep -Ev '^\s*(#|$)' "$case_list")
relative=${cases[$SLURM_ARRAY_TASK_ID]:?Invalid array index}
run_dir=$run_root/$relative
output=$run_dir/output.out

if [[ -f $output ]] \
  && grep -q "SCF run converged" "$output" \
  && grep -q "PROGRAM ENDED AT" "$output"; then
  echo "Already completed: $relative"
  exit 0
fi

if [[ -f $output || -f $run_dir/time.txt ]]; then
  stamp=$(date +%Y%m%dT%H%M%S)
  mkdir -p "$run_dir/incomplete-attempts/$stamp"
  [[ -f $output ]] && mv "$output" "$run_dir/incomplete-attempts/$stamp/"
  [[ -f $run_dir/time.txt ]] && mv "$run_dir/time.txt" "$run_dir/incomplete-attempts/$stamp/"
fi

export GAUXC_SKALA_MODEL=$native_root/models/skala-1.1-rev1.fun
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_DYNAMIC=FALSE
export OMP_STACKSIZE=512M
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none

if [[ $mode == gpu ]]; then
  binary=$native_root/install/cp2k-h100/bin/cp2k.psmp
  export CP2K_DATA_DIR=$native_root/install/cp2k-h100/share/cp2k/data
  export GAUXC_SKALA_CUDA_MODEL=$native_root/models/skala-1.1-rev1-cuda.fun
  export LD_LIBRARY_PATH=$native_root/install/cp2k-h100/lib:$native_root/install/dbcsr-h100/lib:${LD_LIBRARY_PATH:-}
  export OMP_PROC_BIND=spread
  export OMP_PLACES=cores
  launch=(srun --mpi=pmix --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores)
  input_name=input.inp
elif [[ $mode == cpu ]]; then
  binary=$native_root/install/cp2k-cpu/bin/cp2k.psmp
  export CP2K_DATA_DIR=$native_root/cp2k/data
  export CUDA_VISIBLE_DEVICES=-1
  export LD_LIBRARY_PATH=$native_root/install/dbcsr-cpu-serial/lib:$native_root/install/cuda-stubs:${LD_LIBRARY_PATH:-}
  export LD_PRELOAD=$native_root/install/libtorch/lib/libtorch_cpu.so${LD_PRELOAD:+:$LD_PRELOAD}
  export OMP_PROC_BIND=FALSE
  unset OMP_PLACES
  launch=()
  sed 's/NATIVE_GRID_USE_CUDA[[:space:]]*\.TRUE\./NATIVE_GRID_USE_CUDA .FALSE./' \
    "$run_dir/input.inp" >"$run_dir/input.cpu.inp"
  input_name=input.cpu.inp
else
  echo "Unknown PRODUCTION_MODE: $mode" >&2
  exit 2
fi

printf 'host=%s\njob=%s\nstarted=%s\nroute=%s\nmode=%s\nthreads=%s\n' \
  "$(hostname)" "$SLURM_JOB_ID" "$(date --iso-8601=seconds)" "$relative" \
  "$mode" "$SLURM_CPUS_PER_TASK" >"$run_dir/execution.txt"

cd "$run_dir"
/usr/bin/time -v "${launch[@]}" "$binary" -i "$input_name" -o output.out 2>time.txt
grep -q "SCF run converged" output.out
grep -q "PROGRAM ENDED AT" output.out
