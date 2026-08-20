#!/usr/bin/env bash
#SBATCH --job-name=native-skala-species
#SBATCH --partition=gpu-b200-casus
#SBATCH --account=casus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --array=0-479%8
#SBATCH --output=/home/kuehne88/cp2k-native-skala/logs/species-%A_%a.out
#SBATCH --error=/home/kuehne88/cp2k-native-skala/logs/species-%A_%a.err

set -euo pipefail

native_root=/home/kuehne88/cp2k-native-skala
run_root=${native_root}/validation-uzh/diet-production
source "${native_root}/env.sh"

mapfile -t inputs < <(find "${run_root}" -type f -name input.inp | sort)
if ((${#inputs[@]} != 480)); then
  echo "Expected 480 production inputs, found ${#inputs[@]}" >&2
  exit 1
fi

input=${inputs[$SLURM_ARRAY_TASK_ID]}
run_dir=$(dirname "${input}")
output=${run_dir}/output.out

if [[ -f "${output}" ]] && grep -q "SCF run converged" "${output}"; then
  echo "Already converged: ${run_dir}"
  exit 0
fi

export CP2K_DATA_DIR=${native_root}/install/cp2k-b200/share/cp2k/data
export GAUXC_SKALA_MODEL=${native_root}/models/skala-1.1-rev1.fun
export GAUXC_SKALA_CUDA_MODEL=${native_root}/models/skala-1.1-rev1-cuda.fun
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

mkdir -p "${native_root}/logs"
rm -f "${output}" "${run_dir}/time.txt"
cd "${run_dir}"

echo "Starting task ${SLURM_ARRAY_TASK_ID}: ${run_dir}"
/usr/bin/time -v \
  "${native_root}/install/cp2k-b200/bin/cp2k.psmp" \
  -i input.inp -o output.out 2>time.txt

grep -q "SCF run converged" output.out
