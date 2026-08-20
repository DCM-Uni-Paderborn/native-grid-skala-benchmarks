#!/usr/bin/env bash
#SBATCH --job-name=native-skala-diet
#SBATCH --partition=gpu-b200-casus
#SBATCH --account=casus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --gres=gpu:4
#SBATCH --mem=1T
#SBATCH --time=48:00:00
#SBATCH --output=/home/kuehne88/cp2k-native-skala/logs/%x-%j.out
#SBATCH --error=/home/kuehne88/cp2k-native-skala/logs/%x-%j.err

set -euo pipefail

root=/home/kuehne88/cp2k-native-skala/validation-uzh/diet-production
scripts=/home/kuehne88/cp2k-native-skala/validation-uzh

mkdir -p /home/kuehne88/cp2k-native-skala/logs
MAX_B200_JOBS=4 \
MIN_JOB_MEMORY_HEADROOM_BYTES=130000000000 \
bash "${scripts}/run_b200_pool.sh" "${root}"
