#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/native-skala-spark
toolchain=/home/kuehne88/cp2k-gpu-work/cp2k-toolchain-gpu-systemmpi
torch=/home/kuehne88/cp2k-gauxc-skala-toolchain/pytorch-cu130-venv/lib/python3.12/site-packages/torch

set +u
source "$toolchain/setup"
set -u
export LD_LIBRARY_PATH="$toolchain/dbcsr-2.9.1/lib:$torch/lib:${LD_LIBRARY_PATH:-}"
export LD_PRELOAD=/usr/lib/aarch64-linux-gnu/libgfortran.so.5
export CP2K_DATA_DIR=$root/cp2k/data
export GAUXC_SKALA_MODEL=/home/kuehne88/skala-models-2026.8/skala-1.1-rev1.fun
export GAUXC_SKALA_CUDA_MODEL=/home/kuehne88/skala-models-2026.8/skala-1.1-rev1-cuda.fun
export CP2K_BINARY=$root/install-gpu/bin/cp2k.psmp
export PRODUCTION_THREADS=20
export PRODUCTION_LOG=$root/logs/p25-production.log

exec "$root/run_gpu_production_queue.sh" \
  "$root/diet-production-p25" \
  "$root/diet-production-p25/queues/spark-gpu.txt"
