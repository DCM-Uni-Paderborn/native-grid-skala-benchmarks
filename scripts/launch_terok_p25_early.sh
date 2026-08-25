#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/native-skala-terok
production=$root/diet-production-p25

export CP2K_BINARY=$root/build-gauxc-native-host-current/bin/cp2k.psmp
export MATCHED_THREADS=12
export MATCHED_SLOTS=3
export MATCHED_CPU_OFFSET=48
unset MATCHED_SEED_DIR

exec "$root/validation-scripts/run_terok_matched_validation.sh" \
  "$production" "$production/queues/terok-cpu-early.txt"
