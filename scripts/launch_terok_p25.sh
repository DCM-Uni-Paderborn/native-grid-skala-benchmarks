#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/native-skala-terok
production=$root/diet-production-p25

# The padding controls use the same large atom grids. Let them release memory
# before the eight-slot production queue starts.
while tmux ls 2>/dev/null | grep -Eq \
  'padding-p2[245]|padding-high-ascending|native-skala-p25-early'; do
  sleep 30
done

export CP2K_BINARY=$root/build-gauxc-native-host-current/bin/cp2k.psmp
export MATCHED_THREADS=12
export MATCHED_SLOTS=8
export MATCHED_CPU_OFFSET=0
unset MATCHED_SEED_DIR

exec "$root/validation-scripts/run_terok_matched_validation.sh" \
  "$production" "$production/queues/terok-cpu-main.txt"
