#!/usr/bin/env bash
set -eo pipefail
profile=${1:?host}; shift
case "$profile" in
  terok)
    root=/home/kuehne88/work/native-skala-terok
    prefix=/home/kuehne88/micromamba/envs/cp2k-skala-gauxc
    export PATH=$prefix/bin:$PATH
    export GAUXC_SKALA_MODEL=$root/pilot-source/models/skala-1.1-rev1.fun
    export TIME_BINARY=$prefix/bin/time
    ;;
  rosi)
    root=/home/kuehne88/cp2k-native-skala
    source "$root/env.sh"
    export GAUXC_SKALA_MODEL=$root/models/skala-1.1-rev1.fun
    export LD_LIBRARY_PATH=$root/x23-mini/runtime/build-37aedacbb3/src:$root/install/dbcsr-cpu-serial/lib:$root/install/cuda-stubs:${LD_LIBRARY_PATH:-}
    export LD_PRELOAD=$root/install/libtorch/lib/libtorch_cpu.so
    export TIME_BINARY=/usr/bin/time
    ;;
  *) exit 2 ;;
esac
export CP2K_DATA_DIR=$root/cp2k/data
export CUDA_VISIBLE_DEVICES=-1 PYTHONDONTWRITEBYTECODE=1
exec python3 -B "$(dirname "$0")/run_prepared.py" "$@" --host "$profile" --runtime "$root/x23-mini/runtime/runtime-build.json"
