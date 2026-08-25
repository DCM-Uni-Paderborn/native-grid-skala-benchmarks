#!/usr/bin/env bash
set -uo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 DIAGNOSTIC_ROOT CASE_LIST" >&2
  exit 2
fi

root=$1
case_list=$2
binary=${CP2K_BINARY:?Set CP2K_BINARY to the cp2k.psmp executable}

if [[ ! -d $root ]]; then
  echo "Diagnostic root does not exist: $root" >&2
  exit 2
fi
if [[ ! -f $case_list ]]; then
  echo "Case list does not exist: $case_list" >&2
  exit 2
fi

failures=0
while IFS= read -r relative || [[ -n $relative ]]; do
  [[ -z $relative || $relative == \#* ]] && continue
  run_dir=$root/$relative
  input=$run_dir/input.inp
  [[ -f $run_dir/input.restart.inp ]] && input=$run_dir/input.restart.inp
  if [[ ! -f $input ]]; then
    echo "Missing input: $input" >&2
    failures=$((failures + 1))
    continue
  fi

  if [[ -f $run_dir/output.out ]] \
    && grep -q "SCF run converged" "$run_dir/output.out" \
    && grep -q "PROGRAM ENDED AT" "$run_dir/output.out"; then
    echo "Already completed: $relative"
    continue
  fi

  attempt_dir=$run_dir/incomplete-attempts
  if [[ -f $run_dir/output.out || -f $run_dir/time.txt ]]; then
    stamp=$(date +%Y%m%dT%H%M%S)
    mkdir -p "$attempt_dir/$stamp"
    [[ -f $run_dir/output.out ]] && mv "$run_dir/output.out" "$attempt_dir/$stamp/"
    [[ -f $run_dir/time.txt ]] && mv "$run_dir/time.txt" "$attempt_dir/$stamp/"
  fi

  echo "$(date --iso-8601=seconds) starting $relative"
  (
    cd "$run_dir" || exit 1
    /usr/bin/time -v "$binary" -i "$(basename "$input")" -o output.out 2>time.txt
    grep -q "SCF run converged" output.out
    grep -q "PROGRAM ENDED AT" output.out
  ) || {
    echo "$(date --iso-8601=seconds) failed $relative" >&2
    failures=$((failures + 1))
  }
done < "$case_list"

echo "Finished matched diagnostics with $failures failures at $(date --iso-8601=seconds)"
exit "$failures"
