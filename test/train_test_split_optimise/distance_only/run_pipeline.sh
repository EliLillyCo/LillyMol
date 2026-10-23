#!/usr/bin/env bash
set -euo pipefail

exe="$1"
indir="$2"
outdir="$3"
mode="${4:-distance}"

repo="${LILLYMOL_HOME:-$(cd "$indir/../../.." && pwd)}"
export LILLYMOL_HOME="$repo"
installed_bin="$repo/bin/$(uname)"
export PATH="$installed_bin:$PATH"

gfp_nearneighbours="$installed_bin/gfp_nearneighbours_single_file"

"$repo/contrib/bin/gfp_make.sh" "$indir/ttopt.smi" > ttopt.gfp
"$gfp_nearneighbours" -n 7 -S ttopt.nn ttopt.gfp

args=(-Z seed=12345 -f 0.5 -n 1 -S SPLIT -o 20000)
if [[ "$mode" == "activity" ]]; then
  args+=(-A "$indir/activity.txt:buckets=2" -Y 10)
fi
args+=(ttopt.nn)

"$exe" "${args[@]}"
