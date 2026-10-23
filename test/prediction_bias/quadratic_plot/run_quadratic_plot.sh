#!/usr/bin/env bash
set -euo pipefail

exe="$1"
indir="$2"
outdir="$3"

"$exe" -q -A "$indir/activity.txt" -R plot.py "$indir/predicted.txt" > fit.out

if [[ ! -s plot.py ]]; then
  echo "plot.py was not created or is empty" >&2
  exit 1
fi

python3 -m py_compile plot.py
cat fit.out
printf 'plot.py nonempty
'
