#!/usr/bin/env bash

here=$(dirname $0)
fname=$(basename $0)
exec python3 ${here}/${fname%%.sh}.py "$@"
