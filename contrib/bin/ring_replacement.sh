#!/usr/bin/env bash
here=$(dirname $0)
exec python "${0%%.sh}.py" --collections_dir ${here} "$@"
