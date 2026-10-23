#!/usr/bin/env bash

here=$(dirname $(realpath $0))

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname ${here}))
fi

export PYTHONPATH="${LILLYMOL_HOME}/src:${LILLYMOL_HOME}/lib:${PYTHONPATH:-}"

exec python "${LILLYMOL_HOME}/src/Utilities/GFP_Tools/gfp_http_server.py" "$@"
