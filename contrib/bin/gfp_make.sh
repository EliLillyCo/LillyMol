#!/usr/bin/env bash

here=$(dirname $(realpath $0))

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname ${here}))
fi

PATH=${here}:${LILLYMOL_HOME}/bin/${uname}:$PATH

exec perl ${here}/gfp_make.pl "$@"
