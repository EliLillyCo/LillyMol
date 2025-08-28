#!/bin/bash

here=$(dirname $(realpath $0))

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname ${here}))
fi

exec perl ${here}/gfp_make.pl "$@"

