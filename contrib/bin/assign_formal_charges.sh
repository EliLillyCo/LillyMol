#!/bin/bash

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $0)))
fi

charges="${LILLYMOL_HOME}/data/queries/charges/queries"

exe="${LILLYMOL_HOME}/bin/$(uname)/fileconv"
exec ${exe} -i mdlD -i mdlT -g all -E autocreate -i ICTE -N F:${charges} -S - "$@"
