#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  LILLYMOL_HOME="$(dirname $(realpath $0))/../../"
fi

executable="${LILLYMOL_HOME}/bin/Linux/jwdp"

if [[ ! -x ${executable} ]] ; then
  executable='jwdp'
fi

DATA_DIR=${LILLYMOL_HOME}/data
QUERIES=${DATA_DIR}/queries

exec $executable -A I -A C -i ignore_bad_chiral "$@" 
