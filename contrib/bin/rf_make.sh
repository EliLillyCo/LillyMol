#!/usr/bin/env bash
if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $(realpath $0))/../..
fi

export PATH=${LILLYMOL_HOME}/contrib/bin:$PATH

exec python ${LILLYMOL_HOME}/contrib/bin/xgbd/rf_make.py "$@"
