#!/bin/bash
if [[ ! LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $0)/../..
fi

exec python ${LILLYMOL_HOME}/contrib/bin/xgbd/rf_make.py "$@"
