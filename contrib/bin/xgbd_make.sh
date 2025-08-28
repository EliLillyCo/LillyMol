#!/bin/bash
if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $0)/../..
fi

exec python ${LILLYMOL_HOME}/contrib/bin/xgbd/xgbd_make.py "$@"
