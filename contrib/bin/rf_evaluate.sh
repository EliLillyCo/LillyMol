#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $(realpath $0))/../..
fi

export PATH=${LILLYMOL_HOME}/contrib/bin:$PATH

export PYTHONPATH=$(dirname $0):$PYTHONPATH

exec ruby ${LILLYMOL_HOME}/contrib/bin/xgbd/rf_evaluate.rb "$@"
