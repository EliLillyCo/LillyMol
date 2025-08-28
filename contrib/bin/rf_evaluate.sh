#!/bin/bash

if [[ ! LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $0)/../..
fi

export PATH=${LILLYMOL_HOME}/contrib/bin:$PATH

exec ruby ${LILLYMOL_HOME}/contrib/bin/xgbd/rf_evaluate.rb "$@"
