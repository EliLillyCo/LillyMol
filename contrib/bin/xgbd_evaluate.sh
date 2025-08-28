#!/bin/bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $0)/../..
fi

export PATH=${LILLYMOL_HOME}/contrib/bin:$PATH

exec ruby ${LILLYMOL_HOME}/contrib/bin/xgbd/xgbd_evaluate.rb "$@"
