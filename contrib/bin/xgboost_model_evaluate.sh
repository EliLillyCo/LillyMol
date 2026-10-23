#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  LILLYMOL_HOME=$(dirname $0)/../..
fi

export LD_LIBRARY_PATH=${LILLYMOL_HOME}/third_party/lib:${LD_LIBRARY_PATH}
exec ${LILLYMOL_HOME}/bin/$(uname)/xgboost_model_evaluate "$@"
