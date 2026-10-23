#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $0)/../..
fi

script_dir=$(dirname $0)
export PYTHONPATH=${script_dir}:${PYTHONPATH}
exec python ${script_dir}/xgbd/xgbd_make.py "$@"
