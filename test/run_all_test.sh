#!/usr/bin/env bash

export UNAME=$(uname)
export OSTYPE=$(uname)

if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $(realpath $(dirname $0)))
fi

echo "LILLYMOL_HOME ${LILLYMOL_HOME}"
export PATH=${LILLYMOL_HOME}/bin/${OSTYPE}:$PATH

exec ruby run_all_test.rb "$@"
