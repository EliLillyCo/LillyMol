#! /bin/bash

# LILLYMOL_HOME and BUILD_DIR should be exported, use my location if not set.

if [[ ! -v LILLYMOL_HOME ]] ; then
  me=$(readlink -f $0)
  up=$(dirname ${me})  # LillyMol/test
  export LILLYMOL_HOME=$(dirname ${up})
fi

if [[ ! -v BUILD_DIR ]] ; then
  export BUILD_DIR=$(uname)
fi

if [[ ! -v PYTHONPATH ]] ; then
  export PYTHONPATH="${LILLYMOL_HOME}/contrib/python:${LILLYMOL_HOME}/contrib/python/pybase"
elif [[ ! $PYTHONPATH =~ "contrib/python" ]]; then
  export PYTHONPATH="${LILLYMOL_HOME}/contrib/python:${LILLYMOL_HOME}/contrib/python/pybase:${PYTHONPATH}"
fi

for dir in ./*; do
    if [ -d "$dir" ]; then
        current_path= pwd &> /dev/null
        pushd "$dir" &> /dev/null
        if [ -f "run_test.sh" ]
        then
            ./run_test.sh
        fi
        popd $current_path &> /dev/null
    fi
done
