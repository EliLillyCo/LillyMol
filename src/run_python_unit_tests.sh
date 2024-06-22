#!/usr/bin/env bash

# Run all the python unit tests found in the pybind directory.

here=$(dirname $0)

if [[ ! -v PYTHONPATH ]] ; then
  export PYTHONPATH=${here}
fi

if [[ ! -s "${here}/../lib" ]] ; then
  echo "No shared libraries available ${here}, python unit tests not done"
  exit 1
fi

if [[ ! -d "${here}/pybind" ]] ; then
  echo "pybind not found ${here}"
  exit 1
fi

run_python="${here}/../run_python.sh"
if [[ ! -x ${run_python} ]] ; then
  echo "Where is ${run_python}"
  exit 1
fi

for file in ${here}/pybind/*_test.py ; do
  ${run_python} ${file}
done

echo 'Python unit tests complete'
