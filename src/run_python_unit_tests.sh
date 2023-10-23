#!/bin/bash

# Run all the python unit tests found in the pybind directory.

here=$(dirname $0)

if [[ ! -s "${here}/../lib" ]] ; then
  echo "No shared libraries available ${here}, python unit tests not done"
  exit 1
fi

if [[ ! -d "${here}/pybind" ]] ; then
  echo "pybind not found ${here}"
  exit 1
fi

for file in ${here}/pybind/*_test.py ; do
  ${here}/run_python.sh ${file}
done

echo 'Python unit tests complete'
