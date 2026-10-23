#!/usr/bin/env bash

# Run the nanobind-backed LillyMol Python unit tests from the source tree.

here=$(dirname "$0")

if [[ ! -v PYTHONPATH ]] ; then
  export PYTHONPATH=${here}
fi

# If tmpdir is not set, multiple people running absl tests will collide.
# Unique tmpdir for each user.
if [[ ! -v TMPDIR ]] ; then
  export "TMPDIR=/tmp/absl_testing_${USER}"
fi

libdir="${here}/../lib"
if [[ ! -d "${libdir}" ]] ; then
  echo "No shared library directory ${libdir}, python unit tests not done"
  exit 1
fi

if [[ ! -s "${libdir}/lillymol.so" ]] ; then
  echo "No lillymol module ${libdir}/lillymol.so, python unit tests not done"
  exit 1
fi

if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $(realpath ${here}))
fi

if [[ ! -d "${here}/nanobind" ]] ; then
  echo "nanobind not found ${here}"
  exit 1
fi

run_python="${here}/../run_python.sh"
if [[ ! -x ${run_python} ]] ; then
  echo "Where is ${run_python}"
  exit 1
fi

declare -a tests=(
  "${here}/nanobind/lillymol_nb_doc_test.py"
  "${here}/nanobind/lillymol_nb_test.py"
)

if [[ -s "${libdir}/lillymol_gfp_server.so" ]] ; then
  tests+=("${here}/nanobind/gfp_http_server_test.py")
fi

if [[ -s "${libdir}/lillymol_bdb.so" ]] ; then
  tests+=("${here}/nanobind/lillymol_nb_bdb_test.py")
fi

if [[ -s "${libdir}/lillymol_depict.so" ]] ; then
  tests+=("${here}/nanobind/lillymol_nb_depict_test.py")
fi

declare -i failures=0
for file in "${tests[@]}" ; do
  ${run_python} ${file}
  if [[ $? -ne 0 ]] ; then
    echo "${file} failed"
    failures=$(($failures + 1))
  fi
done

echo 'Python unit tests complete'
if [[ ${failures} -gt 0 ]] ; then
  echo "${failures} python tests failed"
  exit 1
fi
