#!/bin/bash
# Create a new test

set -e

# This should be within an existing top level directory
# new_test.sh dicer/new_dicer_test

test_name=$1

if [[ -z "${test_name}" ]] ; then
  echo "Must specify test and case" >&2
  exit 1
fi

mkdir ${test_name}
mkdir ${test_name}/in
mkdir ${test_name}/out
touch ${test_name}/out/stdout

