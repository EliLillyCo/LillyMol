#!/bin/bash

# Run python with PYTHONPATH pointing to the shared libraries generated
# by build_from_src.sh
# Make sure that PATH is set so that the same version of python that
# is in WORKSPACE is used for run-time.
# This will not work otherwise.

# build_from_src.sh must have been run, and one of the last things it
# does is to run copy_shared_libraries.sh, which copies the compiled
# shared libraries out of bazel-bin to ../lib

# Assume that lib, with the shared libraries, is one above where this script is.
my_dir=$(dirname $0)
PYTHONPATH=${my_dir}/../lib:${my_dir}:${PYTHONPATH} python3 "$@"
