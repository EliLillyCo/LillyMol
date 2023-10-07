#!/bin/bash

echo "Builds and installs LillyMol executables"
echo "The assumption is that WORKSPACE and build_deps/install.bzl"
echo "have both been configured."
echo ""
echo "First task is to run C++ unit tests."


# note that we do not check return codes from any of the invocations.

set -x

jobs='8'

# Options that are used by all bazelisk invocations.
# bazel will not work on an NFS mounted file system. So if you are on an NFS
# file system, you must specify a value for --output_user_root that is
# locally mounted.
bazel_options="--output_user_root=/node/scratch/ian"
build_options="--cxxopt=-DGIT_HASH=\"$(git rev-parse --short --verify HEAD)\" --cxxopt=-DTODAY=\"$(date +%Y-%b-%d)\" --jobs=${jobs} -c opt --enable_bzlmod --experimental_cc_shared_library"

echo ${common_options}

# First task is unit tests

bazelisk ${bazel_options} test ${build_options} Foundational/...:all
bazelisk ${bazel_options} test ${build_options} Molecule_Lib:all
bazelisk ${bazel_options} test ${build_options} Molecule_Tools:all
bazelisk ${bazel_options} test ${build_options} Utilities/...:all

# Once the tests pass, then executables can be built

echo "Building tools"
bazelisk ${bazel_options} build ${build_options} Molecule_Tools:all
bazelisk ${bazel_options} build ${build_options} Foundational/iw_tdt:all
bazelisk ${bazel_options} build ${build_options} Utilities/...:all

bazelisk ${bazel_options} build ${build_options} pybind:all

# Now install the targets

echo "Installing tools"
bazelisk ${bazel_options} run ${build_options} Foundational/iw_tdt:install
bazelisk ${bazel_options} run ${build_options} Molecule_Tools:install
bazelisk ${bazel_options} run ${build_options} Utilities/General:install
bazelisk ${bazel_options} run ${build_options} Utilities/GFP_Tools:install
bazelisk ${bazel_options} run ${build_options} Utilities/Distance_Matrix:install
bazelisk ${bazel_options} run ${build_options} Utilities/BerkeleyDB:install

# The python shared libraries. Note that PYTHONPATH will need to be adjusted,
# or copy these to your default PYTHONPATH.
./copy_shared_libraries.sh ../lib64
