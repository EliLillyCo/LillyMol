#!/bin/bash

echo "Builds and installs LillyMol executables"
echo "The assumption is that WORKSPACE and build_deps/install.bzl"
echo "have both been configured."
echo ""
echo "First task is to build and run C++ unit tests."

# note that we do not check return codes from any of the invocations.

# If bazelisk is available use it, otherwise try bazel, otherwise fail.

if [[ ! -z "$(type -p bazelisk)" ]] ; then
  bazel='bazelisk'
elif [[ ! -z "$(type -p bazel)" ]] ; then
  bazel='bazel'
else
  echo "No bazel or bazelisk, build will fail"
  bazel='bazelisk'
fi

# Adjust to local resource availability.
jobs='8'

declare -i inside_lilly
if [[ $(hostname -d) =~ 'lilly.com' ]] ; then
  inside_lilly=1
else
  inside_lilly=0
fi

# Options that are used by all bazelisk invocations.

# bazel will not work on an NFS mounted file system. So if you are on an NFS
# file system, you must specify a value for --output_user_root that is
# locally mounted.
# Note that the bazel cache can get quite large, 1-2GB.

# If inside Lilly, some local scratch storage
if [[ ${inside_lilly} -eq 1 && -d '/node/scratch' ]] ; then
  bazel_options="--output_user_root=/node/scratch/${USER}"
elif [[ $(df -TP ${HOME}) =~ 'nfs' ]] ; then
  echo "Your HOME dir is an NFS mounted file system. bazel will not work."
  echo "Will attempt to use /tmp/ for bazel cache, that will need to be changed."
  bazel_options='--output_user_root=/tmp'
else
  # Even if outside Lilly, you may still need to set this
  bazel_options=""
fi

build_options="--cxxopt=-DGIT_HASH=\"$(git rev-parse --short --verify HEAD)\" --cxxopt=-DTODAY=\"$(date +%Y-%b-%d)\" --jobs=${jobs} -c opt"

# All BerkeleyDB things are now isolated, so this probably makes no difference.
if [[ ! -v BUILD_BDB ]] ; then
  build_options+=' --build_tag_filters=-berkeleydb'
fi

# First task is unit tests

${bazel} ${bazel_options} test ${build_options} Foundational/...:all
${bazel} ${bazel_options} test ${build_options} Molecule_Lib:all
${bazel} ${bazel_options} test ${build_options} Molecule_Tools:all
${bazel} ${bazel_options} test ${build_options} Utilities/...:all

# Currently no tests in these.
if [[ -v BUILD_BDB ]] ; then
  # ${bazel} ${bazel_options} test ${build_options} BerkeleyDB:all
  # ${bazel} ${bazel_options} test ${build_options} Molecule_Tools_Bdb:all
  echo ""
fi

# Once the tests run, then executables can be built.
# donor_acceptor_test frequently fails due to lack of supporting files.

if [[ ! -v BUILD_LIBRARY_ONLY ]] ; then
  echo "Building tools"
  ${bazel} ${bazel_options} build ${build_options} Molecule_Tools:all
  ${bazel} ${bazel_options} build ${build_options} Obsolete:all
  ${bazel} ${bazel_options} build ${build_options} Obsolete/Descriptor_Similarity:all
  ${bazel} ${bazel_options} build ${build_options} Foundational/iw_tdt:all
  ${bazel} ${bazel_options} build ${build_options} Utilities/...:all
fi

if [[ ${inside_lilly} -eq 1 || -v BUILD_VENDOR ]] ; then
  ${bazel} ${bazel_options} build ${build_options} Vendor/...:all
fi

if [[ -v BUILD_BDB ]] ; then
  ${bazel} ${bazel_options} build ${build_options} BerkeleyDB:all
  ${bazel} ${bazel_options} build ${build_options} Molecule_Tools_Bdb:all
fi

# Now install the targets

if [[ ! -v BUILD_LIBRARY_ONLY ]] ; then
  echo "Installing tools"
  ${bazel} ${bazel_options} run ${build_options} Foundational/iw_tdt:install
  ${bazel} ${bazel_options} run ${build_options} Molecule_Tools:install
  ${bazel} ${bazel_options} run ${build_options} Obsolete:install
  ${bazel} ${bazel_options} run ${build_options} Obsolete/Descriptor_Similarity:install
  ${bazel} ${bazel_options} run ${build_options} Utilities/General:install
  ${bazel} ${bazel_options} run ${build_options} Utilities/GFP_Tools:install
  ${bazel} ${bazel_options} run ${build_options} Utilities/GFP_Knn:install
  ${bazel} ${bazel_options} run ${build_options} Utilities/Distance_Matrix:install
fi

if [[ ${inside_lilly} -eq 1 || -v BUILD_VENDOR ]] ; then
  ${bazel} ${bazel_options} run ${build_options} Vendor:install
fi

if [[ -v BUILD_BDB ]] ; then
  ${bazel} ${bazel_options} run ${build_options} BerkeleyDB:install
  ${bazel} ${bazel_options} run ${build_options} Molecule_Tools_Bdb:install
fi

# Python if requested, build, install and test.
# Note that PYTHONPATH will need to be adjusted, or copy the shared
# libraries from LillyMol/lib to your default PYTHONPATH.

if [[ -v BUILD_PYTHON ]] ; then
  ${bazel} ${bazel_options} build ${build_options} pybind:all
  ./copy_shared_libraries.sh ../lib
fi
