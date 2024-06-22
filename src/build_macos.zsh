#!/usr/bin/env zsh

# check if REPO_HOME is set from Makefile
# if directly calling this script, set it
if [[ ! -v REPO_HOME ]] ; then
    echo "REPO_HOME is not set, initializing..."
    REPO_HOME=$(dirname $(dirname  "$(readlink -f "$0")"))
fi

echo "REPO_HOME is: $REPO_HOME"

# If bazelisk is available use it, otherwise try bazel, otherwise fail.
if [[ ! -z "$(type -p bazelisk)" ]] ; then
    bazel='bazelisk'
elif [[ ! -z "$(type -p bazel)" ]] ; then
    bazel='bazel'
else
    echo "No bazel/bazelisk, see README.md" && exit 1
fi

# Step 1: udpate WORKSPACE and install.bzl
# Update WORKSPACE and install.bzl for the current location.
# Must be invoked from the src directory /path/to/LillyMol/src
pushd $REPO_HOME/src

if [[ ! -s 'WORKSPACE' ]] ; then
    echo "Must be invoked in the directory with WORKSPACE" && exit 1
fi

# Only build python if requested
if [[ -v BUILD_PYTHON ]] ; then
    # Use python to update WORKSPACE for python locations.
    if [[ -s 'update_python_in_workspace.py' ]] ; then
        cp WORKSPACE /tmp
        python3 ./update_python_in_workspace.py /tmp/WORKSPACE > WORKSPACE
            if [[ ! -s WORKSPACE ]] ; then
                echo "Updating WORKSPACE failed, restoring orignal, python bindings will not work"
                cp -f /tmp/WORKSPACE WORKSPACE
            fi
        echo "WORKSPACE updated"
    else
        echo "Missing update_python_in_workspace.py, WORKSPACE not updated for python"
    fi
fi

# install.bzl does need to be updated.
echo 'Updating build_deps/install.bzl'
if [[ ! -s 'build_deps/install.bzl' ]] ; then
    echo "build_deps/install.bzl not found" && exit 1
fi

bindir=$(echo $REPO_HOME/bin/$(uname) | sed -e 's/\//\\\//g')

# Make a copy
cp build_deps/install.bzl /tmp/install.bzl.orig

sed -i -e "s/default = \".*\",/default = \"${bindir}\",/" build_deps/install.bzl

# Create bindir if not already present
bindir=$REPO_HOME/bin/$(uname)
if [[ ! -d  ${bindir} ]] ; then
  mkdir -p ${bindir}
fi

# Step 2: build third party dependencies
# Download and built pre-requisites for LillyMol.
# Deliberately no robust error checking.

# Some third party dependencies an be built in parallel.
# Make sure the number of threads is the same as what you
# have available via quota on the machine being used.
# Things will likely go badly if you request 12 cores, but
# are restructed to 1.
THREADS=8

# We expect to be invoked from either the top level directory of
# the repo, or src. We want third party to be at the same level
# as src.
third_party=$REPO_HOME/third_party
echo "third_party in ${third_party}"

if [[ ! -d "${third_party}" ]] ; then
    mkdir -p "${third_party}"
fi

lib=$REPO_HOME/lib
echo "lib in ${lib}"

if [[ ! -d "${lib}" ]] ; then
    mkdir -p "${lib}"
fi


# If you are building with bazel, also need to update this path in WORKSPACE
# If you are building with cmake, also need to update this path in CMakeLists.txt

# The general stragegy here is that if the source directory does not exist, fetch it.
# Then, if some artifact of the build/install is absent, go into that directory, and re/build,
# and re/install.
pushd $third_party

if [[ ! -d ${third_party}/bin ]] ; then
    mkdir -p ${third_party}/bin
fi

# If we clone the repo we must build it, even if the
# file being checked is still present.
declare -i must_build

must_build=0
if [[ ! -s 'f2c.tar.gz' ]] ; then
    wget -O f2c.tar.gz https://www.netlib.org/f2c/src.tgz 
    tar -zxvf f2c.tar.gz
    mv src f2c  # change the non descriptive name
    must_build=1
fi
if [[ ${must_build} -eq 1 || ! -s 'f2c/cds.o' ]] ; then  # check for an arbitrary object file
    cd f2c && make -f makefile.u
    cd ..
fi

must_build=0
if [[ ! -s 'libf2c.zip' ]] ; then
    wget https://www.netlib.org/f2c/libf2c.zip
    mkdir libf2c
    cd libf2c && unzip ../libf2c.zip
    must_build=1
    cd ..
fi
if [[ ${must_build} -eq 1 || ! -s 'libf2c/libf2c.a' ]] ; then
    cd libf2c && make -f makefile.u
    cd ..
fi

# You should examine the BerkeleyDB license terms, it is not necessarily free.
if [[ -v BUILD_BDB ]] ; then
    must_build=0
    bdb_version='18.1.40'
    if [[ ! -s "db-${bdb_version}.tar.gz" ]] ; then
        wget http://download.oracle.com/berkeley-db/db-${bdb_version}.tar.gz
        tar zxf db-${bdb_version}.tar.gz
        must_build=1
    fi
    if [[ ${must_build} -eq 1 || ! -s "${third_party}/BDB/include/db.h" ]] ; then
        (cd "db-${bdb_version}/build_unix" && make realclean)
        (cd "db-${bdb_version}/build_unix" && ../dist/configure --prefix=${third_party}/BDB --enable-cxx --enable-shared=yes --with-repmgr-ssl=no)
        (cd "db-${bdb_version}/build_unix" && make -j 6)
        # fails on some items not installed, but installs what we need.
        (cd "db-${bdb_version}/build_unix" && make -j 4 install)

        echo ""
        echo "Ignore error messages from BerkeleyDB install, it is for components we do not use"
        # Copy shared libraries to our lib folder so python bindings work.
        (cp BDB/lib/lib*.dylib ${lib}) || echo "Did not copy BerkeleyDB shared libraries"
    fi
fi

# Step 3: biuld LillyMol executables
echo "Builds and installs LillyMol executables"
echo "The assumption is that WORKSPACE and build_deps/install.bzl"
echo "have both been configured."
echo ""
echo "First task is to build and run C++ unit tests."

# note that we do not check return codes from any of the invocations.
pushd $REPO_HOME/src
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
build_options="${build_options} --copt=-march=native --cxxopt=-march=native --cxxopt=-mtune=native"

# Seems like splitting out the BerkeleyDB components of the python
# bindings is needlessly complex at this stage. If python is being
# built, also build the BerkeleyDB bindings to avoid that complexity.
if [[ -v BUILD_PYTHON ]] ; then
    BUILD_BDB=1
fi

# All BerkeleyDB things are now isolated, so this probably makes no difference.
if [[ ! -v BUILD_BDB ]] ; then
    build_options="${build_options} --build_tag_filters=-berkeleydb"
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
    ./copy_shared_libraries.sh $REPO_HOME/lib
fi
