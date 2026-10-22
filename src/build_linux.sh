#!/usr/bin/env bash

# check if REPO_HOME is set from Makefile
# if directly calling this script, set it
if [[ ! -v REPO_HOME ]] ; then
    echo "REPO_HOME is not set, initializing..."
    REPO_HOME=$(dirname $(dirname  "$(realpath "$0")"))
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

# Must be invoked from the src directory /path/to/LillyMol/src
pushd $REPO_HOME/src

if [[ ! -s 'MODULE.bazel' ]] ; then
    echo "Must be invoked in the directory with MODULE.bazel" && exit 1
fi

# Homebrew openmp
if [[ $(uname) == 'Darwin' ]] ; then
  grep -q 'homebrew.*openmp' MODULE.bazel
  if [[ $? -ne 0 ]] ; then
    python3 ./add_homebrew_to_module_bazel.py >> MODULE.bazel
  else
    echo "Homebrew openmp already added to MODULE.bazel"
  fi
fi

# Build only the Python bindings and their shared-library dependencies.
# This is intended for container builds that do not need LillyMol executables.
if [[ -v PYTHON_ONLY ]] ; then
  BUILD_PYTHON=1
fi

# Seems like splitting out the BerkeleyDB components of the python
# bindings is needlessly complex at this stage. If python is being
# built, also build the BerkeleyDB bindings to avoid that complexity.

# This must be set before the third_party section so BerkeleyDB can be built.
if [[ -v BUILD_PYTHON ]] ; then
    BUILD_BDB=1
fi

# If building only the python bindings, do not enable BDB - that is just because
# I know that the people using this do not need the BerkeleyDB bindings.
# There is no technical reason for the choice below, change as needed.

if [[ -v PYTHON_ONLY ]] ; then
  unset BUILD_BDB
fi

# Create bindir if not already present
bindir=$REPO_HOME/bin/$(uname)
if [[ ! -d ${bindir} ]] ; then
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

# Shared library components for python go here.
lib=$REPO_HOME/lib
echo "lib in ${lib}"

if [[ ! -d "${lib}" ]] ; then
  mkdir -p "${lib}"
fi

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
    (cd f2c && make -f makefile.u)
fi

must_build=0
if [[ ! -s 'libf2c.zip' ]] ; then
    wget https://www.netlib.org/f2c/libf2c.zip
    mkdir libf2c
    (cd libf2c && unzip ../libf2c.zip)
    must_build=1
fi
if [[ ${must_build} -eq 1 || ! -s 'libf2c/libf2c.a' ]] ; then
    (cd libf2c && make -f makefile.u)
fi

# Planarity support proved problematic in many ways. Originally made part of
# iwdescr and molecule_filter, but it does not really perceive the kinds of
# things we think of as non-planar. Not worth the complexity.

#if [[ ! -d 'dragonbox' ]] ; then
#  git clone https://github.com/jk-jeon/dragonbox
#  # We fix a specific version because third_party/BUILD.bazel needs the full path
#  # of the inlcude directory.
#  # We could do something clever with links...
#  (cd dragonbox && git -c advice.detachedHead=false checkout 1.1.3)
#  (cd dragonbox && mkdir build)
#  (cd dragonbox/build && cmake -DCMAKE_INSTALL_PREFIX=${third_party} ..)
#  (cd dragonbox/build && make install)
#fi

# You should examine the BerkeleyDB license terms, it is not necessarily free.
if [[ -v BUILD_BDB ]] ; then
    must_build=0
    bdb_version='18.1.40'
    tarfile="db-${bdb_version}.tar.gz"
    if [[ ! -s "${tarfile}" ]] ; then
        wget http://download.oracle.com/berkeley-db/${tarfile}
        tar zxf ${tarfile}
        must_build=1
    fi
    if [[ $(uname) == "Linux" ]] ; then
      suffix='so'
    elif [[ $(uname) == "Darwin" ]] ; then
      suffix='dylib'
    fi
    if [[ ${must_build} -eq 1 || ! -s "${third_party}/BDB/include/db.h" ]] ; then
        (cd "db-${bdb_version}/build_unix" && make realclean)
        (cd "db-${bdb_version}/build_unix" && ../dist/configure --prefix=${third_party}/BDB --enable-cxx --enable-shared=yes --with-repmgr-ssl=no)
        (cd "db-${bdb_version}/build_unix" && make -j ${THREADS})
        # fails on some items not installed, but installs what we need.
        (cd "db-${bdb_version}/build_unix" && make -j ${THREADS} install)

        echo ""
        echo "Ignore error messages from BerkeleyDB install, it is for components we do not use"
        # Copy shared libraries to our lib folder so python bindings work.
        (cp BDB/lib/lib*.${suffix} ${lib}) || echo "Did not copy BerkeleyDB shared libraries"
    fi
fi

if [[ -v BUILD_XGBOOST ]] ; then
  if [[ ! -d 'xgboost' ]] ; then
    git clone --recursive https://github.com/dmlc/xgboost
    mkdir xgboost/build
    (cd xgboost/build && cmake -DCMAKE_INSTALL_PREFIX=${third_party} -DCMAKE_INSTALL_LIBDIR=${third_party}/lib ..)
    (cd xgboost/build && make -j${THREADS})
    (cd xgboost/build && make install)
  fi
fi

if [[ -v BUILD_ZEROMQ ]] ; then
  git clone https://github.com/zeromq/libzmq
  (cd libzmq && mkdir build)
  (cd libzmq/build && cmake -DCMAKE_INSTALL_PREFIX=${third_party} -DZMQ_BUILD_TESTS=OFF ..)
  (cd libzmq/build && make -j${THREADS})
  (cd libzmq/build && make install)

  git clone https://github.com/zeromq/cppzmq
  # could not get it to build, but seems like a header only library. This worked.
  cp cppzmq/zmq.hpp ${third_party}/include
fi

if [[ -v BUILD_INCHI ]] ; then
  if [[ ! -d 'InChI' ]] ; then
    git clone https://github.com/IUPAC-InChI/InChI
    (cd InChI && mkdir build)
    (cd InChI/build && cmake -DCMAKE_INSTALL_PREFIX=${third_party} -DCMAKE_BUILD_TYPE=Release ../INCHI-1-SRC/INCHI_API/libinchi/src)
    (cd InChI/build && cmake --build .)
    cp InChI/build/lib/libinchi*${suffix} ${REPO_HOME}/lib
    # The repo does not have an install target. third_party/BUILD.bazel links directly to the files.
    # Ensure canonical form for shell variable.
    BUILD_INCHI=1
  fi
fi

if [[ -v BUILD_NLOPT ]] ; then
  git clone https://github.com/stevengj/nlopt/
  (cd nlopt && mkdir build)
  (cd nlopt/build && cmake -DBUILD_SHARED_LIBS=NO -DCMAKE_INSTALL_PREFIX=${REPO_HOME}/third_party)
  (cd nlopt/build && make -j ${THREADS})
  (cd nlopt/build && make install)
  BUILD_NLOPT=1
fi

if [[ -v BUILD_LIGHTGBM ]] ; then
  git clone --recursive https://github.com/microsoft/LightGBM
  (cd LightGBM && cmake -B build -DBUILD_STATIC_LIB=ON -DCMAKE_INSTALL_PREFIX=${REPO_HOME}/third_party/ -S .)
  (cd LightGBM && cmake --build build -j 4)
  (cd LightGBM && cmake --install build --prefix=${REPO_HOME}/third_party)
  BUILD_LIGHTGBM=1
fi

# Optional 2D coordinate generation support. Build a static library so future
# LillyMol targets can link it without adding a runtime shared-library dependency.
if [[ -v BUILD_SCHRODINGER_2D ]] ; then
  must_build=0
  if [[ ! -d 'coordgenlibs' ]] ; then
    git clone https://github.com/schrodinger/coordgenlibs.git
    must_build=1
  fi
  if [[ ${must_build} -eq 1 || ! -s "${third_party}/lib/libcoordgen.a" ]] ; then
    cmake -S coordgenlibs -B coordgenlibs/build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=${third_party} \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
      -DCOORDGEN_BUILD_EXAMPLE=OFF \
      -DCOORDGEN_BUILD_SHARED_LIBS=OFF \
      -DCOORDGEN_BUILD_TESTS=OFF \
      -DCOORDGEN_RIGOROUS_BUILD=OFF \
      -DCOORDGEN_USE_MAEPARSER=OFF
    cmake --build coordgenlibs/build --parallel ${THREADS}
    cmake --install coordgenlibs/build
  fi
fi

# Step 3: build LillyMol executables
echo "Builds and installs LillyMol executables"
echo ""
echo "First task is to build and run C++ unit tests."

# note that we do not check return codes from any of the invocations.
pushd $REPO_HOME/src

declare -i inside_lilly
if [[ $(hostname -d) =~ 'lilly.com' ]] ; then
    inside_lilly=1
else
    inside_lilly=0
fi

# The script returns --output_user_root=<dir> if a suitable
# directory can be found
bazel_options=$(./output_user_root.sh)

build_options="--cxxopt=-DGIT_HASH=\"$(git rev-parse --short --verify HEAD)\" --cxxopt=-DTODAY=\"$(date +%Y-%b-%d)\" --jobs=${THREADS} -c opt"

# Enable partial builds.
build_options="${build_options} -k"

# Compile for the build machine by default. This has historically made LillyMol
# about 5% faster. Set BUILD_NATIVE=0 if the resulting binaries must run across
# multiple CPU generations or otherwise need a generic architecture.
if [[ ! -v BUILD_NATIVE ]] ; then
  BUILD_NATIVE=1
fi

if [[ ${BUILD_NATIVE} -ne 0 ]] ; then
  build_options="${build_options} --config=native"
fi

if [[ -v BUILD_INCHI ]] ; then
  # Note that this config also sets the variable BUILD_INCHI 
  build_options+=' --config=inchi'
fi

# For things not being built, assemble an array of tag filters.
declare -a build_tag_filters

# All BerkeleyDB things are now isolated, so this probably makes no difference.
if [[ ! -v BUILD_BDB ]] ; then
    build_tag_filters+=('-berkeleydb')
fi

# Building the zeromq based gfp_server is optional and generally
# needs to be configured manually.
if [[ ! -v BUILD_ZEROMQ ]] ; then
    build_tag_filters+=('-gfp_server' '-zeromq')
fi

# Acknowledge that this is not the best way to handle something like this.
if [[ -v BUILD_INCHI ]] ; then
  build_options="${build_options} --define 'BUILD_INCHI=1'"
fi

if [[ ! -v BUILD_XGBOOST ]] ; then
    build_tag_filters+=('-xgboost')
fi

if [[ ! -v BUILD_NLOPT ]] ; then
    build_tag_filters+=('-nlopt')
fi

if [[ ${inside_lilly} -eq 0 ]] ; then
  build_tag_filters+=('-vendor')
fi

if [[ "${#build_tag_filters[@]}" -gt 0 ]] ; then
  tmp=$(IFS=,;printf '%s' "${build_tag_filters[*]}")
  build_options="${build_options} --build_tag_filters=${tmp}"
fi

echo "build_options ${build_options}"

if [[ -v PYTHON_ONLY ]] ; then
    echo "Building Python bindings only"
    python_only_targets=(//nanobind:lillymol //nanobind:lillymol_gfp_server //nanobind:lillymol_soabi)
    if [[ -v BUILD_SCHRODINGER_2D ]] ; then
        python_only_targets+=(//nanobind:lillymol_depict)
    fi
    ${bazel} ${bazel_options} build ${build_options} "${python_only_targets[@]}"
    ./copy_shared_libraries.sh $REPO_HOME/lib
    if [[ -v RUN_PYTHON_TESTS ]] ; then
        ./run_python_unit_tests.sh
    fi
    popd
    popd
    exit 0
fi

# First task is unit tests

set -x

${bazel} ${bazel_options} test ${build_options} Foundational/...:all

${bazel} ${bazel_options} test ${build_options} Molecule_Lib:all
${bazel} ${bazel_options} test ${build_options} Molecule_Tools:all

# Depiction links the coordgen library, which is only built above when asked for.
if [[ -v BUILD_SCHRODINGER_2D ]] ; then
    ${bazel} ${bazel_options} test ${build_options} Depict:all
fi

${bazel} ${bazel_options} test ${build_options} Utilities/...:all

# Currently no tests in these.
if [[ -v BUILD_BDB ]] ; then
    # ${bazel} ${bazel_options} test ${build_options} BerkeleyDB:all
    # ${bazel} ${bazel_options} test ${build_options} Molecule_Tools_Bdb:all
    echo ""
fi

# Once the tests run, then executables can be built.

if [[ -v BUILD_INCHI ]] ; then
  tmp_build_options="${build_options} --config=inchi"
  ${bazel} ${bazel_options} build ${tmp_build_options} Molecule_Tools:fileconv
  ${bazel} ${bazel_options} build ${tmp_build_options} Molecule_Tools:unique_molecules
  echo "Inchi BUILD of fileconv and unique molecules, executables in bazel-bin/Molecule_Tools"
  echo "Copy these files and then rebuild everything else without BUILD_INCHI"
  exit 0
fi

if [[ ! -v BUILD_LIBRARY_ONLY ]] ; then
    echo "Building tools"
    ${bazel} ${bazel_options} build ${build_options} Molecule_Tools:all
    if [[ -v BUILD_SCHRODINGER_2D ]] ; then
        ${bazel} ${bazel_options} build ${build_options} Depict:all
    fi
    ${bazel} ${bazel_options} build ${build_options} Obsolete:all
    ${bazel} ${bazel_options} build ${build_options} Obsolete/Descriptor_Similarity:all
    ${bazel} ${bazel_options} build ${build_options} Foundational/iw_tdt:all
    ${bazel} ${bazel_options} build ${build_options} Utilities/...:all
fi  

# Build the Vendor directory only when explicitly requested or when inside Lilly.
# Executables that link to vendor supplied libraries are tagged with "vendor" and
# those will not be built outside Lilly unless BUILD_VENDOR is set. What will be
# built are pure LillyMol tools that parse the output from Vendor tools.
if [[ ${inside_lilly} -eq 1 || -v BUILD_VENDOR ]] ; then
    echo "Building Vendor target"
    ${bazel} ${bazel_options} build ${build_options} Vendor/...:all
else
  echo "Vendor target not built"
fi

if [[ -v BUILD_BDB ]] ; then
    ${bazel} ${bazel_options} build ${build_options} BerkeleyDB:all
    ${bazel} ${bazel_options} build ${build_options} Molecule_Tools_Bdb:all
fi

if [[ -v BUILD_GO ]] ; then
    ${bazel} ${bazel_options} build ${build_options} go:all
fi

if [[ -v BUILD_XGBOOST ]] ; then
    ${bazel} ${bazel_options} build ${build_options} xgboost:all
fi

# Now install the targets

./install.sh

# uninstall executables not currently being used within Lilly.
if [[ ${inside_lilly} -eq 1 ]] ; then
    ./uninstall.sh
fi

#if [[ -v BUILD_BDB ]] ; then
#    ${bazel} ${bazel_options} run ${build_options} BerkeleyDB:install
#    ${bazel} ${bazel_options} run ${build_options} Molecule_Tools_Bdb:install
#fi

#if [[ -v BUILD_XGBOOST ]] ; then
#    ${bazel} ${bazel_options} run ${build_options} xgboost:install
#fi

# Deal with compressed files and hidden file names
LILLYMOL_HOME=${REPO_HOME} ./uncompress_and_install.sh

# The ruby bindings used by test/run_all_test.rb are needed on every platform,
# not just Macos. Generating them here keeps them matched to the local protoc.
# Note the form of the generated file matters: protoc older than about 3.25
# emits a DescriptorPool DSL that the google-protobuf 4.x gem no longer
# supports, so a file generated long ago will fail to load.
(cd ${REPO_HOME}/test && protoc --ruby_out=. --proto_path=. lillymol_tests.proto)

# If Macos compile some protos.
# There are a great many other protos that could be compiled for the platform.
# bazel query 'kind("proto_library", ...:all)'

if [[ $(uname) == 'Darwin' ]] ; then
  if [[ -v BUILD_XGBOOST ]] ; then
    (cd ${REPO_HOME}/src && protoc --python_out=. xgboost/xgboost_model.proto && cp xgboost/xgboost_model_pb2.py ${REPO_HOME}/contrib/bin/xgbd )
    (cd ${REPO_HOME}/src && protoc --python_out=. xgboost/random_forest_model.proto && cp xgboost/random_forest_model_pb2.py ${REPO_HOME}/contrib/bin/xgbd )
    (cd ${REPO_HOME}/src && protoc --python_out=. Utilities/General/linear_scaling.proto && cp Utilities/General/linear_scaling_pb2.py ${REPO_HOME}/contrib/bin/xgbd )
  fi
fi

# Python if requested, build, install and test.
# Note that PYTHONPATH will need to be adjusted, or copy the shared
# libraries from LillyMol/lib to your default PYTHONPATH.

if [[ -v BUILD_PYTHON ]] ; then
    python_targets=(//nanobind:lillymol //nanobind:lillymol_gfp_server //nanobind:lillymol_soabi)
    if [[ -v BUILD_BDB ]] ; then
        python_targets+=(//nanobind:lillymol_bdb)
    fi
    if [[ -v BUILD_SCHRODINGER_2D ]] ; then
        python_targets+=(//nanobind:lillymol_depict)
    fi
    ${bazel} ${bazel_options} build ${build_options} "${python_targets[@]}"
    ${bazel} ${bazel_options} test ${build_options} //nanobind:all
    if [[ -v BUILD_SCHRODINGER_2D ]] ; then
        ${bazel} ${bazel_options} test ${build_options} //nanobind:lillymol_nb_depict_test
    fi
    ./copy_shared_libraries.sh $REPO_HOME/lib 

    ./run_python_unit_tests.sh
fi
