#!/usr/bin/env bash
# Download and built pre-requisites for LillyMol.
# Deliberately no robust error checking.

# Some third party dependencies an be built in parallel.
# Make sure the number of threads is the same as what you
# have available via quota on the machine being used.
# Things will likely go badly if you request 12 cores, but
# are restructed to 1.
THREADS=8

declare -i inside_lilly
if [[ $(hostname -d) =~ 'lilly.com' ]] ; then
  inside_lilly=1
else
  inside_lilly=0
fi

# We expect to be invoked from either the top level directory of
# the repo, or src. We want third party to be at the same level
# as src.
if [[ ${PWD} =~ LillyMol$ || ${PWD} =~ LillyMolPrivate$ ]] ; then
  third_party='third_party'
elif [[ ${PWD} =~ LillyMol.*/src$ ]] ; then
  third_party='../third_party'
else
  echo "Must be invoked from either LillyMol or LillyMol/src" >&1
  exit 1
fi
echo "third_party in ${third_party}"

# Convert to full path name
third_party=$(readlink -m ${third_party})
if [[ ! -d "${third_party}" ]] ; then
  mkdir -p "${third_party}"
fi

toplevel=$(dirname ${third_party})
if [[ ! -d "${toplevel}/lib" ]] ; then
  mkdir -p "${toplevel}/lib"
fi

uname=$(uname)
if [[ ! -d "${toplevel}/bin/${uname}" ]] ; then
  mkdir -p "${toplevel}/bin/${uname}"
fi

# If you are building with bazel, also need to update this path in WORKSPACE
# If you are building with cmake, also need to update this path in CMakeLists.txt

# The general stragegy here is that if the source directory does not exist, fetch it.
# Then, if some artifact of the build/install is absent, go into that directory, and re/build,
# and re/install.
cd "${third_party}" || echo "Cannot go to ${third_party}" >&2

# If we clone the repo we must build it, even if the
# file being checked is still present.
declare -i must_build
must_build=0

if [[ ! -d ${third_party}/bin ]] ; then
 mkdir -p ${third_party}/bin
fi

# Needed if building with cmake, and protoc is not installed

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
    (cp BDB/lib/lib*.so ${toplevel}/lib) || echo "Did not copy BerkeleyDB shared libraries"
  fi
fi

must_build=0
if [[ ! -s 'pybind11_protobuf/WORKSPACE' ]] ; then
  git clone https://github.com/pybind/pybind11_protobuf
else
  (cd pybind11_protobuf && git pull)
fi

# cilk??

# Not yet in the public release.
# if [[ ! -d "x86-simd-sort" ]] ; then
#   git clone 'https://github.com/intel/x86-simd-sort'
# fi

# Fastfloat get prebuilt single header
# Oct 2022. Not implemented because it seems
# that the existing conversion functions are faster.
# That seems impossible. TODO:ianwatson investigate.
# if [[ ! -d "${third_party}/include" ]] ; then
#   mkdir -p "${third_party}/include"
# fi
# (cd ${third_party}/include && wget https://github.com/fastfloat/fast_float/releases/download/v3.4.0/fast_float.h)

if [[ -v BUILD_GFP_SERVER ]] ; then
  if [[ ! -d libzmq ]] ; then
    git clone https://github.com/zeromq/libzmq
  else
    (cd libzmq && git pull)
  fi
  if [[ ! -d 'libzmq/cmake-build' ]] ; then
    mkdir 'libzmq/cmake-build'
  fi
  (cd libzmq/cmake-build && cmake -DCMAKE_INSTALL_PREFIX=${third_party} -DBUILD_SHARED=0 ..)
  (cd libzmq/cmake-build && make -j 4)
  (cd libzmq/cmake-build && make install)

  if [[ ! -d cppzmq ]] ; then
    git clone https://github.com/zeromq/cppzmq
    mkdir cppzmq/build
  else
    (cd cppzmq && git pull)
  fi
  if [[ ${inside_lilly} -eq 1 ]] ; then
    ${toplevel}/src/make_cppzmq_lly.sh ${third_party}
  else
    (cd cppzmq/build && cmake -DCMAKE_INSTALL_PREFIX=${third_party} ..)
    (cd cppzmq/build && make -j 4 && make install)
  fi
fi
