#!/bin/bash
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

# If you are building with bazel, also need to update this path in WORKSPACE
# If you are building with cmake, also need to update this path in CMakeLists.txt

# We want to update WORKSPACE if the build succeeds, so save our current location
maybe_src=${PWD}

# The general stragegy here is that if the source directory does not exist, fetch it.
# Then, if some artifact of the build/install is absent, go into that directory, and re/build,
# and re/install.
cd "${third_party}" || echo "Cannot go to ${third_party}" >&2

# If we clone the repo we must build it, even if the
# file being checked is still present.
declare -i must_buld
must_buld=0

# Jun 2023 crc32c now obtained from absl
# if [[ ! -d 'crc32c' ]] ; then
#   git clone https://github.com/google/crc32c
#   must_build=1
# fi
# if [[ $must_build -eq 1 || ! -s "${third_party}/include/crc32c/crc32c.h" ]] ; then
#   (cd crc32c && git pull && git submodule update --init --recursive)
#   (cd crc32c && if [[ ! -d build ]] ; then mkdir build ; fi)
#   (cd crc32c/build && cmake -DCMAKE_INSTALL_PREFIX:PATH=${third_party} -DBUILD_SHARED_LIBS=no -DCRC32C_BUILD_TESTS=0 -DCRC32C_BUILD_BENCHMARKS=0 .. && make all install)
# fi

must_build=0
# re2 now from bazel MODULE
#if [[ ! -d 're2' ]] ; then
#  git clone https://github.com/google/re2
#  must_build=1
#fi
#if [[ ${must_build} -eq 1 || ! -s "${third_party}/RE2/include/re2/re2.h" ]] ; then
#  (cd re2 && git pull)
#  (cd re2 && make prefix=${third_party}/RE2)
#  (cd re2 && make prefix=${third_party}/RE2 install)
#fi

# Googletest not handled via MODULE
# if [[ ! -d 'googletest' ]] ; then
#   git clone https://github.com/google/googletest
# fi
# if [[ ! -f "${third_party}/include/gtest/gtest.h" ]] ; then
#   (cd googletest && git pull)
#   (cd googletest && mkdir build)
#   (cd googletest && cd build && cmake -DCMAKE_INSTALL_PREFIX:PATH=${third_party}/GOOGLETEST ..)
#   (cd googletest && cd build && make install)
# fi

# Needed if building with cmake, and protoc is not installed
if [[ ! -d ${third_party}/bin ]] ; then
  mkdir -p ${third_party}/bin
fi

must_build=0
if [[ ! -d 'protobuf' ]] ; then
  git clone https://github.com/protocolbuffers/protobuf
  must_build=1
fi
if [[ ${must_build} -eq 1 || ! -d "${third_party}/include/google/protobuf" ]] ; then
  (cd protobuf && git pull)
  (cd protobuf && git submodule update --init --recursive)
  # if bazel is used, the libraries don't get generated, not sure why...
  (cd protobuf && cmake . -DCMAKE_CXX_STANDARD=14)
  (cd protobuf && cmake --build .)

fi

# oneTBB now available via MODULE
# must_build=0
# if [[ ! -s 'oneTBB' ]] ; then
#   git clone https://github.com/oneapi-src/oneTBB
#   must_build=1
# fi
# if [[ ${must_build} -eq 1 || ! -s "${third_party}/TBB//lib/libtbb.so" ]] ; then
#   (cd oneTBB && git pull)
#   (cd oneTBB && mkdir build)
#   (cd oneTBB && cd build && cmake -DCMAKE_INSTALL_PREFIX:PATH=${third_party}/TBB -DTBB_TEST=OFF -DBUILD_SHARED_LIBS=off ..)
#   (cd oneTBB && cd build && make install)
# fi
# 
# Eigen now available via MODULE
# if [[ ! -s 'eigen-3.4.0.tar.bz2' ]] ; then
#   wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
#   tar -xvjf eigen-3.4.0.tar.bz2
# fi

must_build=0
if [[ ! -d 'highwayhash' ]] ; then
  git clone https://github.com/google/highwayhash
  must_build=1
fi
if [[ ${must_build} -eq 1 || ! -s '${third_party}/highwayhash/lib/libhighwayhash.a' ]] ; then
  (cd highwayhash && git pull)
  (cd highwayhash && make)
fi

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
must_build=0
bdb_version='18.1.40'
if [[ ! -s "db-${bdb_version}.tar.gz" ]] ; then
  wget http://download.oracle.com/berkeley-db/db-${bdb_version}.tar.gz
  tar zxvf db-${bdb_version}.tar.gz
  must_build=1
fi
if [[ ${must_build} -eq 1 || ! -s "${third_party}/BDB/include/db.h" ]] ; then
  (cd "db-${bdb_version}/build_unix" && make realclean)
  (cd "db-${bdb_version}/build_unix" && ../dist/configure --prefix=${third_party}/BDB --enable-cxx --enable-shared=no --with-repmgr-ssl=no)
  (cd "db-${bdb_version}/build_unix" && make)
  # fails on some items not installed, but installs what we need.
  (cd "db-${bdb_version}/build_unix" && make install)
fi

# cilk??

# Not yet in the public release.
if [[ ! -d "x86-simd-sort" ]] ; then
  git clone 'https://github.com/intel/x86-simd-sort'
fi
# Figure out of WORKSPACE is where we can find it.
if [[ -s 'WORKSPACE' ]] ; then
  workspace='WORKSPACE'
elif [[ -s 'src/WORKSPACE' ]] ; then
  workspace='src/WORKSPACE'
else
  echo "No WORKSPACE file found" >&2
  exit 1
fi

# Fastfloat get prebuilt single header
# Oct 2022. Not implemented because it seems
# that the existing conversion functions are faster.
# That seems impossible. TODO:ianwatson investigate.
set -x
if [[ ! -d "${third_party}/include" ]] ; then
  mkdir -p "${third_party}/include"
fi
if [[ ! -s "${third_party}/include/fast_float.h" ]] ; then
  (cd ${third_party}/include && wget https://github.com/fastfloat/fast_float/releases/download/v3.4.0/fast_float.h)
fi

# If you do not want to install MPI dependent tools, omit this AND comment out the mpich
# section in WORKSPACE. Run the build with --build-tag_filters=-mpi which suppresses
# building of tools that depend on MPI.

# Indeed this takes too long to build, and then fails because it needs some
# further system dependencies (Ubuntu 22). The best way to use mpi might be to install it on
# the system if you can. Only one tool depends on it.

# must_build=0
# mpi_version='4.1.2'
# if [[ ! -d ${third_party}/mpich-${mpi_version}/ ]] ; then
#   (cd ${third_party} && wget https://www.mpich.org/static/downloads/4.1.2/mpich-${mpi_version}.tar.gz)
#   (cd ${third_party} && tar zxvf mpich-${mpi_version}.tar.gz)
#   must_build=1
# fi

# if [[ $must_build -eq 1 || ! -s ${third_party}/mpich/src/lib ]] ; then
  (cd ${third_party}/mpich-${mpi_version}/ && ./configure --prefix=${third_party})
#   (cd ${third_party}/mpich-${mpi_version}/ && make)
# fi

# Make an attempt to generate an updated WORKSPACE file.
# First return to our starting point, which hopefully has WORKSPACE

cd ${maybe_src}
echo "Back to ${PWD}"


# OMG there seems to be a bug in sed, and this does not work!
# If I change the third_party to third_partq it works. Amazing. 
# This remains broken until sed works properly. Version 4.8.1 does
# not work, nor does 4.1.5.
# Seems hard to imagine that sed has a bug like this, maybe I am mistaken??!!

# Quote and backslash hell...
third_party=$(echo ${third_party} | sed -e 's/\//\\\//g')
echo "third_party ${third_party}"
tmpworkspace='/tmp/WORKSPACE'
sed --regexp-extended -e "s/path = \"..*\/(..+)\"/path = \"${third_party}\/\\1\"/" ${workspace} > ${tmpworkspace}
echo "Please check ${tmpworkspace} for a possibly ready to use WORKSPACE file" >&2

if [[ -s 'build_deps/install.bzl' ]] ; then
  tmpinstall='/tmp/install.bzl'
  bindir=$(echo ${PWD}/bin/Linux | sed -e 's/\//\\\//g')
  sed -e "s/default = \"..+\"/default =\"${bindir}\"/" build_deps/install.bzl > ${tmpinstall}
  echo "Please check ${tmpinstall} for a possibly ready to use build_deps/install.bzl file" >&2
fi
