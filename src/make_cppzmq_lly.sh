#!/usr/bin/env bash

# This is called by build_third_party if inside Lilly.

# Must be called from the third_party directory.

# cppzmq cannot be built with the default Lilly environment,
# it picks up zeromq instances out of the python virtual environment.
# Set up our own environment to better control things.


export LD_LIBRARY_PATH=$PWD/lib64:/opt/rh/devtoolset-12/root/usr/lib64:/opt/rh/devtoolset-12/root/usr/lib:/usr/lib:/usr/X11R6/lib:/usr/local/lib:/usr/local/lib64:
export PATH=/lrlhps/apps/bazelisk/bazelisk-1.13.2:/opt/rh/devtoolset-12/root/usr/bin:/opt/uge/uge-2022.1.2/bin/lx-amd64:/bin:/usr/bin:/usr/sbin:/usr/local/sbin:
export PKG_CONFIG_PATH=$PWD/lib64/pkgconfig:PKG_CONFIG_PATH

module load cmake/3.28.1

cd cppzmq
if [[ ! -d 'build' ]] ; then
  mkdir build
fi
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PWD} ..
make -j 4
make install
