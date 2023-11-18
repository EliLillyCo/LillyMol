#!/bin/bash

# Update WORKSPACE and install.bzl for the current location.
# Must be invoked from the src directory /path/to/LillyMol/src

if [[ ! -s 'WORKSPACE' ]] ; then
  echo "Must be invoked in the directory with WORKSPACE"
  exit 1
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
  echo "build_deps/install.bzl not found"
  exit 1
fi

bindir=$(echo ${PWD}/../bin/$(uname) | sed -e 's/\//\\\//g')

# Make a copy
tmpinstall='/tmp/install.bzl'
cp build_deps/install.bzl /tmp/install.bzl.orig

sed --in-place --regexp-extended -e "s/default = *\"..+\",/default =\"${bindir}\",/" build_deps/install.bzl > ${tmpinstall}

# Create bindir if not already present
bindir=$(echo ${PWD}/../bin/$(uname))
if [[ ! -d  ${bindir} ]] ; then
  mkdir -p ${bindir}
fi
