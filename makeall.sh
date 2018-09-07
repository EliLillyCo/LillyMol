#!/bin/bash
# $Id$
# 
# Default compilier
export UNAME=Linux-gcc-7.2.1
echo "Start script for build $UNAME"
echo "Clean gcc"
IWPROGRAMMES=${PWD} make veryclean
echo "Copy Include file"
IWPROGRAMMES=${PWD} make copy_include
echo "Building libraries"
IWPROGRAMMES=${PWD} make library
echo "Copying libraries"
IWPROGRAMMES=${PWD} make copy_library
echo "Building executables"
IWPROGRAMMES=${PWD} make exe
echo "Copying executables"
IWPROGRAMMES=${PWD} make copy_exe
echo "End of the build"

echo "All done"
