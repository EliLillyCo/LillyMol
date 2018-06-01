# Copyright 2018 Eli Lilly and Company 
# 
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. 
# You may obtain a copy of the License at  
# 
#     http://www.apache.org/licenses/LICENSE-2.0  
# 
# Unless required by applicable law or agreed to in writing, software 
# distributed under the License is distributed on an "AS IS" BASIS, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and 
# limitations under the License. 
########################################################################

#!/bin/bash
# $Id$
# Build small molecule libraries and executables 
# Some are built with gcc/gfortran
# 

# Required library and package
# 1. GCC >= 6.2.0
# module unload gcc/4.9.1
# module load gcc/6.2.0
# 2. zlib >= 1.2.11 
# module load zlib/1.2.11

echo "Start script for build"
echo "Clean gcc"
make veryclean
echo "Copy Include file"
make copy_include
echo "Building libraries"
make library
echo "Copying libraries"
make copy_library
echo "Building executables"
make exe
echo "Copying executables"
make copy_exe
echo "End of the build"

echo "All done"
