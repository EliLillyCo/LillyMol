#! /bin/bash

if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]
then 
    # undefined BIN_DIR
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME(local path to LillyMol code)"
    echo "Please export BUILD_DIR(the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
else
    BIN_DIR=$LILLYMOL_HOME/bin/$BUILD_DIR
fi

command=$BIN_DIR/normalise
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

name1=stdout
name1_out=out/${name1}

# Support linux and mac 
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    name1_out=out/linux/${name1}
elif [[ "$OSTYPE" == "darwin"* ]]; then
    name1_out=out/osx/${name1}
else
    echo "OS is not supported"
fi

diff_tool=../../fileDiff.sh
$command -i 1 -j $LILLYMOL_HOME/data/descriptor.txt >${name1} 2>err.txt
$diff_tool ${name1} ${name1_out}
ret=$?
if [ $ret -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm err.txt
