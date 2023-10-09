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

command=$BIN_DIR/unique_molecules
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable ${command} is not found"
    exit 1
fi

name1=unique.smi
name1_out=out/unique.smi
name2=duplicate.smi
name2_out=out/duplicate.smi

diff_tool=../../fileDiff.sh
$command -S unique -D duplicate -v -l in/input.smi>>log.txt 2>>err.txt

$diff_tool $name1 $name1_out
ret1=$?

$diff_tool $name2 $name2_out
ret2=$?

if [ $ret1 -eq 1 ] && [ $ret2 -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm duplicate.smi
rm log.txt
rm err.txt
