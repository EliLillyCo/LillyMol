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

command=$BIN_DIR/rxn_standardize
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

name1=output.rxnsmi
name1_out=out/output.rxnsmi
diff_tool=../../fileDiff.sh
#$command -s -c -D x -X igbad -v -C 60 -K -E autocreate -e -o -m -I -b -f gsub in/input.rxnsmi >>output.rxnsmi 2>>err.txt
$command -s -c -X rmdmap -X igbad -v -C  60 -K -E autocreate -o -X rmmr -X  rmiso -D rmdr -f gsub in/input.rxnsmi >>output.rxnsmi 2>>err.txt

$diff_tool $name1 $name1_out
ret1=$?
if [ $ret1 -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi

rm $name1
rm err.txt
