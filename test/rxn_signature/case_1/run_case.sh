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

command=$BIN_DIR/rxn_signature
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

name1=all.sig
name1_out=out/all.sig
name2=Cfile
name2_out=out/Cfile
name3=Ffile
name3_out=out/Ffile
diff_tool=../../fileDiff.sh
$command -v -r 0,1,2 -C Cfile -F Ffile in/all.rxnsmi >all.sig 2>all.log

$diff_tool $name1 $name1_out
ret1=$?
$diff_tool $name2 $name2_out
ret2=$?
$diff_tool $name3 $name3_out
ret3=$?

if [ $ret1 -eq 1 ] && [ $ret2 -eq 1 ] && [ $ret3 -eq 1 ] 
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm $name2
rm $name3
rm all.log
