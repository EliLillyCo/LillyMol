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

command=$BIN_DIR/trxn
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

name1=output.smi
name1_out=out/output.smi
diff_tool=../../fileDiff.sh
$command -v -r in/1.2.1_Aldehyde_reductive_amination_FROM_amines_AND_aldehydes.rxn -Z -z i -M RMX -m RMX -S output in/20180412_amines.smi in/20180412_aldehydes.smi>>log.txt 2>>err.txt

$diff_tool $name1 $name1_out
ret1=$?
if [ $ret1 -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm log.txt
rm err.txt
