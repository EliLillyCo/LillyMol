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

shared_data_dir=$LILLYMOL_HOME/data
command=$BIN_DIR/descriptor_file_to_01_fingerprints
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

name1=out.gfp
name1_out=out/out.gfp

diff_tool=../../fileDiff.sh
$command -F NCFP -S in/pubchem.smi $shared_data_dir/cleaned_descriptor.txt >out.gfp 2>err.txt
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
