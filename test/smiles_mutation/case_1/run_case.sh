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

command=$BIN_DIR/smiles_mutation
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

name1=log.txt
name1_out=out/output.smi

diff_tool=../../fileDiff.sh

$command -N 10000 -n 20 -p 5 -c 15 -C 40 in/pubchem_example.smi >${name1} 2>err.txt

line_count=$(wc -l < "${name1}")
# echo $line_count
# 5000 is an arbitrary number. Beware this test is non deterministic
if [ ${line_count} -ge 5000 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm err.txt
