#!/bin/bash

if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]; then 
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME (local path to LillyMol code)"
    echo "Please export BUILD_DIR (the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
else
    BIN_DIR="$LILLYMOL_HOME/bin/$BUILD_DIR"
fi

test_command=molecules_from_reagents
case=case_1
case_id="Case 1"

test_top="$LILLYMOL_HOME/test"
test_cmd_top="$test_top/$test_command"

diff_tool=../../fileDiff.sh
command="$BIN_DIR/$test_command"

if [ ! -x "$command" ]; then
    echo "$command is not executable or does not exist"
    exit 1
fi

in="$test_cmd_top/$case/in/seeds.smi"
out=seeds_found.smi
cmp_out="$test_cmd_top/$case/out/seeds_found.smi"

echo "Testing: $command"

queries_dir="$LILLYMOL_HOME/contrib/data/queries"

$command -M NN -v -R $test_cmd_top/$case/in/set1.smi \
    -R $test_cmd_top/$case/in/set2.smi "$in" 1>> "$out" 2>> err.txt
$diff_tool "$out" "$cmp_out"
ret=$?

if [ $ret -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi


rm -f "$out"
rm -f err.txt
rm -f NN0.smi
rm -f NN1.smi
