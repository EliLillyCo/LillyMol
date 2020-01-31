#!/bin/bash

if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]; then 
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME (local path to LillyMol code)"
    echo "Please export BUILD_DIR (the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
else
    BIN_DIR="$LILLYMOL_HOME/contrib/script/py/mmp/"
fi

test_command=getMMPfromSMI
case=case_4
case_id="Case 4"

test_top="$LILLYMOL_HOME/test"
test_cmd_top="$test_top/$test_command"

diff_tool=../../fileDiff.sh

command="$BIN_DIR/$test_command.py"

if [ ! -x "$command" ]; then
    echo "$command is not executable or does not exist"
    exit 1
fi

in="$test_cmd_top/$case/in/test_data_02.smi"
out=test_data_04.pairs
gold_out="$test_cmd_top/$case/out/test_data_04.pairs"

echo "Testing: $command"

$command -i "$in" -o "$out" -c DOUBLE -t 2>>err.txt
$diff_tool "$out" "$gold_out"
ret=$?

if [ $ret -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi

rm -f "$out"
rm -f err.txt
