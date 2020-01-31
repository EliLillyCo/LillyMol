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
case=case_1
case_id="Case 1"

test_top="$LILLYMOL_HOME/test"
test_cmd_top="$test_top/$test_command"

diff_tool=../../fileDiff.sh

command="$BIN_DIR/$test_command.py"

if [ ! -x "$command" ]; then
    echo "$command is not executable or does not exist"
    exit 1
fi

in="$test_cmd_top/$case/in/test_data_01.smi"
out=test_data_01.pairs
gold_err="$test_cmd_top/$case/out/err.txt"

echo "Testing: $command"

$command -i "$in" -o "$out" 2>>err.txt
# take first line only to avoid path
head -n 1 err.txt > err.head
$diff_tool err.head $gold_err
ret=$?

# This is a negative pass fail - we expect no output
if [ $ret -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi

rm -f $out
rm -f err.txt
rm -f err.head
