#!/bin/bash

if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]; then 
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME (local path to LillyMol code)"
    echo "Please export BUILD_DIR (the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/path/to/LillyMolPrivate"
    echo "Example: export BUILD_DIR=Linux" 
    exit 1
else
    BIN_DIR="$LILLYMOL_HOME/bin/$BUILD_DIR"
fi

test_command=iwdescr
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

in="$test_cmd_top/$case/in/test.smi"
out=test.out
cmp_out="$test_cmd_top/$case/out/test.out"

# Support linux and mac 
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    cmp_out=out/linux/test.out
elif [[ "$OSTYPE" == "darwin"* ]]; then
    cmp_out=out/osx/test.out
else
    echo "OS is not supported"
fi


echo "Testing: $command"

queries_dir="$LILLYMOL_HOME/data/queries"

$command \
   -N F:${queries_dir}/charges/queries \
   -H a=F:${queries_dir}/hbonds/acceptor \
   -H d=${queries_dir}/hbonds/donor.qry -H label \
   -l -g all -A D -u 0 -b 5 -O complex -O dm -B quiet -E autocreate \
   "$in" > "$out" 2>> err.txt
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
