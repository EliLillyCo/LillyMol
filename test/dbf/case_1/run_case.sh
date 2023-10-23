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

test_command=dbf
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

in="$test_cmd_top/$case/in/test.sdf"
out=test.out
cmp_out="$test_cmd_top/$case/out/test.out"

echo "Testing: $command"

queries_dir="$LILLYMOL_HOME/data/queries"

$command \
   -N F:${queries_dir}/charges/queries \
   -H a=F:${queries_dir}/hbonds/acceptor \
   -H d=${queries_dir}/hbonds/donor.qry \
   -H label -s '[*+]' -s '[*-]' -s '[1*,2*]' -s '[2*,3*]' \
   -u 0 -p 2d -p 3d -r -A D -A I -A C \
   -i ignore_bad_chiral -i sdf "$in" > "$out" 2>> err.txt
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
