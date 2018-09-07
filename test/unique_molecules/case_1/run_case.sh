#! /bin/bash

if [[ -z "$BIN_DIR" ]]
then
# undefined BIN_DIR
    BIN_DIR="../../../bin/Linux-gcc-7.2.1"
fi
command=$BIN_DIR/unique_molecules
case_id="Case 1"
echo "Testing:  $command"
name1=unique.smi
name1_out=out/unique.smi
name2=duplicate.smi
name2_out=out/duplicate.smi

diff_tool=../../fileDiff.sh
$command -S unique in/input.smi -D duplicate -v -l>>log.txt 2>>err.txt

$diff_tool $name1 $name1_out
ret1=$?

$diff_tool $name2 $name2_out
ret2=$?

if [ $ret1 -eq 1 ] && [ $ret2 -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm duplicate.smi
rm log.txt
rm err.txt
