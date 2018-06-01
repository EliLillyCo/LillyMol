#! /bin/bash

name1=unique.smi
name1_out=out/unique.smi
name2=duplicate.smi
name2_out=out/duplicate.smi

diff_tool=../../fileDiff.sh
../../../bin/Linux-gcc-6.2.0/unique_molecules -S unique in/input.smi -D duplicate -v -l>>log.txt 2>>err.txt

$diff_tool $name1 $name1_out
ret1=$?

$diff_tool $name2 $name2_out
ret2=$?

if [ $ret1 -eq 1 ] && [ $ret2 -eq 1 ]
then
        echo "TEST PASS"
else
        echo "TEST FAIL"
fi
rm $name1
rm duplicate.smi
rm log.txt
rm err.txt
