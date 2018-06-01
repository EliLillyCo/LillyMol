#! /bin/bash

name1=log.txt
name1_out=out/log.txt
name2=err.txt
name2_out=out/err.txt

diff_tool=../../fileDiff.sh
../../../bin/Linux-gcc-6.2.0/unique_rows -c 1 -c 2 in/input.dat >>log.txt 2>>err.txt
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

rm log.txt
rm err.txt
