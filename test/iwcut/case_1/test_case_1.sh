#! /bin/bash

name1=log.txt
name1_out=out/log.txt
diff_tool=../../fileDiff.sh
../../../bin/Linux-gcc-6.2.0/iwcut -f 5,3 in/input.txt >>log.txt 2>>err.txt
$diff_tool $name1 $name1_out
ret=$?
if [ $ret -eq 1 ]
then
        echo "TEST PASS"
else
        echo "TEST FAIL"
fi
rm $name1
rm err.txt
