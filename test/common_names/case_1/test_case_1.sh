#! /bin/bash

name1=output.smi
name1_out=out/output.smi
diff_tool=../../fileDiff.sh
../../../bin/Linux-gcc-6.2.0/common_names in/input1.smi in/input2.smi -S output -s 10000 -D + -v >>log.txt 2>>err.txt
$diff_tool $name1 $name1_out
ret=$?
if [ $ret == 1 ]
then
        echo "TEST PASS"
else
        echo "TEST FAIL"
fi
rm $name1
rm log.txt
rm err.txt
