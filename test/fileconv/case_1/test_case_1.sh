#! /bin/bash

name1=selection.smi
name1_out=out/selection.smi
diff_tool=../../fileDiff.sh
../../../bin/Linux-gcc-6.2.0/fileconv -F 6 -c 4 -C 14 -v -i smi in/list.smi -S selection >>log.txt 2>>err.txt
$diff_tool $name1 $name1_out
ret=$?
if [ $ret -eq 1 ]
then
        echo "TEST PASS"
else
        echo "TEST FAIL"
fi
rm $name1
rm log.txt
rm err.txt
