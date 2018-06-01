#! /bin/bash

name1=output.rxnsmi
name1_out=out/output.rxnsmi
diff_tool=../../fileDiff.sh
../../../bin/Linux-gcc-6.2.0/rxn_standardize -s -c -D x -X igbad -v -C 60 -K -E autocreate -e -o -m -I -b -f gsub  in/input.rxnsmi > output.rxnsmi 2>>err.txt

$diff_tool $name1 $name1_out
ret1=$?
if [ $ret1 -eq 1 ]
then
        echo "TEST PASS"
else
        echo "TEST FAIL"
fi

rm $name1
rm err.txt
