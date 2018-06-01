#! /bin/bash

name1=log.txt
name1_out=out/log.txt
name2=err.txt
name2_out=out/err.txt
diff_tool=../../fileDiff.sh
../../../bin/Linux-gcc-6.2.0/retrosynthesis -Y all -X kg -X kekule -X ersfrm -a 2 -q f -v -R 1 -I in/CentroidRxnSmi_1 -P UST:AZUCORS -M ncon -M ring -M unsat -M arom in/10Cmpds.smi>>log.txt 2>>err.txt
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
rm $name1 $name2

