#! /bin/bash

name1=output.smi
name1_out=out/output.smi
diff_tool=../../fileDiff.sh

../../../bin/Linux-gcc-6.2.0/trxn -v -r in/1.2.1_Aldehyde_reductive_amination_FROM_amines_AND_aldehydes.rxn -Z -z i -M RMX -m RMX -S output in/20180412_amines.smi in/20180412_aldehydes.smi>>log.txt 2>>err.txt

$diff_tool $name1 $name1_out
ret1=$?
if [ $ret1 -eq 1 ]
then
        echo "TEST PASS"
else
        echo "TEST FAIL"
fi
rm $name1
rm log.txt
rm err.txt
