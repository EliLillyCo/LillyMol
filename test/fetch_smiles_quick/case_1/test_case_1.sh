#! /bin/bash

name1=notInRecord
name1_out=out/notInRecord

name2=notInIdentifier
name2_out=out/notInIdentifier

diff_tool=../../fileDiff.sh
../../../bin/Linux-gcc-6.2.0/fetch_smiles_quick -c 1 -C 2 -X notInRecord -Y notInIdentifier in/record.w  in/identifier.smi >>log.txt 2>>err.txt
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
rm log.txt
rm err.txt
