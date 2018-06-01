#! /bin/bash

name1=hits.smi
name1_out=out/hits.smi
name2=nonhits.smi
name2_out=out/nonhits.smi

diff_tool=../../fileDiff.sh
../../../bin/Linux-gcc-6.2.0/tsubstructure -s 'C(C)(=O)C' -m hits.smi -n nonhits.smi in/list.smi >>log.txt 2>>err.txt
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
