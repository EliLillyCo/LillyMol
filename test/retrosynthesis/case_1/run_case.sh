#! /bin/bash

if [[ -z "$BIN_DIR" ]]
then
# undefined BIN_DIR
    BIN_DIR="../../../bin/Linux-gcc-7.2.1"
fi
command=$BIN_DIR/retrosynthesis
case_id="Case 1"
echo "Testing:  $command"

name1=log.txt
name1_out=out/log.txt
name2=err.txt
name2_out=out/err.txt
diff_tool=../../fileDiff.sh
$command -Y all -X kg -X kekule -X ersfrm -a 2 -q f -v -R 1 -I in/CentroidRxnSmi_1 -P UST:AZUCORS in/10Cmpds.smi>>log.txt 2>>err.txt
$diff_tool $name1 $name1_out
ret1=$?

$diff_tool $name2 $name2_out
ret2=$?

if [ $ret1 -eq 1 ] && [ $ret2 -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1 $name2

