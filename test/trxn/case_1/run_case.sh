#! /bin/bash

if [[ -z "$BIN_DIR" ]]
then
# undefined BIN_DIR
    BIN_DIR="../../../bin/Linux-gcc-7.2.1"
fi
command=$BIN_DIR/trxn
case_id="Case 1"

echo "Testing:  $command"

name1=output.smi
name1_out=out/output.smi
diff_tool=../../fileDiff.sh
$command -v -r in/1.2.1_Aldehyde_reductive_amination_FROM_amines_AND_aldehydes.rxn -Z -z i -M RMX -m RMX -S output in/20180412_amines.smi in/20180412_aldehydes.smi>>log.txt 2>>err.txt

$diff_tool $name1 $name1_out
ret1=$?
if [ $ret1 -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm log.txt
rm err.txt
