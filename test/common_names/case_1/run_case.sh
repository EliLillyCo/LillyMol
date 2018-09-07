#! /bin/bash
if [[ -z "$BIN_DIR" ]]
then
# undefined BIN_DIR
    BIN_DIR="../../../bin/Linux-gcc-7.2.1"
fi
command=$BIN_DIR/common_names
case_id="Case 1"
echo "Testing:  $command"
name1=output.smi
name1_out=out/output.smi
diff_tool=../../fileDiff.sh
$command in/input1.smi in/input2.smi -S output -s 10000 -D + -v >>log.txt 2>>err.txt

$diff_tool $name1 $name1_out
ret=$?
if [ $ret == 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm log.txt
rm err.txt
