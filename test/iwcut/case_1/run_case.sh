#! /bin/bash
if [[ -z "$BIN_DIR" ]]
then
# undefined BIN_DIR
    BIN_DIR="../../../bin/Linux-gcc-7.2.1"
fi
command=$BIN_DIR/iwcut
case_id="Case 1"
echo "Testing:  $command"

name1=log.txt
name1_out=out/log.txt
diff_tool=../../fileDiff.sh
$command -f 5,3 in/input.txt >>log.txt 2>>err.txt
$diff_tool $name1 $name1_out
ret=$?
if [ $ret -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm err.txt
