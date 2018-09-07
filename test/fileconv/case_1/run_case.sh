#! /bin/bash

if [[ -z "$BIN_DIR" ]]
then
# undefined BIN_DIR
    BIN_DIR="../../../bin/Linux-gcc-7.2.1"
fi
command=$BIN_DIR/fileconv
case_id="Case 1"
echo "Testing:  $command"

name1=selection.smi
name1_out=out/selection.smi
diff_tool=../../fileDiff.sh
$command -F 6 -c 4 -C 14 -v -i smi in/list.smi -S selection >>log.txt 2>>err.txt
$diff_tool $name1 $name1_out
ret=$?
if [ $ret -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm log.txt
rm err.txt
