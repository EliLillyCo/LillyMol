#! /bin/bash

if [[ -z "$BIN_DIR" ]]
then
# undefined BIN_DIR
    BIN_DIR="../../../bin/Linux-gcc-7.2.1"
fi
command=$BIN_DIR/rxn_standardize
case_id="Case 1"
echo "Testing:  $command"
name1=output.rxnsmi
name1_out=out/output.rxnsmi
diff_tool=../../fileDiff.sh
$command -s -c -D x -X igbad -v -C 60 -K -E autocreate -e -o -m -I -b -f gsub  in/input.rxnsmi > output.rxnsmi 2>>err.txt

$diff_tool $name1 $name1_out
ret1=$?
if [ $ret1 -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi

rm $name1
rm err.txt
