#! /bin/bash

if [[ -z "$BIN_DIR" ]]
then
# undefined BIN_DIR
    BIN_DIR="../../../bin/Linux-gcc-7.2.1"
fi
command=$BIN_DIR/rxn_signature
case_id="Case 1"
echo "Testing:  $command"

name1=all.sig
name1_out=out/all.sig
name2=all.log
name2_out=out/all.log
name3=Cfile
name3_out=out/Cfile
name4=Ffile
name4_out=out/Ffile
diff_tool=../../fileDiff.sh
$command -v -r 0,1,2 -C Cfile -F Ffile in/all.rxnsmi >all.sig 2>all.log

$diff_tool $name1 $name1_out
ret1=$?
$diff_tool $name2 $name2_out
ret2=$?
$diff_tool $name3 $name3_out
ret3=$?
$diff_tool $name4 $name4_out
ret4=$?

if [ $ret1 -eq 1 ] && [ $ret2 -eq 1 ] && [ $ret3 -eq 1 ] && [ $ret4 -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi
rm $name1
rm all.log
rm Cfile
rm Ffile
