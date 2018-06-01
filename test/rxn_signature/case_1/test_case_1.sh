#! /bin/bash

name1=all.sig
name1_out=out/all.sig
name2=all.log
name2_out=out/all.log
name3=Cfile
name3_out=out/Cfile
name4=Ffile
name4_out=out/Ffile
diff_tool=../../fileDiff.sh

../../../bin/Linux-gcc-6.2.0/rxn_signature -v -r 0,1,2 -C Cfile -F Ffile in/all.rxnsmi >all.sig 2>all.log

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
        echo "TEST PASS"
else
        echo "TEST FAIL"
fi
rm $name1
rm all.log
rm Cfile
rm Ffile
