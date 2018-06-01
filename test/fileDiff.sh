#! /bin/bash

in_name=$1
out_name=$2
ret=0
#echo name1
#echo name1_out
result=$(diff -y -W 72 $in_name $out_name)
if [ $? -eq 0 ]
then
        ret=1
fi
#echo $ret
exit $ret
