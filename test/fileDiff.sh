#! /bin/bash

in_name=$1
out_name=$2

if [[ -z "${in_name}" || -z "${out_name}" ]] ; then
  echo "one or more file(s) missing '${in_name}' '${out_name}'" >&1
  exit 1
fi

ret=0
#echo name1
#echo name1_out

if [[ ! -f ${in_name} ]] ; then
  echo "${in_name} not found" >&2
  exit 2
fi

if [[ ! -f ${out_name} ]] ; then
  echo "${out_name} not found" >&2
  exit 2
fi

result=$(diff -y -W 72 $in_name $out_name)

if [ $? -eq 0 ]
then
    #return 1 for file match
    ret=1
else
    line_count_1=$(wc -l < $in_name)
    #echo $line_count_1
    line_count_2=$(wc -l < $out_name)
    #echo $line_count_2
    if [ $line_count_1 = $line_count_2 ]
    then
        #return 2 for line count match
        #echo "line count match"
        ret=2
    fi
fi
#echo $ret
exit $ret
