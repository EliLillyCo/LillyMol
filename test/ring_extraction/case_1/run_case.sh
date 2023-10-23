#! /bin/bash

if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]
then 
    # undefined BIN_DIR
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME(local path to LillyMol code)"
    echo "Please export BUILD_DIR(the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
else
    BIN_DIR=$LILLYMOL_HOME/bin/$BUILD_DIR
fi

command=$BIN_DIR/ring_extraction
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

stdout='stdout'
stderr='stderr'
diff_tool=../../fileDiff.sh

pid=${$}
$command -X Ar:Al -R 7 -k -c -v -S /tmp/${pid}ring in/rings.smi >${stdout} 2>${stderr}

status='PASS'
for file in /tmp/${pid}ring_*smi ; do
  diff -q $file out/$(basename ${file/${pid}/})
  if [[ $? -ne 0 ]] ; then
    status='FAIL'
  fi
done
echo "$case_id : TEST ${status}"

rm ${stdout} ${stderr}
rm /tmp/${pid}ring_*smi
