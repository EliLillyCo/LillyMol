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

shared_data_dir=$LILLYMOL_HOME/data
command=$BIN_DIR/gfp_to_descriptors
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

stdout='stdout'
stderr='stderr'

golden='out/out.txt'
# Support linux and mac 
if [[ "$OSTYPE" == "linux-gnu" ]]; then
  golden=out/linux/out.txt
elif [[ "$OSTYPE" == "darwin"* ]]; then
  golden=out/osx/out.txt
else
  echo "OS is not supported"
fi

diff_tool=../../fileDiff.sh
same_bits=../../same_bits.py
${command} -b -n 50 -f -F FPDSC ${shared_data_dir}/pubchem.gfp >${stdout} 2>${stderr}
$diff_tool ${stdout} ${golden}

if [ $? -eq 1 ] ; then
  echo "$case_id : TEST PASS"
  rm ${stdout} ${stderr}
  exit
fi

$same_bits ${stdout} ${golden}
if [[ $? -eq 0 ]] ; then
  echo "$case_id : TEST PASS"
  rm ${stdout} ${stderr}
  exit
fi

echo "$case_id : TEST FAIL"

rm ${stdout} ${stderr}
