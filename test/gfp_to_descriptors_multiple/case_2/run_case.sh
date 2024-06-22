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

command=$BIN_DIR/gfp_to_descriptors_multiple
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

stdout='stdout'
stderr='stderr'

input="${LILLYMOL_HOME}/data/chembl10.gfp"
golden='out/chembl10.txt.gz'
unzipped='unzipped'
gunzip --stdout ${golden} > ${unzipped}
golden=${unzipped}

diff_tool='../../same_bits.py'
${command} ${input} > ${stdout} 2> ${stderr}
$diff_tool ${stdout} ${unzipped}

if [[ $? -eq 0 ]] ; then
  echo "${case_id} : TEST PASS"
else
  echo "${case_id} : TEST FAIL"
  head ${stderr}
fi

rm ${stdout} ${stderr} ${unzipped}
