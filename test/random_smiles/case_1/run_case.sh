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

command=$BIN_DIR/random_smiles
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "${command}" ]
then
  echo "Executable ${command} not found"
  exit 1
fi

diff_tool='../../fileDiff.sh'

golden='out/output.smi'
stdout='stdout'
stderr='stderr'

${command} -n 5 -a -e -A D 'in/pubchem_example.smi' > ${stdout} 2> ${stderr}

${diff_tool} ${stdout} ${golden}

if [ $? == 2 ]
then
  echo "${case_id} : TEST PASS"
else
  echo "${case_id} : TEST FAIL"
fi
rm ${stdout} ${stderr}
