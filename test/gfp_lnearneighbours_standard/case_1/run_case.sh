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

command=$BIN_DIR/gfp_lnearneighbours_standard
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
  echo "Executable ${command} not found"
  exit 1
fi

diff_tool='../../fileDiff.sh'

golden='out/log.txt'
stdout='stdout'
stderr='stderr'

${command} -p "${LILLYMOL_HOME}/data/pubchem_example.gfp" -T 0.5 "${LILLYMOL_HOME}/data/pubchem_example.gfp" > ${stdout} 2> ${stderr}
${diff_tool} ${stdout} ${golden}

if [ $? -eq 1 ]
then
  echo "$case_id : TEST PASS"
else
  echo "$case_id : TEST FAIL"
  cat ${stderr}
fi

rm ${stdout} ${stderr}
