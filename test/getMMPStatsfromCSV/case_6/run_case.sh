#!/bin/bash

if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]; then 
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME (local path to LillyMol code)"
    echo "Please export BUILD_DIR (the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
else
    BIN_DIR="${LILLYMOL_HOME}/contrib/python/mmp"
fi

test_command='getMMPStatsfromCSV'
case='case_6'
case_id="Case 6"

test_top="${LILLYMOL_HOME}/test"
test_cmd_top="${test_top}/${test_command}"

diff_tool='../../fileDiff.sh'

command="${BIN_DIR}/$test_command.py"

if [ ! -x "${command}" ]; then
  echo "${command} is not executable or does not exist"
  exit 1
fi

in="${test_cmd_top}/${case}/in/dalkeJCM2018_ChEMBL_hERGsubset.csv"
gold_out="${test_cmd_top}/${case}/out/dalke_c6.pairs.sum"
stdout='dalke_c6.pairs'
stderr='stderr'

echo "Testing: ${command}"

${command} -i "${in}" -o "${stdout}" -s SMILES -n CHEMBL_ID -a hERG_pIC50 -A MEAN_DIFF_INVLOG -p BASIC 2> ${stderr}
${diff_tool} ${stdout}.sum ${gold_out}

if [ $? -eq 1 ]
then
  echo "${case_id} : TEST PASS"
else
  echo "${case_id} : TEST FAIL"
  diff ${stdout} ${gold_out} | head -n 50
  cat ${stderr}
fi

rm -f ${stdout}.sum ${stdout}.pairs ${stderr}
