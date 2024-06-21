#!/bin/bash

if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]; then 
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME (local path to LillyMol code)"
    echo "Please export BUILD_DIR (the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
else
    BIN_DIR="$LILLYMOL_HOME/contrib/python/mmp"
fi

test_command=getMMPSeriesSuggestfromCSV
case=case_2
case_id="Case 2"

test_top="$LILLYMOL_HOME/test"
test_cmd_top="$test_top/$test_command"

diff_tool=../../fileDiff.sh

command="$BIN_DIR/$test_command.py"

if [ ! -x "$command" ]; then
    echo "$command is not executable or does not exist"
    exit 1
fi

in="${test_cmd_top}/${case}/in/dalkeJCM2018_ChEMBL_hERG_series.csv"
gold_out="${test_cmd_top}/${case}/out/dalke_case2.series"
stdout='dalke_case2.series'
stderr='stderr'

echo "Testing: $command"

${command} -i "${in}" -o "${stdout}" -s SMILES -n CHEMBL_ID -a hERG_pIC50 -p 5 -q 9 -y 2> ${stderr}
$diff_tool "${stdout}" "${gold_out}"

if [ $? -eq 1 ]
then
    echo "${case_id} : TEST PASS"
else
    echo "${case_id} : TEST FAIL"
    diff ${stdout} ${gold_out}
    cat ${stderr}
fi

rm -f ${stdout} ${stderr}
