#!/bin/bash

if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]; then 
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME (local path to LillyMol code)"
    echo "Please export BUILD_DIR (the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
else
    BIN_DIR="$LILLYMOL_HOME/bin/$BUILD_DIR"
fi

test_command=random_molecular_permutations
case=case_1
case_id="Case 1"

test_top="$LILLYMOL_HOME/test"
test_cmd_top="$test_top/$test_command"

diff_tool='../../fileDiff.sh'
command="$BIN_DIR/$test_command"

if [ ! -x "${command}" ] ; then
  echo "Executable ${command} not found"
fi

in="${test_cmd_top}/${case}/in/test.smi"
golden="${test_cmd_top}/${case}/out/test.out"

stdout='stdout'
stderr='stderr'

echo "Testing: ${command}"

# take 10 copies, do 100 mutations, max ring size 8, 10-40 atoms,
# max num rings 6, remove duplicates, use seed for deterministic test
${command} -x 10 -p 100 -c 10 -c 40 -R 8 -Y 7 -U all -s 1234567 "$in" > "${stdout}" 2> ${stderr}
${diff_tool} ${stdout} ${golden}

if [ $? -eq 1 ]
then
  echo "$case_id : TEST PASS"
else
  echo "$case_id : TEST FAIL"
  cat ${stderr}
fi

rm ${stdout} ${stderr}
