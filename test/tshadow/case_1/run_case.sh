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

command=$BIN_DIR/tshadow
case_id="Case 1"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

stdout='stdout'
stderr='stderr'

# Support linux and mac 
if [[ "${UNAME}" == "Linux" ]]; then
    golden=out/${UNAME}/stdout
elif [[ "${UNAME}" == "darwin"* ]]; then
    golden=out/osx/stdout
else
    echo "${UNAME} is not supported"
    exit 1
fi

diff_tool=../../fileDiff.sh
$command -L -G in/in.sdf > stdout 2>stderr
$diff_tool ${stdout} ${golden}
ret1=$?

if [ $ret1 -eq 1 ]
then
        echo "$case_id : TEST PASS"
else
        echo "$case_id : TEST FAIL"
fi

rm ${stdout}
rm ${stderr}

