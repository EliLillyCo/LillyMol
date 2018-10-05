#! /bin/bash
# LILLYMOL_HOME and BUILD_DIR are required to be exported before running test
if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]
then 
    # undefined BIN_DIR
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME(local path to LillyMol code)"
    echo "Please export BUILD_DIR(the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
fi

for dir in ./*; do
    if [ -d "$dir" ]; then
        current_path= pwd &> /dev/null
        pushd "$dir" &> /dev/null
        if [ -f "run_test.sh" ]
        then
            ./run_test.sh
        fi
        popd $current_path &> /dev/null
    fi
done
