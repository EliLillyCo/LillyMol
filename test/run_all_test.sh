#! /bin/bash
# BIN_DIR is required to be exported before running test
# export BIN_DIR=../../../../../bin/Linux-gcc-7.2.1
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
