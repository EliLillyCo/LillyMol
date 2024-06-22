#! /bin/bash
echo "$(basename $(dirname $(readlink -f $0))) test disabled"
exit
for dir in ./*; do
    if [ -d "$dir" ]; then
        current_path= pwd &> /dev/null
        pushd "$dir" &> /dev/null
        if [ -f "run_case.sh" ]
        then
            ./run_case.sh
        fi
        popd $current_path &> /dev/null
    fi
done
