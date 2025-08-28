#! /bin/bash
declare -i failures=0
for dir in ./*; do
    if [ -d "$dir" ]; then
        if [ -x "${dir}/run_case.sh" ]
        then
            cd ${dir} && ./run_case.sh
            echo "Return code $?"
            if [[ $? -ne 0 ]] ; then
              failures=$((failures+1))
            fi
        fi
    fi
done

exit ${failures}
