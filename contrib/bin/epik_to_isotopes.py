#!/bin/bash
${C3TK_HOME}/utils/telemetry.sh $0

if [[ ! -v C3TK_HOME || ! -d "${C3TK_HOME}"  ]]; then
    echo $(basename $0)": C3TK_HOME not set: did you forget to run 'module load gc3tk'?" >&2 && exit 1
fi

program=$(basename $0)
exec python ${C3TK_HOME}/bin/py/pytk/${program%%.sh}.py "$@"
