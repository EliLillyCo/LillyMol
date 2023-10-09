#!/bin/bash

# The directory in which the queries are found
dir=$(dirname ${0%%.sh})
program="${dir}/../../../bin/Linux/pubchem_fingerprints"

query_dir="${dir}/../../data/pubchem_fingerprints"
options="-q section4=${query_dir}//section4.smt "`
        `"-q section5=${query_dir}/section5.smt "`
        `"-q section6=${query_dir}/section6.smt "`
        `"-q section7=${query_dir}/section7.smt"

    exec ${program} ${options} "$@"
