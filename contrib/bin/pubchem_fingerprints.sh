#!/bin/bash

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $(readlink -e $0))))
fi

# The directory in which the queries are found

query_dir="${LILLYMOL_HOME}/data/pubchem_fingerprints"

options="-q section4=${query_dir}//section4.smt "`
        `"-q section5=${query_dir}/section5.smt "`
        `"-q section6=${query_dir}/section6.smt "`
        `"-q section7=${query_dir}/section7.smt"

exe=${LILLYMOL_HOME}/bin/$(uname)/pubchem_fingerprints
exec ${exe} ${options} "$@"
