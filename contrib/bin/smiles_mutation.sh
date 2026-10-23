#!/usr/bin/env bash
here=$(dirname $(realpath $0))

if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $(dirname ${here}))
fi

libfile="${LILLYMOL_HOME}/data/smiles_mutation_library"

exec ${LILLYMOL_HOME}/bin/$(uname)/smiles_mutation -S ${libfile} "$@"
