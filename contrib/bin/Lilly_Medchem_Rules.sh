#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  here=$(realpath $0)
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname ${here})))
fi

here=$(realpath $(dirname $0))
exec ruby ${here}/Lilly_Medchem_Rules.rb "$@"
