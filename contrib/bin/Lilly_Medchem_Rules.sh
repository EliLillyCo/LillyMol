#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  here=$(readlink -f $0)
  echo ${here}
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname ${here})))
fi

exec ruby $(dirname ${here})/Lilly_Medchem_Rules.rb "$@"
