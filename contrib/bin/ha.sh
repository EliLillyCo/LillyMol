#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  here=$(basename $0)
  export LILLYMOL_HOME=${here}/../..
fi

charges=${LILLYMOL_HOME}/data/queries/charges/queries
queries=${LILLYMOL_HOME}/data/queries/heteroatoms/unique_queries_tnass
exec "${LILLYMOL_HOME}/bin/$(uname)/tnass" -a -a -q F:${queries} $@
