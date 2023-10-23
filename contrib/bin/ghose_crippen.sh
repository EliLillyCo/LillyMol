#!/bin/bash

# Set up queries for ghose_crippen logp calculation
# Note: the results are not very good. We don't use
# this, perhaps the program will be useful for some
# other additive contribution model.

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $0)))
fi

set -x
QUERIES=${LILLYMOL_HOME}/data/wildman_crippen.dat
CHARGES=${LILLYMOL_HOME}/data/queries/charges/queries

exe=${LILLYMOL_HOME}/bin/$(uname)/ghose_crippen
exec ${exe} -F ${QUERIES} -N F:${CHARGES} "$@"
