#!/bin/bash

# Set up queries for ghose_crippen logp calculation
# Note: the results are not very good. We don't use
# this, perhaps the program will be useful for some
# other additive contribution model.

if [ ! -v LILLYMOL_HOME ] ; then
  echo "Must set shell variable LILLYMOL_HOME" >&2
  exit 1
fi

set -x
QUERIES=${LILLYMOL_HOME}/data/wildman_crippen.dat
CHARGES=${LILLYMOL_HOME}/contrib/data/queries/charges/queries

exec ${LILLYMOL_HOME}/bin/Linux/ghose_crippen -F ${QUERIES} -N F:${CHARGES} "$@"
